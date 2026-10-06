"""Diagnose saved WGBS v4 errors using truth only for evaluation, never inference."""
import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd
import polars as pl
import pyarrow.parquet as pq
from epykit.dmc import _load_sample_chrom


def distances(pos, regions):
    if len(regions)==0:return np.full(len(pos),np.inf)
    regions=regions.sort_values('start');starts=regions.start.to_numpy();ends=regions.end.to_numpy()
    i=np.searchsorted(starts,pos,side='right')-1
    left=np.where(i>=0,np.maximum(pos-ends[np.maximum(i,0)],0),np.inf)
    right=np.where(i+1<len(starts),np.maximum(starts[np.minimum(i+1,len(starts)-1)]-pos,0),np.inf)
    return np.minimum(left,right)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--benchmark',type=Path,default=Path('/scratch/wgbs_benchmark_v4'))
    ap.add_argument('--out',type=Path,default=Path('benchmark/results/wgbs_calibration'));args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    summaries=[];errors=[];regions=[]
    for mode in ['signal','null']:
        dataset=args.benchmark/'data'/f'genome_autosomes_{mode}'
        result=args.benchmark/'results/genome_baseline/simulated'/dataset.name
        sheet=pd.read_csv(dataset/'samples.tsv',sep='\t')
        truth_regions=pd.read_csv(dataset/'dmr_truth.tsv',sep='\t')
        # Only significant DSS sites are retained from the large TSV.
        dss=pd.concat([c[c.qvalue<=.05] for c in pd.read_csv(result/'DSS/dml.tsv',sep='\t',usecols=['chrom','pos','qvalue'],chunksize=750000)],ignore_index=True)
        for path in sorted((dataset/'truth_parts').glob('*.parquet')):
            chrom=path.stem
            truth=pq.read_table(path,columns=['pos','control_mu','tau','region_id']).to_pandas().set_index('pos')
            raw=pl.read_parquet(result/'epykit/store/.cache/dmc/lr'/f'chrom={chrom}.parquet',columns=['pos','pvalue','qvalue','meth_diff'])
            pos=raw['pos'].to_numpy();f=raw.to_pandas().set_index('pos')
            f=f.join(truth,how='left',validate='one_to_one');assert f.region_id.notna().all()
            cov=np.zeros(len(pos),dtype=float)
            store=result/'epykit/store/.cache/raw'
            for sample in sheet.sample_id:
                _,coverage=_load_sample_chrom(store,chrom,sample,raw.select('pos'));cov+=coverage
            f['mean_coverage']=cov/len(sheet)
            f['positive']=f.region_id>0;f['called']=f.qvalue<=.05
            f['true_positive']=f.called&f.positive;f['false_positive']=f.called&~f.positive
            f['mu_bin']=pd.cut(f.control_mu,[0,.05,.2,.8,.95,1],include_lowest=True).astype(str)
            f['tau_bin']=pd.cut(f.tau,[0,.02,.05,.1,.2,.5,1,5],include_lowest=True).astype(str)
            f['coverage_bin']=pd.cut(f.mean_coverage,[0,5,10,20,40,80,np.inf],include_lowest=True).astype(str)
            f['raw_p_lt_05']=f.pvalue<.05
            for dim in ['mu_bin','tau_bin','coverage_bin']:
                g=f.groupby(dim,observed=True).agg(sites=('called','size'),positive_sites=('positive','sum'),calls=('called','sum'),tp=('true_positive','sum'),fp=('false_positive','sum'),raw_p_lt_05=('raw_p_lt_05','sum')).reset_index().rename(columns={dim:'stratum'})
                g['dataset']=mode;g['dimension']=dim;summaries.append(g)
            for tool,called_pos in [('epykit',f.index[f.called].to_numpy()),('DSS',dss.loc[dss.chrom==chrom,'pos'].to_numpy())]:
                annot=truth.reindex(called_pos);assert annot.region_id.notna().all()
                positive=annot.region_id.to_numpy()>0
                dist=distances(called_pos,truth_regions[truth_regions.chrom==chrom])
                false_dist=dist[~positive]
                for label,select in [('0-100bp',false_dist<=100),('101-500bp',(false_dist>100)&(false_dist<=500)),('501-1000bp',(false_dist>500)&(false_dist<=1000)),('>1000bp',false_dist>1000)]:
                    errors.append(dict(dataset=mode,tool=tool,distance=label,false_positive_sites=int(select.sum()),true_positive_sites=int(positive.sum()) if label=='0-100bp' else 0))
        calls=pd.read_csv(result/'epykit/dmr.tsv',sep='\t')
        regions.append(dict(dataset=mode,calls=len(calls),mean_effect_below_10pct=int((calls.mean_meth_diff.abs()<.1).sum()),median_cpgs=float(calls.n_cpgs.median()),median_length_bp=float((calls.end-calls.start+1).median())))
        print(f'Diagnosed {mode}',flush=True)
    grouped=pd.concat(summaries).groupby(['dataset','dimension','stratum'],as_index=False).sum(numeric_only=True)
    grouped['recall']=grouped.tp/grouped.positive_sites.replace(0,np.nan)
    grouped['false_positive_site_rate']=grouped.fp/(grouped.sites-grouped.positive_sites).replace(0,np.nan)
    grouped.to_csv(args.out/'cpg_strata.csv',index=False)
    error=pd.DataFrame(errors).groupby(['dataset','tool','distance'],as_index=False).sum(numeric_only=True)
    error.to_csv(args.out/'false_positive_distances.csv',index=False)
    pd.DataFrame(regions).to_csv(args.out/'region_diagnostics.csv',index=False)
    print(error.to_string(index=False))

if __name__=='__main__':main()
