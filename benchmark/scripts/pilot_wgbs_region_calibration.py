"""Development chr22 test of complete-pipeline label randomization on frozen counts."""
import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import polars as pl
import pyarrow.parquet as pq
import epykit as ep
from epykit.dmr import (call_dmr_chain_merge,apply_region_qfilter,_chain_merge_perm_survivors,
    _permutation_assignment,_is_self_or_mirror_perm,_run_permutations,_aggregate_region_perm_results)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--benchmark',type=Path,default=Path('/scratch/wgbs_benchmark_v4'))
    ap.add_argument('--out',type=Path,default=Path('benchmark/results/wgbs_calibration/chr22_pilot'))
    ap.add_argument('--permutations',type=int,default=99);ap.add_argument('--jobs',type=int,default=4)
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(args.benchmark/'src'))
    from wgbs_v3.score import score
    rows=[]
    for mode in ['signal','null']:
        original=args.benchmark/'data'/f'genome_autosomes_{mode}'
        prior=args.benchmark/'results/genome_baseline/simulated'/original.name/'epykit'
        target=args.out/mode;target.mkdir()
        sheet=pd.read_csv(original/'samples.tsv',sep='\t')
        obs=pl.from_pandas(sheet.assign(treatment=(sheet.group=='treatment').astype(int)))
        store=target/'raw';store.mkdir()
        for sample in sheet.sample_id:
            dest=store/f'sample={sample}/chrom=chr22';dest.mkdir(parents=True)
            shutil.copy2(prior/'store/.cache/raw'/f'sample={sample}/chrom=chr22/part-0.parquet',dest/'part-0.parquet')
        md=ep.MethylData(obs=obs,store=str(store.resolve()),assembly='hg38',analysis_root=str(target.resolve()))
        ep.pp.set_unite_type(md,type='intersect')
        ep.tl.dmc(md,test='lr',dispersion='eb',reference='adaptive',backend='sequential',smoothing=False,chromosomes=['chr22'])
        cm=dict(alpha=.05,min_abs_meth_diff=.1,min_cpgs=5,minlen_bp=50,dis_merge_bp=500,pct_sig=.5,use_q_for_sig=False)
        observed=apply_region_qfilter(call_dmr_chain_merge(md.dmc,**cm),.05)
        print(mode,'observed candidates',len(observed),flush=True)
        dmc=dict(test='lr',unite=True,min_samples_treatment=0,min_samples_control=0,dispersion='eb',reference='adaptive',smoothing=False,backend='sequential')
        def run_one(i):
            treatment,control=_permutation_assignment(i,seed=20261006,samples_treatment=md.treatment_ids,samples_control=md.control_ids,empirical_strata=None)
            is_self=_is_self_or_mirror_perm(perm_treatment=treatment,observed_treatment=md.treatment_ids,observed_control=md.control_ids)
            scores=_chain_merge_perm_survivors(methylstore_path=str(store.resolve()),samples_treatment=treatment,samples_control=control,chromosomes=['chr22'],dmc_kwargs=dmc,dmc_fdr_method='fdr_bh',chain_merge_kwargs=cm,min_mean_qvalue=.05)
            return is_self,scores
        scans=_run_permutations(run_one,args.permutations,args.jobs)
        (target/'permutation_scores.json').write_text(json.dumps([dict(self_or_mirror=flag,pvalues=None if arr is None else arr.tolist()) for flag,arr in scans])+'\n')
        # Evaluation-only chr22 view. Inference never reads latent truth.
        view=target/'evaluation';view.mkdir()
        for name in ['truth_parts','eligibility_parts']:
            (view/name).mkdir();shutil.copy2(original/name/'chr22.parquet',view/name/'chr22.parquet')
        tr=pd.read_csv(original/'dmr_truth.tsv',sep='\t');tr=tr[tr.chrom=='chr22'];tr.to_csv(view/'dmr_truth.tsv',sep='\t',index=False)
        sheet.to_csv(view/'samples.tsv',sep='\t',index=False)
        truth_frame=pq.read_table(view/'truth_parts/chr22.parquet',columns=['is_dml']).to_pandas()
        eligible_frame=pq.read_table(view/'eligibility_parts/chr22.parquet',columns=['eligible']).to_pandas()
        metadata=json.loads((original/'manifest.json').read_text());metadata.update(chromosomes=['chr22'],n_sites=len(truth_frame),n_regions=len(tr),n_positive_cpg=int(truth_frame.is_dml.sum()),n_common_eligible=int(eligible_frame.eligible.sum()),evaluation_view_of=str(original),evaluation_note='Chr22-only development view; sample sheet identifies original full-genome count sources.')
        (view/'manifest.json').write_text(json.dumps(metadata)+'\n')
        for method in ['asymptotic','region','max_t']:
            frame=observed
            if method!='asymptotic':
                p,q,fdp=_aggregate_region_perm_results(observed_pvalues=observed['combined_pvalue'].to_numpy(),results=scans,n_perm=args.permutations,fdr_method=method)
                frame=frame.with_columns(pl.Series('empirical_pvalue',p),pl.Series('empirical_qvalue',q)).filter(pl.col('empirical_qvalue')<=.05)
            destination=target/method;destination.mkdir()
            frame.with_columns((pl.col('end')-1).alias('end')).write_csv(destination/'dmr.tsv',separator='\t')
            metrics=score(view,None,destination/'dmr.tsv',destination/'score.json')['dmr']
            rows.append(dict(dataset=mode,method=method,**metrics))
        print(mode,'complete',flush=True)
    pd.DataFrame(rows).to_csv(args.out/'comparison.csv',index=False)
    print(pd.DataFrame(rows)[['dataset','method','called_regions','matched_regions','region_precision','region_recall','region_f1']].to_string(index=False))

if __name__=='__main__':main()
