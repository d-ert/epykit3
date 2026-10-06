"""Compare existing CpG variance options on the development null chromosome 22."""
import json
from pathlib import Path
import numpy as np
import polars as pl
from epykit.dmc import process_chromosomes_dmc
ROOT=Path('/scratch/wgbs_benchmark_v4')
OUT=Path('benchmark/results/wgbs_calibration/depth_check');OUT.mkdir(exist_ok=False)
base=ROOT/'results/genome_baseline/simulated/genome_autosomes_null/epykit'
reference=pl.read_parquet(base/'store/.cache/dmc/lr/chrom=chr22.parquet',columns=['pos','pvalue'])
rows=[]
for dispersion in ['eb','bb_eb']:
 result=process_chromosomes_dmc(methylstore_path=str(base/'store/.cache/raw'),
    samples_treatment=[f'treatment_{i}' for i in range(1,6)],samples_control=[f'control_{i}' for i in range(1,6)],
    chromosomes=['chr22'],out_dir=str(OUT/dispersion),test='lr',unite=True,
    dispersion=dispersion,reference='adaptive',smoothing=False,backend='sequential',return_store=False)
 assert result['pos'].to_list()==reference['pos'].to_list()
 p=result['pvalue'].to_numpy()
 if dispersion=='eb':np.testing.assert_allclose(p,reference['pvalue'].to_numpy(),equal_nan=True)
 rows.append(dict(dispersion=dispersion,sites=len(p),**{str(t):float(np.mean(p<=t)) for t in [.05,.01,.001,.0001,.00001]}))
(OUT/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
