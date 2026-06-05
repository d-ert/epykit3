# Chain-merge replication of GSE263850 (paper-matched parameters)

Per-CpG: `ep.tl.dmc(test='lr', dispersion='site', smoothing=True, smoothing_span_bp=500)`

DMR aggregation: `ep.tl.dmr(method='chain_merge', alpha=1e-5, min_abs_meth_diff=0, minlen_bp=50, min_cpgs=3, dis_merge_bp=100, pct_sig=0.5)`

Annotation: HOMER-equivalent (UCSC refGene) via `ep.tl.annotate(refgene=...)`


## Pipeline

- CpGs tested: **21,993,377**
- Significant DMCs (q < 0.05): **17,326**
- DMRs called (chain_merge, alpha=1e-5, BH q < 0.05): **852**
- DMC runtime: 72.5 s
- DMR runtime: 0.0 s


## DMR morphology

- Hyper / hypo: **673 / 179** (79.0% hyper)
- Median length: **125 bp** (IQR 89-185 bp)
- Mean length: 152 bp
- Range: 50-1851 bp


## 100 kb DMR-gene linkage

- (DMR, gene) pairs: **1,501**
- Unique genes within 100 kb of any DMR: **1,290**
- DMRs with >= 1 linked gene: **576 / 852**
- Genes per DMR: median 2, mean 2.6, max 45
