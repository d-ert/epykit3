# DSS replication of GSE263850 (paper-matched parameters)

Per-CpG: `DSS::DMLfit.multiFactor(smoothing=TRUE)` -> `DMLtest.multiFactor(coef = groupWT)`

DMR calling: `DSS::callDMR(p.threshold=1e-5, delta=0, minlen=50, minCG=3, dis.merge=100, pct.sig=0.5)`

Note: callDMR on multifactor output returns `chr, start, end, length, nCG, areaStat` — no `diff.Methy`. Per-DMR per-group mean methylation (and hence diff.Methy / dmr_type) is derived here directly from the 6 per-CpG BEDs.


## Pipeline

- DMRs called: **922**
- Hyper / hypo: **688 / 234** (74.6% hyper)
- Median length: **241 bp** (IQR 159-357 bp)
- Range: 51-2679 bp


## 100 kb DMR-gene linkage

- (DMR, gene) pairs: **1,649**
- Unique genes within 100 kb: **1,467**
- DMRs with >= 1 linked gene: **617 / 922**
- Genes per DMR: median 2, mean 2.7, max 70


## Resources

- Resume wall time: **341.9 s**
- Initial run (DMLfit smoothing+test) wall: **~2,044 s** (see resume_log.txt header for details).
- Peak RSS across the whole run: **14320 MB**
- Mean RSS: 8550 MB
- Peak USS: 14315 MB
- Peak CPU% (1 core = 100%; logical cores: 24): 980.0%
- Mean CPU%: 114.0%
- Cumulative R-process CPU time observed: **427 s**
- Peak threads in R tree: 13
- Total samples collected: 1544
