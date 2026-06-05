# DSS replication of GSE263850 (paper-matched parameters)

Clean reproduction of the differential methylation analysis in
Farhangdoost et al. (2024) on GSE263850 using DSS — matching the paper's
exact `DMLfit.multiFactor(smoothing=TRUE)` + `callDMR` call with
parameters `p.threshold=1e-5, delta=0, minlen=50, minCG=3, dis.merge=100,
pct.sig=0.5`. Annotation uses the same UCSC refGene catalog HOMER ships,
reimplemented in R for portability.

Note: callDMR on multifactor DMLtest output returns
`[chr, start, end, length, nCG, areaStat]` — no `diff.Methy`. Per-DMR
per-group mean methylation (and hence `diff.Methy` / `dmr_type`) is
derived here by averaging per-sample (sumM / sumT) within each DMR and
then averaging within each group. This matches the paper's Table 5
`diff.meth_mean` column derivation.

## Files

| File | Description |
|---|---|
| `dmr_dss.tsv` / `dmr_dss.csv` | Annotated DMR table (chrom, start, end, length, n_cpgs, areaStat, dmr_id, meanMethy_treatment, meanMethy_control, diff.Methy, dmr_type, feature_type, feature_gene, nearest_tss_gene, nearest_tss_distance). |
| `dmr_dss_raw.tsv` | Raw `callDMR` output before any per-DMR mean-methylation derivation. |
| `dmr_gene_links_100kb.csv` | Long-form (DMR, gene) pairs where the gene's canonical TSS is within 100 kb of the DMR midpoint. |
| `dmltest_per_cpg.tsv.gz` | `DMLtest.multiFactor` output for every CpG (cached; used by the resume script). |
| `step_timings.tsv` / `step_timings_resume.tsv` | Per-step R-side wall / CPU / R-mem-peak. |
| `resources.csv` | Per-1s OS-level samples: RSS, USS, CPU%, threads. Two phases marked: `initial` (DMLfit+DMLtest) and `resume` (callDMR + downstream). |
| `resources.json` | Aggregated peak/mean/totals + host info. |
| `parameters.json` | DSS<->paper parameter mapping. |
| `summary.md` | Headline numbers + paper-coord comparison + resource summary. |
| `run_log.txt` | Concatenated log of both phases. |
| `dss_session_info.txt` | `sessionInfo()` for reproducibility. |
