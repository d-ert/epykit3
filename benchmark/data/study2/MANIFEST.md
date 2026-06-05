# Study 2 — head-to-head DMC results (auditability bundle)

These files back **Table 1** (DMC × coverage) and **Table 2** (DMC × replicate
count) of [`paper.md`](../../paper/paper.md) §3.2 and the corresponding
[`REPORT.md`](../../paper/report/REPORT.md) §2. They were added so the
head-to-head numbers — in particular the n = 2 headline — are reproducible from
the committed bundle, not only from the analysis host.

Both tools were run on the same simulated inputs (the held-out Piao
re-implementation, coverage/replicate grid; see Methods §2.1, §2.7) under the
same evaluation harness.

## Files

| Path | Content |
|---|---|
| `methylkit_results/dmc_cov{5,10,15,20,25}.tsv.gz` | methylKit `calculateDiffMeth` per-CpG output, coverage sweep (3 v 3). Columns: `chr, start, end, strand, pvalue, qvalue, meth.diff` (1-based; `meth.diff` on the −100..+100 percent scale). |
| `methylkit_results/dmc_rep{2,4,6,8,10}.tsv.gz` | Same, replicate-count sweep (n = 2…10 total) at 10× coverage. `dmc_rep2.tsv.gz` is the n = 2 headline. |
| `methylkit_results/timings.tsv` | Per-scenario wall-clock / CPU / RSS (backs §3.2.4 Table 4). |
| `epykit_results/dmc_{scenario}_lr.parquet` | epykit `lr` (default quasi-binomial) per-CpG DMC output for the matching scenario. Columns: `chrom, pos, n_case, n_control, mean_beta_case, mean_beta_control, pvalue, log2_odds_ratio, meth_diff, meth_diff_ci_*, qvalue, reject` (0-based `pos`; `meth_diff` fractional). |
| `ground_truth/dmc_truth.parquet` | Per-CpG truth: `chrom, pos, true_meth_diff, is_dmc, direction, meth_diff_bin`. 100,000 CpGs, 19,999 true DMCs (20.0 %). Shared with Study 1 (same simulator). |

Gzip is used for the methylKit TSVs (≈ 6× smaller; ~16 MB vs ~96 MB raw);
`zcat`/polars read them natively. The ~1.2 GB DMR-simulation per-CpG TSVs
(`dmr_cov*_dmc.tsv`) are **not** bundled — they are large and reproducible from
the committed simulator; Table 3 (DMR detection) is scored from them on the run
host.

## Significance cut (both tools)

`qvalue < 0.05` AND `|meth_diff| ≥ 0.25` (= `|meth.diff| ≥ 25` on methylKit's
percent scale).

## Reproduce the n = 2 headline

```python
import polars as pl
t = pl.read_parquet("ground_truth/dmc_truth.parquet")
n_true = int(t["is_dmc"].cast(bool).sum())                       # 19999

e = pl.read_parquet("epykit_results/dmc_rep2_lr.parquet")
e = e.with_columns(((pl.col("qvalue") < 0.05) & (pl.col("meth_diff").abs() >= 0.25)).alias("sig"))
j = e.join(t.select(["pos", "is_dmc"]), on="pos", how="inner")
print("epykit lr n=2 n_sig:", int(e["sig"].sum()))               # 11283
print("epykit lr n=2 TPR :", j.filter(pl.col("sig") & pl.col("is_dmc")).height / n_true)  # 0.564

m = pl.read_csv("methylkit_results/dmc_rep2.tsv.gz", separator="\t")
m = m.with_columns(((pl.col("qvalue") < 0.05) & (pl.col("meth.diff").abs() >= 25)).alias("sig"))
print("methylKit n=2 n_sig:", int(m["sig"].sum()))               # 6030
```

These reproduce Table 2's n = 2 row (epykit TPR 0.564 / n_sig 11,283;
methylKit n_sig 6,030) exactly.
