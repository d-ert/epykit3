# Paper-ready data export

Curated TSV + markdown export of every comparison that ships with the epykit Phase 4 paper. All parquets have been converted to TSV for easy inspection in Excel / R / Python; original .parquet sources stay in `benchmark/data/` for byte-exact reproducibility.

Generated 2026-06-01 via the converter script — see `06_methodology/` for documentation of how each comparison was produced.

## Contents

### `01_headline_piao/` — Piao-as-distributed (Studies 1+2)
The headline benchmark cell sweep that the paper's main Tables 1 and 2 cite.

| File | Rows × Cols | What |
|---|---:|---|
| `eval_summary_post_phase3.tsv` | 563 × 25 | TPR / FPR / F1 / AUROC per (tool, scenario, parameter_value, threshold, threshold_kind, meth_diff_bin) for 15 tools (epykit_lr, epykit_lrplus, epykit_welch_t, epykit_fisher, 4 epykit_dmr_* methods, methylkit, methylkit_tuned, dss, biseq, methylsig, radmeth, fisher) on Piao at 15 cells (5 dmc_coverage + 5 dmc_replicate + 5 dmr_coverage). Includes Wilson CIs on TPR/FPR. |
| `timings_post_phase3_epykit.tsv` | 57 × 9 | Per-engine, per-cell wallclock for all epykit engines on the headline sweep, fresh same-machine-state regeneration. |
| `methylkit_piao_cov10.timing.tsv` | 3 × 3 | Phase breakdown (read / unite / diffmeth) for methylKit on the headline cell (Piao cov=10). |
| `dss_piao_cov10.timing.tsv` | 3 × 3 | Same for DSS with `smoothing=FALSE`. |
| `dss_piao_cov10_smoothed.timing.tsv` | 3 × 3 | Same for DSS with `smoothing=TRUE`. |

### `02_multiseed_simulator/` — held-out intrinsic-truth simulator (20 seeds)
Median + IQR aggregations for the supplementary tables on simulator variance.

| File | Rows × Cols | What |
|---|---:|---|
| `eval_per_seed_epykit.tsv` | 80 × 30 | epykit (4 engines) accuracy + wallclock per seed (20 seeds × 4 engines) at cov=10. |
| `eval_seed_iqr_epykit.tsv` | 4 × 14 | Median + IQR across 20 seeds for the 4 epykit engines. |
| `eval_per_seed_external.tsv` | 540 × 18 | methylKit + DSS (smoothing TRUE and FALSE) accuracy on the full scoring grid × 20 seeds. |
| `eval_seed_iqr_external.tsv` | 3 × 15 | Median + IQR across 20 seeds for methylKit + DSS variants. |
| `timings_per_seed_external.tsv` | 60 × 3 | methylKit + DSS wallclock per seed (sum of phases). |
| `timings_iqr_external.tsv` | 3 × 7 | Median + IQR of methylKit + DSS wallclock across 20 seeds. |
| `timings_per_phase_per_seed_external.tsv` | 180 × 5 | Per-phase breakdown (read / fit / test) for methylKit + DSS, every seed. |
| `single_seed_2026000_intrinsic_scoring_grid.tsv` | 27 × 18 | Full scoring grid (p-value thresholds × q-value × per-bin) for methylKit + DSS on seed=2026000 only — the original Task 5 single-seed result. |
| `parallel_column_summary.md` | — | Human-readable summary of the seven-tool comparison at headline cell of seed=2026000 with FDR breach annotations. |

### `03_gse_real_data/` — GSE263850 cross-tool agreement (Study 3)
Real-cohort concordance metrics that the paper cites for the "external validation" section.

| File | Rows × Cols | What |
|---|---:|---|
| `dmr_iou.tsv` | 3 × 12 | Pairwise IoU of DMR sets across epykit ∩ DSS ∩ methylKit on GSE263850. |
| `per_dmr_stat_concordance.tsv` | 2149 × 17 | Per-DMR statistical concordance — pairwise meth_diff and -log10(qvalue) correlations across the same three tools. |

### `04_null_calibration/` — observed FDR under label shuffling (Phase 4 Task 7)
The Table S-Calib that the paper's Methods §3.Y cites.

| File | Rows × Cols | What |
|---|---:|---|
| `summary.tsv` | 12 × 11 | 12 of 13 (engine, dataset, scenario) cells — observed FDR median + IQR + Wilson CI from k=10–20 label-shuffles per cell. `fisher@gse263850` deliberately omitted (compute infeasible). |
| `MANIFEST.txt` | — | Sweep date, engine versions, per-cell wallclock summary, methodology paragraph. |

### `05_bug_fix_audit/` — pre/post-fix per-cell delta (Phase 4 Task 8)
The Table S-Fix that the paper's Limitations §10.5 cites.

| File | Rows × Cols | What |
|---|---:|---|
| `bug_fix_deltas.tsv` | 619 × 7 | Per-cell delta between pre-Phase-1 baseline and post-Phase-3 freeze, attributed to specific P0/P1 fix commits via `Affects:` trailers. |
| `bug_fix_deltas.md` | — | Markdown formatting of the same data for paper insertion. |
| `commits.json` | — | The hand-curated commit manifest used by `bug_fix_audit.py` (includes `_original_body` and `_note` audit fields for every entry). |

### `06_methodology/` — write-ups
Reference documents the paper's Methods and Discussion cite.

| File | What |
|---|---|
| `timing-comparison.md` | The comprehensive seven-tool timing + accuracy write-up. Includes single-cell Piao numbers, multi-seed simulator IQRs, per-tool deep dives, FDR breach analysis, Pareto frontier, and reproduction commands. |
| `headline-engine-selection.md` | Rationale for which engines appear in the headline benchmark vs which are reserved for the null-calibration table only (the `epykit_glm` exclusion explanation). |

## Reading the files in Excel / R / Python

All TSVs are tab-separated, UTF-8, with headers in row 1. No quoting on string fields. Numeric columns may carry NaN values where the metric is undefined for that row (typically AUROC on per-bin rows where the binomial event count is too small).

```python
# Python / polars
import polars as pl
df = pl.read_csv("01_headline_piao/eval_summary_post_phase3.tsv", separator="\t")

# Python / pandas
import pandas as pd
df = pd.read_csv("01_headline_piao/eval_summary_post_phase3.tsv", sep="\t")

# R
df <- read.table("01_headline_piao/eval_summary_post_phase3.tsv",
                 sep="\t", header=TRUE, na.strings=c("", "NA"))
```

## Provenance

Every TSV in this directory was derived from a `.parquet` (or pre-existing `.tsv`) under `benchmark/data/` at commit reachable from the `p0-fixes` branch on `origin/`. The parquets are the authoritative source — these TSVs are convenience copies for inspection and Excel work.

If you regenerate any source parquet (e.g. by re-running `run_epykit_study1.py` or `run_external_simulator_sweep.py`), re-run the converter script to refresh this directory.
