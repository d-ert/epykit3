# all_data_tsv -- every benchmark data table, as TSV

One browsable folder with the underlying data for **all four studies**.
Every `.parquet` under `benchmark/data/` is converted to `.tsv` here and
every `.csv` is copied as `.tsv`, so you can open any number behind any
figure or report table in Excel.

Regenerate: `uv run python benchmark/scripts/export_all_data_tsv.py`

## Layout

| Folder | Study | Key files |
|---|---|---|
| `study1_panel/` | Study 1 -- epykit vs 8 published baselines on the Piao 2021 simulated grid | `eval_summary_all_tools.tsv` (TPR/FPR/F1/AUROC per tool x coverage x replicate x effect-bin), `epykit_engine_timings.tsv` |
| `study1b_simulator/` | Held-out 21-seed simulator (intrinsic + threshold truth) | `epykit_seed_iqr.tsv`, `intrinsic_truth_iqr.tsv`, `dual_truth_iqr.tsv` (M3), per-seed variants, timings |
| `study2_head_to_head/` | Study 2 -- epykit vs methylKit, same grid | **shares `study1_panel/eval_summary_all_tools.tsv`** -- filter `tool` to `epykit_lr` / `methylkit` / `methylkit_tuned`. No separate file. |
| `study3_real/` | Study 3 -- GSE263850 real WGBS | DMR call sets (chain_merge, DSS, tile), gene links, annotation, concordance, panel-E; `figure_source_tsv/` holds the exact data plotted by F1-F10 |
| `dis_merge_sweep/` | Study 3 chain_merge dis.merge sweep | `sweep_summary.tsv`, per-setting DMR call sets (100-500 bp), recall/precision vs DSS-922 |
| `null_calibration/` | Null p-value calibration | `summary_all_engines.tsv`, GSE263850 lr per-partition + p-value histogram |
| `sep_threshold/` | lr+ separation-fallback ROC | simulator ROC + GSE263850 prevalence (0 candidate sites) |
| `bug_fix_audit/` | Pre/post bug-fix per-cell deltas | `bug_fix_deltas.tsv` |

## Which file feeds which figure / table

**Study 1 figures** (`figures/study1_simulated_allPackages/F1-F9`) and the
REPORT §1 tables: all from `study1_panel/eval_summary_all_tools.tsv`
(+ `epykit_engine_timings.tsv` for runtime). Columns: `tool, scenario,
parameter_value, meth_diff_bin, threshold_kind, threshold, tpr, fpr, f1,
auroc` plus Wilson/bootstrap CI columns.

**Study 2 figures** (`figures/study2_simulated_headToHead/F1-F9`) and
REPORT §2 tables: same `eval_summary_all_tools.tsv`, filtered to epykit
vs methylKit; runtimes in `study1b_simulator/external_timings_*` and
`epykit_engine_timings_per_seed.tsv`.

**Study 3 figures** (`figures/study3_real_GSE263850/three_way/F1-F10`):
see `study3_real/figure_source_tsv/` -- one TSV per figure with the exact
plotted values, plus its own README.

**Abstract / §3.4 held-out simulator** (21-seed median+IQR):
`study1b_simulator/intrinsic_truth_iqr.tsv` (+ `dual_truth_iqr.tsv` for
the intrinsic-vs-threshold M3 panel).

**§3.5 null calibration**: `null_calibration/gse263850_lr_summary.tsv`
+ `gse263850_lr_pvalue_histogram_50bins.tsv`.

## External inputs (NOT in this folder -- they live outside the repo)
- Paper Supp Table 5: `epykit2/GSE263850_RAW/Paper resources/DMR_total_list.xlsx`
- methylKit-tile real-data DMRs + step_benchmarks: `methylkit_realResults/...`

## See also
- `benchmark/paper_data/` -- the curated, paper-section-organised TSV
  mirror the manuscript cites (01_headline_piao ... 06_methodology).
- This folder is the *complete* dump; paper_data is the *curated* subset.
