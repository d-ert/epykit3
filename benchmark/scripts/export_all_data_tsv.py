"""Export EVERY committed benchmark data table (all studies) to TSV in one
browsable folder, so any number behind any figure or table can be
inspected manually in Excel.

Output root: benchmark/all_data_tsv/
  study1_panel/            Study 1 -- 8-tool panel on the Piao simulator
  study1b_simulator/       Multi-seed held-out simulator (epykit + external)
  study2_head_to_head/     Study 2 -- epykit vs methylKit (shares study1 grid)
  study3_real/             Study 3 -- GSE263850 real WGBS (+ figure source data)
  null_calibration/        Null p-value calibration
  sep_threshold/           sep_threshold ROC (lr+ separation fallback)
  dis_merge_sweep/         chain_merge dis.merge sweep (Study 3)
  bug_fix_audit/           Pre/post bug-fix deltas
  README.md                master map: file -> what it is -> which fig/table

Every .parquet under benchmark/data/ is converted to .tsv; every .csv is
copied as .tsv. One very large file (lr_pvalues.parquet, ~22M rows) is
summarised to a histogram instead of dumped in full.
"""

from __future__ import annotations
import sys
import io
import shutil
from pathlib import Path
import numpy as np
import polars as pl

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

BENCH = Path(__file__).resolve().parents[1]
DATA = BENCH / "data"
OUT = BENCH / "all_data_tsv"
THREE_WAY_TSV = BENCH / "figures" / "study3_real_GSE263850" / "three_way" / "source_data_tsv"

# (source under data/, destination subfolder, output filename)
# parquet -> tsv conversions and csv copies, grouped by study.
JOBS = [
    # ---- Study 1: 8-tool panel on the Piao simulator -------------------
    ("study1/eval_summary_post_phase3.parquet", "study1_panel", "eval_summary_all_tools.tsv"),
    ("study1/timings_table.csv", "study1_panel", "epykit_engine_timings.tsv"),

    # ---- Study 1b: multi-seed held-out simulator ----------------------
    ("study1b_simulator/eval_per_seed.parquet", "study1b_simulator", "epykit_per_seed.tsv"),
    ("study1b_simulator/eval_seed_iqr.parquet", "study1b_simulator", "epykit_seed_iqr.tsv"),
    ("study1b_simulator/eval_external_timings_per_seed.parquet", "study1b_simulator", "external_timings_per_seed.tsv"),
    ("study1b_simulator/eval_external_timings_iqr.parquet", "study1b_simulator", "external_timings_iqr.tsv"),
    ("study1b_simulator/eval_frozen_grid.parquet", "study1b_simulator", "frozen_grid_control.tsv"),
    ("study1b_simulator/eval_simulator_intrinsic.parquet", "study1b_simulator", "intrinsic_truth.tsv"),
    ("study1b_simulator/eval_simulator_intrinsic_iqr.parquet", "study1b_simulator", "intrinsic_truth_iqr.tsv"),
    ("study1b_simulator/eval_simulator_intrinsic_per_seed.parquet", "study1b_simulator", "intrinsic_truth_per_seed.tsv"),
    ("study1b_simulator/eval_simulator_intrinsic_truth_both_iqr.parquet", "study1b_simulator", "dual_truth_iqr.tsv"),
    ("study1b_simulator/eval_simulator_intrinsic_truth_both_per_seed.parquet", "study1b_simulator", "dual_truth_per_seed.tsv"),
    ("study1b_simulator/timings_simulator.parquet", "study1b_simulator", "epykit_engine_timings_per_seed.tsv"),

    # ---- Study 2: head-to-head (shares the study1 eval grid) ----------
    # The head-to-head numbers live in study1/eval_summary (tool column
    # includes methylkit / methylkit_tuned). A pointer file is written
    # in the README; no separate parquet exists.

    # ---- Study 3: real WGBS GSE263850 ---------------------------------
    ("study3/chain_merge/dmr_chain_merge.csv", "study3_real", "epykit_chain_merge_100_dmrs.tsv"),
    ("study3/chain_merge/dmr_gene_links_100kb.csv", "study3_real", "epykit_chain_merge_100_gene_links_100kb.tsv"),
    ("study3/dss/dmr_dss.csv", "study3_real", "dss_from_scratch_dmrs.tsv"),
    ("study3/dss/dmr_gene_links_100kb.csv", "study3_real", "dss_gene_links_100kb.tsv"),
    ("study3/dss/resources.csv", "study3_real", "dss_resource_samples.tsv"),
    ("study3/dmr_significant_lenient.csv", "study3_real", "epykit_tile_dmrs_lenient.tsv"),
    ("study3/comparisons/annotation_distribution.csv", "study3_real", "annotation_distribution.tsv"),
    ("study3/comparisons/per_dmr_stat_concordance.csv", "study3_real", "per_dmr_concordance.tsv"),
    ("study3/comparisons/epykit_vs_dss/panel_e_capture_dss.csv", "study3_real", "panel_e_capture_dss.tsv"),
    ("study3/comparisons/epykit_vs_dss/heatmap_gene_hits_dss.csv", "study3_real", "heatmap_gene_hits_dss.tsv"),
    # NOTE: these two use a DIFFERENT (all-q<0.05, 36,811-DMR) chain_merge
    # set, not the 852 paper-faithful headline set -- segregated under
    # _other_views/ so the main folder stays headline-consistent.
    ("study3/comparisons_post_phase3/dmr_iou.parquet", "study3_real/_other_views", "dmr_iou_ALLSIG_chainmerge36811.tsv"),
    ("study3/comparisons_post_phase3/per_dmr_stat_concordance.parquet", "study3_real/_other_views", "per_dmr_concordance_phase3_allpairs.tsv"),
    ("study3/samplesheet.csv", "study3_real", "samplesheet.tsv"),

    # ---- dis.merge sweep (Study 3) ------------------------------------
    ("multi_thread_and_chain_sweep/chain_merge_dis_merge_sweep/sweep_summary.csv", "dis_merge_sweep", "sweep_summary.tsv"),
    ("multi_thread_and_chain_sweep/chain_merge_dis_merge_sweep/dis_merge_vs_dss_sensitivity.csv", "dis_merge_sweep", "dis_merge_vs_dss_sensitivity.tsv"),
    ("multi_thread_and_chain_sweep/chain_merge_dis_merge_sweep/dis_merge_100/dmr.csv", "dis_merge_sweep", "dmrs_dis_merge_100.tsv"),
    ("multi_thread_and_chain_sweep/chain_merge_dis_merge_sweep/dis_merge_150/dmr.csv", "dis_merge_sweep", "dmrs_dis_merge_150.tsv"),
    ("multi_thread_and_chain_sweep/chain_merge_dis_merge_sweep/dis_merge_200/dmr.csv", "dis_merge_sweep", "dmrs_dis_merge_200.tsv"),
    ("multi_thread_and_chain_sweep/chain_merge_dis_merge_sweep/dis_merge_250/dmr.csv", "dis_merge_sweep", "dmrs_dis_merge_250.tsv"),
    ("multi_thread_and_chain_sweep/chain_merge_dis_merge_sweep/dis_merge_250/dmr_gene_links_100kb.csv", "dis_merge_sweep", "dmrs_dis_merge_250_gene_links_100kb.tsv"),
    ("multi_thread_and_chain_sweep/chain_merge_dis_merge_sweep/dis_merge_500/dmr.csv", "dis_merge_sweep", "dmrs_dis_merge_500.tsv"),

    # ---- Null calibration ---------------------------------------------
    ("null_calibration/summary.parquet", "null_calibration", "summary_all_engines.tsv"),
    ("null_calibration/gse263850/lr_summary.parquet", "null_calibration", "gse263850_lr_summary.tsv"),
    ("null_calibration/gse263850/lr.parquet", "null_calibration", "gse263850_lr_per_partition.tsv"),

    # ---- sep_threshold ROC --------------------------------------------
    ("sep_threshold_roc/sep_threshold_roc_simulator.csv", "sep_threshold", "roc_simulator.tsv"),
    ("sep_threshold_roc/sep_threshold_roc_summary.csv", "sep_threshold", "roc_summary.tsv"),
    ("sep_threshold_roc/sep_prevalence_gse263850.csv", "sep_threshold", "prevalence_gse263850.tsv"),

    # ---- Bug-fix audit ------------------------------------------------
    ("audit/bug_fix_deltas.parquet", "bug_fix_audit", "bug_fix_deltas.tsv"),
]


def convert(src: Path, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.suffix == ".parquet":
        df = pl.read_parquet(src)
    else:
        df = pl.read_csv(src, infer_schema_length=10000)
    df.write_csv(dest, separator="\t")
    return df.height, df.width


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Remove files that were renamed/moved in a prior layout so they don't
    # linger (overwrite-in-place; safe to re-run with the folder open).
    for stale in ("study3_real/dmr_iou_three_way.tsv",
                  "study3_real/per_dmr_concordance_phase3.tsv"):
        p = OUT / stale
        if p.exists():
            p.unlink()
    print(f"Exporting to {OUT}")
    n_ok = 0
    for rel, subdir, name in JOBS:
        src = DATA / rel
        if not src.exists():
            print(f"  SKIP (missing): {rel}")
            continue
        rows, cols = convert(src, OUT / subdir / name)
        print(f"  {subdir}/{name}: {rows}x{cols}")
        n_ok += 1

    # Giant null-calibration p-value vector -> histogram instead of full dump.
    pvals = DATA / "null_calibration" / "gse263850" / "lr_pvalues.parquet"
    if pvals.exists():
        s = pl.read_parquet(pvals)
        col = s.columns[0] if "p" not in s.columns else "p"
        vals = s[s.columns[0]].to_numpy()
        # find the p-value column robustly
        for c in s.columns:
            if "p" in c.lower():
                vals = s[c].to_numpy(); col = c; break
        hist, edges = np.histogram(vals[~np.isnan(vals)], bins=50, range=(0, 1))
        hdf = pl.DataFrame({
            "bin_lo": edges[:-1], "bin_hi": edges[1:], "count": hist,
        })
        hdf.write_csv(OUT / "null_calibration" / "gse263850_lr_pvalue_histogram_50bins.tsv", separator="\t")
        print(f"  null_calibration/gse263850_lr_pvalue_histogram_50bins.tsv "
              f"(histogram of {len(vals):,} per-CpG null p-values, column '{col}')")

    # Copy the per-figure source TSVs for the Study 3 three-way figures.
    if THREE_WAY_TSV.exists():
        dst = OUT / "study3_real" / "figure_source_tsv"
        shutil.copytree(THREE_WAY_TSV, dst, dirs_exist_ok=True)
        n = len(list(dst.glob("*.tsv")))
        print(f"  study3_real/figure_source_tsv/  ({n} per-figure TSVs copied)")

    readme = """# all_data_tsv -- every benchmark data table, as TSV

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

## A note on `study3_real/_other_views/`
Two files there answer a DIFFERENT question than the headline paper and
are kept separate so they don't confuse:
- `dmr_iou_ALLSIG_chainmerge36811.tsv` -- three-way IoU using epykit's
  **all-significant** (q<0.05) DMR set (**36,811** regions), NOT the 852
  paper-faithful chain_merge set the paper reports on.
- `per_dmr_concordance_phase3_allpairs.tsv` -- an earlier Phase-4 all-tool
  concordance computation (2,149 rows). The headline F9 concordance is
  `study3_real/per_dmr_concordance.tsv` (1,352 ek-vs-DSS matched pairs).
Everything else in `study3_real/` is the headline 852/922/1,139 data.

## External inputs (NOT in this folder -- they live outside the repo)
- Paper Supp Table 5: `epykit2/GSE263850_RAW/Paper resources/DMR_total_list.xlsx`
- methylKit-tile real-data DMRs + step_benchmarks: `methylkit_realResults/...`

## See also
- `benchmark/paper_data/` -- the curated, paper-section-organised TSV
  mirror the manuscript cites (01_headline_piao ... 06_methodology).
- This folder is the *complete* dump; paper_data is the *curated* subset.
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")
    total = len(list(OUT.rglob("*.tsv")))
    print(f"\nDone. {total} TSVs + README under {OUT}")


if __name__ == "__main__":
    main()
