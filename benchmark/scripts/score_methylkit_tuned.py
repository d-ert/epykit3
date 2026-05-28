"""score_methylkit_tuned.py -- Phase 4 Task 4 Step 2.

Score the methylKit Stouffer-tuned per-scenario DMC TSVs against the
Study 1 ground truth and APPEND the resulting rows (under tool name
``methylkit_tuned``) to ``benchmark/data/study1/eval_summary_post_phase3.parquet``.

Step 1 (commit f2cb369) produced the 10 tuned TSVs under
``benchmark/data/study2/methylkit_tuned/`` by applying
``methylkit_stouffer_combine.R`` to the existing per-scenario methylKit
DMC TSVs. This step closes the loop by re-scoring those tuned calls so
reviewers can see the tuning delta alongside the un-tuned ``methylkit``
rows (which are *transcribed published baselines* from Piao et al. 2021;
new rows are recomputed from the actual per-CpG calls and therefore
also carry tp/fp/tn/fn counts so Wilson CIs are populated).

Schema parity:
    Existing ``methylkit`` rows (n=30) are 10 cells x 3 meth_diff_bins
    ({0.2-0.4, 0.4-0.6, 0.6-0.8}) at threshold_kind='qvalue',
    threshold=0.05, test=null, with only ``tpr``/``fpr`` populated and
    everything else null. The ``methylkit_tuned`` rows we emit use the
    SAME (scenario, parameter, parameter_value, meth_diff_bin,
    threshold_kind, threshold, test) layout so they pair 1:1 with the
    untuned rows -- but we additionally fill the integer counts and
    derived per-bin metrics (precision/f1) because we have the per-CpG
    data. AUROC is left null at the per-bin level (AUROC is defined over
    the full ROC, not a meth_diff stratum -- matches Task 2's per-bin
    pattern for epykit rows).

Scoring rule:
    The TSVs carry both raw (``pvalue``/``qvalue``) and combined
    (``pvalue_combined``/``qvalue_combined``) columns. We score using
    the combined values (this is the entire point of tuning), falling
    back to raw per-CpG when combined is null (which can happen for
    isolated CpGs with no neighbours within --max-gap-bp). This
    mirrors epykit's lr+ scoring convention in ``_join_with_truth``.

Usage:
    uv run python benchmark/scripts/score_methylkit_tuned.py
    uv run python benchmark/scripts/score_methylkit_tuned.py --skip-ci
    uv run python benchmark/scripts/score_methylkit_tuned.py \
        --eval-summary <alt> --tuned-dir <alt> --truth <alt>

Exit non-zero on any per-cell failure.
"""

from __future__ import annotations

import argparse
import logging
import re
import subprocess
import sys
from pathlib import Path

import polars as pl


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent  # benchmark/
EVAL_SUMMARY = ROOT / "data" / "study1" / "eval_summary_post_phase3.parquet"
TUNED_DIR = ROOT / "data" / "study2" / "methylkit_tuned"
TRUTH_PATH = ROOT / "data" / "study1" / "ground_truth" / "dmc_truth.parquet"
MANIFEST_PATH = ROOT / "data" / "study1" / "epykit_post_phase3" / "MANIFEST.txt"

# Shared scoring contract.
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from _epykit_scoring import Q_THRESHOLD, split_for_ci  # noqa: E402


# methylkit-row schema constants (mirrors existing untuned rows).
METHYLKIT_BINS = ("0.2-0.4", "0.4-0.6", "0.6-0.8")
TUNED_TOOL = "methylkit_tuned"

# Filename -> (scenario, parameter, parameter_value).
FNAME_RE = re.compile(r"^dmc_(cov|rep)(\d+)_tuned\.tsv$")


logger = logging.getLogger("score_methylkit_tuned")


def parse_cell(fname: str) -> tuple[str, str, int]:
    """Map ``dmc_cov10_tuned.tsv`` -> (scenario, parameter, parameter_value).

    Mirrors the convention used for the existing methylkit rows in
    ``eval_summary_post_phase3.parquet``:

        cov<N>  -> scenario='dmc_coverage',  parameter='coverage', value=N
        rep<N>  -> scenario='dmc_replicate', parameter='n_total',  value=N
    """
    m = FNAME_RE.match(fname)
    if not m:
        raise ValueError(f"unrecognised tuned-TSV filename: {fname}")
    kind, num = m.group(1), int(m.group(2))
    if kind == "cov":
        return ("dmc_coverage", "coverage", num)
    return ("dmc_replicate", "n_total", num)


def load_tuned_calls(tsv: Path) -> pl.DataFrame:
    """Read a tuned methylKit TSV and project to (chrom, pos, pvalue, qvalue).

    Uses ``pvalue_combined``/``qvalue_combined`` when present and non-null,
    falling back to the raw ``pvalue``/``qvalue`` per row (isolated CpGs
    with no neighbour Stouffer combine retain null combined values, per
    the R script's docstring).
    """
    df = pl.read_csv(tsv, separator="\t")
    required = {"chr", "start", "pvalue", "qvalue",
                "pvalue_combined", "qvalue_combined"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{tsv.name}: missing columns {sorted(missing)}")
    # Cell-by-cell fallback: combined where non-null, raw otherwise.
    return df.select(
        pl.col("chr").alias("chrom"),
        pl.col("start").alias("pos"),
        pl.coalesce(
            pl.col("pvalue_combined").cast(pl.Float64),
            pl.col("pvalue").cast(pl.Float64),
        ).alias("pvalue"),
        pl.coalesce(
            pl.col("qvalue_combined").cast(pl.Float64),
            pl.col("qvalue").cast(pl.Float64),
        ).alias("qvalue"),
    )


def score_one_cell(
    tsv: Path, truth: pl.DataFrame,
) -> list[dict]:
    """Score one tuned TSV at q < ``Q_THRESHOLD`` per meth_diff bin.

    Emits one row per (cell, bin) for bins in ``METHYLKIT_BINS`` --
    matches the existing methylkit-row layout (no 'all' row, no
    pvalue-threshold sweep, no 0.8-1.0 bin). Rows carry integer counts
    + per-bin precision/f1; AUROC is null at the per-bin level.
    """
    scenario, parameter, value = parse_cell(tsv.name)
    calls = load_tuned_calls(tsv)
    joined = truth.join(calls, on=["chrom", "pos"], how="left")

    rows: list[dict] = []
    sig_expr = pl.col("qvalue").fill_null(1.0) < Q_THRESHOLD
    sub = joined.with_columns(sig=sig_expr)
    fp_total = sub.filter(pl.col("sig") & ~pl.col("is_dmc")).height
    tn_total = sub.filter(~pl.col("sig") & ~pl.col("is_dmc")).height

    for bin_label in METHYLKIT_BINS:
        in_bin = pl.col("is_dmc") & (pl.col("meth_diff_bin") == bin_label)
        tp = sub.filter(pl.col("sig") & in_bin).height
        fn = sub.filter(~pl.col("sig") & in_bin).height
        n_pos = tp + fn
        # FPR uses the genome-wide negative pool (same convention as the
        # legacy per-bin TPR rows in _epykit_scoring.score_dmc_parquet
        # and the published methylkit baselines).
        n_neg = fp_total + tn_total
        denom_f1 = 2 * tp + fp_total + fn
        rows.append({
            "tool": TUNED_TOOL,
            "scenario": scenario,
            "parameter": parameter,
            "parameter_value": value,
            "test": None,
            "meth_diff_bin": bin_label,
            "threshold_kind": "qvalue",
            "threshold": Q_THRESHOLD,
            "tp": tp, "fp": fp_total, "tn": tn_total, "fn": fn,
            "tpr": tp / n_pos if n_pos else 0.0,
            "fpr": fp_total / n_neg if n_neg else 0.0,
            "precision": tp / (tp + fp_total) if (tp + fp_total) else 0.0,
            "f1": (2 * tp / denom_f1) if denom_f1 else 0.0,
            "auroc": None,
        })
    return rows


def append_ci_for_new_rows(new_df: pl.DataFrame) -> pl.DataFrame:
    """Run ``evaluate.py --ci-only`` on the new rows in isolation.

    Same pattern as ``run_epykit_simulator.add_ci_via_evaluate``: split
    rows with counts (all of ours, in this case) from rows without,
    write a temp parquet, invoke the sealed evaluate.py, and read back.
    """
    has_counts, no_counts = split_for_ci(new_df)
    if has_counts.height == 0:
        return new_df

    tmp = (TUNED_DIR / "_ci_tmp_methylkit_tuned.parquet")
    has_counts.write_parquet(tmp)
    proc = subprocess.run(
        [sys.executable, str(HERE / "evaluate.py"),
         "--ci-only", "--eval-summary", str(tmp)],
        check=False, capture_output=True, text=True,
    )
    if proc.returncode != 0:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(
            f"evaluate.py --ci-only failed on methylkit_tuned rows:\n"
            f"{proc.stderr}"
        )
    with_ci = pl.read_parquet(tmp)
    tmp.unlink(missing_ok=True)

    if no_counts.height:
        for c in ("tpr_ci_lo", "tpr_ci_hi", "fpr_ci_lo", "fpr_ci_hi",
                  "auroc_ci_lo", "auroc_ci_hi", "f1_ci_lo", "f1_ci_hi"):
            if c not in no_counts.columns:
                no_counts = no_counts.with_columns(
                    pl.lit(float("nan")).alias(c).cast(pl.Float64)
                )
        no_counts = no_counts.select(with_ci.columns)
        return pl.concat([with_ci, no_counts], how="vertical_relaxed")
    return with_ci


def append_to_eval_summary(
    new_rows: list[dict],
    eval_summary_path: Path,
    *,
    skip_ci: bool = False,
) -> pl.DataFrame:
    """Append new methylkit_tuned rows to the existing eval summary.

    * Reads the existing parquet (must already contain Phase-3 rows).
    * Refuses to clobber: errors out if any ``methylkit_tuned`` row is
      already present (use ``--force`` to re-run; we keep the rerun-
      safety explicit because the eval summary is downstream-consumed).
    * Aligns schemas via vertical_relaxed concat.
    """
    existing = pl.read_parquet(eval_summary_path)
    if (existing["tool"] == TUNED_TOOL).any():
        raise RuntimeError(
            f"{eval_summary_path.name} already contains "
            f"{TUNED_TOOL} rows; refusing to clobber."
        )

    new_df = pl.DataFrame(new_rows)
    if not skip_ci:
        new_df = append_ci_for_new_rows(new_df)

    # Schema alignment: union columns, fill missing with null on both
    # sides, vertical_relaxed for dtype promotion.
    all_cols = sorted(set(existing.columns) | set(new_df.columns))
    for c in all_cols:
        if c not in existing.columns:
            existing = existing.with_columns(pl.lit(None).alias(c))
        if c not in new_df.columns:
            new_df = new_df.with_columns(pl.lit(None).alias(c))
    # Preserve the existing parquet's column order so a column-order
    # diff doesn't show up downstream.
    out_cols = existing.columns
    new_df = new_df.select(out_cols)
    combined = pl.concat(
        [existing, new_df], how="vertical_relaxed",
    )
    return combined


def update_manifest(n_added: int, manifest_path: Path) -> None:
    """Append a 'methylkit_tuned' section to the post-Phase-3 MANIFEST."""
    section = [
        "",
        "Methylkit-tuned rows added -- Phase 4 Task 4 Step 2",
        "----------------------------------------------------",
        "Tool name in eval_summary_post_phase3.parquet: methylkit_tuned",
        f"Rows added: {n_added}  (10 cells x 3 meth_diff_bins)",
        "Source TSVs: benchmark/data/study2/methylkit_tuned/dmc_*_tuned.tsv",
        f"Threshold: qvalue < {Q_THRESHOLD}, bins {list(METHYLKIT_BINS)}",
        "Scoring column: qvalue_combined (coalesce -> qvalue per row)",
        "Truth: benchmark/data/study1/ground_truth/dmc_truth.parquet",
        "Script: benchmark/scripts/score_methylkit_tuned.py",
        "Pair: row pairs 1:1 with the existing 'methylkit' rows (same",
        "  scenario/parameter/parameter_value/meth_diff_bin keys), with",
        "  integer tp/fp/tn/fn populated so Wilson CIs are real (the",
        "  untuned methylkit rows carry transcribed Piao-paper TPR/FPR",
        "  with null counts, hence NaN CIs).",
        "",
    ]
    text = manifest_path.read_text(encoding="utf-8") + "\n".join(section)
    manifest_path.write_text(text, encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eval-summary", type=Path, default=EVAL_SUMMARY)
    parser.add_argument("--tuned-dir", type=Path, default=TUNED_DIR)
    parser.add_argument("--truth", type=Path, default=TRUTH_PATH)
    parser.add_argument("--manifest", type=Path, default=MANIFEST_PATH)
    parser.add_argument(
        "--skip-ci", action="store_true",
        help="Don't run evaluate.py --ci-only (CI cols will be null).",
    )
    parser.add_argument(
        "--skip-manifest", action="store_true",
        help="Don't update MANIFEST.txt.",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )

    truth = pl.read_parquet(args.truth)

    tsvs = sorted(args.tuned_dir.glob("dmc_*_tuned.tsv"))
    if not tsvs:
        logger.error("no tuned TSVs found in %s", args.tuned_dir)
        return 2
    logger.info("scoring %d tuned TSVs", len(tsvs))

    all_rows: list[dict] = []
    for tsv in tsvs:
        rows = score_one_cell(tsv, truth)
        logger.info("[%s] +%d rows", tsv.name, len(rows))
        all_rows.extend(rows)

    combined = append_to_eval_summary(
        all_rows, args.eval_summary, skip_ci=args.skip_ci,
    )
    combined.write_parquet(args.eval_summary)
    logger.info(
        "wrote %s (%d rows; +%d methylkit_tuned)",
        args.eval_summary.name, combined.height, len(all_rows),
    )

    if not args.skip_manifest:
        update_manifest(len(all_rows), args.manifest)
        logger.info("appended methylkit_tuned section to %s",
                    args.manifest.name)

    return 0


if __name__ == "__main__":
    sys.exit(main())
