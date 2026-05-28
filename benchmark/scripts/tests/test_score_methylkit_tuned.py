"""Tests for score_methylkit_tuned.py (Phase 4 Task 4 Step 2).

Constructs a tiny synthetic tuned methylKit TSV + truth parquet,
exercises the scorer end-to-end, and verifies the row schema, the
``methylkit_tuned`` tool tag, and the ``pvalue_combined`` selection
rule (i.e. tuning actually changes scoring vs. raw pvalues).
"""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl
import pytest

SCRIPTS = Path(__file__).resolve().parents[1]
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import score_methylkit_tuned as smkt  # noqa: E402


# ---------------------------------------------------------------------------
# Synthetic fixtures
# ---------------------------------------------------------------------------


def _make_truth(tmp_path: Path) -> pl.DataFrame:
    """Tiny truth: 6 CpGs on chr1, 3 DMC (one per bin), 3 not-DMC."""
    return pl.DataFrame({
        "chrom": ["chr1"] * 6,
        "pos":   [100, 200, 300, 400, 500, 600],
        "is_dmc": [True, True, True, False, False, False],
        "meth_diff_bin": [
            "0.2-0.4", "0.4-0.6", "0.6-0.8",
            "none", "none", "none",
        ],
    })


def _make_tuned_tsv(
    tmp_path: Path,
    *,
    raw_pvals: list[float],
    combined_pvals: list[float] | None,
    name: str = "dmc_cov10_tuned.tsv",
) -> Path:
    """Write a 6-row TSV mirroring the methylkit_stouffer_combine.R schema."""
    positions = [100, 200, 300, 400, 500, 600]
    raw_q = raw_pvals  # ok for tiny test -- not actually BH-adjusted
    comb_q = combined_pvals
    rows: dict[str, list] = {
        "chr": ["chr1"] * 6,
        "start": positions,
        "end": positions,
        "strand": ["+"] * 6,
        "pvalue": raw_pvals,
        "qvalue": raw_q,
        "meth.diff": [30.0, 40.0, 60.0, 5.0, 5.0, 5.0],
        "pvalue_combined": (
            combined_pvals if combined_pvals is not None else [None] * 6
        ),
        "qvalue_combined": (
            comb_q if comb_q is not None else [None] * 6
        ),
    }
    out = tmp_path / name
    pl.DataFrame(rows).write_csv(out, separator="\t")
    return out


def _make_existing_eval(tmp_path: Path) -> Path:
    """Write a stub eval_summary parquet carrying the 25-col schema."""
    schema = {
        "auroc": pl.Float64, "f1": pl.Float64,
        "fn": pl.Int64, "fp": pl.Int64, "fpr": pl.Float64,
        "meth_diff_bin": pl.Utf8, "parameter": pl.Utf8,
        "parameter_value": pl.Int64, "precision": pl.Float64,
        "scenario": pl.Utf8, "test": pl.Utf8,
        "threshold": pl.Float64, "threshold_kind": pl.Utf8,
        "tn": pl.Int64, "tool": pl.Utf8, "tp": pl.Int64, "tpr": pl.Float64,
        "tpr_ci_lo": pl.Float64, "tpr_ci_hi": pl.Float64,
        "fpr_ci_lo": pl.Float64, "fpr_ci_hi": pl.Float64,
        "auroc_ci_lo": pl.Float64, "auroc_ci_hi": pl.Float64,
        "f1_ci_lo": pl.Float64, "f1_ci_hi": pl.Float64,
    }
    df = pl.DataFrame(
        {
            "auroc": [None], "f1": [None], "fn": [None], "fp": [None],
            "fpr": [0.023], "meth_diff_bin": ["0.2-0.4"],
            "parameter": ["coverage"], "parameter_value": [10],
            "precision": [None], "scenario": ["dmc_coverage"],
            "test": [None], "threshold": [0.05], "threshold_kind": ["qvalue"],
            "tn": [None], "tool": ["methylkit"], "tp": [None],
            "tpr": [0.963],
            "tpr_ci_lo": [float("nan")], "tpr_ci_hi": [float("nan")],
            "fpr_ci_lo": [float("nan")], "fpr_ci_hi": [float("nan")],
            "auroc_ci_lo": [float("nan")], "auroc_ci_hi": [float("nan")],
            "f1_ci_lo": [float("nan")], "f1_ci_hi": [float("nan")],
        },
        schema=schema,
    )
    out = tmp_path / "eval_summary_post_phase3.parquet"
    df.write_parquet(out)
    return out


# ---------------------------------------------------------------------------
# parse_cell
# ---------------------------------------------------------------------------


def test_parse_cell_coverage_and_replicate():
    assert smkt.parse_cell("dmc_cov10_tuned.tsv") == (
        "dmc_coverage", "coverage", 10,
    )
    assert smkt.parse_cell("dmc_rep4_tuned.tsv") == (
        "dmc_replicate", "n_total", 4,
    )


def test_parse_cell_rejects_unknown():
    with pytest.raises(ValueError):
        smkt.parse_cell("not_a_tuned_file.tsv")


# ---------------------------------------------------------------------------
# score_one_cell -- schema + tool tag
# ---------------------------------------------------------------------------


def test_score_one_cell_row_schema(tmp_path):
    truth = _make_truth(tmp_path)
    # Combined p-values clearly significant; raw not so much.
    tsv = _make_tuned_tsv(
        tmp_path,
        raw_pvals=[0.04, 0.04, 0.04, 0.04, 0.04, 0.04],
        combined_pvals=[1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9],
    )
    rows = smkt.score_one_cell(tsv, truth)

    # 3 bins per cell, all tagged methylkit_tuned.
    assert len(rows) == 3
    assert {r["tool"] for r in rows} == {"methylkit_tuned"}
    assert {r["meth_diff_bin"] for r in rows} == set(smkt.METHYLKIT_BINS)
    for r in rows:
        assert r["scenario"] == "dmc_coverage"
        assert r["parameter"] == "coverage"
        assert r["parameter_value"] == 10
        assert r["threshold_kind"] == "qvalue"
        assert r["threshold"] == 0.05
        assert r["test"] is None
        assert r["auroc"] is None
        # Integer counts populated -> Wilson CIs will be real downstream.
        assert isinstance(r["tp"], int)
        assert isinstance(r["fp"], int)
        assert isinstance(r["tn"], int)
        assert isinstance(r["fn"], int)


# ---------------------------------------------------------------------------
# pvalue_combined vs raw -- the regression check
# ---------------------------------------------------------------------------


def test_combined_drives_scoring_when_present(tmp_path):
    """When combined q-values are significant and raw are not, all 3
    truth DMCs are called -- TPR = 1.0 per bin (one DMC per bin)."""
    truth = _make_truth(tmp_path)
    tsv = _make_tuned_tsv(
        tmp_path,
        raw_pvals=[0.5, 0.5, 0.5, 0.5, 0.5, 0.5],          # all NS
        combined_pvals=[1e-9, 1e-9, 1e-9, 0.5, 0.5, 0.5],  # DMCs sig
    )
    rows = smkt.score_one_cell(tsv, truth)
    tprs = {r["meth_diff_bin"]: r["tpr"] for r in rows}
    assert tprs == {"0.2-0.4": 1.0, "0.4-0.6": 1.0, "0.6-0.8": 1.0}


def test_falls_back_to_raw_when_combined_null(tmp_path):
    """When combined columns are entirely null, the scorer must fall
    back to raw pvalues (isolated-CpG case in the R script). Same calls
    as raw -> identical TPR to a 'raw scored' run."""
    truth = _make_truth(tmp_path)
    # Tighten so the BH-q is still small (we use raw_p == q in fixture).
    tsv = _make_tuned_tsv(
        tmp_path,
        raw_pvals=[1e-6, 1e-6, 1e-6, 0.5, 0.5, 0.5],
        combined_pvals=None,
    )
    rows = smkt.score_one_cell(tsv, truth)
    tprs = {r["meth_diff_bin"]: r["tpr"] for r in rows}
    assert tprs == {"0.2-0.4": 1.0, "0.4-0.6": 1.0, "0.6-0.8": 1.0}


def test_combined_and_raw_disagree_changes_metrics(tmp_path):
    """Regression check: if combined p-values shift the calls relative
    to raw, scoring MUST move with them. This is the entire point of
    Task 4 Step 2."""
    truth = _make_truth(tmp_path)
    # Raw: all sig -> TPR=1, but inflates FP via the 'not-DMC' CpGs.
    # Combined: only DMC positions sig -> TPR=1 with FP=0.
    tsv = _make_tuned_tsv(
        tmp_path,
        raw_pvals=[1e-6] * 6,
        combined_pvals=[1e-9, 1e-9, 1e-9, 0.5, 0.5, 0.5],
    )
    rows = smkt.score_one_cell(tsv, truth)
    for r in rows:
        assert r["fp"] == 0, f"combined should suppress non-DMC calls; got fp={r['fp']}"
        assert r["tp"] == 1


# ---------------------------------------------------------------------------
# append_to_eval_summary -- refuses to clobber + preserves existing rows
# ---------------------------------------------------------------------------


def test_append_refuses_clobber(tmp_path):
    """If methylkit_tuned rows already exist, the appender must error
    rather than silently duplicate."""
    existing = _make_existing_eval(tmp_path)
    # Hand-add a methylkit_tuned row so the gate trips.
    df = pl.read_parquet(existing)
    df = df.with_columns(
        pl.when(pl.col("tool") == "methylkit").then(
            pl.lit("methylkit_tuned")
        ).otherwise(pl.col("tool")).alias("tool")
    )
    df.write_parquet(existing)

    with pytest.raises(RuntimeError, match="already contains"):
        smkt.append_to_eval_summary(
            [{"tool": "methylkit_tuned", "scenario": "x", "parameter": "y",
              "parameter_value": 1, "test": None, "meth_diff_bin": "0.2-0.4",
              "threshold_kind": "qvalue", "threshold": 0.05,
              "tp": 1, "fp": 0, "tn": 0, "fn": 0,
              "tpr": 1.0, "fpr": 0.0, "precision": 1.0, "f1": 1.0,
              "auroc": None}],
            existing, skip_ci=True,
        )


def test_append_preserves_methylkit_rows(tmp_path):
    """Existing methylkit rows must be byte-identical after append --
    only methylkit_tuned rows are added."""
    existing_path = _make_existing_eval(tmp_path)
    before = pl.read_parquet(existing_path).filter(
        pl.col("tool") == "methylkit"
    )

    combined = smkt.append_to_eval_summary(
        [{"tool": "methylkit_tuned", "scenario": "dmc_coverage",
          "parameter": "coverage", "parameter_value": 10, "test": None,
          "meth_diff_bin": "0.2-0.4", "threshold_kind": "qvalue",
          "threshold": 0.05,
          "tp": 5, "fp": 1, "tn": 9, "fn": 1,
          "tpr": 5/6, "fpr": 1/10,
          "precision": 5/6, "f1": 0.83, "auroc": None}],
        existing_path, skip_ci=True,
    )
    after = combined.filter(pl.col("tool") == "methylkit")
    # Same shape + same values, column-by-column.
    assert before.shape == after.shape
    assert before.equals(after.select(before.columns))
    # New row landed.
    assert (combined["tool"] == "methylkit_tuned").sum() == 1
