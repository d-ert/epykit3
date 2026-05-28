"""Tests for the pure scoring helpers in ``_epykit_scoring.py``.

The module is pure (no I/O paths, no global state), so these tests
construct small in-memory frames and assert exact numeric outputs. The
goal is to lock in the Phase 4 contract -- the lr+ pvalue_combined
selection rule and the AUROC all-equal-scores guard -- so future
refactors cannot silently change either.
"""

from __future__ import annotations

import polars as pl
import pytest


# conftest.py wires benchmark/scripts/ onto sys.path
from _epykit_scoring import (
    P_THRESHOLDS,
    Q_THRESHOLD,
    STALE_EPYKIT_TOOLS,
    _auroc,
    _confusion,
    _dmc_kwargs,
    _join_with_truth,
    score_dmc_parquet,
    split_for_ci,
)


# ---------------------------------------------------------------------------
# _dmc_kwargs: lr+ contract single source of truth
# ---------------------------------------------------------------------------


def test_dmc_kwargs_lrplus_returns_lr_backend():
    backend, kw = _dmc_kwargs("lr+", allow_n1=False)
    assert backend == "lr"
    # The power-stack four: each is README/CLAUDE.md-validated.
    assert kw["test"] == "lr"
    assert kw["allow_n1"] is False
    assert kw["neighbour_combine"] is True
    assert kw["neighbour_bp"] == 500
    assert kw["sep_fallback"] is True
    assert kw["sep_threshold"] == 0.9
    assert kw["fdr_method"] == "fdr_tsbh"
    assert kw["dispersion"] == "eb"


def test_dmc_kwargs_plain_lr_no_power_stack():
    backend, kw = _dmc_kwargs("lr", allow_n1=False)
    assert backend == "lr"
    assert kw == {"test": "lr", "allow_n1": False}
    # No power-stack opts on plain lr
    assert "neighbour_combine" not in kw
    assert "fdr_method" not in kw


def test_dmc_kwargs_propagates_allow_n1():
    _, kw_n1 = _dmc_kwargs("fisher", allow_n1=True)
    assert kw_n1["allow_n1"] is True


# ---------------------------------------------------------------------------
# _auroc edge cases (Phase 4 Task 3 hardening)
# ---------------------------------------------------------------------------


def test_auroc_all_equal_scores_returns_nan():
    """When an engine emits a constant p-value the ranking is degenerate.

    Reporting 0.5 from the tied-rank formula would silently misreport;
    NaN forces the consumer to acknowledge the degeneracy.
    """
    joined = pl.DataFrame({
        "is_dmc": [True, False, True, False, True, False],
        "pvalue": [0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
    })
    assert _auroc(joined) != _auroc(joined) or (
        # NaN != NaN; the canonical NaN check is `x != x` or math.isnan
        str(_auroc(joined)) == "nan"
    )


def test_auroc_no_positives_returns_nan():
    joined = pl.DataFrame({
        "is_dmc": [False, False, False],
        "pvalue": [0.1, 0.5, 0.9],
    })
    out = _auroc(joined)
    assert str(out) == "nan"


def test_auroc_no_negatives_returns_nan():
    joined = pl.DataFrame({
        "is_dmc": [True, True, True],
        "pvalue": [0.1, 0.5, 0.9],
    })
    out = _auroc(joined)
    assert str(out) == "nan"


def test_auroc_perfect_separation_is_one():
    """All positives have lower p than all negatives -> AUROC = 1.0."""
    joined = pl.DataFrame({
        "is_dmc": [True, True, True, False, False, False],
        "pvalue": [0.01, 0.02, 0.03, 0.7, 0.8, 0.9],
    })
    assert _auroc(joined) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# _join_with_truth: lr+ pvalue_combined selection rule
# ---------------------------------------------------------------------------


def test_join_with_truth_prefers_combined_columns_when_present():
    """When pvalue_combined is present (lr+), it MUST become the scoring
    pvalue. Otherwise lr+ would silently collapse to plain lr.
    """
    truth = pl.DataFrame({
        "chrom": ["chr1", "chr1"], "pos": [100, 200],
        "is_dmc": [True, False],
    })
    epy = pl.DataFrame({
        "chrom": ["chr1", "chr1"], "pos": [100, 200],
        "pvalue": [0.5, 0.5],           # raw -- uninformative
        "qvalue": [0.6, 0.6],
        "pvalue_combined": [0.001, 0.9],  # combined -- informative
        "qvalue_combined": [0.002, 0.95],
    })
    joined = _join_with_truth(epy, truth)
    # The 'pvalue' column in the joined frame is the COMBINED value
    assert joined.filter(pl.col("pos") == 100)["pvalue"][0] == pytest.approx(0.001)
    assert joined.filter(pl.col("pos") == 200)["pvalue"][0] == pytest.approx(0.9)


def test_join_with_truth_falls_back_to_raw_when_no_combined():
    truth = pl.DataFrame({
        "chrom": ["chr1"], "pos": [100], "is_dmc": [True],
    })
    epy = pl.DataFrame({
        "chrom": ["chr1"], "pos": [100],
        "pvalue": [0.01], "qvalue": [0.02],
    })
    joined = _join_with_truth(epy, truth)
    assert joined["pvalue"][0] == pytest.approx(0.01)


# ---------------------------------------------------------------------------
# split_for_ci
# ---------------------------------------------------------------------------


def test_split_for_ci_separates_null_counts():
    df = pl.DataFrame({
        "tool": ["a", "b", "c"],
        "tp": [10, None, 5],
        "fp": [1, 2, None],
        "tn": [100, 200, 300],
        "fn": [1, 2, 3],
    })
    has, no = split_for_ci(df)
    assert has.height == 1
    assert no.height == 2
    assert has["tool"][0] == "a"
    assert set(no["tool"].to_list()) == {"b", "c"}


# ---------------------------------------------------------------------------
# score_dmc_parquet end-to-end smoke
# ---------------------------------------------------------------------------


def test_score_dmc_parquet_full_threshold_grid(tmp_path):
    """One synthetic cell scored end-to-end -- assert row count + schema."""
    epy = pl.DataFrame({
        "chrom": ["chr1"] * 10,
        "pos": list(range(100, 1100, 100)),
        "pvalue": [0.001] * 5 + [0.5] * 5,
        "qvalue": [0.01] * 5 + [0.6] * 5,
    })
    truth = pl.DataFrame({
        "chrom": ["chr1"] * 10,
        "pos": list(range(100, 1100, 100)),
        "is_dmc": [True] * 5 + [False] * 5,
        "meth_diff_bin": ["0.4-0.6"] * 10,
    })
    pq = tmp_path / "cell.parquet"
    epy.write_parquet(pq)

    rows = score_dmc_parquet(
        pq, truth, tool="epykit_lr", scenario="test",
        parameter="coverage", parameter_value=10, test="lr",
    )

    # 4 p-thresholds + 1 q-threshold + 4 meth_diff_bins = 9 rows
    assert len(rows) == len(P_THRESHOLDS) + 1 + 4

    # All threshold rows have valid tp/fp/tn/fn
    for r in rows:
        assert r["tp"] is not None
        assert r["fp"] is not None
        assert r["tn"] is not None
        assert r["fn"] is not None

    # The q@0.05 row: 5 TP, 0 FP, 5 TN, 0 FN -> perfect
    q_row = [
        r for r in rows
        if r["threshold_kind"] == "qvalue" and r["meth_diff_bin"] == "all"
    ][0]
    assert q_row["tp"] == 5
    assert q_row["fp"] == 0
    assert q_row["tn"] == 5
    assert q_row["fn"] == 0
    assert q_row["tpr"] == 1.0
    assert q_row["fpr"] == 0.0


# ---------------------------------------------------------------------------
# Stale-engine filter contract
# ---------------------------------------------------------------------------


def test_stale_epykit_tools_includes_bb_lr():
    """bb_lr was removed in Phase 3; rows MUST be filtered on reassembly."""
    assert "epykit_bb_lr" in STALE_EPYKIT_TOOLS
