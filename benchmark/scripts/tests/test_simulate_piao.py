"""Tests for the Piao 2021 binomial simulator re-implementation."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest


def test_simulate_dmc_truth_schema_matches_dmc_truth_parquet(tmp_path):
    """The simulator's truth parquet must match the columns + dtypes of
    evaluate.py's expected `dmc_truth.parquet`."""
    from simulate_piao import simulate_dmc

    result = simulate_dmc(
        n_cpgs=2000,
        n_per_group=3,
        coverage=10,
        seed=42,
        out_dir=tmp_path,
    )
    truth = pl.read_parquet(result["truth"])
    assert set(truth.columns) == {
        "chrom", "pos", "mean_beta_treat", "mean_beta_ctrl",
        "true_meth_diff", "is_dmc", "direction", "meth_diff_bin",
    }, f"Schema mismatch: {truth.columns}"
    # Dtype check on the columns that downstream code casts.
    assert truth["chrom"].dtype == pl.Utf8
    assert truth["pos"].dtype == pl.Int64
    assert truth["true_meth_diff"].dtype == pl.Float64
    assert truth["is_dmc"].dtype == pl.Boolean


def test_simulate_dmc_is_deterministic_with_seed(tmp_path):
    """Two calls with the same seed must produce bit-identical truth + reads."""
    from simulate_piao import simulate_dmc

    a = simulate_dmc(n_cpgs=1000, n_per_group=3, coverage=10, seed=7, out_dir=tmp_path / "a")
    b = simulate_dmc(n_cpgs=1000, n_per_group=3, coverage=10, seed=7, out_dir=tmp_path / "b")
    truth_a = pl.read_parquet(a["truth"])
    truth_b = pl.read_parquet(b["truth"])
    # Same is_dmc vector and same true_meth_diff vector.
    np.testing.assert_array_equal(truth_a["is_dmc"].to_numpy(), truth_b["is_dmc"].to_numpy())
    np.testing.assert_array_equal(truth_a["true_meth_diff"].to_numpy(), truth_b["true_meth_diff"].to_numpy())
    # Same per-sample AMP files.
    for i in range(1, 7):
        f_a = tmp_path / "a" / f"amp.coverage=10.sample{i}.txt"
        f_b = tmp_path / "b" / f"amp.coverage=10.sample{i}.txt"
        assert f_a.read_bytes() == f_b.read_bytes(), f"sample{i} differs"


def test_simulate_dmc_marginals_match_design(tmp_path):
    """~20% true DMCs with |meth_diff| in [0.2, 1.0] and 50/50 direction split."""
    from simulate_piao import simulate_dmc

    result = simulate_dmc(
        n_cpgs=10000,
        n_per_group=3,
        coverage=10,
        seed=1,
        out_dir=tmp_path,
    )
    truth = pl.read_parquet(result["truth"])

    n_dmc = int(truth["is_dmc"].sum())
    frac = n_dmc / len(truth)
    assert 0.18 <= frac <= 0.22, f"DMC fraction {frac:.3f} outside design 0.20 ± 0.02"

    # Among true DMCs, |true_meth_diff| ~ U(0.2, 1.0): mean ≈ 0.6, min >= 0.2, max <= 1.0.
    dmc_only = truth.filter(pl.col("is_dmc"))
    abs_diff = dmc_only["true_meth_diff"].abs().to_numpy()
    assert abs_diff.min() >= 0.20, f"min |meth_diff| {abs_diff.min():.3f} below 0.2"
    assert abs_diff.max() <= 1.00, f"max |meth_diff| {abs_diff.max():.3f} above 1.0"
    assert 0.55 <= abs_diff.mean() <= 0.65, f"mean |meth_diff| {abs_diff.mean():.3f} outside U(0.2,1.0) expectation 0.6"

    # Direction split is ~50/50 among true DMCs.
    n_hyper = int((dmc_only["direction"] == "hyper").sum())
    n_hypo = int((dmc_only["direction"] == "hypo").sum())
    assert abs(n_hyper - n_hypo) / n_dmc < 0.05, (
        f"direction split unbalanced: {n_hyper} hyper vs {n_hypo} hypo (ratio "
        f"{n_hyper / n_dmc:.3f})"
    )
