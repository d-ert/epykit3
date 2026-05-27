"""Tests for run_null_calibration.py.

The real null-calibration runner dispatches to epykit engines via
`ep.tl.dmc`, which is heavy and integration-test-only. These tests use
a mock engine that returns prebuilt q-value arrays so we can verify the
shuffle loop, the FDR computation, and the Wilson CI integration without
running real DMC calls.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest


def test_compute_observed_fdr_matches_definition():
    """Observed FDR at nominal q<0.05 = (# q<0.05) / total. On a fixed
    array of qvalues, the function must compute this exactly."""
    from run_null_calibration import compute_observed_fdr

    qvals = np.array([0.001, 0.01, 0.04, 0.05, 0.1, 0.5, 0.9])
    # 3 sites with q<0.05 out of 7.
    out = compute_observed_fdr(qvals, q_thresh=0.05)
    assert out["n_called"] == 3
    assert out["n_total"] == 7
    assert abs(out["observed_fdr"] - 3 / 7) < 1e-12


def test_compute_observed_fdr_ignores_nan():
    """NaN q-values are excluded from the denominator."""
    from run_null_calibration import compute_observed_fdr

    qvals = np.array([0.001, 0.01, np.nan, np.nan, 0.5])
    out = compute_observed_fdr(qvals, q_thresh=0.05)
    assert out["n_called"] == 2
    assert out["n_total"] == 3
    assert abs(out["observed_fdr"] - 2 / 3) < 1e-12


def test_run_shuffles_returns_one_row_per_shuffle_with_ci(tmp_path):
    """The runner returns a frame with k rows per (engine, scenario),
    each carrying observed_fdr and tpr_ci_lo/tpr_ci_hi (here repurposed
    as Wilson CI bounds on the observed FDR proportion)."""
    from run_null_calibration import run_null_calibration

    # Mock engine: returns deterministic q-values per shuffle seed.
    def mock_engine(samples_treatment, samples_control, seed=0, n_sites=200):
        rng = np.random.default_rng(seed)
        # Under-null engine: q-values uniform-ish on [0, 1].
        return rng.uniform(0, 1, size=n_sites)

    samples = [f"s{i}" for i in range(1, 7)]
    out = run_null_calibration(
        engine_fn=mock_engine,
        engine_name="mock_lr",
        scenario_name="cov10_3v3",
        samples=samples,
        n_per_group=3,
        k_shuffles=10,
        q_thresh=0.05,
        seed=42,
    )
    assert isinstance(out, pl.DataFrame)
    assert set(out.columns) >= {
        "engine", "scenario", "k_shuffle", "observed_fdr",
        "n_called", "n_total",
        "observed_fdr_ci_lo", "observed_fdr_ci_hi",
    }
    assert len(out) == 10
    # All entries reference the same engine + scenario.
    assert out["engine"].unique().to_list() == ["mock_lr"]
    assert out["scenario"].unique().to_list() == ["cov10_3v3"]
    # Wilson CI bounds bracket observed_fdr.
    assert (out["observed_fdr_ci_lo"] <= out["observed_fdr"]).all()
    assert (out["observed_fdr_ci_hi"] >= out["observed_fdr"]).all()


def test_run_shuffles_is_deterministic_with_seed(tmp_path):
    """Two runs with same seed give identical observed_fdr column."""
    from run_null_calibration import run_null_calibration

    def mock_engine(samples_treatment, samples_control, seed=0, n_sites=200):
        rng = np.random.default_rng(seed)
        return rng.uniform(0, 1, size=n_sites)

    samples = [f"s{i}" for i in range(1, 7)]
    a = run_null_calibration(
        engine_fn=mock_engine, engine_name="mock", scenario_name="s1",
        samples=samples, n_per_group=3, k_shuffles=5, q_thresh=0.05, seed=7,
    )
    b = run_null_calibration(
        engine_fn=mock_engine, engine_name="mock", scenario_name="s1",
        samples=samples, n_per_group=3, k_shuffles=5, q_thresh=0.05, seed=7,
    )
    np.testing.assert_array_equal(
        a["observed_fdr"].to_numpy(), b["observed_fdr"].to_numpy(),
    )
