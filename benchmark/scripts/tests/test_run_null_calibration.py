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
    """The runner returns ``(per_shuffle_fdr, pvalue_samples)``. The
    first frame has k rows per (engine, scenario) with observed_fdr +
    Wilson CI bounds on the FDR proportion.

    Mock engines that only return q-values (the legacy single-array
    contract) are unpacked to ``(empty, qvals)`` so the per-shuffle
    FDR frame populates correctly while ``pvalue_samples`` stays empty.
    """
    from run_null_calibration import run_null_calibration

    def mock_engine(samples_treatment, samples_control, seed=0, n_sites=200):
        rng = np.random.default_rng(seed)
        return rng.uniform(0, 1, size=n_sites)

    samples = [f"s{i}" for i in range(1, 7)]
    per_shuffle, pvalue_samples = run_null_calibration(
        engine_fn=mock_engine,
        engine_name="mock_lr",
        scenario_name="cov10_3v3",
        samples=samples,
        n_per_group=3,
        k_shuffles=10,
        q_thresh=0.05,
        seed=42,
    )
    assert isinstance(per_shuffle, pl.DataFrame)
    assert isinstance(pvalue_samples, pl.DataFrame)
    assert set(per_shuffle.columns) >= {
        "engine", "scenario", "k_shuffle", "observed_fdr",
        "n_called", "n_total",
        "observed_fdr_ci_lo", "observed_fdr_ci_hi",
    }
    assert len(per_shuffle) == 10
    assert per_shuffle["engine"].unique().to_list() == ["mock_lr"]
    assert per_shuffle["scenario"].unique().to_list() == ["cov10_3v3"]
    assert (per_shuffle["observed_fdr_ci_lo"] <= per_shuffle["observed_fdr"]).all()
    assert (per_shuffle["observed_fdr_ci_hi"] >= per_shuffle["observed_fdr"]).all()
    # Legacy single-array mock engine yields no p-values for Q-Q.
    assert pvalue_samples.is_empty()


def test_run_shuffles_is_deterministic_with_seed(tmp_path):
    """Two runs with same seed give identical observed_fdr column."""
    from run_null_calibration import run_null_calibration

    def mock_engine(samples_treatment, samples_control, seed=0, n_sites=200):
        rng = np.random.default_rng(seed)
        return rng.uniform(0, 1, size=n_sites)

    samples = [f"s{i}" for i in range(1, 7)]
    a, _ = run_null_calibration(
        engine_fn=mock_engine, engine_name="mock", scenario_name="s1",
        samples=samples, n_per_group=3, k_shuffles=5, q_thresh=0.05, seed=7,
    )
    b, _ = run_null_calibration(
        engine_fn=mock_engine, engine_name="mock", scenario_name="s1",
        samples=samples, n_per_group=3, k_shuffles=5, q_thresh=0.05, seed=7,
    )
    np.testing.assert_array_equal(
        a["observed_fdr"].to_numpy(), b["observed_fdr"].to_numpy(),
    )


def test_run_shuffles_collects_pvalues_when_engine_returns_tuple():
    """When the engine returns ``(pvalues, qvalues)``, the runner
    populates the pvalue_samples frame for the first ``qq_shuffles``
    iterations."""
    from run_null_calibration import run_null_calibration

    def mock_engine_with_pvals(
        samples_treatment, samples_control, seed=0, n_sites=200,
    ):
        rng = np.random.default_rng(seed)
        pvals = rng.uniform(0, 1, size=n_sites)
        qvals = np.clip(pvals * 2.0, 0.0, 1.0)
        return pvals, qvals

    samples = [f"s{i}" for i in range(1, 7)]
    per_shuffle, pvalue_samples = run_null_calibration(
        engine_fn=mock_engine_with_pvals,
        engine_name="mock_lr",
        scenario_name="cov10_3v3",
        samples=samples,
        n_per_group=3,
        k_shuffles=15,
        q_thresh=0.05,
        seed=42,
        qq_shuffles=3,
        qq_samples_per_shuffle=100,
    )
    assert per_shuffle.height == 15
    assert not pvalue_samples.is_empty()
    # Only the first 3 shuffles contribute p-values.
    assert pvalue_samples["k_shuffle"].n_unique() == 3
    assert pvalue_samples["k_shuffle"].max() == 2
    assert pvalue_samples.height <= 3 * 100


def test_summarize_observed_fdr_emits_median_iqr():
    """``summarize_observed_fdr`` reduces per-shuffle rows to a single
    median + IQR row."""
    from run_null_calibration import (
        run_null_calibration,
        summarize_observed_fdr,
    )

    def mock_engine(samples_treatment, samples_control, seed=0, n_sites=200):
        rng = np.random.default_rng(seed)
        return rng.uniform(0, 1, size=n_sites)

    samples = [f"s{i}" for i in range(1, 7)]
    per_shuffle, _ = run_null_calibration(
        engine_fn=mock_engine, engine_name="mock",
        scenario_name="cov10_3v3", samples=samples,
        n_per_group=3, k_shuffles=20, q_thresh=0.05, seed=0,
    )
    summary = summarize_observed_fdr(per_shuffle)
    assert summary.height == 1
    row = summary.row(0, named=True)
    assert row["k_shuffles"] == 20
    assert row["q25_observed_fdr"] <= row["median_observed_fdr"]
    assert row["median_observed_fdr"] <= row["q75_observed_fdr"]


def test_ks_against_uniform_under_truly_uniform():
    """Under a true Uniform(0, 1) draw with large n the K-S test should
    not reject."""
    from run_null_calibration import ks_against_uniform

    rng = np.random.default_rng(42)
    p = rng.uniform(0, 1, size=20_000)
    res = ks_against_uniform(p)
    assert res["n"] == 20_000
    # K-S should NOT reject under the null at alpha = 0.01.
    assert res["pvalue"] > 0.01
    # D-statistic is small for n=20k uniform draws.
    assert res["statistic"] < 0.02


def test_ks_against_uniform_under_conservative_test():
    """A test that produces p-values stochastically larger than uniform
    (e.g., all p-values shifted toward 1) should be rejected by K-S."""
    from run_null_calibration import ks_against_uniform

    rng = np.random.default_rng(42)
    # Shifted Beta(2, 1): CDF F(x) = x^2, so stochastically larger than U(0,1).
    p = rng.beta(2.0, 1.0, size=20_000)
    res = ks_against_uniform(p)
    assert res["n"] == 20_000
    assert res["pvalue"] < 1e-6
