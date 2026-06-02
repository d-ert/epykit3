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


def _piao_sample_path() -> Path | None:
    """Return the path to a Piao coverage=10 sample if it exists locally."""
    candidates = [
        # User-confirmed location (Phase 2 Task 1 archive):
        Path("_legacy_benchmark/deneme2/raw_sim_data/simulated_datasets/dmc_simulation/coverage/amp.coverage=10.sample1.txt"),
        Path("benchmark/data/study1/raw_sim_data/simulated_datasets/dmc_simulation/coverage/amp.coverage=10.sample1.txt"),
        Path("D:/Coding/Projeler/methyl_lib/benchmarkin_merges/epykit_vs_allPackages(simulated_approxData)/raw_sim_data/simulated_datasets/dmc_simulation/coverage/amp.coverage=10.sample1.txt"),
    ]
    for p in candidates:
        if p.exists():
            return p
    return None


def test_simulator_marginals_match_piao_within_tolerance(tmp_path):
    """The simulator's per-CpG count_M distribution at coverage=10 should
    have the same first two moments as Piao's distributed sample within
    Monte Carlo noise on ~100K CpGs. This is the 'did I rebuild the right
    simulator?' check.

    Tolerance is loose because we don't know Piao's exact baseline model:
    mean within ±10%, std within ±20%. Failure here means investigate the
    `_draw_baseline_beta` mixture parameters.
    """
    from simulate_piao import simulate_dmc

    piao = _piao_sample_path()
    if piao is None:
        pytest.skip("Piao raw data not available locally; skipping marginal match")

    # Read Piao's count_M distribution.
    piao_df = pl.read_csv(piao, separator="\t")
    piao_count_M = (piao_df["coverage"].cast(pl.Float64) * piao_df["freqC"] / 100.0).round().cast(pl.Int64)
    piao_mean = float(piao_count_M.mean())
    piao_std = float(piao_count_M.std())

    # Match Piao's CpG count (100K for DMC sim) and coverage (10).
    res = simulate_dmc(n_cpgs=len(piao_count_M), n_per_group=3, coverage=10,
                       seed=12345, out_dir=tmp_path)
    sim_amp = res["amp_files"][0]
    sim_df = pl.read_csv(sim_amp, separator="\t")
    sim_count_M = (sim_df["coverage"].cast(pl.Float64) * sim_df["freqC"] / 100.0).round().cast(pl.Int64)
    sim_mean = float(sim_count_M.mean())
    sim_std = float(sim_count_M.std())

    # Tolerances are loose: we don't claim Piao's exact baseline model.
    rel_mean_err = abs(sim_mean - piao_mean) / piao_mean
    rel_std_err = abs(sim_std - piao_std) / piao_std
    assert rel_mean_err < 0.10, (
        f"simulator count_M mean {sim_mean:.2f} vs Piao {piao_mean:.2f}: "
        f"rel error {rel_mean_err:.3f} > 10%"
    )
    assert rel_std_err < 0.20, (
        f"simulator count_M std {sim_std:.2f} vs Piao {piao_std:.2f}: "
        f"rel error {rel_std_err:.3f} > 20%"
    )


def test_simulator_truth_dmc_count_close_to_piao_design(tmp_path):
    """Piao's design has exactly 20,000 / 100,000 = 20% true DMCs.
    The simulator with default dmc_fraction=0.2 should land at 20% ± 0.5%
    on 100K CpGs (~50 std error)."""
    from simulate_piao import simulate_dmc

    res = simulate_dmc(n_cpgs=100_000, n_per_group=3, coverage=10,
                       seed=2026, out_dir=tmp_path)
    truth = pl.read_parquet(res["truth"])
    n_dmc = int(truth["is_dmc"].sum())
    assert 19_500 <= n_dmc <= 20_500, (
        f"n_dmc = {n_dmc:,}; design is 20,000 (20% of 100,000). "
        f"Outside ±0.5% tolerance suggests a bug in _assign_dmcs."
    )
