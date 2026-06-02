"""clinical / cohort QC pack."""

from __future__ import annotations

import warnings

import numpy as np
import polars as pl
import pytest

import epykit as ep
from epykit.qc import (
    _classify_sex_from_values,
    contamination_estimate,
    power,
    sample_correlation,
    sex_check,
)


def test_sex_check_runs_with_or_without_chrx(synth_md_filtered):
    """sex_check must not crash even when chrX is absent from the fixture."""
    md = synth_md_filtered
    samples = md.obs.get_column("sample_id").to_list()
    result = sex_check(md.store, samples)
    assert isinstance(result, pl.DataFrame)
    assert {"sample_id", "mean_chrx_beta", "inferred_sex", "mismatch"} <= set(
        result.columns
    )


def test_contamination_score_in_range(synth_md_filtered):
    md = synth_md_filtered
    samples = md.obs.get_column("sample_id").to_list()
    score = contamination_estimate(md.store, samples[0])
    if score == score:  # not NaN
        assert 0.0 <= score <= 1.0


def test_sample_correlation_matrix_shape(synth_md_filtered):
    md = synth_md_filtered
    samples = md.obs.get_column("sample_id").to_list()
    corr = sample_correlation(md.store, samples, method="spearman")
    assert isinstance(corr, pl.DataFrame)
    if len(corr) > 0:
        # Diagonal == 1
        diag = corr.filter(pl.col("sample_a") == pl.col("sample_b"))
        diag_vals = diag.get_column("correlation").to_numpy()
        assert (abs(diag_vals - 1.0) < 1e-9).all()


def test_power_calc_increases_with_n():
    """Power monotone-increasing with n_per_group at fixed effect."""
    p_small = power(meth_diff=0.20, coverage=15, n_per_group=3)
    p_med   = power(meth_diff=0.20, coverage=15, n_per_group=10)
    p_big   = power(meth_diff=0.20, coverage=15, n_per_group=30)
    assert p_small < p_med < p_big
    assert 0.0 <= p_small <= 1.0
    assert 0.0 <= p_big   <= 1.0


def test_power_solves_for_n():
    """When `power=` is passed, returns the smallest n hitting the target."""
    n_needed = power(meth_diff=0.10, coverage=20, power=0.80)
    assert isinstance(n_needed, int)
    assert n_needed >= 2


def test_sex_check_unimodal_cohort_falls_back_to_threshold():
    """P1-9: on a synthetic all-female cohort (unimodal chrX-beta
    distribution), the dip-test must trigger and the function must fall
    back to the fixed 0.25 threshold with a UserWarning."""
    rng = np.random.default_rng(0)
    # All female: chrX beta clustered near 0.45 (unimodal, well above 0.25).
    beta_values = rng.normal(0.45, 0.02, size=8)
    sample_ids = [f"S{i}" for i in range(8)]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = _classify_sex_from_values(sample_ids, beta_values)

    # All should be called female (beta > 0.25).
    assert all(v == "female" for v in result.values()), (
        f"Expected all-female assignment; got {result}"
    )
    user_warns = [w for w in caught if issubclass(w.category, UserWarning)]
    assert user_warns, "Expected UserWarning for single-sex unimodal cohort"
    assert any(
        "single-sex" in str(w.message).lower() or "unimodal" in str(w.message).lower()
        for w in user_warns
    ), (
        "Warning should mention single-sex or unimodal; "
        f"got: {[str(w.message) for w in user_warns]}"
    )


def test_sex_check_bimodal_cohort_uses_clustering():
    """P1-9: bimodal cohort (mixed male+female) must NOT trigger the dip-test
    fallback — largest-gap clustering must be used instead."""
    rng = np.random.default_rng(42)
    # Male: beta ~0.07; Female: beta ~0.45
    male_betas   = rng.normal(0.07, 0.01, size=5)
    female_betas = rng.normal(0.45, 0.02, size=5)
    beta_values  = np.concatenate([male_betas, female_betas])
    sample_ids   = [f"M{i}" for i in range(5)] + [f"F{i}" for i in range(5)]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = _classify_sex_from_values(sample_ids, beta_values)

    # Males (first 5) -> "male", Females (last 5) -> "female"
    males   = [result[f"M{i}"] for i in range(5)]
    females = [result[f"F{i}"] for i in range(5)]
    assert all(s == "male"   for s in males),   f"Expected all-male; got {males}"
    assert all(s == "female" for s in females), f"Expected all-female; got {females}"
    # No unimodal warning should be emitted
    unimodal_warns = [
        w for w in caught
        if issubclass(w.category, UserWarning)
        and ("single-sex" in str(w.message).lower() or "unimodal" in str(w.message).lower())
    ]
    assert not unimodal_warns, (
        f"Unexpected unimodal warning for bimodal cohort: {unimodal_warns}"
    )


def test_tl_qc_opt_in_flags(synth_md_filtered):
    md = synth_md_filtered
    ep.tl.qc(
        md,
        run_sex_check=True,
        run_contamination=True,
        run_sample_correlation=True,
    )
    # New obs columns from opt-in metrics
    expected = {"contamination_score", "min_pairwise_corr"}
    cols = set(md.obs.columns)
    assert expected & cols
    assert "qc_sample_correlation" in md.uns
