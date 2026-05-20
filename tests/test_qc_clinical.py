"""clinical / cohort QC pack."""

from __future__ import annotations

import polars as pl
import pytest

import epykit as ep
from epykit.qc import (
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
