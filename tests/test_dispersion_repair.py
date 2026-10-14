"""Parameter-recovery and null-reference checks independent of EB internals."""

import numpy as np
import polars as pl
import pytest
from scipy.stats import chi2, f

from epykit._glm import compute_dispersion_phi, reference_pvalues
from epykit.dmc import _smyth_eb_moderation, combine_neighbour_pvalues


def test_heterogeneous_df_prior_is_invariant_to_common_dispersion_scale():
    rng = np.random.default_rng(491)
    df = np.tile([2.0, 8.0], 50000)
    true_phi = 5.0 * 10.0 / rng.chisquare(10.0, len(df))
    chi = true_phi * rng.chisquare(df)
    a = _smyth_eb_moderation(chi, df, np.ones(len(df), bool), min_dispersion=.01)
    b = _smyth_eb_moderation(3.0 * chi, df, np.ones(len(df), bool), min_dispersion=.01)
    assert b[2] == pytest.approx(a[2], rel=.002)
    assert b[3] == pytest.approx(3.0 * a[3], rel=.002)


@pytest.mark.parametrize("prior_df", [10.0, 30.0])
def test_inverse_chi_square_prior_parameters_are_recovered(prior_df):
    rng = np.random.default_rng(9302)
    df = np.full(200000, 4.0)
    phi = 5.0 * prior_df / rng.chisquare(prior_df, len(df))
    chi = phi * rng.chisquare(df)
    _, _, estimated_df, estimated_scale = _smyth_eb_moderation(chi, df, np.ones(len(df), bool))
    assert estimated_df == pytest.approx(prior_df, rel=.06)
    assert estimated_scale == pytest.approx(5.0, rel=.025)


def test_high_dispersion_is_not_forced_to_twenty_and_null_reference_is_calibrated():
    rng = np.random.default_rng(491)
    df = np.full(100000, 4.0)
    phi, df_mod, _, prior = _smyth_eb_moderation(
        50.0 * rng.chisquare(df), df, np.ones(len(df), bool)
    )
    numerator = 50.0 * rng.chisquare(1, len(df))
    p = reference_pvalues(numerator / phi, phi, df_mod, reference="F")
    assert prior == pytest.approx(50, rel=.025)
    assert np.isfinite(p).all()
    assert .04 < np.mean(p < .05) < .06


@pytest.mark.parametrize("reference", ["adaptive", "F"])
def test_infinite_prior_df_uses_limiting_chi_square_reference(reference):
    stat = np.array([1.0, 4.0, 10.0])
    got = reference_pvalues(stat, np.full(3, 5.0), np.full(3, np.inf), reference)
    assert got == pytest.approx(chi2.sf(stat, 1))


def test_raw_residual_is_pooled_before_posterior_dispersion_floor():
    rng = np.random.default_rng(121)
    df = np.full(10000, 4.0)
    chi = 5 * 10 / rng.chisquare(10, len(df)) * rng.chisquare(df)
    chi[0] = .1
    phi, df_mod, d0, prior = _smyth_eb_moderation(chi, df, np.ones(len(df), bool))
    want = max(1.0, (chi[0] + d0 * prior) / df_mod[0])
    assert phi[0] == pytest.approx(want)


def test_glm_dispersion_uses_same_uncapped_eb_prior():
    rng = np.random.default_rng(200)
    df = np.full(100000, 4.0)
    _, scale, _ = compute_dispersion_phi(50 * rng.chisquare(df), df, dispersion="eb")
    assert scale == pytest.approx(50, rel=.03)


def test_isolated_untested_cpg_does_not_count_as_a_contributor():
    data = pl.DataFrame({"chrom": ["chr1", "chr1"], "pos": [100, 1000],
                         "pvalue": [float("nan"), .01], "meth_diff": [.1, .1]})
    result = combine_neighbour_pvalues(data)
    assert result["pvalue_combined_n_neighbours"].to_list() == [0, 1]
    assert np.isnan(result["pvalue_combined"][0])


def test_small_invalid_residual_sets_return_finite_conservative_dispersion():
    phi, df, _, _ = _smyth_eb_moderation(
        np.array([np.nan, -1.0, 0.0, 8.0]), np.array([0.0, 2.0, 2.0, 4.0]),
        np.ones(4, bool), min_disp_sites=100,
    )
    assert np.isfinite(phi).all()
    assert (phi >= 1).all()
    assert np.isfinite(df).all()
