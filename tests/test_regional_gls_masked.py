"""Exact observed-submatrix GLS oracles for incomplete biological samples."""

import numpy as np
import pytest
from scipy.linalg import block_diag
from scipy.optimize import minimize_scalar
from statsmodels.regression.linear_model import GLS

from epykit._regional_gls_native import evaluate_covariance, fit_regions, profile_covariance


def masked_fixture(seed=891):
    rng = np.random.default_rng(seed)
    p = np.cumsum(rng.integers(5, 50, 8)).astype(float)
    h = np.sqrt(rng.integers(5, 50, 8)).astype(float)
    g = np.r_[np.zeros(3), np.ones(3)]
    mask = rng.random((6, 8)) > .3
    mask[[0, 3], :] = True
    y = rng.normal(size=(6, 8))
    y[~mask] = 0
    y[0, 0] = 0  # An observed zero residual is NOT missing.
    return p, y, h, g, mask


def dense_masked_oracle(p, y, h, g, mask, length, c):
    k = len(p)
    car = np.exp(-abs(p[:, None]-p[None, :])/length) if length else np.eye(k)
    covariance = (1-c)*car + c*np.outer(h, h)/np.median(h*h)
    responses, designs, blocks = [], [], []
    for s in range(len(g)):
        indices = np.flatnonzero(mask[s])
        responses.extend(y[s, indices])
        designs.append(np.column_stack([np.eye(k)[indices]*h[indices, None], g[s]*h[indices]]))
        blocks.append(covariance[np.ix_(indices, indices)])
    x = np.vstack(designs)
    sigma = block_diag(*blocks)
    fit = GLS(np.array(responses), x, sigma=sigma).fit()
    info = x.T @ np.linalg.solve(sigma, x)
    objective = (fit.df_resid*np.log(fit.ssr/fit.df_resid)
                 + np.linalg.slogdet(sigma)[1] + np.linalg.slogdet(info)[1])
    return objective, fit.params[-1], fit.bse[-1], fit.df_resid


@pytest.mark.parametrize("length,c", [(0., 0.), (10., .3), (120., .3), (10000., .98)])
def test_masked_fixed_covariance_matches_observed_dense_gls(length, c):
    p, y, h, g, mask = masked_fixture()
    got = evaluate_covariance(p, y, h, g, length, c, obs_mask=mask)
    want = dense_masked_oracle(p, y, h, g, mask, length, c)
    assert got == pytest.approx(want[:3], rel=2e-8, abs=2e-8)


def test_absent_sample_cells_do_not_attenuate_effect_or_invent_df():
    p, h = np.array([0., 20., 40.]), np.ones(3)
    g = np.array([0., 0., 0., 1., 1., 1.])
    y = np.tile(np.array([-.3, -.1, 0., .1, .3, 0.])[:, None], (1, 3))
    mask = np.ones_like(y, dtype=bool)
    mask[[2, 5], :] = False
    got = evaluate_covariance(p, y, h, g, 0., 0., obs_mask=mask)
    assert got[1] == pytest.approx(.4)
    assert got[2] == pytest.approx(.07071067811865474)
    assert dense_masked_oracle(p, y, h, g, mask, 0., 0.)[3] == 8


def test_observed_zero_cells_retain_their_design_information():
    p, y, h, g, mask = masked_fixture()
    got = evaluate_covariance(p, y, h, g, 80., .1, obs_mask=mask)
    want = dense_masked_oracle(p, y, h, g, mask, 80., .1)
    assert got == pytest.approx(want[:3], rel=1e-9)


def test_unidentifiable_masked_treatment_is_rejected_explicitly():
    p, y, h, g, mask = masked_fixture()
    mask[:] = False
    mask[:3, :4] = True
    mask[3:, 4:] = True
    with pytest.raises(ValueError, match="identifiable"):
        evaluate_covariance(p, y, h, g, 0., 0., obs_mask=mask)


def test_masked_donor_profile_matches_numerical_observed_likelihood():
    p, y, h, g, mask = masked_fixture()
    got = profile_covariance(p, y, h, g, 120., obs_mask=mask)
    def objective(c):
        return dense_masked_oracle(p, y, h, g, mask, 120., c)[0]
    numeric = minimize_scalar(objective, bounds=(0., .98), method="bounded",
                              options={"xatol": 1e-11})
    best = min(numeric.fun, objective(0.), objective(.98))
    assert got[0] == pytest.approx(best, abs=2e-7)
    assert (got[0], got[2], got[3]) == pytest.approx(
        dense_masked_oracle(p, y, h, g, mask, 120., got[1])[:3], abs=2e-7
    )


def test_masked_donor_profile_does_not_miss_a_minimum_near_zero():
    rng = np.random.default_rng(985)
    p, h = np.arange(8, dtype=float)*20, np.ones(8)
    g = np.r_[np.zeros(3), np.ones(3)]
    for _ in range(7):
        mask = rng.random((6, 8)) > .2
        mask[[0, 3]] = True
        y = rng.normal(size=(6, 8)) + .5*rng.normal(size=(6, 1))
        y[~mask] = 0
    got = profile_covariance(p, y, h, g, 0., obs_mask=mask)
    assert got[1] == pytest.approx(.04138004, abs=2e-6)
    assert got[0] == pytest.approx(16.449197738, abs=2e-7)


def test_masked_batch_parallel_and_mirror_match():
    p, y, h, g, mask = masked_fixture()
    bounds = [[0, len(p)]]
    a = fit_regions(p, y, h, g, bounds, obs_mask=mask, workers=1)
    b = fit_regions(p, y, h, 1-g, bounds, obs_mask=mask, workers=2)
    assert a[:, 2] == pytest.approx(-b[:, 2], abs=1e-7)
    assert a[:, [0, 1, 3, 4]] == pytest.approx(b[:, [0, 1, 3, 4]], rel=1e-6)


def test_complete_mask_retains_complete_model_effect_and_standard_error():
    p, y, h, g, _ = masked_fixture()
    y -= y.mean(axis=0)
    mask = np.ones_like(y, bool)
    old = evaluate_covariance(p, y, h, g, 120., .3)
    new = evaluate_covariance(p, y, h, g, 120., .3, obs_mask=mask)
    assert new[1:] == pytest.approx(old[1:], abs=1e-9)


def test_count_level_region_fit_uses_only_positive_coverage_cells():
    from epykit._regional_gls import fit_region

    p, _, h, g, mask = masked_fixture()
    rng = np.random.default_rng(892)
    n = np.full((len(p), len(g)), 30.)
    n[~mask.T] = 0
    m = rng.binomial(n.astype(int), .4).astype(float)
    y = np.arcsin(2*(m+.05)/(np.maximum(n, 1)+.1)-1).T
    h = np.sqrt(np.full(len(p), 30.))
    want = dense_masked_oracle(p, y*h, h, g, mask, 120., .3)
    got = fit_region(p, m, n, g, covariance=(120., .3))
    assert got["effect_transformed"] == pytest.approx(want[1], abs=1e-9)
    assert got["standard_error"] == pytest.approx(want[2], abs=1e-9)
    assert got["residual_df"] == want[3]
