"""Independent GLM oracles and isolation of unidentifiable sites."""

import numpy as np
import pytest
import statsmodels.api as sm

from epykit import _glm


def fixture():
    rng = np.random.default_rng(3841)
    x = np.column_stack([np.ones(10), np.repeat([0., 1.], 5), np.linspace(-1, 1, 10)])
    n = rng.integers(40, 80, (12, 10))
    m = rng.binomial(n, 1 / (1 + np.exp(-(x @ [.2, .5, -.4]))))
    return x, m, n


def test_glm_valid_sites_match_statsmodels_with_joint_covariance():
    x, m, n = fixture()
    beta, se, dev, pearson, count, covariance = _glm.irls_binomial_batch(m, n, x, return_cov=True)
    for i in range(len(m)):
        ref = sm.GLM(np.column_stack([m[i], n[i]-m[i]]), x, family=sm.families.Binomial()).fit(tol=1e-12)
        np.testing.assert_allclose(beta[i], ref.params, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(se[i], ref.bse, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(covariance[i], ref.cov_params(), rtol=1e-6, atol=1e-8)
        assert dev[i] == pytest.approx(ref.deviance, rel=1e-6, abs=1e-8)
        assert pearson[i] == pytest.approx(ref.pearson_chi2, rel=1e-6, abs=1e-8)
    assert count.tolist() == [10]*len(m)


def test_covered_design_rank_matters_even_with_enough_observations():
    # Covered rows have three observations but treatment is unidentifiable.
    x = np.column_stack([np.ones(6), [0., 0., 0., 1., 1., 1.]])
    n = np.array([[30, 30, 30, 0, 0, 0], [30]*6, [0]*6])
    m = np.array([[10, 12, 15, 0, 0, 0], [10, 12, 15, 16, 19, 17], [0]*6])
    result = _glm.irls_binomial_batch(m, n, x, return_cov=True)
    for array in [*result[:4], result[5]]:
        assert np.isnan(array[[0, 2]]).all()
        assert np.isfinite(array[1]).all()
    assert result[4].tolist() == [3, 6, 0]


def test_one_singular_system_does_not_scalar_solve_every_valid_site(monkeypatch):
    x = np.column_stack([np.ones(4), [0., 0., 1., 1.]])
    w = np.ones((101, 4))
    w[-1, 2:] = 0
    z = np.tile([-.2, .1, .2, .3], (101, 1))
    calls = []
    original = np.linalg.solve

    def counted(a, b):
        calls.append(a.shape[0] if a.ndim == 3 else 1)
        return original(a, b)

    monkeypatch.setattr(np.linalg, 'solve', counted)
    actual = _glm._solve_weighted_lsq(x, z, w)
    np.testing.assert_allclose(actual[:-1], np.tile([-.05, .3], (100, 1)), atol=1e-12)
    assert np.isnan(actual[-1]).all()
    assert calls.count(1) < 20


def test_stacked_solver_isolates_singular_and_nonfinite_systems():
    a = np.tile(np.eye(2), (103, 1, 1))
    a[17] = 0
    a[32] = np.nan
    b = np.broadcast_to(np.eye(2), a.shape).copy()
    got = _glm._solve_stacked_systems(a, b)
    assert np.isnan(got[[17, 32]]).all()
    np.testing.assert_array_equal(got[np.setdiff1d(np.arange(103), [17, 32])],
                                  np.tile(np.eye(2), (101, 1, 1)))


def test_estimability_keeps_full_rank_ill_conditioned_design():
    x = np.column_stack([np.ones(6), np.arange(6)*1e-6])
    mask = np.ones((3, 6), bool)
    mask[1] = False
    mask[2, 1:] = False
    assert _glm._estimable_sites(x, mask).tolist() == [True, False, False]
