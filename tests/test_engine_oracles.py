"""Oracle tests for the DMC engines' *inferential* output (M9).

The engines' shapes and aggregate power/FDR are tested elsewhere; this module
pins the actual p-values / SEs against independent references (scipy /
statsmodels), which is what a methods reviewer asks for. These are fast
(direct function calls, no methylstore), so they run in the default PR tier
on every platform -- the calibration/recovery suites that need a store stay
`slow`-marked.

Oracles
-------
* welch_t  -> scipy.stats.ttest_ind(equal_var=False) on the per-replicate betas
* lr       -> the 2x2 G-test (likelihood-ratio), scipy.stats.chi2_contingency
              with lambda_="log-likelihood", at the dispersion-clamped (phi=1)
              regime where the adaptive reference is chi2(1)
* glm      -> statsmodels GLM(Binomial): coefficient SE and the deviance-LR
              statistic of full-vs-reduced
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy import stats as sp_stats

from epykit.dmc import (
    _beta_binom_mom_from_welford,
    _score_finalize,
    _welford_init,
    _welford_update,
)

# --- welch_t vs scipy.stats.ttest_ind(equal_var=False) ----------------------

def test_welch_t_matches_scipy_ttest_ind():
    # Per-replicate (meth, cov); betas = meth/cov.
    case = [(7, 10), (8, 10), (9, 10)]      # betas 0.7, 0.8, 0.9
    ctrl = [(2, 10), (3, 10), (1, 10)]      # betas 0.2, 0.3, 0.1

    mean_c, M2_c, nv_c = _welford_init(1)
    mean_k, M2_k, nv_k = _welford_init(1)
    for m, c in case:
        _welford_update(mean_c, M2_c, nv_c,
                        np.array([m], dtype=np.int32), np.array([c], dtype=np.int32))
    for m, c in ctrl:
        _welford_update(mean_k, M2_k, nv_k,
                        np.array([m], dtype=np.int32), np.array([c], dtype=np.int32))

    p_epykit = _beta_binom_mom_from_welford(mean_c, M2_c, nv_c, mean_k, M2_k, nv_k)[0][0]

    betas_c = [m / c for m, c in case]
    betas_k = [m / c for m, c in ctrl]
    p_scipy = sp_stats.ttest_ind(betas_c, betas_k, equal_var=False).pvalue

    assert np.isclose(p_epykit, p_scipy, rtol=1e-9, atol=1e-12)


# --- lr quasi-binomial LRT vs the 2x2 G-test --------------------------------

def _score_p(meth_cov_case, meth_cov_ctrl):
    """Build _score_finalize accumulators from per-replicate counts and return
    the lr p-value at site 0."""
    sn_c = np.array([sum(c for _, c in meth_cov_case)], dtype=np.float64)
    sm_c = np.array([sum(m for m, _ in meth_cov_case)], dtype=np.float64)
    s2_c = np.array([sum(m * m / c for m, c in meth_cov_case)], dtype=np.float64)
    nv_c = np.array([len(meth_cov_case)], dtype=np.int32)
    sn_k = np.array([sum(c for _, c in meth_cov_ctrl)], dtype=np.float64)
    sm_k = np.array([sum(m for m, _ in meth_cov_ctrl)], dtype=np.float64)
    s2_k = np.array([sum(m * m / c for m, c in meth_cov_ctrl)], dtype=np.float64)
    nv_k = np.array([len(meth_cov_ctrl)], dtype=np.int32)
    out = _score_finalize(
        sn_c, sm_c, s2_c, nv_c, sn_k, sm_k, s2_k, nv_k,
        dispersion="site", statistic="lr", reference="adaptive",
    )
    pvals, _log2or, _pc, _pk, _phi_hat, phi_eff, _df_phi = out
    return float(pvals[0]), float(phi_eff[0])


def test_lr_matches_g_test_at_unit_dispersion():
    # Identical replicates within each group -> zero between-replicate
    # variance -> per-site dispersion clamps to phi=1 -> adaptive reference is
    # chi2(1), so the lr statistic is exactly the 2x2 G-test (LRT).
    case = [(8, 10), (8, 10), (8, 10)]   # pooled 24 / 30
    ctrl = [(2, 10), (2, 10), (2, 10)]   # pooled  6 / 30
    p_epykit, phi_eff = _score_p(case, ctrl)
    assert phi_eff == pytest.approx(1.0), "fixture must clamp dispersion to 1"

    table = [[24, 6], [6, 24]]  # [meth, unmeth] x [case, ctrl]
    _g, p_gtest, _dof, _exp = sp_stats.chi2_contingency(
        table, lambda_="log-likelihood", correction=False
    )
    assert np.isclose(p_epykit, p_gtest, rtol=1e-9, atol=1e-12)


def test_lr_null_table_is_nonsignificant():
    # Equal proportions -> LRT ~ 0 -> p ~ 1.
    case = [(5, 10), (5, 10), (5, 10)]
    ctrl = [(5, 10), (5, 10), (5, 10)]
    p_epykit, _phi = _score_p(case, ctrl)
    assert p_epykit == pytest.approx(1.0, abs=1e-9)


# --- GLM SE + deviance-LR vs statsmodels GLM(Binomial) ----------------------

def test_glm_se_and_deviance_match_statsmodels():
    sm = pytest.importorskip("statsmodels.api")
    from epykit._glm import irls_binomial_batch

    meth = np.array([7, 8, 9, 2, 3, 1], dtype=np.int32)
    cov  = np.array([10, 10, 10, 10, 10, 10], dtype=np.int32)
    # Design: intercept + treatment (coef index 1).
    X = np.array(
        [[1, 1], [1, 1], [1, 1], [1, 0], [1, 0], [1, 0]], dtype=np.float64
    )
    X_red = X[:, [0]]

    beta, se, dev_full, _pearson, _n_eff = irls_binomial_batch(
        meth[None, :], cov[None, :], X
    )
    _b2, _s2, dev_red, _p2, _n2 = irls_binomial_batch(meth[None, :], cov[None, :], X_red)

    # statsmodels reference: endog = [successes, failures].
    endog = np.column_stack([meth, cov - meth]).astype(float)
    full = sm.GLM(endog, X, family=sm.families.Binomial()).fit()
    red  = sm.GLM(endog, X_red, family=sm.families.Binomial()).fit()

    # Coefficient + its Wald SE on the treatment term.
    assert np.isclose(beta[0, 1], full.params[1], rtol=1e-4, atol=1e-5)
    assert np.isclose(se[0, 1], full.bse[1], rtol=1e-3, atol=1e-5)

    # Deviance-LR statistic (full vs reduced) -> chi2(1) p-value.
    lr_epykit = float(dev_red[0] - dev_full[0])
    lr_sm = float(red.deviance - full.deviance)
    assert np.isclose(lr_epykit, lr_sm, rtol=1e-3, atol=1e-4)
    p_epykit = float(sp_stats.chi2.sf(lr_epykit, df=1))
    p_sm = float(sp_stats.chi2.sf(lr_sm, df=1))
    assert np.isclose(p_epykit, p_sm, rtol=1e-3, atol=1e-6)
