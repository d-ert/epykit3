"""Regression test for GLM separation handling (B3-followup).

When the binomial GLM in ``_glm.irls_binomial_batch`` hits a site where one
covered sample's linear predictor saturates the eta clip bound (logistic
"separation"), the per-sample Pearson denominator ``n * mu * (1 - mu)``
collapses to ``n * _PROP_CLIP``. Without the saturation guard the per-site
Pearson chi-square blows up to ~10^7 at those samples and drives the
chrom-pooled dispersion estimate from O(1) to O(10^6) -- exactly the bug
reported on real data (LR phi=2.5, GLM phi=1.35M on the same tiles).

These tests construct a tiny mixed batch of well-fitted and separated
sites and confirm:

1. Well-fitted sites produce sane finite Pearson values (single digits).
2. Separated sites are NaN'd in deviance, pearson, beta, and se_beta.
3. ``compute_dispersion_phi`` on the mixed batch yields O(1) phi, not the
   millions you get if separated sites contribute.
"""

from __future__ import annotations

import numpy as np

from epykit._glm import irls_binomial_batch, compute_dispersion_phi


def _make_design_treatment_donor() -> np.ndarray:
    """Tiny 6-sample design: intercept + treatment + donor (3 params)."""
    # Sample order:    0    1    2    3    4    5
    # treatment:       0    0    0    1    1    1
    # donor=d2:        0    1    0    0    1    0
    # (donor=d1 is the reference level)
    X = np.array(
        [
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 0.0],
        ],
        dtype=np.float64,
    )
    return X


def test_well_fitted_site_produces_sane_pearson():
    """Site with moderate methylation in both groups -- no separation."""
    X = _make_design_treatment_donor()
    # 50% methylation in controls, 60% in treated. n=20 reads per sample.
    meth = np.array([[10, 10, 10, 12, 12, 12]], dtype=np.int32)
    cov = np.full((1, 6), 20, dtype=np.int32)
    _beta, _se, deviance, pearson, n_eff = irls_binomial_batch(meth, cov, X)
    assert np.isfinite(deviance[0])
    assert np.isfinite(pearson[0])
    assert pearson[0] < 50.0, (
        f"well-fitted Pearson should be small; got {pearson[0]}"
    )
    assert n_eff[0] == 6


def test_separated_site_is_nan():
    """Perfectly methylated treatment, perfectly unmethylated control.

    With ``~ treatment + donor``, the treatment coefficient goes to +infinity;
    IRLS saturates eta at the clip bound. The site must be NaN'd.
    """
    X = _make_design_treatment_donor()
    # control (samples 0-2): 0/20 methylated; treated (3-5): 20/20.
    meth = np.array([[0, 0, 0, 20, 20, 20]], dtype=np.int32)
    cov = np.full((1, 6), 20, dtype=np.int32)
    beta, se, deviance, pearson, n_eff = irls_binomial_batch(meth, cov, X)
    assert np.isnan(deviance[0]), "separated site must NaN deviance"
    assert np.isnan(pearson[0]), "separated site must NaN Pearson"
    assert np.all(np.isnan(beta[0])), "separated site must NaN beta"
    assert np.all(np.isnan(se[0])), "separated site must NaN se_beta"


def test_dispersion_pool_excludes_separated():
    """Chrom-pooled phi must remain O(1) when one site separates."""
    X = _make_design_treatment_donor()
    rng = np.random.default_rng(0)
    n_good = 200
    # Good sites: binomial draws around p=0.5 with small treatment effect.
    p_ctrl, p_treat = 0.45, 0.55
    cov = np.full((n_good + 1, 6), 25, dtype=np.int32)
    meth = np.zeros_like(cov)
    for i in range(n_good):
        meth[i, :3] = rng.binomial(25, p_ctrl, size=3)
        meth[i, 3:] = rng.binomial(25, p_treat, size=3)
    # One separated site appended at the end.
    meth[-1] = [0, 0, 0, 25, 25, 25]

    _beta, _se, _dev, pearson, n_eff = irls_binomial_batch(meth, cov, X)
    df_per_site = n_eff.astype(np.float64) - float(X.shape[1])

    _phi_eff, phi_hat = compute_dispersion_phi(
        pearson_per_site=pearson,
        df_per_site=df_per_site,
        dispersion="chrom",
        chrom_name="test",
    )
    # Without the saturation guard the separated site alone would push
    # phi_hat well into the thousands. With the guard, the pool only
    # contains the n_good well-fitted sites and phi should sit near 1-3.
    assert 0.5 < phi_hat < 10.0, (
        f"chrom-pooled phi should be O(1); got {phi_hat}. The separation "
        f"guard in irls_binomial_batch is probably no longer NaN-ing "
        f"saturated sites."
    )


def test_dispersion_handles_all_separated():
    """If every site separates, dispersion falls back to min_dispersion=1.0."""
    X = _make_design_treatment_donor()
    meth = np.tile(np.array([[0, 0, 0, 20, 20, 20]], dtype=np.int32), (50, 1))
    cov = np.full((50, 6), 20, dtype=np.int32)
    _beta, _se, _dev, pearson, n_eff = irls_binomial_batch(meth, cov, X)
    df_per_site = n_eff.astype(np.float64) - float(X.shape[1])
    _phi_eff, phi_hat = compute_dispersion_phi(
        pearson_per_site=pearson,
        df_per_site=df_per_site,
        dispersion="chrom",
        chrom_name="test_all_separated",
    )
    # No usable sites -> fallback to min_dispersion (1.0).
    assert phi_hat == 1.0
