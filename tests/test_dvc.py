"""P1-7: DVC per-site variance test uses Brown-Forsythe (median-centred
Levene), not Bartlett. Verified against scipy.stats.levene(center='median')."""
from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import levene


def test_dvc_per_site_stat_matches_brown_forsythe():
    """Brown-Forsythe F-stat must match scipy.stats.levene(center='median')
    on a synthetic bimodal beta dataset to 1e-6 rtol."""
    from epykit.dvc import _per_site_variance_test

    rng = np.random.default_rng(0)
    n_sites = 50
    n_per_group = 6
    out_f = np.empty(n_sites)
    ref_f = np.empty(n_sites)
    for i in range(n_sites):
        a = rng.beta(0.5, 0.5, size=n_per_group)
        b = rng.beta(0.5, 0.5, size=n_per_group)
        result = _per_site_variance_test(a, b)
        # result is (f_stat, p_val)
        out_f[i] = result[0]
        ref_f[i] = levene(a, b, center="median").statistic
    np.testing.assert_allclose(
        out_f, ref_f, rtol=1e-6, atol=1e-9,
        err_msg="Brown-Forsythe F-stat must match scipy levene(center='median')",
    )


def test_dvc_per_site_nan_on_insufficient_data():
    """_per_site_variance_test returns (nan, nan) when either group < 2 obs."""
    from epykit.dvc import _per_site_variance_test

    f, p = _per_site_variance_test(np.array([0.5]), np.array([0.2, 0.8]))
    assert np.isnan(f) and np.isnan(p), "single-obs group_a should give nan"

    f, p = _per_site_variance_test(np.array([0.3, 0.7]), np.array([0.9]))
    assert np.isnan(f) and np.isnan(p), "single-obs group_b should give nan"
