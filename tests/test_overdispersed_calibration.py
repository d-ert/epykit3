"""Two-sided null calibration of the lr engine under realistic overdispersion
(M10).

The existing null-calibration suite bounds FPR only from above and runs on the
under-dispersed (phi~0.4) synthetic generator. Real WGBS is phi~1.5-5. This
drives lr's _score_finalize directly on a beta-binomial null (no methylstore,
so it runs in the default PR tier) at phi~4 and pins two things:

  1. dispersion="eb" -- the tl.dmc default -- is well-calibrated (type-I error
     ~ nominal, two-sided) even at n=3 per group. This is the load-bearing
     guarantee: it justifies eb being the default.
  2. dispersion="site" (the noisy raw per-site estimator) is anti-conservative
     at small n, which is exactly what eb's shrinkage toward the chromosome
     pool corrects. Documented and bounded so a regression in either direction
     is caught.
"""
from __future__ import annotations

import numpy as np

from epykit.dmc import _score_finalize


def _betabinom_group(p, n_rep, cov, kappa, rng):
    """Accumulators (sn, sm, sm2n, nv) for one group under a beta-binomial.

    Replicate rates ~ Beta(p*kappa, (1-p)*kappa) give between-replicate
    overdispersion (small kappa -> large phi); meth ~ Binomial(cov, rate).
    """
    n_sites = p.shape[0]
    sn = np.zeros(n_sites, dtype=np.float64)
    sm = np.zeros(n_sites, dtype=np.float64)
    s2 = np.zeros(n_sites, dtype=np.float64)
    nv = np.zeros(n_sites, dtype=np.int32)
    a = p * kappa
    b = (1.0 - p) * kappa
    c = float(cov)
    for _ in range(n_rep):
        rate = rng.beta(a, b)
        m = rng.binomial(cov, rate).astype(np.float64)
        sn += c
        sm += m
        s2 += m * m / c
        nv += 1
    return sn, sm, s2, nv


def _fpr(dispersion, *, n_rep=3, cov=30, kappa=8.0, n_sites=8000, seed=7):
    # kappa=8 at cov=30 implies ICC ~1/9 -> phi ~ 4 (realistic WGBS).
    rng = np.random.default_rng(seed)
    p = rng.uniform(0.2, 0.8, size=n_sites)  # H0: same per-site mean both groups
    gc = _betabinom_group(p, n_rep, cov, kappa, rng)
    kc = _betabinom_group(p, n_rep, cov, kappa, rng)
    out = _score_finalize(
        *gc, *kc, dispersion=dispersion, statistic="lr", reference="adaptive",
    )
    pvals, phi_eff = out[0], out[5]
    finite = np.isfinite(pvals)
    return {
        "fpr05": float(np.mean(pvals[finite] < 0.05)),
        "fpr01": float(np.mean(pvals[finite] < 0.01)),
        "med_phi": float(np.median(phi_eff[finite])),
        "n_finite": int(finite.sum()),
    }


def test_eb_default_is_calibrated_under_overdispersion():
    eb = _fpr("eb", n_rep=3)
    # The fixture must actually be overdispersed, else the test is vacuous.
    assert eb["med_phi"] > 2.0, f"fixture not overdispersed (median phi={eb['med_phi']:.2f})"
    # Default tl.dmc dispersion is ~ nominal at n=3, phi~4 (two-sided).
    assert 0.025 <= eb["fpr05"] <= 0.075, f"eb FPR@.05={eb['fpr05']:.4f}"
    assert 0.003 <= eb["fpr01"] <= 0.020, f"eb FPR@.01={eb['fpr01']:.4f}"


def test_site_is_more_liberal_than_eb_at_small_n():
    eb = _fpr("eb", n_rep=3)
    site = _fpr("site", n_rep=3)
    # The raw per-site phi (4 residual df) is noisy at n=3, so 'site' is
    # anti-conservative; eb shrinkage toward the chrom pool is the better
    # default. Pin the relationship + an upper bound (not the exact value).
    assert site["fpr05"] > eb["fpr05"]
    assert site["fpr05"] <= 0.16, f"site FPR@.05={site['fpr05']:.4f} (wildly off)"


def test_eb_calibrated_across_replicate_counts():
    for n_rep in (2, 6, 10):
        eb = _fpr("eb", n_rep=n_rep)
        assert 0.02 <= eb["fpr05"] <= 0.085, (
            f"eb FPR@.05 at n={n_rep} is {eb['fpr05']:.4f} under phi~4"
        )
