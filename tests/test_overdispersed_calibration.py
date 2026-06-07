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
import pytest

from epykit.dmc import _score_finalize, _beta_binom_mom_from_welford


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


# ---------------------------------------------------------------------------
# welch_t calibration across realistic overdispersion (D18)
#
# welch_t derives its variance from between-replicate scatter rather than a
# dispersion model, so it should adapt to overdispersion automatically. This
# pins that it is not anti-conservative across phi 1.5-5 (the real-WGBS range
# CLAUDE.md names as the calibration risk) -- the lr tests above cover the
# count engine; welch_t had no overdispersed FPR test before.
# ---------------------------------------------------------------------------

def _kappa_for_phi(phi: float, cov: int) -> float:
    """Beta-binomial concentration kappa giving target dispersion phi.

    phi ~ 1 + (cov - 1) * ICC and ICC = 1 / (kappa + 1), so
    kappa = (cov - 1) / (phi - 1) - 1.
    """
    return (cov - 1) / (phi - 1.0) - 1.0


def _welford_group(p, n_rep, cov, kappa, rng):
    """Per-replicate Welford accumulators (mean, M2, n_valid) of beta under a
    beta-binomial null with concentration ``kappa``."""
    n_sites = p.shape[0]
    a = p * kappa
    b = (1.0 - p) * kappa
    betas = np.empty((n_rep, n_sites), dtype=np.float64)
    for r in range(n_rep):
        rate = rng.beta(a, b)
        betas[r] = rng.binomial(cov, rate).astype(np.float64) / cov
    mean = betas.mean(axis=0)
    M2 = ((betas - mean) ** 2).sum(axis=0)
    nv = np.full(n_sites, n_rep, dtype=np.int32)
    return mean, M2, nv


def _welch_fpr(phi, *, n_rep=3, cov=30, n_sites=8000, seed=11):
    rng = np.random.default_rng(seed)
    p = rng.uniform(0.2, 0.8, size=n_sites)  # H0: both groups same per-site mean
    kappa = _kappa_for_phi(phi, cov)
    mc, M2c, nc = _welford_group(p, n_rep, cov, kappa, rng)
    mk, M2k, nk = _welford_group(p, n_rep, cov, kappa, rng)
    pvals, _ = _beta_binom_mom_from_welford(mc, M2c, nc, mk, M2k, nk)
    finite = np.isfinite(pvals)
    return float(np.mean(pvals[finite] < 0.05)), int(finite.sum())


@pytest.mark.parametrize("phi", [1.5, 2.0, 3.0, 5.0])
def test_welch_t_not_anticonservative_across_overdispersion(phi):
    fpr05, n_finite = _welch_fpr(phi, n_rep=3)
    assert n_finite > 2000, f"too few finite p-values ({n_finite}) at phi={phi}"
    # welch_t adapts variance from replicates -> should stay near nominal.
    # Upper bound catches a genuinely anti-conservative engine; the
    # generous ceiling reflects that welch_t is a documented weak fallback
    # (t(~4) at n=3, mild boundary skew). Do NOT loosen to make a future
    # regression pass -- an anti-conservative welch_t is a finding (D18).
    assert fpr05 <= 0.10, (
        f"welch_t is anti-conservative at phi={phi}: FPR@0.05={fpr05:.4f} (> 0.10)"
    )


# ---------------------------------------------------------------------------
# Real-engine (store-backed) overdispersed null for glm + welch_t (D18).
#
# test_null_calibration.py runs glm/welch_t under the null but at the synth
# generator's modest default replicate_sd (~phi 1.1). This drives the *real*
# engines through a methylstore at an elevated replicate_sd so the count
# model genuinely sees overdispersion, closing the "no overdispersed FPR
# test for glm/welch_t" gap. Slow (builds a store + runs the engine).
# ---------------------------------------------------------------------------

@pytest.mark.slow
@pytest.mark.parametrize("name,kwargs", [
    ("glm", dict(formula="~ group", contrast="group")),
    ("welch_t", dict(test="welch_t")),
])
def test_engine_not_anticonservative_overdispersed_store(tmp_path, name, kwargs):
    import epykit as ep
    from tests.fixtures.synth import SimConfig, generate

    cfg = SimConfig(
        n_per_group=4,
        n_scattered_dmcs=0, n_dmrs=0,
        dmc_effect=0.0, dmr_effect=0.0,
        replicate_sd=0.10,           # ~3x the default -> real overdispersion
        seed=20260608,
    )
    res = generate(cfg, tmp_path / "od")
    md = ep.read_bismark(
        res["samplesheet"], treatment_group="treatment", control_group="control",
        store_dir=str(tmp_path / "od_store"),
    )
    ep.pp.filter_coverage(md, lo_count=5, hi_perc=99.9)
    ep.pp.set_unite_type(md, type="intersect")
    ep.tl.dmc(md, min_samples_treatment=0, min_samples_control=0, **kwargs)

    p = md.dmc["pvalue"].to_numpy()
    p = p[np.isfinite(p)]
    assert p.size > 1000, f"{name}: too few finite p-values ({p.size})"
    fpr05 = float((p < 0.05).mean())
    fpr01 = float((p < 0.01).mean())
    # Upper bound only: a conservative engine is fine (and is what the paper
    # claims); anti-conservative is the bug. Don't loosen on a future
    # regression -- it's a finding (D18).
    assert fpr05 <= 0.085, (
        f"{name} anti-conservative under overdispersed null: "
        f"FPR@0.05={fpr05:.4f} (>0.085), FPR@0.01={fpr01:.4f}"
    )
