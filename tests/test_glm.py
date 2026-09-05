"""P1-11 GLM half: coef_treatment_log2.
P1-4: reference_level kwarg for patsy Treatment coding.
P1-5: NaN-mask non-converged IRLS sites + log fraction."""
from __future__ import annotations

import logging
import warnings

import numpy as np
import polars as pl
import pytest

import epykit as ep


def test_glm_emits_coef_treatment_log2(synth_md_filtered):
    md = synth_md_filtered
    md.obs = md.obs.with_columns(
        (pl.col("group") == "treatment").cast(int).alias("treatment")
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ep.tl.dmc(md, test="glm", formula="~ treatment")
    df = md.dmc
    assert "coef_treatment_log2" in df.columns, (
        f"missing coef_treatment_log2; got {df.columns}"
    )
    assert "log2_odds_ratio" in df.columns
    legacy = df["log2_odds_ratio"].to_numpy()
    assert np.isnan(legacy).all()
    fut = [w for w in caught if issubclass(w.category, FutureWarning)]
    assert fut and "log2_odds_ratio" in str(fut[0].message)


def test_reference_level_respected(synth_md_filtered, tmp_path):
    """Passing reference_level= sets patsy Treatment coding reference.
    The coefficient sign flips when the reference is swapped."""
    import shutil

    md = synth_md_filtered
    # Add a factor with two levels 'A' and 'B' mapped from the existing group.
    md.obs = md.obs.with_columns(
        pl.when(pl.col("group") == "treatment")
        .then(pl.lit("A"))
        .otherwise(pl.lit("B"))
        .alias("group_ab")
    )
    # Ensure a treatment column (1 for A, 0 for B) is present for the GLM.
    md.obs = md.obs.with_columns(
        (pl.col("group_ab") == "A").cast(int).alias("treatment")
    )

    # Default reference: alphabetical -> 'A' is reference (first alphabetically).
    # With reference='A': coefficient for 'B' vs 'A'.
    ep.tl.dmc(md, test="glm", formula="~ group_ab", contrast="group_ab")
    coef_default = md.dmc["coef_treatment_log2"].to_numpy().copy()

    # Clear the DMC cache so the second call recomputes with a different design.
    from pathlib import Path
    dmc_cache = Path(md.store).parent / "dmc"
    if dmc_cache.exists():
        shutil.rmtree(dmc_cache)

    # Explicit reference='B': coefficient for 'A' vs 'B' -> sign should flip.
    ep.tl.dmc(md, test="glm", formula="~ group_ab", contrast="group_ab",
              reference_level="B")
    coef_swapped = md.dmc["coef_treatment_log2"].to_numpy().copy()

    finite = np.isfinite(coef_default) & np.isfinite(coef_swapped)
    assert finite.sum() > 100, "too few finite sites for a meaningful test"
    np.testing.assert_allclose(
        coef_default[finite], -coef_swapped[finite], atol=1e-9,
        err_msg=(
            "reference_level should flip the sign of the coefficient; "
            "default and swapped should sum to ~0"
        ),
    )


def test_glm_reports_covariate_adjusted_effect_and_ci(synth_md_filtered):
    """M2: the GLM (single-coef contrast) reports a covariate-ADJUSTED effect
    size and a delta-method CI that brackets it -- not the raw marginal
    difference of group means sitting next to an adjusted p-value.

    The distinguishing signal is that ``meth_diff`` no longer equals
    ``mean_beta_case - mean_beta_control`` (which is exactly what the buggy
    pre-M2 code emitted): it is the delta-method central estimate from the
    fitted logit coefficient, and the reported CI brackets *that* value.
    """
    md = synth_md_filtered
    groups = md.obs.get_column("group").to_list()
    md.obs = md.obs.with_columns(
        (pl.col("group") == "treatment").cast(int).alias("treatment")
    )
    # Confound a batch covariate with group (most treatment -> b1, most
    # control -> b2) but cross one sample per group so it isn't collinear.
    seen_t = seen_c = 0
    batch: list[str] = []
    for g in groups:
        if g == "treatment":
            seen_t += 1
            batch.append("b2" if seen_t == 1 else "b1")
        else:
            seen_c += 1
            batch.append("b1" if seen_c == 1 else "b2")
    md.obs = md.obs.with_columns(pl.Series("batch", batch))

    ep.tl.dmc(md, test="glm", formula="~ treatment + batch", contrast="treatment")
    df = md.dmc

    for col in (
        "meth_diff", "meth_diff_ci_lo", "meth_diff_ci_hi",
        "mean_beta_case", "mean_beta_control", "coef_se",
    ):
        assert col in df.columns, f"missing {col}: {df.columns}"

    eff = df["meth_diff"].to_numpy()
    lo = df["meth_diff_ci_lo"].to_numpy()
    hi = df["meth_diff_ci_hi"].to_numpy()
    mbc = df["mean_beta_case"].to_numpy()
    mbk = df["mean_beta_control"].to_numpy()
    raw_marginal = mbc - mbk
    finite = (
        np.isfinite(eff) & np.isfinite(lo) & np.isfinite(hi)
        & np.isfinite(raw_marginal)
    )
    assert finite.sum() > 100, "too few finite sites for a meaningful test"

    # (1) THE fix: the reported effect is the ADJUSTED estimate, not the raw
    # marginal difference of group means. Pre-M2 these were identical (diff
    # == 0); post-M2 they differ at a substantial fraction of sites.
    differs = np.abs(eff[finite] - raw_marginal[finite]) > 1e-3
    assert differs.mean() > 0.3, (
        f"adjusted meth_diff matches the raw marginal at "
        f"{(1 - differs.mean()) * 100:.0f}% of sites -- the delta-method "
        "effect is not wired in"
    )

    # (2) The CI brackets the reported (adjusted) effect and stays in [-1, 1].
    assert np.all(lo[finite] <= eff[finite] + 1e-6)
    assert np.all(eff[finite] <= hi[finite] + 1e-6)
    assert np.all(lo[finite] >= -1.0 - 1e-6) and np.all(hi[finite] <= 1.0 + 1e-6)

    # (3) mean_beta_case / mean_beta_control remain the RAW marginal group
    # means (still proper proportions in [0, 1]); the adjustment lives only in
    # meth_diff / coef_*.
    mb_finite = np.isfinite(mbc) & np.isfinite(mbk)
    assert np.all((mbc[mb_finite] >= -1e-9) & (mbc[mb_finite] <= 1 + 1e-9))
    assert np.all((mbk[mb_finite] >= -1e-9) & (mbk[mb_finite] <= 1 + 1e-9))


def test_glm_meth_diff_ci_widens_with_dispersion():
    """The GLM CI is phi-scaled: delta_method_meth_diff_ci on a dispersion-
    inflated SE is wider than on the raw binomial SE (the consistency that
    makes the CI agree with the phi-corrected p-value)."""
    from epykit._glm import delta_method_meth_diff_ci

    coef = np.full(20, 0.4)
    coef_se = np.full(20, 0.1)
    ref_eta = np.full(20, 0.0)  # control mean = 0.5
    phi = 4.0

    _, hi_binom = delta_method_meth_diff_ci(coef, coef_se, ref_eta=ref_eta)
    _, hi_disp = delta_method_meth_diff_ci(
        coef, np.sqrt(phi) * coef_se, ref_eta=ref_eta,
    )
    # sqrt(4)=2x the SE -> a wider upper bound (same centre).
    assert np.all(hi_disp > hi_binom)


def test_nonconverged_irls_sites_are_nan(synth_md_filtered, caplog):
    """P1-5: non-converged IRLS sites must have NaN Wald statistics
    and a WARNING must be logged when the fraction exceeds 1%."""
    import epykit._glm as _glm_mod

    md = synth_md_filtered
    md.obs = md.obs.with_columns(
        (pl.col("group") == "treatment").cast(int).alias("treatment")
    )

    # Monkeypatch irls_dispatch so every call injects max_iter=1.
    # With only 1 IRLS iteration the vast majority of sites won't converge.
    _real_dispatch = _glm_mod.irls_dispatch

    def _dispatch_max1(meth, cov, X, *, backend="cpu", **kwargs):
        kwargs["max_iter"] = 1
        return _real_dispatch(meth, cov, X, backend=backend, **kwargs)

    _glm_mod.irls_dispatch = _dispatch_max1
    try:
        caplog.set_level(logging.WARNING, logger="epykit._glm")
        ep.tl.dmc(md, test="glm", formula="~ treatment")
    finally:
        _glm_mod.irls_dispatch = _real_dispatch

    df = md.dmc
    pvalues = df["pvalue"].to_numpy(allow_copy=True).astype(float)
    n_nan = int(np.isnan(pvalues).sum())
    n_total = df.height
    assert n_nan > 0, (
        f"Expected NaN p-values for non-converged sites; got 0/{n_total}. "
        "P1-5 fix may not be in place."
    )

    # If more than 1% of sites are non-converged, a WARNING must be logged.
    if n_nan / n_total > 0.01:
        warning_records = [r for r in caplog.records if r.levelname == "WARNING"]
        assert any("converg" in r.message.lower() for r in warning_records), (
            f"Expected WARNING about non-convergence; got: "
            f"{[r.message for r in warning_records]}"
        )
