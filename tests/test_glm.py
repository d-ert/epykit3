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
