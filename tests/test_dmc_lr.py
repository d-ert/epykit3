"""P1-11: log2_odds_ratio renamed per backend.

lr / fisher: log2_odds_ratio_pooled (semantics unchanged, name clearer).
glm:         coef_treatment_log2  (was always the logit coefficient, not log2(OR)).

Both backends emit a transitional log2_odds_ratio column NaN-filled
with a FutureWarning once per tl.dmc call.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

import epykit as ep


def test_lr_emits_log2_odds_ratio_pooled(synth_md_filtered):
    md = synth_md_filtered
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ep.tl.dmc(md, test="lr")
    df = md.dmc
    assert "log2_odds_ratio_pooled" in df.columns, (
        f"missing log2_odds_ratio_pooled; got {df.columns}"
    )
    assert "log2_odds_ratio" in df.columns, (
        "transitional log2_odds_ratio column must be present (NaN-filled)"
    )
    legacy = df["log2_odds_ratio"].to_numpy()
    assert np.isnan(legacy).all(), (
        f"transitional column must be all-NaN; got {legacy[:5]}"
    )
    fut = [w for w in caught if issubclass(w.category, FutureWarning)]
    assert fut and "log2_odds_ratio" in str(fut[0].message), (
        f"expected FutureWarning mentioning log2_odds_ratio; got {[str(w.message) for w in fut]}"
    )


def test_fisher_emits_log2_odds_ratio_pooled(synth_md_filtered):
    md = synth_md_filtered
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ep.tl.dmc(md, test="fisher")
    df = md.dmc
    assert "log2_odds_ratio_pooled" in df.columns
    assert "log2_odds_ratio" in df.columns
    assert np.isnan(df["log2_odds_ratio"].to_numpy()).all()
    fut = [w for w in caught if issubclass(w.category, FutureWarning)]
    assert fut and "log2_odds_ratio" in str(fut[0].message), (
        f"expected FutureWarning for fisher backend; got {[str(w.message) for w in fut]}"
    )


def test_lr_meth_diff_ci_is_asymmetric_near_boundary(synth_md_filtered):
    """P1-3: lr emits Newcombe (asymmetric) CIs, not Wald (symmetric).

    At sites where mean_beta is near 0 or 1, Newcombe CIs are meaningfully
    asymmetric: |hi - point| != |point - lo| by more than 1e-4. With the old
    Wald CI both half-widths are identical (symmetric) to floating-point
    precision (~1e-7). This test checks that the *meaningful* asymmetry
    threshold (> 1e-4) is exceeded at >30% of near-boundary sites.
    """
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    df = md.dmc

    # Basic shape checks.
    assert "meth_diff_ci_lo" in df.columns
    assert "meth_diff_ci_hi" in df.columns

    lo = df["meth_diff_ci_lo"].to_numpy()
    hi = df["meth_diff_ci_hi"].to_numpy()
    point = df["meth_diff"].to_numpy()
    valid = np.isfinite(lo) & np.isfinite(hi) & np.isfinite(point)

    # CIs must bracket the point estimate.
    assert ((lo[valid] <= point[valid] + 1e-4) & (point[valid] - 1e-4 <= hi[valid])).all(), (
        "CI does not bracket meth_diff"
    )

    # CIs must be within [-1, 1].
    assert (lo[valid] >= -1.0 - 1e-4).all() and (hi[valid] <= 1.0 + 1e-4).all(), (
        "CI outside [-1, 1]"
    )

    # Near-boundary sites (mean_beta_treat or mean_beta_ctrl < 0.1 or > 0.9):
    # Newcombe is meaningfully asymmetric there: the two half-widths differ
    # by > 1e-4 (not just float noise).  Wald produces symmetric CIs where
    # both half-widths agree to ~1e-7.
    mean_t = df["mean_beta_case"].to_numpy()
    mean_c = df["mean_beta_control"].to_numpy()
    boundary = (mean_t < 0.1) | (mean_t > 0.9) | (mean_c < 0.1) | (mean_c > 0.9)
    boundary_valid = boundary & valid
    if boundary_valid.sum() > 10:
        width_hi = hi[boundary_valid] - point[boundary_valid]
        width_lo = point[boundary_valid] - lo[boundary_valid]
        # Meaningful asymmetry: at least 1e-4 difference in the two half-widths.
        # Wald produces < 1e-7 difference (floating-point noise only).
        asymmetric = np.abs(width_hi - width_lo) > 1e-4
        assert asymmetric.mean() > 0.3, (
            f"Expected >30% of {boundary_valid.sum()} near-boundary sites "
            f"to have meaningfully asymmetric CIs (|width_hi - width_lo| > 1e-4); "
            f"got {asymmetric.mean():.1%}. This suggests Wald CI is still being "
            f"used instead of Newcombe."
        )
