"""Tests for Wilson CI / bootstrap CI helpers."""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest


def test_add_wilson_ci_for_tpr_matches_scipy_reference():
    """Wilson 95% CI on a known proportion. Reference: scipy directly."""
    from wilson_bootstrap_ci import add_wilson_ci

    # 90 successes out of 100 -> p_hat = 0.9
    df = pl.DataFrame({
        "tp": [90], "fp": [10], "fn": [10], "tn": [90],
        "tpr": [0.9], "fpr": [10 / 100],
    })
    out = add_wilson_ci(df, rate="tpr", k_col="tp", n_col_expr=lambda d: d["tp"] + d["fn"])
    assert "tpr_ci_lo" in out.columns
    assert "tpr_ci_hi" in out.columns

    # Scipy reference for p_hat=0.9, n=100:
    from scipy.stats import binomtest
    ref = binomtest(90, 100).proportion_ci(method="wilson", confidence_level=0.95)
    assert abs(out["tpr_ci_lo"][0] - ref.low) < 1e-10
    assert abs(out["tpr_ci_hi"][0] - ref.high) < 1e-10


def test_add_wilson_ci_handles_zero_count_edges():
    """k=0 and k=n must not crash and must produce sensible intervals."""
    from wilson_bootstrap_ci import add_wilson_ci

    df = pl.DataFrame({
        "tp": [0, 100], "fp": [0, 0], "fn": [100, 0], "tn": [100, 100],
        "tpr": [0.0, 1.0], "fpr": [0.0, 0.0],
    })
    out = add_wilson_ci(df, rate="tpr", k_col="tp", n_col_expr=lambda d: d["tp"] + d["fn"])
    # k=0/n=100 Wilson lo = 0, hi ~ 0.037
    assert out["tpr_ci_lo"][0] == 0.0
    assert 0.02 < out["tpr_ci_hi"][0] < 0.05
    # k=100/n=100 Wilson lo ~ 0.963, hi = 1.0
    assert 0.95 < out["tpr_ci_lo"][1] < 0.98
    assert out["tpr_ci_hi"][1] == 1.0


def test_add_wilson_ci_zero_denominator_returns_nan():
    """If tp+fn = 0 (no positives in the truth set), CI is NaN."""
    from wilson_bootstrap_ci import add_wilson_ci

    df = pl.DataFrame({
        "tp": [0], "fp": [5], "fn": [0], "tn": [100],
        "tpr": [0.0], "fpr": [5 / 105],
    })
    out = add_wilson_ci(df, rate="tpr", k_col="tp", n_col_expr=lambda d: d["tp"] + d["fn"])
    assert np.isnan(out["tpr_ci_lo"][0])
    assert np.isnan(out["tpr_ci_hi"][0])
