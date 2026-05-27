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


def test_bootstrap_auroc_ci_contains_point_estimate():
    """Bootstrap CI for AUROC must contain the original point estimate."""
    from wilson_bootstrap_ci import bootstrap_auroc_ci

    rng = np.random.default_rng(42)
    n = 1000
    is_dmc = rng.random(n) < 0.2
    # Pvalues correlated with is_dmc but noisy.
    pvalues = np.where(is_dmc, rng.beta(0.5, 5.0, n), rng.beta(5.0, 0.5, n))
    point = _auroc_reference(is_dmc, pvalues)

    lo, hi = bootstrap_auroc_ci(
        is_dmc=is_dmc, pvalues=pvalues, B=200, seed=42, confidence=0.95,
    )
    assert lo < point < hi, f"CI [{lo:.4f}, {hi:.4f}] does not contain point {point:.4f}"
    assert (hi - lo) < 0.10, f"CI [{lo:.4f}, {hi:.4f}] too wide (>0.1) for n=1000"


def test_bootstrap_auroc_ci_is_deterministic_with_seed():
    """Same seed -> same CI bounds."""
    from wilson_bootstrap_ci import bootstrap_auroc_ci

    rng = np.random.default_rng(7)
    n = 500
    is_dmc = rng.random(n) < 0.2
    pvalues = np.where(is_dmc, rng.beta(0.5, 5.0, n), rng.beta(5.0, 0.5, n))

    lo_a, hi_a = bootstrap_auroc_ci(is_dmc=is_dmc, pvalues=pvalues, B=100, seed=99)
    lo_b, hi_b = bootstrap_auroc_ci(is_dmc=is_dmc, pvalues=pvalues, B=100, seed=99)
    assert lo_a == lo_b
    assert hi_a == hi_b


def test_bootstrap_f1_ci_contains_point_estimate():
    """Same shape for F1 at a fixed q-threshold."""
    from wilson_bootstrap_ci import bootstrap_f1_ci

    rng = np.random.default_rng(1)
    n = 1000
    is_dmc = rng.random(n) < 0.2
    qvalues = np.where(is_dmc, rng.beta(0.5, 5.0, n), rng.beta(5.0, 0.5, n))

    # Compute point F1 at q < 0.05.
    pred = qvalues < 0.05
    tp = int((pred & is_dmc).sum())
    fp = int((pred & ~is_dmc).sum())
    fn = int((~pred & is_dmc).sum())
    point_f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else 0.0

    lo, hi = bootstrap_f1_ci(
        is_dmc=is_dmc, qvalues=qvalues, threshold=0.05, B=200, seed=1,
    )
    assert lo < point_f1 < hi, f"CI [{lo:.4f}, {hi:.4f}] does not contain point {point_f1:.4f}"


# --- helper for tests -------------------------------------------------------


def _auroc_reference(is_dmc: np.ndarray, pvalues: np.ndarray) -> float:
    """Reference AUROC via Mann-Whitney U with average-rank tie handling."""
    score = 1.0 - pvalues
    n_pos = int(is_dmc.sum())
    n_neg = len(is_dmc) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    from scipy.stats import rankdata
    ranks = rankdata(score, method="average")
    sum_ranks_pos = ranks[is_dmc].sum()
    u = sum_ranks_pos - n_pos * (n_pos + 1) / 2
    return u / (n_pos * n_neg)
