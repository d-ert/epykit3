"""Wilson 95% CIs on per-cell TPR/FPR and bootstrap CIs on AUROC/F1.

Reads the existing eval_summary.parquet schema:
  tool, scenario, parameter, parameter_value, test, meth_diff_bin,
  threshold_kind, threshold, tp, fp, tn, fn, tpr, fpr, precision, f1, auroc

Adds:
  tpr_ci_lo, tpr_ci_hi (Wilson, 95%)
  fpr_ci_lo, fpr_ci_hi (Wilson, 95%)
  auroc_ci_lo, auroc_ci_hi (bootstrap, 95%, B=1000) -- Task 5
  f1_ci_lo, f1_ci_hi (bootstrap, 95%, B=1000) -- Task 5

This module is pure-Python; it does not call any epykit engine. Re-runs
of the engines are not needed -- CIs operate on the counts already in
eval_summary.parquet.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
import polars as pl
from scipy.stats import binomtest


def _wilson_single(k: int, n: int, confidence: float = 0.95) -> tuple[float, float]:
    """Return Wilson lo/hi for k successes / n trials. NaN/NaN if n == 0."""
    if n == 0:
        return (float("nan"), float("nan"))
    res = binomtest(k, n).proportion_ci(method="wilson", confidence_level=confidence)
    return (float(res.low), float(res.high))


def add_wilson_ci(
    df: pl.DataFrame,
    rate: str,
    k_col: str,
    n_col_expr: Callable[[pl.DataFrame], pl.Series],
    confidence: float = 0.95,
) -> pl.DataFrame:
    """Add `<rate>_ci_lo` and `<rate>_ci_hi` columns via Wilson interval.

    Parameters
    ----------
    df : input DataFrame with `<rate>`, `<k_col>` columns.
    rate : the rate column name (e.g. "tpr", "fpr").
    k_col : the integer-count column for successes (e.g. "tp" for TPR).
    n_col_expr : callable mapping the input DataFrame to a polars Series
        of trial counts (denominators). For TPR: lambda d: d["tp"] + d["fn"].
        For FPR: lambda d: d["fp"] + d["tn"].
    confidence : confidence level, default 0.95.
    """
    k = df[k_col].to_numpy().astype(np.int64)
    n = n_col_expr(df).to_numpy().astype(np.int64)
    lo = np.empty(len(k), dtype=np.float64)
    hi = np.empty(len(k), dtype=np.float64)
    for i in range(len(k)):
        lo[i], hi[i] = _wilson_single(int(k[i]), int(n[i]), confidence)
    return df.with_columns([
        pl.Series(f"{rate}_ci_lo", lo),
        pl.Series(f"{rate}_ci_hi", hi),
    ])


def add_wilson_ci_for_tpr_fpr(
    df: pl.DataFrame, confidence: float = 0.95,
) -> pl.DataFrame:
    """Convenience: add both tpr and fpr Wilson CIs to an eval_summary frame."""
    df = add_wilson_ci(
        df, rate="tpr", k_col="tp",
        n_col_expr=lambda d: d["tp"] + d["fn"],
        confidence=confidence,
    )
    df = add_wilson_ci(
        df, rate="fpr", k_col="fp",
        n_col_expr=lambda d: d["fp"] + d["tn"],
        confidence=confidence,
    )
    return df


def main(argv: list[str] | None = None) -> None:
    """CLI: `python wilson_bootstrap_ci.py --eval eval_summary.parquet --out eval_summary_with_ci.parquet`"""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eval", required=True, type=str, help="Path to eval_summary.parquet")
    parser.add_argument("--out", required=True, type=str, help="Output parquet path")
    args = parser.parse_args(argv)

    df = pl.read_parquet(args.eval)
    df = add_wilson_ci_for_tpr_fpr(df)
    # Bootstrap CIs for AUROC/F1 land in Task 5; this CLI stub leaves them.
    df.write_parquet(args.out)
    print(f"wrote {args.out} with TPR/FPR Wilson CIs added")


# --- Bootstrap CIs for AUROC and F1 ----------------------------------------


def _auroc_mwu(is_dmc: np.ndarray, score: np.ndarray) -> float:
    """AUROC via Mann-Whitney U with proper average-rank tie handling
    (matches scipy.stats.rankdata(method='average') and evaluate.py)."""
    n_pos = int(is_dmc.sum())
    n_neg = len(is_dmc) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    from scipy.stats import rankdata
    ranks = rankdata(score, method="average")
    sum_ranks_pos = ranks[is_dmc].sum()
    u = sum_ranks_pos - n_pos * (n_pos + 1) / 2.0
    return float(u / (n_pos * n_neg))


def bootstrap_auroc_ci(
    is_dmc: np.ndarray,
    pvalues: np.ndarray,
    B: int = 1000,
    seed: int = 0,
    confidence: float = 0.95,
) -> tuple[float, float]:
    """Bootstrap AUROC CI by resampling CpGs with replacement.

    Returns (lo, hi) at the given confidence level (two-sided percentile).
    Uses ``np.nanpercentile`` because a bootstrap draw where all resampled
    CpGs share one class returns NaN from ``_auroc_mwu``; those draws are skipped.
    """
    rng = np.random.default_rng(seed)
    n = len(is_dmc)
    score = 1.0 - np.asarray(pvalues, dtype=np.float64)
    is_dmc = np.asarray(is_dmc, dtype=bool)

    boot = np.empty(B, dtype=np.float64)
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        boot[b] = _auroc_mwu(is_dmc[idx], score[idx])
    alpha = (1.0 - confidence) / 2.0
    lo = float(np.nanpercentile(boot, 100 * alpha))
    hi = float(np.nanpercentile(boot, 100 * (1.0 - alpha)))
    return (lo, hi)


def bootstrap_f1_ci(
    is_dmc: np.ndarray,
    qvalues: np.ndarray,
    threshold: float = 0.05,
    B: int = 1000,
    seed: int = 0,
    confidence: float = 0.95,
) -> tuple[float, float]:
    """Bootstrap F1 CI at a fixed q-threshold, resampling CpGs with replacement."""
    rng = np.random.default_rng(seed)
    n = len(is_dmc)
    is_dmc = np.asarray(is_dmc, dtype=bool)
    pred = np.asarray(qvalues, dtype=np.float64) < threshold

    boot = np.empty(B, dtype=np.float64)
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        tp = int((pred[idx] & is_dmc[idx]).sum())
        fp = int((pred[idx] & ~is_dmc[idx]).sum())
        fn = int((~pred[idx] & is_dmc[idx]).sum())
        denom = 2 * tp + fp + fn
        boot[b] = (2 * tp / denom) if denom > 0 else 0.0
    alpha = (1.0 - confidence) / 2.0
    lo = float(np.nanpercentile(boot, 100 * alpha))
    hi = float(np.nanpercentile(boot, 100 * (1.0 - alpha)))
    return (lo, hi)


if __name__ == "__main__":
    main()
