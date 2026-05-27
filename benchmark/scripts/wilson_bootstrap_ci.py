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


if __name__ == "__main__":
    main()
