"""evaluate.py — CI finalisation for eval_summary.parquet.

Phase 3 scope: --ci-only mode that adds Wilson CIs on TPR/FPR and
bootstrap CIs on AUROC/F1 (NaN when per-CpG cache unavailable) to an
existing eval_summary.parquet in place.

Usage:
    python evaluate.py --ci-only --eval-summary <path>
    python evaluate.py --ci-only --eval-summary <path> --per-cpg-dir <dir>
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl

# Make sibling scripts importable.
_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from wilson_bootstrap_ci import (  # noqa: E402
    add_wilson_ci_for_tpr_fpr,
    bootstrap_auroc_ci,
    bootstrap_f1_ci,
)


def _add_bootstrap_ci(df: pl.DataFrame, per_cpg_dir: Path | None) -> pl.DataFrame:
    """Add bootstrap CIs for AUROC and F1.

    Requires per-CpG joined parquets at per_cpg_dir/<tool>_<scenario>.parquet.
    Falls back to NaN when the cache is unavailable (Phase 3 seed; Phase 4
    populates the cache during the main eval pass).
    """
    n = df.height
    auroc_lo = np.full(n, np.nan)
    auroc_hi = np.full(n, np.nan)
    f1_lo = np.full(n, np.nan)
    f1_hi = np.full(n, np.nan)

    if per_cpg_dir is not None and per_cpg_dir.exists():
        for i, row in enumerate(df.iter_rows(named=True)):
            cache = per_cpg_dir / f"{row['tool']}_{row['scenario']}.parquet"
            if not cache.exists():
                continue
            j = pl.read_parquet(str(cache))
            if "is_dmc" not in j.columns:
                continue
            is_dmc = j["is_dmc"].to_numpy().astype(bool)
            pvals = j["pvalue"].to_numpy() if "pvalue" in j.columns else None
            qvals = j["qvalue"].to_numpy() if "qvalue" in j.columns else pvals
            if pvals is None:
                continue
            seed = abs(
                hash((str(row.get("tool", "")), str(row.get("scenario", "")), 0.05))
            ) % (2**32)
            try:
                auroc_lo[i], auroc_hi[i] = bootstrap_auroc_ci(
                    is_dmc=is_dmc, pvalues=pvals, B=1000, seed=seed,
                )
                if qvals is not None:
                    f1_lo[i], f1_hi[i] = bootstrap_f1_ci(
                        is_dmc=is_dmc, qvalues=qvals, threshold=0.05, B=1000, seed=seed,
                    )
            except Exception:
                pass

    return df.with_columns([
        pl.Series("auroc_ci_lo", auroc_lo),
        pl.Series("auroc_ci_hi", auroc_hi),
        pl.Series("f1_ci_lo", f1_lo),
        pl.Series("f1_ci_hi", f1_hi),
    ])


def add_ci_columns(
    df: pl.DataFrame,
    per_cpg_dir: Path | None = None,
) -> pl.DataFrame:
    """Add all CI columns to an eval_summary DataFrame."""
    df = add_wilson_ci_for_tpr_fpr(df)
    df = _add_bootstrap_ci(df, per_cpg_dir)
    return df


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ci-only", action="store_true",
        help="Append CI columns to an existing eval_summary.parquet in place.",
    )
    parser.add_argument(
        "--eval-summary", required=True, type=Path,
        help="Path to eval_summary.parquet.",
    )
    parser.add_argument(
        "--per-cpg-dir", type=Path, default=None,
        help="Optional dir of per-CpG joined parquets for bootstrap CIs.",
    )
    args = parser.parse_args(argv)

    if not args.ci_only:
        parser.error("Only --ci-only mode is implemented in 0.7.5.")

    df = pl.read_parquet(str(args.eval_summary))
    out = add_ci_columns(df, per_cpg_dir=args.per_cpg_dir)
    out.write_parquet(str(args.eval_summary))
    print(f"updated {args.eval_summary} with CI columns ({out.height} rows)")


if __name__ == "__main__":
    main()
