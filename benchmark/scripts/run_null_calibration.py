"""Null-calibration runner: empirical FDR + p-value Q-Q under label shuffles.

For each (engine, scenario), randomly shuffles treatment/control labels
over the same samples, re-runs the engine, and records:

  - Per-shuffle observed FDR at the nominal threshold (n_called / n_total).
  - A sampled subset of per-CpG p-values from a configurable number of
    representative shuffles, so a Q-Q plot vs Uniform(0, 1) and a
    Kolmogorov-Smirnov test can be computed across all sampled draws.

The pre-1.0 version used k=20 shuffles and within-shuffle Wilson CIs on
the binomial proportion. That is insufficient for a Nature/Genome
Biology submission: across-shuffle variance is not bounded by a 20-row
sample, and the calibration-vs-conservatism question (is the test
calibrated, or merely conservative?) cannot be answered without a
distributional check on the raw p-values. This module now defaults to
``k_shuffles = 1000`` and emits a Q-Q + KS bundle alongside the FDR
table.

Outputs:
  ``--out`` per-shuffle FDR parquet (columns:
      engine, scenario, k_shuffle, n_called, n_total, observed_fdr,
      observed_fdr_ci_lo, observed_fdr_ci_hi).
  ``--summary-out`` across-shuffle summary parquet with median + IQR
      of observed_fdr.
  ``--pvalue-out`` sampled p-value parquet for Q-Q (columns:
      engine, scenario, k_shuffle, pvalue).
  ``--ks-out`` KS test result JSON (statistic, p-value, n).

See ``docs/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md``
sec 2.2 for the upstream design intent. The Q-Q + KS additions land
under Phase 1.2 of the GB resubmission plan.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Callable, Optional

import numpy as np
import polars as pl
from scipy import stats as sp_stats

from wilson_bootstrap_ci import _wilson_single

logger = logging.getLogger(__name__)

# Per-shuffle p-value sample budget. Storing every CpG p-value from
# every shuffle would dominate output size; sampling a fixed number of
# p-values per sampled shuffle keeps the Q-Q + KS bundle bounded.
DEFAULT_QQ_SAMPLES_PER_SHUFFLE: int = 50_000
DEFAULT_QQ_SHUFFLES: int = 10


def compute_observed_fdr(qvalues: np.ndarray, q_thresh: float = 0.05) -> dict:
    """Per-shuffle stats: n_called, n_total, observed_fdr."""
    q = np.asarray(qvalues, dtype=np.float64)
    finite = np.isfinite(q)
    q_clean = q[finite]
    n_called = int((q_clean < q_thresh).sum())
    n_total = int(len(q_clean))
    if n_total == 0:
        return {"n_called": 0, "n_total": 0, "observed_fdr": float("nan")}
    return {
        "n_called": n_called,
        "n_total": n_total,
        "observed_fdr": float(n_called / n_total),
    }


def _sample_pvalues(
    pvalues: np.ndarray,
    n_samples: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Down-sample p-values without replacement; return all of them
    if fewer than n_samples are available."""
    p = np.asarray(pvalues, dtype=np.float64)
    finite = p[np.isfinite(p)]
    if finite.size <= n_samples:
        return finite
    idx = rng.choice(finite.size, size=n_samples, replace=False)
    return finite[idx]


def ks_against_uniform(pvalues: np.ndarray) -> dict:
    """Kolmogorov-Smirnov test of pvalues vs Uniform(0, 1).

    Returns dict with ``statistic``, ``pvalue``, ``n``. Under the null
    (calibrated test), the K-S p-value is approximately uniform on
    [0, 1]; under conservatism, p-values are stochastically larger than
    uniform and K-S detects this with high power at n >> 1e3.
    """
    p = np.asarray(pvalues, dtype=np.float64)
    finite = p[np.isfinite(p)]
    if finite.size == 0:
        return {"statistic": float("nan"), "pvalue": float("nan"), "n": 0}
    res = sp_stats.kstest(finite, "uniform", args=(0.0, 1.0))
    return {
        "statistic": float(res.statistic),
        "pvalue": float(res.pvalue),
        "n": int(finite.size),
    }


def summarize_observed_fdr(per_shuffle: pl.DataFrame) -> pl.DataFrame:
    """Compute median + IQR of observed_fdr across the shuffles."""
    if per_shuffle.is_empty():
        return pl.DataFrame()
    finite = per_shuffle.filter(pl.col("observed_fdr").is_finite())
    if finite.is_empty():
        return pl.DataFrame()
    vals = finite["observed_fdr"].to_numpy()
    return pl.DataFrame({
        "engine": [finite["engine"][0]],
        "scenario": [finite["scenario"][0]],
        "k_shuffles": [finite.height],
        "median_observed_fdr": [float(np.median(vals))],
        "q25_observed_fdr": [float(np.quantile(vals, 0.25))],
        "q75_observed_fdr": [float(np.quantile(vals, 0.75))],
        "min_observed_fdr": [float(np.min(vals))],
        "max_observed_fdr": [float(np.max(vals))],
    })


def run_null_calibration(
    engine_fn: Callable,
    engine_name: str,
    scenario_name: str,
    samples: list[str],
    n_per_group: int,
    k_shuffles: int = 1000,
    q_thresh: float = 0.05,
    seed: int = 0,
    qq_shuffles: int = DEFAULT_QQ_SHUFFLES,
    qq_samples_per_shuffle: int = DEFAULT_QQ_SAMPLES_PER_SHUFFLE,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Run k label-shuffles and aggregate observed FDR + sampled p-values.

    Parameters
    ----------
    engine_fn : Callable
        ``(samples_treatment, samples_control, seed) -> (pvalues, qvalues)``.
        Each array is a 1d numpy array of per-CpG values from one
        shuffle. The legacy contract returned only qvalues; closures in
        ``_null_engines.py`` are updated to return both.
    k_shuffles : int
        Total number of label shuffles to run. 1000 by default for
        usable across-shuffle variance estimates.
    qq_shuffles : int
        Number of shuffles whose p-values are sampled for the Q-Q + KS
        bundle. The first ``qq_shuffles`` shuffles are sampled; the
        rest record only the per-shuffle FDR summary.
    qq_samples_per_shuffle : int
        Cap on the number of p-values retained per Q-Q shuffle to keep
        output size bounded.

    Returns
    -------
    per_shuffle : pl.DataFrame
        One row per shuffle with FDR statistics.
    pvalue_samples : pl.DataFrame
        Long-form ``(engine, scenario, k_shuffle, pvalue)`` rows for
        the sampled p-values.
    """
    rng = np.random.default_rng(seed)
    n = len(samples)
    if n_per_group * 2 != n:
        raise ValueError(
            f"need {n_per_group * 2} samples for {n_per_group}v{n_per_group}, got {n}"
        )

    fdr_rows: list[dict] = []
    pvalue_rows: list[dict] = []

    for k in range(k_shuffles):
        shuffled = rng.permutation(samples).tolist()
        treat = shuffled[:n_per_group]
        ctrl = shuffled[n_per_group:]

        out = engine_fn(samples_treatment=treat, samples_control=ctrl, seed=seed + k + 1)
        pvals, qvals = _unpack_engine_output(out)
        stats_dict = compute_observed_fdr(qvals, q_thresh=q_thresh)
        lo, hi = _wilson_single(stats_dict["n_called"], stats_dict["n_total"])
        fdr_rows.append({
            "engine": engine_name,
            "scenario": scenario_name,
            "k_shuffle": k,
            "n_called": stats_dict["n_called"],
            "n_total": stats_dict["n_total"],
            "observed_fdr": stats_dict["observed_fdr"],
            "observed_fdr_ci_lo": lo,
            "observed_fdr_ci_hi": hi,
        })

        if k < qq_shuffles and pvals.size > 0:
            sampled = _sample_pvalues(pvals, qq_samples_per_shuffle, rng)
            for pv in sampled:
                pvalue_rows.append({
                    "engine": engine_name,
                    "scenario": scenario_name,
                    "k_shuffle": k,
                    "pvalue": float(pv),
                })

        if (k + 1) % 100 == 0:
            logger.info(
                "null calibration %s/%s: %d/%d shuffles done",
                engine_name, scenario_name, k + 1, k_shuffles,
            )

    return pl.DataFrame(fdr_rows), pl.DataFrame(pvalue_rows)


def _unpack_engine_output(out) -> tuple[np.ndarray, np.ndarray]:
    """Accept either the new ``(pvalues, qvalues)`` tuple contract or
    the legacy single-array ``qvalues`` contract from ``_null_engines.py``.

    The legacy contract degrades the Q-Q bundle: with no p-values
    available, the sampled p-value parquet is empty and the K-S test
    reports n=0. The summary FDR still computes correctly.
    """
    if isinstance(out, tuple) and len(out) == 2:
        pvals = np.asarray(out[0], dtype=np.float64)
        qvals = np.asarray(out[1], dtype=np.float64)
        return pvals, qvals
    arr = np.asarray(out, dtype=np.float64)
    # Legacy single-array path. Cannot compute Q-Q.
    return np.array([], dtype=np.float64), arr


def main(argv: Optional[list[str]] = None) -> None:
    """CLI: real-engine label-shuffle calibration.

    python run_null_calibration.py \\
        --engine lr \\
        --methylstore <path-to-store-root> \\
        --scenario cov10_3v3 \\
        --k-shuffles 1000 \\
        --qq-shuffles 10 \\
        --seed 0 \\
        --out      benchmark/data/null_calibration/cov10_3v3/lr.parquet \\
        --summary-out benchmark/data/null_calibration/cov10_3v3/lr_summary.parquet \\
        --pvalue-out  benchmark/data/null_calibration/cov10_3v3/lr_pvalues.parquet \\
        --ks-out      benchmark/data/null_calibration/cov10_3v3/lr_ks.json
    """
    import argparse
    import sys

    _scripts_dir = Path(__file__).resolve().parent
    if str(_scripts_dir) not in sys.path:
        sys.path.insert(0, str(_scripts_dir))

    from _null_engines import ENGINE_REGISTRY
    from epykit import MethylData

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", required=True,
                        choices=sorted(ENGINE_REGISTRY), help="Engine name")
    parser.add_argument("--methylstore", required=True, type=Path,
                        help="Path to the methylstore directory (contains .cache/)")
    parser.add_argument("--scenario", required=True, type=str)
    parser.add_argument("--k-shuffles", type=int, default=1000,
                        help="Total label shuffles (default: 1000)")
    parser.add_argument("--qq-shuffles", type=int, default=DEFAULT_QQ_SHUFFLES,
                        help="Number of shuffles to sample p-values from "
                             "for Q-Q + KS (default: 10)")
    parser.add_argument("--qq-samples-per-shuffle", type=int,
                        default=DEFAULT_QQ_SAMPLES_PER_SHUFFLE,
                        help="Per-shuffle p-value sample cap "
                             "(default: 50000)")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True, type=Path,
                        help="Per-shuffle FDR parquet")
    parser.add_argument("--summary-out", type=Path, default=None,
                        help="Across-shuffle median+IQR summary parquet")
    parser.add_argument("--pvalue-out", type=Path, default=None,
                        help="Sampled p-value parquet (for Q-Q plotting)")
    parser.add_argument("--ks-out", type=Path, default=None,
                        help="K-S test result JSON")
    parser.add_argument("--q-thresh", type=float, default=0.05)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    md = MethylData.load(str(args.methylstore))
    closure = ENGINE_REGISTRY[args.engine](md)

    samples = list(md.treatment_ids) + list(md.control_ids)
    n_per_group = len(md.treatment_ids)

    per_shuffle, pvalue_samples = run_null_calibration(
        engine_fn=closure,
        engine_name=args.engine,
        scenario_name=args.scenario,
        samples=samples,
        n_per_group=n_per_group,
        k_shuffles=args.k_shuffles,
        q_thresh=args.q_thresh,
        seed=args.seed,
        qq_shuffles=args.qq_shuffles,
        qq_samples_per_shuffle=args.qq_samples_per_shuffle,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    per_shuffle.write_parquet(str(args.out))
    logger.info("wrote per-shuffle FDR: %s (%d rows)", args.out, per_shuffle.height)

    if args.summary_out is not None:
        summary = summarize_observed_fdr(per_shuffle)
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        summary.write_parquet(str(args.summary_out))
        logger.info("wrote summary: %s", args.summary_out)

    if args.pvalue_out is not None and not pvalue_samples.is_empty():
        args.pvalue_out.parent.mkdir(parents=True, exist_ok=True)
        pvalue_samples.write_parquet(str(args.pvalue_out))
        logger.info("wrote sampled pvalues: %s (%d rows)",
                    args.pvalue_out, pvalue_samples.height)

    if args.ks_out is not None and not pvalue_samples.is_empty():
        ks = ks_against_uniform(pvalue_samples["pvalue"].to_numpy())
        args.ks_out.parent.mkdir(parents=True, exist_ok=True)
        args.ks_out.write_text(json.dumps(ks, indent=2) + "\n")
        logger.info("wrote KS test: %s (D=%.4f, p=%.4g, n=%d)",
                    args.ks_out, ks["statistic"], ks["pvalue"], ks["n"])

    # CLI tool: stdout result line is intentional (this is a script,
    # not library code; mirrors the convention from epykit.cli).
    print(f"wrote {args.out} ({per_shuffle.height} rows, engine={args.engine})")


if __name__ == "__main__":
    main()
