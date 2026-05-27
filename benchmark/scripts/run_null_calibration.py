"""Null-calibration runner: empirical FDR under label-shuffled data.

For each (engine, scenario), randomly shuffles treatment/control labels
over the same samples, re-runs the engine, and records the observed
proportion of sites called significant at the nominal threshold. With
no true DMCs in the shuffled design, observed FDR at nominal q < 0.05
should be ~ 0.05 if the test is well-calibrated, OR much lower if the
test is conservative on the input data's noise regime.

See `docs/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md`
§2.2.

This module exposes:

- `compute_observed_fdr(qvalues, q_thresh)`: the pure-arithmetic kernel.
- `run_null_calibration(...)`: orchestrates k shuffles, calls the engine
  closure once per shuffle, returns a polars DataFrame with one row per
  shuffle and Wilson CI bounds on observed FDR.

The engine closure has signature
    engine_fn(samples_treatment, samples_control, seed=int) -> np.ndarray
returning per-site q-values for that shuffle. Real callers wrap epykit's
``ep.tl.dmc`` (or equivalent) in such a closure; tests use a fake.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
import polars as pl

from wilson_bootstrap_ci import _wilson_single


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


def run_null_calibration(
    engine_fn: Callable,
    engine_name: str,
    scenario_name: str,
    samples: list[str],
    n_per_group: int,
    k_shuffles: int = 20,
    q_thresh: float = 0.05,
    seed: int = 0,
) -> pl.DataFrame:
    """Run k label-shuffles and aggregate observed FDR per shuffle.

    Returns a DataFrame with columns:
      engine, scenario, k_shuffle, n_called, n_total, observed_fdr,
      observed_fdr_ci_lo, observed_fdr_ci_hi.

    Wilson CIs treat each shuffle's observed FDR as a binomial proportion
    of `n_called` / `n_total`. They quantify the within-shuffle estimation
    uncertainty; for across-shuffle variability, compute the median + IQR
    over rows externally.
    """
    rng = np.random.default_rng(seed)
    n = len(samples)
    if n_per_group * 2 != n:
        raise ValueError(
            f"need {n_per_group * 2} samples for {n_per_group}v{n_per_group}, got {n}"
        )

    rows = []
    for k in range(k_shuffles):
        # Local RNG per shuffle so a single global seed reproduces the run.
        shuffled = rng.permutation(samples).tolist()
        treat = shuffled[:n_per_group]
        ctrl = shuffled[n_per_group:]

        qvals = engine_fn(samples_treatment=treat, samples_control=ctrl, seed=seed + k + 1)
        stats = compute_observed_fdr(qvals, q_thresh=q_thresh)
        lo, hi = _wilson_single(stats["n_called"], stats["n_total"])
        rows.append({
            "engine": engine_name,
            "scenario": scenario_name,
            "k_shuffle": k,
            "n_called": stats["n_called"],
            "n_total": stats["n_total"],
            "observed_fdr": stats["observed_fdr"],
            "observed_fdr_ci_lo": lo,
            "observed_fdr_ci_hi": hi,
        })
    return pl.DataFrame(rows)


def main(argv: list[str] | None = None) -> None:
    """CLI stub. Real callers wire ``ep.tl.dmc`` here; this stub demonstrates
    the integration surface with a deterministic noise engine."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", default="noise", help="Engine label for reporting")
    parser.add_argument("--scenario", default="demo", help="Scenario label")
    parser.add_argument("--k-shuffles", type=int, default=20)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True, type=str, help="Output parquet path")
    args = parser.parse_args(argv)

    def noise_engine(samples_treatment, samples_control, seed=0):
        rng = np.random.default_rng(seed)
        return rng.uniform(0, 1, size=1000)

    df = run_null_calibration(
        engine_fn=noise_engine,
        engine_name=args.engine,
        scenario_name=args.scenario,
        samples=["s1", "s2", "s3", "s4", "s5", "s6"],
        n_per_group=3,
        k_shuffles=args.k_shuffles,
        seed=args.seed,
    )
    df.write_parquet(args.out)
    print(f"wrote {args.out} with {len(df)} shuffle rows")


if __name__ == "__main__":
    main()
