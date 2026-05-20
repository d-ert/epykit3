"""Partially methylated domains (PMD) caller -- single-sample, megabase-scale.

Two-state HMM over coverage-weighted, smoothed beta per chromosome:
  - State 0: PMD (low / partially methylated)
  - State 1: non-PMD (highly methylated euchromatin)

Per-sample, NOT per-group. Output lands in ``md.uns["pmd"]`` with one
row per called PMD: ``(sample_id, chrom, start, end, length_bp, mean_beta, n_cpgs)``.

Why HMM over hard thresholding?
  Naive threshold-on-smoothed-beta (beta < 0.7 -> PMD) gives jagged
  boundaries and breaks under coverage drop-outs. The 2-state HMM with
  a sticky transition prior smooths the assignment without
  oversmoothing the underlying signal.

This caller routes through the shared
:func:`epykit._compute.run_chrom_pipeline`, so it benefits from the
0.4 distributed backend automatically -- call with ``backend="dask"`` to
parallelise across chromosomes.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import polars as pl

from ._compute import run_chrom_pipeline
from ._hmm import runs_of_state, segment
from .dmc import _detect_chromosomes, _load_sample_chrom

logger = logging.getLogger(__name__)


_PMD_SCHEMA = {
    "sample_id": pl.Utf8,
    "chrom":     pl.Utf8,
    "start":     pl.Int32,
    "end":       pl.Int32,
    "length_bp": pl.Int32,
    "mean_beta": pl.Float32,
    "n_cpgs":    pl.Int32,
}


def _gaussian_smooth_beta(
    positions: np.ndarray,
    beta: np.ndarray,
    coverage: np.ndarray,
    bandwidth_bp: float,
) -> np.ndarray:
    """Coverage-weighted Gaussian smoother.

    Vectorised with an O(n_cpgs * window) inner loop using positional
    binary search to avoid the full O(n^2) kernel. For PMD work
    (bandwidth ~10 kb, CpG density ~1 per 100 bp) the effective window
    is ~100 CpGs -- tractable.
    """
    n = len(positions)
    out = np.full(n, np.nan, dtype=np.float64)
    if n == 0:
        return out
    # cap the kernel at 4 sigma
    cap = 4.0 * bandwidth_bp
    for i in range(n):
        p_i = positions[i]
        # Find sites within cap of i
        lo = np.searchsorted(positions, p_i - cap, side="left")
        hi = np.searchsorted(positions, p_i + cap, side="right")
        if hi <= lo:
            continue
        d = positions[lo:hi] - p_i
        w = np.exp(-0.5 * (d / bandwidth_bp) ** 2) * coverage[lo:hi]
        w_sum = w.sum()
        if w_sum > 0:
            out[i] = float((w * beta[lo:hi]).sum() / w_sum)
    return out


def call_pmd_one_sample(
    store: Path,
    sample: str,
    *,
    chromosomes: Optional[list[str]] = None,
    bandwidth_bp: float = 10_000,
    beta_threshold: float = 0.55,
    min_pmd_bp: int = 100_000,
    self_loop: float = 0.999,
    backend: str = "sequential",
    n_workers: Optional[int] = None,
) -> pl.DataFrame:
    """Call PMDs for a single sample across all detected chromosomes."""
    store = Path(store)
    if chromosomes is None:
        chromosomes = _detect_chromosomes(store)

    state_means = np.array([beta_threshold * 0.5, min(beta_threshold + 0.25, 0.95)])

    def _pmd_chrom_handler(chrom: str) -> Optional[pl.DataFrame]:
        # Build a canonical sorted (pos, strand) frame from the sample's
        # partition. PMDs are per-sample so we don't need an intersect.
        part = store / f"sample={sample}" / f"chrom={chrom}" / "part-0.parquet"
        if not part.exists():
            return None
        df = pl.read_parquet(str(part), columns=["pos", "N_meth", "coverage"]).sort("pos")
        if df.height < 10:
            return None
        positions = df["pos"].to_numpy().astype(np.int32)
        cov = df["coverage"].to_numpy().astype(np.float64)
        n_meth = df["N_meth"].to_numpy().astype(np.float64)
        beta = np.where(cov > 0, n_meth / np.maximum(cov, 1.0), np.nan)

        beta_smooth = _gaussian_smooth_beta(
            positions.astype(np.float64), beta, cov, bandwidth_bp=bandwidth_bp,
        )
        viterbi = segment(
            beta_smooth, n_states=2, state_means=state_means, self_loop=self_loop,
        )
        # State 0 is the PMD (low beta) state.
        runs = runs_of_state(viterbi, target_state=0, positions=positions)
        rows: list[dict[str, object]] = []
        for run_start, run_end, run_len_sites in runs:
            length_bp = run_end - run_start
            if length_bp < min_pmd_bp:
                continue
            # Recover mean_beta over the run from the actual sites.
            mask = (positions >= run_start) & (positions < run_end)
            if not mask.any():
                continue
            sel = beta[mask]
            sel_cov = cov[mask]
            valid = np.isfinite(sel)
            if valid.sum() == 0:
                continue
            mean_beta = float(np.average(sel[valid], weights=np.maximum(sel_cov[valid], 1.0)))
            rows.append({
                "chrom": chrom,
                "start": int(run_start),
                "end": int(run_end),
                "length_bp": int(length_bp),
                "mean_beta": float(mean_beta),
                "n_cpgs": int(run_len_sites),
            })
        if not rows:
            return None
        return pl.DataFrame(
            rows,
            schema={k: v for k, v in _PMD_SCHEMA.items() if k != "sample_id"},
        )

    parts: list[pl.DataFrame] = []
    for chrom, chrom_result in run_chrom_pipeline(
        chromosomes, _pmd_chrom_handler,
        backend=backend, n_workers=n_workers, label=f"PMD[{sample}]",
    ):
        parts.append(chrom_result)
    if not parts:
        return pl.DataFrame(schema={k: v for k, v in _PMD_SCHEMA.items() if k != "sample_id"})
    return pl.concat(parts, how="vertical_relaxed")


def pmd(
    md,
    *,
    samples: Optional[list[str]] = None,
    bandwidth_bp: float = 10_000,
    beta_threshold: float = 0.55,
    min_pmd_bp: int = 100_000,
    chromosomes: Optional[list[str]] = None,
    backend: str = "sequential",
    n_workers: Optional[int] = None,
) -> None:
    """Call PMDs across all (or specified) samples; store in ``md.uns["pmd"]``."""
    md_samples = md.obs.get_column("sample_id").to_list()
    if samples is None:
        samples = md_samples
    bad = [s for s in samples if s not in md_samples]
    if bad:
        raise ValueError(f"unknown samples: {bad[:5]}")

    parts: list[pl.DataFrame] = []
    for sample in samples:
        df = call_pmd_one_sample(
            md.store, sample,
            chromosomes=chromosomes,
            bandwidth_bp=bandwidth_bp,
            beta_threshold=beta_threshold,
            min_pmd_bp=min_pmd_bp,
            backend=backend, n_workers=n_workers,
        )
        if df.height > 0:
            df = df.with_columns(pl.lit(sample).alias("sample_id"))
            parts.append(df)

    if parts:
        combined = pl.concat(parts, how="vertical_relaxed").sort(["sample_id", "chrom", "start"])
    else:
        combined = pl.DataFrame(schema=_PMD_SCHEMA)

    md.uns["pmd"] = combined
    md.uns["pmd_params"] = {
        "n_samples": len(samples),
        "n_pmds": int(combined.height),
        "bandwidth_bp": bandwidth_bp,
        "beta_threshold": beta_threshold,
        "min_pmd_bp": min_pmd_bp,
    }


__all__ = ["call_pmd_one_sample", "pmd"]
