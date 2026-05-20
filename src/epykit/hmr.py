"""HMR / LMR caller -- single-sample hypo- and low-methylated regions.

MethylSeekR-style. Two-state HMM on per-CpG beta (not smoothed --
HMR/LMR are short-range, ~hundreds bp to few kb, and smoothing washes
them out).

  - HMR (hypo-methylated region): contiguous run of low-beta CpGs at
    average ``beta < hmr_threshold``. Includes CpG islands, promoters,
    enhancers.
  - LMR (low-methylated region): contiguous low-beta run NOT in a CpG-
    island context, i.e. mostly distal regulatory elements. LMRs are
    a subset of "hypo-state" runs filtered on CpG density:
    LMR = hypo-state run with low CpG density (< lmr_max_density).

Results land in ``md.uns["hmr"]`` and ``md.uns["lmr"]`` -- two parallel
frames with ``(sample_id, chrom, start, end, n_cpgs, mean_beta, length_bp)``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import polars as pl

from ._compute import run_chrom_pipeline
from ._hmm import runs_of_state, segment
from .dmc import _detect_chromosomes

logger = logging.getLogger(__name__)


_HMR_SCHEMA = {
    "sample_id": pl.Utf8,
    "chrom":     pl.Utf8,
    "start":     pl.Int32,
    "end":       pl.Int32,
    "length_bp": pl.Int32,
    "n_cpgs":    pl.Int32,
    "mean_beta": pl.Float32,
    "kind":      pl.Utf8,    # "HMR" or "LMR"
}


def _segment_chrom(
    store: Path, sample: str, chrom: str,
    *, hmr_threshold: float, self_loop: float,
) -> Optional[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Load one chrom's beta and return (positions, beta, viterbi)."""
    part = store / f"sample={sample}" / f"chrom={chrom}" / "part-0.parquet"
    if not part.exists():
        return None
    df = pl.read_parquet(str(part), columns=["pos", "N_meth", "coverage"]).sort("pos")
    if df.height < 5:
        return None
    positions = df["pos"].to_numpy().astype(np.int32)
    cov = df["coverage"].to_numpy().astype(np.float64)
    n_meth = df["N_meth"].to_numpy().astype(np.float64)
    beta = np.where(cov > 0, n_meth / np.maximum(cov, 1.0), np.nan)

    # 2-state HMM: state 0 = hypo (mean ~= hmr_threshold * 0.5),
    # state 1 = hyper (mean ~= average of hmr_threshold and 1).
    state_means = np.array([
        max(hmr_threshold * 0.5, 0.05),
        min((hmr_threshold + 1.0) / 2.0, 0.95),
    ])
    viterbi = segment(beta, n_states=2, state_means=state_means, self_loop=self_loop)
    return positions, beta, viterbi


def call_hmr_one_sample(
    store: Path,
    sample: str,
    *,
    chromosomes: Optional[list[str]] = None,
    hmr_threshold: float = 0.30,
    lmr_max_density: float = 0.020,   # CpGs per bp; below this -> LMR
    min_cpgs: int = 4,
    self_loop: float = 0.85,
    backend: str = "sequential",
    n_workers: Optional[int] = None,
) -> pl.DataFrame:
    """Call HMRs (and tag LMR subset) for one sample.

    The default ``self_loop=0.85`` is deliberately looser than the PMD
    caller's 0.999 -- HMR/LMR are short-range (hundreds bp to a few kb)
    so a too-sticky chain misses them between filler high-beta CpGs.
    """
    store = Path(store)
    if chromosomes is None:
        chromosomes = _detect_chromosomes(store)

    def _hmr_chrom_handler(chrom: str) -> Optional[pl.DataFrame]:
        pkg = _segment_chrom(
            store, sample, chrom,
            hmr_threshold=hmr_threshold, self_loop=self_loop,
        )
        if pkg is None:
            return None
        positions, beta, viterbi = pkg
        runs = runs_of_state(viterbi, target_state=0, positions=positions)
        rows: list[dict[str, object]] = []
        for run_start, run_end, run_len_sites in runs:
            if run_len_sites < min_cpgs:
                continue
            length_bp = run_end - run_start
            mask = (positions >= run_start) & (positions < run_end)
            if not mask.any():
                continue
            sel_beta = beta[mask]
            valid = np.isfinite(sel_beta)
            if valid.sum() == 0:
                continue
            mean_beta = float(sel_beta[valid].mean())
            # HMR if mean beta < hmr_threshold; skip otherwise.
            if mean_beta >= hmr_threshold:
                continue
            # LMR if CpG density (n_cpgs / length_bp) is low.
            density = run_len_sites / max(length_bp, 1)
            kind = "LMR" if density < lmr_max_density else "HMR"
            rows.append({
                "chrom": chrom,
                "start": int(run_start),
                "end": int(run_end),
                "length_bp": int(length_bp),
                "n_cpgs": int(run_len_sites),
                "mean_beta": float(mean_beta),
                "kind": kind,
            })
        if not rows:
            return None
        return pl.DataFrame(
            rows,
            schema={k: v for k, v in _HMR_SCHEMA.items() if k != "sample_id"},
        )

    parts: list[pl.DataFrame] = []
    for chrom, chrom_result in run_chrom_pipeline(
        chromosomes, _hmr_chrom_handler,
        backend=backend, n_workers=n_workers, label=f"HMR[{sample}]",
    ):
        parts.append(chrom_result)
    if not parts:
        return pl.DataFrame(schema={k: v for k, v in _HMR_SCHEMA.items() if k != "sample_id"})
    return pl.concat(parts, how="vertical_relaxed")


def hmr(
    md,
    *,
    samples: Optional[list[str]] = None,
    hmr_threshold: float = 0.30,
    lmr_max_density: float = 0.020,
    min_cpgs: int = 4,
    chromosomes: Optional[list[str]] = None,
    backend: str = "sequential",
    n_workers: Optional[int] = None,
) -> None:
    """Call HMR + LMR across samples; store ``md.uns["hmr"]`` and ``md.uns["lmr"]``."""
    md_samples = md.obs.get_column("sample_id").to_list()
    if samples is None:
        samples = md_samples
    bad = [s for s in samples if s not in md_samples]
    if bad:
        raise ValueError(f"unknown samples: {bad[:5]}")

    parts: list[pl.DataFrame] = []
    for sample in samples:
        df = call_hmr_one_sample(
            md.store, sample,
            chromosomes=chromosomes,
            hmr_threshold=hmr_threshold,
            lmr_max_density=lmr_max_density,
            min_cpgs=min_cpgs,
            backend=backend, n_workers=n_workers,
        )
        if df.height:
            df = df.with_columns(pl.lit(sample).alias("sample_id"))
            parts.append(df)

    if parts:
        combined = pl.concat(parts, how="vertical_relaxed")
    else:
        combined = pl.DataFrame(schema=_HMR_SCHEMA)
    combined = combined.sort(["sample_id", "chrom", "start"])

    md.uns["hmr"] = combined.filter(pl.col("kind") == "HMR")
    md.uns["lmr"] = combined.filter(pl.col("kind") == "LMR")
    md.uns["hmr_params"] = {
        "n_samples": len(samples),
        "n_hmr": int(md.uns["hmr"].height),
        "n_lmr": int(md.uns["lmr"].height),
        "hmr_threshold": hmr_threshold,
        "lmr_max_density": lmr_max_density,
        "min_cpgs": min_cpgs,
    }


__all__ = ["call_hmr_one_sample", "hmr"]
