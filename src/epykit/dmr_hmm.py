"""HMM-based DMR caller.

Three-state HMM on the per-CpG ``meth_diff`` signal:
  - State 0: hypo (meth_diff < 0)
  - State 1: no change
  - State 2: hyper (meth_diff > 0)

Models spatial correlation between adjacent CpGs explicitly -- fixes a
known weakness of the tile / sliding-window engines (they treat
adjacent CpGs as independent). Operates on the existing per-CpG DMC
table (``md.varm["dmc_lr"]`` etc.); does NOT need a re-run of the DMC
pass.

Schema parity with ``call_dmr_tile_based``: the result frame uses the
same column names so ``pl.dmr_boxplot``, ``dmrs_to_bed``, and the
report harness work unchanged.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import polars as pl

from ._hmm import runs_of_state, segment
from .dmr import _DMR_TILE_SCHEMA

logger = logging.getLogger(__name__)


def _state_means_for_meth_diff(meth_diff: np.ndarray) -> np.ndarray:
    """Pick 3-state Gaussian means for the meth_diff signal.

    The states are {hypo, neutral, hyper}. Means are placed at a fixed
    offset rather than estimated from the data so a chromosome with no
    real signal still has well-defined targets (avoids the degenerate
    'fit three states to noise' behaviour).
    """
    return np.array([-0.20, 0.00, 0.20])


def call_dmr_hmm(
    dmc_results: pl.DataFrame,
    *,
    self_loop: float = 0.95,
    min_cpgs: int = 5,
    min_abs_meth_diff: float = 0.10,
    alpha: float = 0.05,
) -> pl.DataFrame:
    """Run HMM segmentation on a DMC table and emit a DMR frame.

    Parameters
    ----------
    dmc_results
        Any DMC table with at least ``chrom``, ``pos``, ``meth_diff``,
        and (optionally) ``pvalue`` / ``qvalue`` columns. Schema-
        compatible with the output of ``ep.tl.dmc``.
    self_loop
        Sticky-chain transition prior. Higher -> broader regions.
    min_cpgs, min_abs_meth_diff, alpha
        Per-region filters applied AFTER segmentation:
        - region n_cpgs >= min_cpgs
        - |mean meth_diff| >= min_abs_meth_diff
        - if a ``qvalue`` column exists, mean(qvalue) < alpha is
          required for the region to be called significant.

    Returns
    -------
    pl.DataFrame
        DMR frame with the same schema as ``call_dmr_tile_based``.
    """
    if dmc_results.height == 0:
        return pl.DataFrame(schema=_DMR_TILE_SCHEMA)
    required = {"chrom", "pos", "meth_diff"}
    missing = required - set(dmc_results.columns)
    if missing:
        raise ValueError(f"dmc_results missing required columns: {sorted(missing)}")

    state_means = _state_means_for_meth_diff(
        dmc_results["meth_diff"].to_numpy()
    )

    out_rows: list[dict[str, object]] = []
    for chrom_grp in dmc_results.partition_by("chrom", maintain_order=True):
        chrom = chrom_grp["chrom"][0]
        chrom_sorted = chrom_grp.sort("pos")
        positions = chrom_sorted["pos"].to_numpy().astype(np.int32)
        meth_diff = chrom_sorted["meth_diff"].to_numpy().astype(np.float64)

        viterbi = segment(
            meth_diff, n_states=3,
            state_means=state_means,
            self_loop=self_loop,
            emission="gaussian",
            emission_sd=0.10,
        )

        for state_idx, label in ((0, "hypo"), (2, "hyper")):
            runs = runs_of_state(viterbi, target_state=state_idx, positions=positions)
            for run_start, run_end, n_cpgs_run in runs:
                if n_cpgs_run < min_cpgs:
                    continue
                mask = (positions >= run_start) & (positions < run_end)
                if not mask.any():
                    continue
                run_md = meth_diff[mask]
                valid = np.isfinite(run_md)
                if valid.sum() == 0:
                    continue
                mean_md = float(run_md[valid].mean())
                if abs(mean_md) < min_abs_meth_diff:
                    continue
                # If the DMC table carries q-values, require region mean < alpha.
                if "qvalue" in chrom_sorted.columns:
                    run_q = chrom_sorted.filter(
                        (pl.col("pos") >= run_start) & (pl.col("pos") < run_end)
                    )["qvalue"].to_numpy()
                    run_q_valid = run_q[np.isfinite(run_q)]
                    if run_q_valid.size and float(run_q_valid.mean()) >= alpha:
                        continue
                # Build the row matching the tile DMR schema.
                out_rows.append({
                    "chrom":            chrom,
                    "start":            int(run_start),
                    "end":              int(run_end),
                    "n_cpgs":           int(n_cpgs_run),
                    "n_case":           0,    # unknown from DMC table alone
                    "n_control":        0,
                    "mean_beta_case":   float("nan"),
                    "mean_beta_control": float("nan"),
                    "meth_diff":        float(mean_md),
                    "log2_odds_ratio":  float("nan"),
                    "pvalue":           float("nan"),
                    "qvalue":           float("nan"),
                    "dmr_type":         label,
                })

    if not out_rows:
        return pl.DataFrame(schema=_DMR_TILE_SCHEMA)

    # Build the frame respecting the tile schema's types.
    return pl.DataFrame(
        out_rows,
        schema={
            "chrom":            pl.Utf8,
            "start":            pl.Int32,
            "end":              pl.Int32,
            "n_cpgs":           pl.Int32,
            "n_case":           pl.Int32,
            "n_control":        pl.Int32,
            "mean_beta_case":   pl.Float32,
            "mean_beta_control": pl.Float32,
            "meth_diff":        pl.Float32,
            "log2_odds_ratio":  pl.Float64,
            "pvalue":           pl.Float64,
            "qvalue":           pl.Float64,
            "dmr_type":         pl.Utf8,
        },
    ).sort(["chrom", "start"])


__all__ = ["call_dmr_hmm"]
