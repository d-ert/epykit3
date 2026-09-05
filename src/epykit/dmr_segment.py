"""Rule-based segmentation DMR caller (renamed from dmr_hmm in 0.7.5).

Three-state decoder on the per-CpG ``meth_diff`` signal with **fixed**
state means and emission SDs -- not a fitted HMM. The name
``dmr_segment`` reflects this honestly.

P2-4 fix: per-segment p-values are Stouffer-combined from constituent
CpG p-values and BH-corrected per chromosome. The pre-0.7.5
implementation emitted NaN p/q-values for every segment.
"""

from __future__ import annotations

import logging

import numpy as np
import polars as pl

from ._hmm import runs_of_state, segment
from .dmr import _DMR_TILE_SCHEMA, _stouffer_combine_signed

logger = logging.getLogger(__name__)


def _state_means_for_meth_diff(meth_diff: np.ndarray) -> np.ndarray:
    """Fixed 3-state targets {hypo, neutral, hyper} at -0.20 / 0.0 / +0.20."""
    return np.array([-0.20, 0.00, 0.20])


def _bh_per_chrom(pvals: np.ndarray, chroms: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjustment applied per chromosome."""
    out = np.full_like(pvals, np.nan, dtype=np.float64)
    for chrom in np.unique(chroms):
        mask = chroms == chrom
        p = pvals[mask]
        finite = np.isfinite(p)
        if not finite.any():
            continue
        p_finite = p[finite]
        order = np.argsort(p_finite)
        ranks = np.empty_like(order)
        ranks[order] = np.arange(1, p_finite.size + 1)
        q = p_finite * p_finite.size / ranks
        q_sorted = q[order]
        q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
        q_back = np.empty_like(q_sorted)
        q_back[order] = q_sorted
        q_full = np.full_like(p, np.nan)
        q_full[finite] = np.clip(q_back, 0.0, 1.0)
        out[mask] = q_full
    return out


def call_dmr_rule_segment(
    dmc_results: pl.DataFrame,
    *,
    self_loop: float = 0.95,
    min_cpgs: int = 5,
    min_abs_meth_diff: float = 0.10,
    alpha: float = 0.05,
) -> pl.DataFrame:
    """Rule-based segmentation DMR caller.

    Parameters
    ----------
    dmc_results
        DMC table with at least ``chrom``, ``pos``, ``meth_diff`` columns.
        If a ``pvalue`` column is present, per-segment p-values are
        Stouffer-combined from constituent CpG p-values and BH-corrected
        per chromosome. Without ``pvalue``, p/q-values are NaN.
    self_loop, min_cpgs, min_abs_meth_diff, alpha
        Same semantics as the old ``call_dmr_hmm``.

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
    has_pvalue = "pvalue" in dmc_results.columns

    state_means = _state_means_for_meth_diff(dmc_results["meth_diff"].to_numpy())

    out_rows: list[dict[str, object]] = []
    for chrom_grp in dmc_results.partition_by("chrom", maintain_order=True):
        chrom = chrom_grp["chrom"][0]
        chrom_sorted = chrom_grp.sort("pos")
        positions = chrom_sorted["pos"].to_numpy().astype(np.int64)
        meth_diff = chrom_sorted["meth_diff"].to_numpy().astype(np.float64)
        pvals_per_cpg = chrom_sorted["pvalue"].to_numpy().astype(np.float64) if has_pvalue else None

        viterbi = segment(
            meth_diff,
            n_states=3,
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
                # Signed Stouffer combine (D1): the unsigned two-sided
                # variant added |z| regardless of direction, so a region's
                # combined p shrank toward 0 as it grew even when per-CpG
                # effects were mixed -- anti-conservative. The signed combine
                # (shared with the tile/sliding-window callers) cancels
                # opposing directions so only directionally-coherent regions
                # get small p-values.
                seg_p = (
                    _stouffer_combine_signed(pvals_per_cpg[mask], run_md)
                    if pvals_per_cpg is not None
                    else float("nan")
                )
                out_rows.append(
                    {
                        "chrom": str(chrom),
                        "start": int(run_start),
                        "end": int(run_end),
                        "n_cpgs": int(n_cpgs_run),
                        "n_case": 0,
                        "n_control": 0,
                        "mean_beta_case": float("nan"),
                        "mean_beta_control": float("nan"),
                        "meth_diff": float(mean_md),
                        "log2_odds_ratio": float("nan"),
                        "pvalue": float(seg_p),
                        "qvalue": float("nan"),  # filled by BH below
                        "dmr_type": label,
                    }
                )

    if not out_rows:
        return pl.DataFrame(schema=_DMR_TILE_SCHEMA)

    df = pl.DataFrame(
        out_rows,
        schema={
            "chrom": pl.Utf8,
            "start": pl.Int32,
            "end": pl.Int32,
            "n_cpgs": pl.Int32,
            "n_case": pl.Int32,
            "n_control": pl.Int32,
            "mean_beta_case": pl.Float32,
            "mean_beta_control": pl.Float32,
            "meth_diff": pl.Float32,
            "log2_odds_ratio": pl.Float64,
            "pvalue": pl.Float64,
            "qvalue": pl.Float64,
            "dmr_type": pl.Utf8,
        },
    ).sort(["chrom", "start"])

    # BH per chromosome.
    if has_pvalue:
        qvals = _bh_per_chrom(df["pvalue"].to_numpy(), df["chrom"].to_numpy())
        df = df.with_columns(pl.Series("qvalue", qvals))
        # Significance gate on per-segment q (replaces old mean(qvalue)<alpha hack).
        df = df.filter(pl.col("qvalue") < alpha)

    return df


__all__ = ["call_dmr_rule_segment"]
