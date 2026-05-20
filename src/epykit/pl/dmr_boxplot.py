"""Per-sample beta strip-plot panels for the top DMRs.

For each of the top ``n`` DMRs (ranked by qvalue, pvalue, or any column
on the DMR table), this function queries the per-sample beta within the
DMR's coordinates and renders a strip/box plot grouped by
``md.obs[group_by]``.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from .._style import PALETTE
from ._utils import _save_fig


_MAX_PANELS = 25  # hard cap; see plan Sec.6 implementation notes


def dmr_boxplot(
    md,
    *,
    top_n: int = 10,
    by: str = "qvalue",
    group_by: str = "group",
    figsize=None,
    save: str | None = None,
):
    """Strip-plot of per-sample beta across the top ``n`` DMRs.

    Parameters
    ----------
    top_n : int
        Number of DMRs to draw, ranked by ``by`` ascending. Capped at 25.
    by : str
        DMR-table column used for ranking. Defaults to ``"qvalue"`` (the
        tile-method output); falls back to ``"pvalue"`` if absent.
    group_by : str
        ``md.obs`` column used to colour samples.
    figsize, save : matplotlib plumbing.
    """
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    dmr_df = md.uns.get("dmr")
    if not isinstance(dmr_df, pl.DataFrame) or len(dmr_df) == 0:
        raise ValueError(
            "No DMR table on md.uns['dmr']. Run ep.tl.dmr(md) first."
        )
    rank_col = by if by in dmr_df.columns else "pvalue"
    if rank_col not in dmr_df.columns:
        rank_col = dmr_df.columns[0]
    n = min(int(top_n), _MAX_PANELS, len(dmr_df))
    if n < int(top_n):
        import warnings
        warnings.warn(
            f"top_n capped at {n} (asked for {top_n}); maximum is "
            f"{_MAX_PANELS} or len(dmr) = {len(dmr_df)}.",
            UserWarning, stacklevel=2,
        )
    top = dmr_df.sort(rank_col, descending=False).head(n)

    cols = min(n, 5)
    rows = (n + cols - 1) // cols
    if figsize is None:
        figsize = (3.4 * cols, 3.0 * rows)
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(rows, cols, figure=fig, hspace=0.55, wspace=0.35)
    axes: list = []

    obs = md.obs
    groups_per_sample = (
        obs.get_column(group_by).to_list()
        if group_by in obs.columns else
        ["sample"] * len(obs)
    )
    samples_obs = obs.get_column("sample_id").to_list()
    rng = np.random.default_rng(0)

    for i, row in enumerate(top.iter_rows(named=True)):
        chrom = row["chrom"]
        start = int(row.get("start", row.get("pos", 0)))
        end = int(row.get("end", start + 1))
        ax = fig.add_subplot(gs[i // cols, i % cols])
        axes.append(ax)
        beta_df = md.region_beta(chrom, start, end)
        # Align beta_df ordering to obs sample order
        beta_lookup = {
            r["sample"]: r["mean_beta"]
            for r in beta_df.iter_rows(named=True)
        }
        unique_groups = sorted(set(groups_per_sample))
        x_positions = {g: pos for pos, g in enumerate(unique_groups)}
        for sample, g in zip(samples_obs, groups_per_sample):
            y = beta_lookup.get(sample, float("nan"))
            if not np.isfinite(y):
                continue
            x = x_positions[g] + rng.uniform(-0.15, 0.15)
            color = (
                PALETTE.get("treatment") if g == "treatment"
                else PALETTE.get("control") if g == "control"
                else PALETTE.get("neutral", "#888")
            )
            ax.scatter(x, y, color=color, s=30, alpha=0.85)
        ax.set_xticks(list(x_positions.values()))
        ax.set_xticklabels(unique_groups, fontsize=8)
        ax.set_ylim(0.0, 1.0)
        ax.set_ylabel("beta" if i % cols == 0 else "")
        q = row.get(rank_col)
        q_str = f"{q:.2g}" if isinstance(q, (int, float)) and np.isfinite(q) else "--"
        ax.set_title(
            f"{chrom}:{start:,}-{end:,}\n{rank_col}={q_str}",
            fontsize=8,
        )

    if save:
        _save_fig(md, fig, save)
    return fig, axes


__all__ = ["dmr_boxplot"]
