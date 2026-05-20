from __future__ import annotations

import numpy as np

from .._style import PALETTE
from ._compute import (
    compute_volcano_data,
    compute_ma_data,
    compute_manhattan_data,
)
from ._utils import _get_ax, _save_fig
from ..methyldata import MethylData


def volcano(
    md: MethylData,
    *,
    alpha: float = 0.05,
    min_abs_diff: float = 0.1,
    ax=None,
    figsize=(6, 5),
    save: str | None = None,
):
    data = compute_volcano_data(md, alpha=alpha, min_abs_diff=min_abs_diff)
    diff = data.meth_diff
    y = data.neg_log_p
    ns = ~data.sig

    fig, ax = _get_ax(ax, figsize)
    ax.scatter(diff[ns], y[ns], s=4, color=PALETTE["neutral"], alpha=0.4, rasterized=True)
    ax.scatter(diff[data.hypo], y[data.hypo], s=4, color=PALETTE["hypo"], alpha=0.7, rasterized=True)
    ax.scatter(diff[data.hyper], y[data.hyper], s=4, color=PALETTE["hyper"], alpha=0.7, rasterized=True)

    ax.axhline(-np.log10(alpha), color="grey", lw=0.8, ls="--")
    ax.axvline(min_abs_diff, color="grey", lw=0.8, ls="--")
    ax.axvline(-min_abs_diff, color="grey", lw=0.8, ls="--")

    n_hyper = int(data.hyper.sum())
    n_hypo = int(data.hypo.sum())
    ax.set_title(f"DMC volcano  |  hyper={n_hyper:,}  hypo={n_hypo:,}")
    ax.set_xlabel("Methylation difference (treatment - control)")
    ax.set_ylabel(f"-log_1_0({data.p_col})")

    if save:
        _save_fig(md, fig, save)
    return fig, ax


def ma_plot(
    md: MethylData,
    *,
    alpha: float = 0.05,
    min_abs_diff: float = 0.1,
    ax=None,
    figsize=(7, 5),
    save: str | None = None,
):
    """MA plot: mean beta vs methylation difference."""
    data = compute_ma_data(md, alpha=alpha, min_abs_diff=min_abs_diff)
    ns = ~data.sig

    fig, ax = _get_ax(ax, figsize)
    ax.scatter(data.mean_beta[ns], data.meth_diff[ns], s=4, color=PALETTE["neutral"], alpha=0.4, rasterized=True)
    ax.scatter(data.mean_beta[data.hypo], data.meth_diff[data.hypo], s=4, color=PALETTE["hypo"], alpha=0.7, rasterized=True)
    ax.scatter(data.mean_beta[data.hyper], data.meth_diff[data.hyper], s=4, color=PALETTE["hyper"], alpha=0.7, rasterized=True)

    ax.axhline(0, color="black", lw=1)
    ax.axhline(min_abs_diff, color="grey", lw=0.8, ls="--", alpha=0.5)
    ax.axhline(-min_abs_diff, color="grey", lw=0.8, ls="--", alpha=0.5)

    n_hyper = int(data.hyper.sum())
    n_hypo = int(data.hypo.sum())
    ax.set_title(f"MA plot  |  hyper={n_hyper:,}  hypo={n_hypo:,}")
    ax.set_xlabel("Mean methylation")
    ax.set_ylabel("Methylation difference (treatment - control)")

    if save:
        _save_fig(md, fig, save)
    return fig, ax


def manhattan(
    md: MethylData,
    *,
    alpha: float = 0.05,
    ax=None,
    figsize=(14, 4),
    save: str | None = None,
):
    """Manhattan plot: genome-wide significance."""
    data = compute_manhattan_data(md, alpha=alpha)

    fig, ax = _get_ax(ax, figsize)
    colors = [PALETTE["hypo"], PALETTE["hyper"]]
    for i, block in enumerate(data.chrom_blocks):
        ax.scatter(block["x"], block["y"], s=3, color=colors[i % 2], alpha=0.6, rasterized=True)
    ax.axhline(data.alpha_line_y, color="red", lw=1, ls="--", label=f"alpha={alpha}")
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(f"-log_1_0({data.p_col})")
    ax.set_title("Manhattan plot")
    ax.legend()
    ax.set_xticks(data.tick_pos)
    ax.set_xticklabels(data.tick_label, fontsize=8)

    if save:
        _save_fig(md, fig, save)
    return fig, ax


__all__ = ["volcano", "ma_plot", "manhattan"]
