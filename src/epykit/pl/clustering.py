from __future__ import annotations

import numpy as np

from .._style import PALETTE
from ._compute import compute_pca
from ._utils import _get_ax, _save_fig
from ..methyldata import MethylData


def pca(
    md: MethylData,
    *,
    n_sites: int = 10000,
    ax=None,
    figsize=(6, 5),
    save: str | None = None,
):
    """PCA of per-sample methylation profiles.

    Thin renderer over :func:`epykit.pl._compute.compute_pca`. Single-pass
    over the methylstore with hash-based site subsampling; see the compute
    function's docstring for the memory model.

    Parameters
    ----------
    md : MethylData
        Methylation data object with filtered store.
    n_sites : int
        Number of sites to sample for PCA. Default 10000.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, create new figure.
    figsize : tuple
        Figure size (width, height) in inches.
    save : str, optional
        Path to save figure. If None, don't save.
    """
    res = compute_pca(md, n_sites=n_sites)

    fig, ax = _get_ax(ax, figsize)
    unique_groups = sorted(set(res.groups), key=lambda g: (g is None, str(g)))
    palette_cycle = [
        PALETTE.get("control"), PALETTE.get("treatment"),
        PALETTE.get("hyper"), PALETTE.get("hypo"), PALETTE.get("neutral"),
    ]
    for i, g in enumerate(unique_groups):
        mask = np.array([gg == g for gg in res.groups])
        if not mask.any():
            continue
        color = palette_cycle[i % len(palette_cycle)] or PALETTE.get("neutral")
        ax.scatter(
            res.coords[mask, 0], res.coords[mask, 1],
            s=100, alpha=0.6, label=str(g), color=color,
        )
    ax.set_xlabel(f"PC1 ({res.explained_var[0]:.1%})")
    ax.set_ylabel(f"PC2 ({res.explained_var[1]:.1%})")
    ax.set_title(
        f"PCA of methylation profiles  |  n_sites={res.n_sites_used:,}"
    )
    if res.group_col and res.group_col != "all":
        ax.legend(title=res.group_col)
    ax.grid(True, alpha=0.3)

    if save:
        _save_fig(md, fig, save)
    return fig, ax


__all__ = ["pca"]
