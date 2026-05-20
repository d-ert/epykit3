"""UMAP embedding of per-sample methylation profiles ."""

from __future__ import annotations

import numpy as np

from .._style import PALETTE
from ._utils import _get_ax, _save_fig, build_sample_site_matrix


def umap(
    md,
    *,
    n_sites: int = 10_000,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    n_components: int = 2,
    metric: str = "euclidean",
    color: str | None = None,
    seed: int = 42,
    ax=None,
    figsize=(6, 5),
    save: str | None = None,
):
    """UMAP scatter of per-sample methylation profiles.

    Lazy-imports ``umap-learn`` from the optional ``[viz]`` extra. Falls
    back to a clear ``ImportError`` when the dependency is missing.

    Parameters
    ----------
    md : MethylData
    n_sites : int
        Number of common CpGs to sample for the embedding basis.
    n_neighbors, min_dist, n_components, metric, seed
        UMAP hyperparameters (see ``umap.UMAP``). Defaults match the
        scanpy / single-cell convention.
    color : str, optional
        Name of an ``md.obs`` column to colour points by. Defaults to
        ``"group"`` when present, then ``"treatment"``.
    ax, figsize, save : matplotlib plumbing -- see :mod:`epykit.pl`.
    """
    try:
        import umap as umap_lib  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "umap-learn is required for pl.umap. "
            "Install with: pip install 'epykit[viz]'"
        ) from exc

    matrix, samples = build_sample_site_matrix(md, n_sites=n_sites)
    reducer = umap_lib.UMAP(
        n_neighbors=min(int(n_neighbors), max(matrix.shape[0] - 1, 2)),
        min_dist=min_dist,
        n_components=n_components,
        metric=metric,
        random_state=seed,
    )
    coords = reducer.fit_transform(matrix)

    if color is None:
        if "group" in md.obs.columns:
            color = "group"
        elif "treatment" in md.obs.columns:
            color = "treatment"
    groups = (
        md.obs.get_column(color).to_list()
        if color and color in md.obs.columns
        else ["sample"] * len(samples)
    )

    fig, ax = _get_ax(ax, figsize)
    unique_groups = sorted(set(groups), key=lambda g: (g is None, str(g)))
    palette = (
        [PALETTE.get("control"), PALETTE.get("treatment"),
         PALETTE.get("hyper"), PALETTE.get("hypo"), PALETTE.get("neutral")]
    )
    for i, g in enumerate(unique_groups):
        mask = np.array([row == g for row in groups])
        if mask.any():
            color_val = palette[i % len(palette)]
            ax.scatter(
                coords[mask, 0], coords[mask, 1], s=80, alpha=0.7,
                label=str(g), color=color_val,
            )
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(f"UMAP (n_neighbors={n_neighbors}, min_dist={min_dist})")
    if color:
        ax.legend(title=color)
    ax.grid(True, alpha=0.3)

    if save:
        _save_fig(md, fig, save)
    return fig, ax


__all__ = ["umap"]
