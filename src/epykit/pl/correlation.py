"""Clustered sample correlation heatmap ."""

from __future__ import annotations

import numpy as np
import polars as pl

from ._utils import _get_ax, _save_fig


def sample_correlation(
    md,
    *,
    method: str = "spearman",
    cluster: bool = True,
    ax=None,
    figsize=(7, 6),
    save: str | None = None,
):
    """Plot the all-vs-all sample correlation matrix as a heatmap.

    Reads ``md.uns["qc_sample_correlation"]`` (long-form DataFrame
    produced by :func:`epykit.qc.sample_correlation`) and renders it as
    a heatmap with optional hierarchical clustering.

    When the cached correlation table is absent, falls back to computing
    it inline via :func:`epykit.qc.sample_correlation` with the
    requested ``method``.
    """
    corr_df: pl.DataFrame
    if "qc_sample_correlation" in md.uns:
        corr_df = md.uns["qc_sample_correlation"]
    else:
        from ..qc import sample_correlation as _samp_corr
        samples = md.obs.get_column("sample_id").to_list()
        corr_df = _samp_corr(md.store, samples, method=method)
        md.uns["qc_sample_correlation"] = corr_df

    samples = sorted(set(
        corr_df.get_column("sample_a").to_list()
        + corr_df.get_column("sample_b").to_list()
    ))
    n = len(samples)
    if n == 0:
        raise ValueError("Empty correlation table")
    idx = {s: i for i, s in enumerate(samples)}
    mat = np.full((n, n), np.nan, dtype=np.float64)
    for row in corr_df.iter_rows(named=True):
        i = idx[row["sample_a"]]
        j = idx[row["sample_b"]]
        mat[i, j] = row["correlation"]

    order = list(range(n))
    if cluster and n >= 3:
        from scipy.cluster import hierarchy
        from scipy.spatial.distance import squareform
        # Convert correlation to distance (1 - r). Fill NaN with mean.
        dist = 1.0 - np.nan_to_num(mat, nan=float(np.nanmean(mat)))
        np.fill_diagonal(dist, 0.0)
        try:
            condensed = squareform((dist + dist.T) / 2.0, checks=False)
            linkage = hierarchy.linkage(condensed, method="average")
            order = hierarchy.leaves_list(linkage).tolist()
        except Exception:
            order = list(range(n))

    mat_ord = mat[np.ix_(order, order)]
    labels = [samples[i] for i in order]

    fig, ax = _get_ax(ax, figsize)
    im = ax.imshow(mat_ord, vmin=-1.0, vmax=1.0, cmap="RdBu_r", aspect="equal")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=90, fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)
    fig.colorbar(im, ax=ax, label=f"{method} rho")
    ax.set_title("Sample correlation")
    if save:
        _save_fig(md, fig, save)
    return fig, ax


__all__ = ["sample_correlation"]
