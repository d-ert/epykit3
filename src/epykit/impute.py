"""kNN methylation beta imputation.

Useful for filling sparse coverage gaps before PCA / UMAP / clustering /
regression-on-beta analyses where missing values either propagate as NaN
through the downstream math or force every analysis to do its own
imputation. Two entry points:

* :func:`impute_knn_beta` -- pure-numpy: ``(n_samples, n_sites)`` beta
  matrix in, imputed matrix out. Per-chromosome inverse-distance kNN.
  No heavy dependencies.
* :func:`impute_knn_anndata` -- operates on an :class:`anndata.AnnData`
  with sample x site layout (the orientation epykit's ``to_anndata``
  emits). Imputes one chromosome at a time using ``adata.var['chrom']``
  / ``adata.var['pos']`` so the same algorithm works on a stacked
  multi-chromosome matrix.

The model is intentionally simple: missing beta at (sample s, site j)
becomes the inverse-distance-weighted mean of beta values at sample s
across the k nearest *covered* CpGs within ``max_distance_bp`` of
position j on the same chromosome. No cross-sample sharing (so it
won't smooth across treatment / control groups by accident); no
spatial Gaussian beyond inverse-distance weighting; no probabilistic
imputation (use mice / ashr-style methods if you need uncertainty
propagation).
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


def impute_knn_beta(
    positions: np.ndarray,
    beta_matrix: np.ndarray,
    *,
    k: int = 5,
    max_distance_bp: Optional[int] = 10_000,
) -> np.ndarray:
    """Fill NaN entries of ``beta_matrix`` by k-nearest-neighbour kNN.

    Parameters
    ----------
    positions : (n_sites,) int
        Genomic position of each site. Must be sorted ascending and
        contain CpGs from a *single* chromosome -- call this once per
        chromosome and concatenate the result, or use
        :func:`impute_knn_anndata`.
    beta_matrix : (n_samples, n_sites) float
        Per-sample beta values; missing entries are NaN. Modified copy is
        returned; the input is not mutated.
    k : int
        Number of nearest covered CpGs used in the inverse-distance
        weighted average. Default 5.
    max_distance_bp : int or None
        If set, drop candidate neighbours farther than this many bp
        from the missing site. ``None`` disables the cap (use the
        nearest k regardless of distance). Default 10 kb (typical
        WGBS CpG-correlation decay scale).

    Returns
    -------
    np.ndarray
        Copy of ``beta_matrix`` with NaN entries replaced where at
        least one neighbour was found. Entries where every candidate
        was filtered out by ``max_distance_bp`` remain NaN.
    """
    if positions.ndim != 1:
        raise ValueError("positions must be 1-D")
    if beta_matrix.ndim != 2 or beta_matrix.shape[1] != positions.shape[0]:
        raise ValueError(
            f"beta_matrix shape {beta_matrix.shape} doesn't align with "
            f"positions length {positions.shape[0]}"
        )
    if not np.all(np.diff(positions) >= 0):
        raise ValueError("positions must be sorted ascending")
    if k < 1:
        raise ValueError(f"k must be >=1, got {k}")

    n_samples, n_sites = beta_matrix.shape
    out = beta_matrix.astype(np.float64, copy=True)

    for s in range(n_samples):
        row = out[s]
        covered_mask = ~np.isnan(row)
        if not covered_mask.any():
            continue  # this sample has nothing to impute *from*
        covered_idx = np.where(covered_mask)[0]
        covered_pos = positions[covered_idx]
        covered_beta = row[covered_idx]
        missing_idx = np.where(~covered_mask)[0]
        if missing_idx.size == 0:
            continue

        # For each missing site, locate the 2k window of nearest
        # covered sites via binary search, then take the actual top-k
        # by absolute distance.
        for j in missing_idx:
            pos_j = positions[j]
            # search_pos is the insertion index in covered_pos.
            search_pos = np.searchsorted(covered_pos, pos_j)
            lo = max(0, search_pos - k)
            hi = min(covered_pos.size, search_pos + k)
            if hi <= lo:
                continue
            cand_pos = covered_pos[lo:hi]
            cand_beta = covered_beta[lo:hi]
            distances = np.abs(cand_pos - pos_j)
            if max_distance_bp is not None:
                mask = distances <= max_distance_bp
                if not mask.any():
                    continue
                cand_beta = cand_beta[mask]
                distances = distances[mask]
            # Top-k by smallest distance.
            if distances.size > k:
                order = np.argpartition(distances, k)[:k]
                cand_beta = cand_beta[order]
                distances = distances[order]
            # Inverse-distance weights (+1 to avoid /0 at exact-overlap
            # positions, which can happen on multi-context inputs).
            weights = 1.0 / (distances + 1.0)
            out[s, j] = float(np.sum(weights * cand_beta) / weights.sum())

    return out


def impute_knn_anndata(
    adata,
    *,
    k: int = 5,
    max_distance_bp: Optional[int] = 10_000,
    layer: Optional[str] = None,
    inplace: bool = False,
):
    """Impute missing beta values in an :class:`anndata.AnnData`.

    The AnnData must follow epykit's convention: rows = samples, cols =
    sites, with ``adata.var['chrom']`` and ``adata.var['pos']`` set. The
    matrix to impute is ``adata.X`` by default, or
    ``adata.layers[layer]`` when ``layer`` is given.

    Parameters
    ----------
    adata : anndata.AnnData
        Methylation AnnData from :func:`epykit.anndata_io.to_anndata`.
    k, max_distance_bp
        Forwarded to :func:`impute_knn_beta`.
    layer : str, optional
        Which layer to impute. Default ``None`` -> ``adata.X``.
    inplace : bool
        If True, write the imputed matrix back into the source layer
        and return ``adata``. If False (default), return the imputed
        matrix as a fresh ndarray and leave ``adata`` untouched.

    Returns
    -------
    np.ndarray or anndata.AnnData
        Imputed matrix (or modified ``adata`` when ``inplace=True``).
    """
    import numpy as np

    for required in ("chrom", "pos"):
        if required not in adata.var.columns:
            raise ValueError(
                f"impute_knn_anndata: adata.var['{required}'] is missing. "
                "Build the AnnData via ep.to_anndata(md) so the site axis "
                "carries chrom / pos."
            )
    if layer is None:
        beta = np.asarray(adata.X, dtype=np.float64).copy()
    else:
        if layer not in adata.layers:
            raise ValueError(f"layer {layer!r} not in adata.layers")
        beta = np.asarray(adata.layers[layer], dtype=np.float64).copy()

    chroms = adata.var["chrom"].to_numpy()
    positions_all = adata.var["pos"].to_numpy().astype(np.int64)

    # Impute one chromosome at a time so cross-chromosome distances
    # can't accidentally pull in a neighbour from a different contig.
    for chrom in np.unique(chroms):
        col_mask = chroms == chrom
        cols = np.where(col_mask)[0]
        sub_positions = positions_all[cols]
        order = np.argsort(sub_positions, kind="stable")
        sorted_cols = cols[order]
        sub_beta = beta[:, sorted_cols]
        imputed = impute_knn_beta(
            sub_positions[order], sub_beta,
            k=k, max_distance_bp=max_distance_bp,
        )
        beta[:, sorted_cols] = imputed

    if inplace:
        if layer is None:
            adata.X = beta
        else:
            adata.layers[layer] = beta
        return adata
    return beta


__all__ = ["impute_knn_beta", "impute_knn_anndata"]
