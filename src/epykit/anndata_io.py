"""AnnData export for ecosystem interop (scanpy / muon / multi-omics).

A MethylData object is conceptually shaped exactly like an AnnData:

    obs       = md.obs                 (n_samples x covariate-columns)
    var       = (chrom, pos)           (n_sites x 2)
    X         = beta  or  N_meth  or  coverage   (n_samples x n_sites)
    layers    = {"beta", "coverage", "N_meth", "N_unmeth"}

Memory strategy
---------------
A dense ``(n_samples x n_sites)`` matrix on real WGBS is enormous -- 8
samples x 28 M CpGs x 4 bytes (float32) ~= 900 MB per layer. The naive
"pivot the long-form DataFrame" approach holds the source rows, the
pivot, and the densified matrix simultaneously, which is typically 3-4x
the final array size and OOMs on real datasets.

To stay in the same ballpark as the final array, ``to_anndata`` runs
streamed, per-sample, per-chromosome:

1. Scan only ``(chrom, pos)`` lazily across the store to build the
   sorted site index and a ``(chrom, pos) -> row_idx`` lookup. This step
   uses a few hundred MB at WGBS scale.
2. Allocate ``X`` as a single ``float32`` array of shape
   ``(n_samples, n_sites)``.
3. For each sample, scan only that sample's partition for the requested
   layer column, look the site index up, and fill the matching row of
   ``X``. No long-form intermediate is ever materialised.

Additional layers (``populate_layers=True``) cost one extra dense array
per layer plus one extra streaming pass; default is **only the requested
layer**, so a user who just wants beta does not pay 4x memory.

A ``pp.unite`` is required before calling: with a united site set every
sample shares the same axis, the dense layout is sane, and the streaming
fill cannot leave holes.
"""

from __future__ import annotations

import logging
from pathlib import Path

import polars as pl

from .methyldata import MethylData

logger = logging.getLogger(__name__)


_VALID_LAYERS = {"beta", "coverage", "N_meth", "N_unmeth"}


def _site_index(store: str, *, unite_type: str):
    """Build the sorted (chrom, pos) site index plus per-chromosome position
    arrays for sorted-search lookup.

    Runs **one chromosome at a time** and reduces per-sample position
    arrays in numpy. Peak transient memory is roughly
    ``n_samples x max_chrom_sites x 8 B`` (~= 240 MB for chr1 on 6 samples),
    far below the multi-GB dedup-hashtable spike of a single-shot
    ``scan_parquet(...).unique()`` across the whole store.

    Parameters
    ----------
    store : str
        Methylstore root (Hive-partitioned by sample then chrom).
    unite_type : {"union", "intersect"}
        Per-chromosome reduction across samples. ``pp.unite`` is a
        metadata-only no-op (it doesn't rewrite the store), so the join
        has to happen here at export time.

    Returns
    -------
    var : pl.DataFrame
        Sorted (chrom, pos) DataFrame, length ``n_sites``.
    chrom_index : dict[str, tuple[int, np.ndarray]]
        ``{chrom: (start_idx, positions_int64)}`` -- the slice
        ``var[start_idx:start_idx + len(positions_int64)]`` is the
        chromosome's contiguous block. Empty chromosomes are dropped.
    """
    import numpy as np
    from functools import reduce

    if unite_type not in {"union", "intersect"}:
        raise ValueError(
            f"unite_type must be 'union' or 'intersect'; got {unite_type!r}"
        )

    store_p = Path(store)
    sample_dirs = sorted(d for d in store_p.glob("sample=*") if d.is_dir())
    if not sample_dirs:
        raise FileNotFoundError(f"No sample=* partitions under {store!r}")

    # Union of chromosome names seen across any sample. For 'intersect' we
    # still walk this union and let the per-chrom reduction short-circuit
    # to empty when a sample is missing the chrom -- matches DMC's
    # per-chrom join semantics.
    chrom_set: set[str] = set()
    for s_dir in sample_dirs:
        for d in s_dir.glob("chrom=*"):
            if d.is_dir():
                chrom_set.add(d.name.split("=", 1)[1])
    # Lex sort matches polars' Utf8 .sort("chrom") (chr1, chr10, chr11, ...,
    # chr2, ...) -- preserved for compatibility with prior output ordering.
    chroms_sorted = sorted(chrom_set)

    logger.info(
        "_site_index: %d sample(s), %d chromosome(s), unite=%s",
        len(sample_dirs), len(chroms_sorted), unite_type,
    )
    per_chrom_log = logger.debug if len(chroms_sorted) > 50 else logger.info

    chrom_index: dict[str, tuple[int, "np.ndarray"]] = {}
    chrom_arrays: list[tuple[str, "np.ndarray"]] = []
    running_start = 0

    for chrom in chroms_sorted:
        per_sample_pos: list[np.ndarray] = []
        empty = False
        for s_dir in sample_dirs:
            chrom_dir = s_dir / f"chrom={chrom}"
            if not chrom_dir.exists():
                if unite_type == "intersect":
                    empty = True
                    break
                continue
            arr = (
                pl.scan_parquet(f"{chrom_dir}/part-*.parquet")
                .select(pl.col("pos").cast(pl.Int64))
                .unique()
                .sort("pos")
                .collect()["pos"]
                .to_numpy(zero_copy_only=False)
            )
            if arr.dtype != np.int64:
                raise AssertionError(
                    f"_site_index: expected int64 positions on "
                    f"{s_dir.name}/{chrom}, got {arr.dtype}"
                )
            if len(arr) == 0:
                if unite_type == "intersect":
                    empty = True
                    break
                continue
            per_sample_pos.append(arr)

        if empty or not per_sample_pos:
            positions = np.empty(0, dtype=np.int64)
        elif unite_type == "union":
            positions = np.unique(np.concatenate(per_sample_pos))
        else:  # intersect
            positions = reduce(np.intersect1d, per_sample_pos)

        n = int(positions.size)
        per_chrom_log(
            "_site_index: chrom=%s n_sites=%d (start_idx=%d)",
            chrom, n, running_start,
        )
        if n == 0:
            continue
        chrom_index[chrom] = (running_start, positions)
        chrom_arrays.append((chrom, positions))
        running_start += n

    if not chrom_arrays:
        var = pl.DataFrame(
            {
                "chrom": pl.Series("chrom", [], dtype=pl.Utf8),
                "pos": pl.Series("pos", [], dtype=pl.Int64),
            }
        )
    else:
        all_pos = np.concatenate([a for _, a in chrom_arrays])
        chrom_col = pl.concat(
            [
                pl.Series("chrom", [c] * len(a), dtype=pl.Utf8)
                for c, a in chrom_arrays
            ]
        )
        var = pl.DataFrame(
            {
                "chrom": chrom_col,
                "pos": pl.Series("pos", all_pos, dtype=pl.Int64),
            }
        )

    return var, chrom_index


def _value_lazyframe(
    lf: "pl.LazyFrame",
    chrom: str,
    layer: str,
    *,
    pl_value_dtype: "pl.DataType",
) -> "pl.LazyFrame":
    """Return a lazy frame with columns (pos: Int64, value: pl_value_dtype) for one chromosome.

    The value dtype is pinned in polars so the downstream
    ``df["value"].to_numpy()`` cannot pick an unexpected numpy dtype.
    """
    lf = lf.filter(pl.col("chrom") == chrom)
    pos = pl.col("pos").cast(pl.Int64)
    if layer == "beta":
        return (
            lf.select(["pos", "N_meth", "coverage"])
            .filter(pl.col("coverage") > 0)
            .with_columns(
                (pl.col("N_meth").cast(pl_value_dtype) / pl.col("coverage").cast(pl_value_dtype))
                .alias("value")
            )
            .select([pos, pl.col("value").cast(pl_value_dtype)])
        )
    if layer == "coverage":
        return lf.select([pos, pl.col("coverage").cast(pl_value_dtype).alias("value")])
    if layer == "N_meth":
        return lf.select([pos, pl.col("N_meth").cast(pl_value_dtype).alias("value")])
    if layer == "N_unmeth":
        return (
            lf.select(["pos", "N_meth", "coverage"])
            .with_columns((pl.col("coverage") - pl.col("N_meth")).alias("value"))
            .select([pos, pl.col("value").cast(pl_value_dtype)])
        )
    raise ValueError(f"unknown layer {layer!r}")  # pragma: no cover


def _fill_layer(
    store: str,
    samples: list[str],
    n_sites: int,
    layer: str,
    chrom_index: dict[str, tuple[int, "object"]],
    dtype,
):
    """Stream sample x chromosome, filling a (n_samples x n_sites) array.

    For each (sample, chromosome) we read just that partition file, then
    use ``np.searchsorted`` against the chromosome's pre-built position
    array to compute the destination column indices in O(n_rows * log n_chrom).
    No Python-side dict is ever built -- the only persistent index
    structures are the per-chromosome int64 position arrays.
    """
    import numpy as np

    np_dtype = np.dtype(dtype)
    # Map numpy float dtype -> polars float dtype so the value column is
    # pinned end-to-end. Anything other than float32/float64 is rejected
    # by to_anndata() above, so this mapping is exhaustive.
    if np_dtype == np.float32:
        pl_value_dtype = pl.Float32
    elif np_dtype == np.float64:
        pl_value_dtype = pl.Float64
    else:
        raise ValueError(f"unsupported dtype {np_dtype!r}; use float32 or float64")

    out = np.full((len(samples), n_sites), np.nan, dtype=np_dtype)
    store_p = Path(store)

    for s_idx, sample in enumerate(samples):
        sample_dir = store_p / f"sample={sample}"
        if not sample_dir.exists():
            continue
        for chrom, (chrom_start, var_positions) in chrom_index.items():
            chrom_part = sample_dir / f"chrom={chrom}"
            if not chrom_part.exists():
                continue
            lf = pl.scan_parquet(f"{chrom_part}/part-*.parquet")
            df = _value_lazyframe(lf, chrom, layer, pl_value_dtype=pl_value_dtype).collect()
            if len(df) == 0:
                continue

            sample_positions = df["pos"].to_numpy(zero_copy_only=False)
            if sample_positions.dtype != np.int64:
                raise AssertionError(
                    f"_fill_layer: expected int64 sample_positions on {chrom!r}, "
                    f"got {sample_positions.dtype}"
                )
            values = df["value"].to_numpy(zero_copy_only=False)
            if values.dtype != np_dtype:
                raise AssertionError(
                    f"_fill_layer: expected {np_dtype} values on {chrom!r}, "
                    f"got {values.dtype}"
                )

            # Look up the column index of each sample CpG inside the
            # chromosome's var slice via a single sorted-search call.
            local_idx = np.searchsorted(var_positions, sample_positions)
            # Guard against positions that fall off the end of the array
            # (would happen on a *union* store where a sample contributes
            # a CpG that no other sample has -- searchsorted returns
            # len(arr)). Bounds-check first, then equality-check.
            in_range = local_idx < len(var_positions)
            local_idx_safe = np.where(in_range, local_idx, 0)
            hit = in_range & (var_positions[local_idx_safe] == sample_positions)
            if not hit.any():
                continue
            global_idx = chrom_start + local_idx_safe[hit]
            out[s_idx, global_idx] = values[hit]
    return out


def to_anndata(
    md: MethylData,
    *,
    layer: str = "beta",
    populate_layers: bool = False,
    dtype: str = "float32",
):
    """Materialise a MethylData as an AnnData object.

    Memory-conscious by default: only the layer you asked for is dense.
    On a typical 8-sample x 28 M-CpG run this produces a ~900 MB
    ``adata.X`` and nothing else, vs. ~3.5 GB of intermediate state under
    the previous pivot-based implementation. Pass
    ``populate_layers=True`` to also fill ``adata.layers["coverage"]``,
    etc., paying one extra dense array per layer.

    Parameters
    ----------
    md : MethylData
        Methylation data. Must have been ``pp.unite``'d so every sample
        shares the same site set.
    layer : {"beta", "coverage", "N_meth", "N_unmeth"}
        Which value to put in ``adata.X``. Default "beta".
    populate_layers : bool
        If True, also fill ``adata.layers`` with every layer in
        :data:`_VALID_LAYERS`. Default **False** -- the old default of
        True densified four matrices and was the most common cause of
        OOMs on real WGBS. Layer matrices share the site index so an
        extra pass per layer is cheap on disk but doubles dense RAM.
    dtype : str
        NumPy dtype for the dense matrices. Default ``"float32"``
        (4 bytes / cell). Use ``"float64"`` if you want to feed the
        result into algorithms that demand it; halves the per-sample row
        fill speed and doubles RAM.
    Returns
    -------
    anndata.AnnData
        Shape ``(n_samples, n_sites)`` matching ``md.obs`` x the united
        site axis.
    """
    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError(
            "anndata is required for to_anndata(). "
            "Install it with: pip install 'epykit[anndata]'"
        ) from exc

    if not md._united:
        raise ValueError(
            "to_anndata() requires a united methylstore so all samples share "
            "the same site axis. Run ep.pp.unite(md, type='intersect') first "
            "(or type='union' if you want every site that appears in any "
            "sample -- note the densification cost)."
        )

    if layer not in _VALID_LAYERS:
        raise ValueError(
            f"layer must be one of {sorted(_VALID_LAYERS)}; got {layer!r}"
        )

    import numpy as np
    import pandas as pd

    np_dtype = np.dtype(dtype)
    samples = md.obs.get_column("sample_id").to_list()

    unite_type = md.uns["unite"]["type"]
    logger.info(
        "to_anndata: building site index from %s (unite=%s)",
        md.store, unite_type,
    )
    var, chrom_index = _site_index(md.store, unite_type=unite_type)
    n_sites = len(var)
    dense_gib = len(samples) * n_sites * np_dtype.itemsize / (1024 ** 3)
    logger.info(
        "to_anndata: %d sample(s) x %d site(s) -- estimated dense size %.2f GiB per layer",
        len(samples), n_sites, dense_gib,
    )
    # Loud warning if the user is about to allocate something huge. 4 GiB
    # crosses the threshold where even 32 GiB workstations start swapping
    # once anndata's own copies are factored in.
    if dense_gib > 4.0:
        logger.warning(
            "to_anndata: dense X is %.1f GiB. Consider filtering to a smaller "
            "site set first (e.g. ep.pp.filter_coverage + ep.pp.unite with "
            "type='intersect', or ep.pp.aggregate_regions to a BED of "
            "promoters/peaks) before exporting.",
            dense_gib,
        )

    X = _fill_layer(md.store, samples, n_sites, layer, chrom_index, np_dtype)

    obs_pd = md.obs.to_pandas().set_index("sample_id")

    # Build var via Arrow -> pandas in one shot (no Python-level 42M-element
    # lists). The "{chrom}:{pos}" index is built vectorised on pandas
    # Series rather than with a Python list comprehension, which on real
    # WGBS removes several GB of transient Python object overhead.
    var_pd = var.to_pandas()
    var_pd.index = (var_pd["chrom"].astype(str) + ":" + var_pd["pos"].astype(str)).values

    adata = ad.AnnData(X=X, obs=obs_pd, var=var_pd)

    if populate_layers:
        for extra in _VALID_LAYERS:
            if extra == layer:
                continue
            logger.info("to_anndata: filling layer %r", extra)
            adata.layers[extra] = _fill_layer(
                md.store, samples, n_sites, extra, chrom_index, np_dtype,
            )

    adata.uns["epykit_assembly"] = md.assembly
    adata.uns["epykit_context"] = md.context
    adata.uns["epykit_state"] = list(md.state)

    logger.info(
        "to_anndata: built AnnData shape=%s layers=%s",
        adata.shape, list(adata.layers.keys()),
    )
    return adata


def to_mudata(
    md,
    *,
    layer: str = "beta",
    other_modalities: dict | None = None,
):
    """Return a ``MuData`` with methylation as ``'meth'`` modality.

    Parameters
    ----------
    md : MethylData
        Must be ``pp.unite``-d first (same precondition as ``to_anndata``).
    layer : str
        Which methylation matrix to embed as the methylation modality's
        ``X``. Forwarded to :func:`to_anndata`.
    other_modalities : dict[str, AnnData], optional
        Additional modalities to bundle into the MuData.
    """
    try:
        import mudata as md_lib  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "mudata is required for to_mudata. "
            "Install with: pip install 'epykit[anndata]'"
        ) from exc

    adata = to_anndata(md, layer=layer)
    modalities = {"meth": adata}
    if other_modalities:
        modalities.update(other_modalities)
    return md_lib.MuData(modalities)


__all__ = ["to_anndata", "to_mudata"]