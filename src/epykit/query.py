"""Tabix-style random-access queries on the partitioned Parquet methylstore.

This module exposes :func:`query_region`, :func:`query_regions`, and
:func:`query_sites` -- three single-purpose entry points for fetching
methylation data at specific genomic loci without scanning the whole
store. New in 0.4.0.

Why no extra dependency?
------------------------
The store is already hive-partitioned by ``sample=*/chrom=*``. Polars
lazy I/O (``pl.scan_parquet``) reads parquet row-group statistics --
``min(pos)`` / ``max(pos)`` per row group, written automatically when
``convert_sample`` ran -- and prunes row groups whose statistics don't
overlap the query before reading any data. That gives tabix-equivalent
random access on the existing store, with no separate index file.

The chrom partition is selected by direct file-path construction
(O(1)); the pos filter inside that partition uses predicate pushdown
(no separate index needed; polars + pyarrow handle it transparently).
Stores written with a small ``row_group_size`` (e.g. 10 000) get the
fastest single-locus queries; stores written with the default
(1 000 000) still benefit from the chrom-partition filter and merely
read more rows per query.

Examples
--------
::

    import epykit as ep
    df = ep.query.query_region(md.store, "chr7", 140_453_000, 140_500_000)
    # -> pl.DataFrame with one row per (sample_id, pos) in the region.

    regions = pl.DataFrame({"chrom": ["chr1", "chr2"],
                             "start": [1_000_000, 2_000_000],
                             "end":   [1_100_000, 2_100_000]})
    df = ep.query.query_regions(md.store, regions, samples=["s1", "s2"])

    sites = pl.DataFrame({"chrom": ["chr1", "chr1"], "pos": [12345, 67890]})
    df = ep.query.query_sites(md.store, sites)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Optional

import polars as pl

logger = logging.getLogger(__name__)


_DEFAULT_COLUMNS = ["pos", "strand", "N_meth", "coverage"]


def _list_samples(store: Path) -> list[str]:
    return sorted(
        d.name.removeprefix("sample=")
        for d in store.glob("sample=*")
        if d.is_dir()
    )


def _resolve_samples(store: Path, samples: Optional[Iterable[str]]) -> list[str]:
    if samples is None:
        return _list_samples(store)
    return list(samples)


def query_region(
    store: str | Path,
    chrom: str,
    start: int,
    end: int,
    *,
    samples: Optional[Iterable[str]] = None,
) -> pl.DataFrame:
    """Fetch methylation data for one genomic region across samples.

    Parameters
    ----------
    store
        Methylstore root (the same path you pass to :func:`load`).
    chrom
        Chromosome name (must match the partition key, e.g. ``"chr1"``).
    start, end
        Half-open BED-style interval. ``start`` is inclusive, ``end``
        is exclusive -- matches the convention used everywhere else in
        epykit (see the ``convert.py`` module docstring on coordinates).
    samples
        Sample ids to include. ``None`` returns every sample present in
        the store.

    Returns
    -------
    pl.DataFrame
        Columns ``(sample_id, chrom, pos, strand, N_meth, coverage,
        beta)``. ``beta`` is computed on the fly as ``N_meth / coverage``
        (NaN where ``coverage == 0``). Empty frame (zero rows) when no
        rows overlap the region.
    """
    if end <= start:
        raise ValueError(f"end ({end}) must be > start ({start})")
    store_p = Path(store)
    samples_resolved = _resolve_samples(store_p, samples)
    if not samples_resolved:
        return _empty_frame()

    parts: list[pl.DataFrame] = []
    for sample in samples_resolved:
        part_file = store_p / f"sample={sample}" / f"chrom={chrom}" / "part-0.parquet"
        if not part_file.exists():
            continue
        # scan_parquet enables row-group pruning via min/max statistics.
        # The filter (start <= pos < end) is pushed down so partitions
        # entirely outside the range are skipped without reading data.
        df = (
            pl.scan_parquet(str(part_file))
            .filter((pl.col("pos") >= start) & (pl.col("pos") < end))
            .select(_DEFAULT_COLUMNS)
            .collect()
        )
        if df.height == 0:
            continue
        df = df.with_columns(
            pl.lit(sample).alias("sample_id"),
            pl.lit(chrom).alias("chrom"),
        )
        parts.append(df)

    if not parts:
        return _empty_frame()
    out = pl.concat(parts, how="vertical_relaxed")
    out = out.with_columns(
        (pl.col("N_meth") / pl.col("coverage").cast(pl.Float64))
        .fill_nan(None).alias("beta")
    )
    return out.select(["sample_id", "chrom", "pos", "strand", "N_meth", "coverage", "beta"]).sort(
        ["sample_id", "pos"]
    )


def query_regions(
    store: str | Path,
    regions: pl.DataFrame,
    *,
    samples: Optional[Iterable[str]] = None,
) -> pl.DataFrame:
    """Batched region query -- concatenates results across multiple regions.

    ``regions`` must have at least columns ``chrom``, ``start``, ``end``
    (any extra columns are ignored). Returns the same schema as
    :func:`query_region` plus a ``region_id`` integer column indicating
    the row index of the source region.
    """
    required = {"chrom", "start", "end"}
    missing = required - set(regions.columns)
    if missing:
        raise ValueError(f"regions missing required columns: {sorted(missing)}")

    parts: list[pl.DataFrame] = []
    for region_id, row in enumerate(regions.iter_rows(named=True)):
        df = query_region(
            store, row["chrom"], int(row["start"]), int(row["end"]),
            samples=samples,
        )
        if df.height:
            parts.append(df.with_columns(pl.lit(region_id).cast(pl.Int32).alias("region_id")))

    if not parts:
        return _empty_frame().with_columns(pl.lit(None).cast(pl.Int32).alias("region_id"))
    return pl.concat(parts, how="vertical_relaxed")


def query_sites(
    store: str | Path,
    sites: pl.DataFrame,
    *,
    samples: Optional[Iterable[str]] = None,
) -> pl.DataFrame:
    """Exact-position queries -- return rows at specified (chrom, pos) sites.

    Useful for epigenetic clock CpGs, validation against reference panels,
    or any analysis that targets a fixed set of CpGs rather than a
    continuous region.

    ``sites`` must have columns ``chrom`` and ``pos``. The result is the
    same long-form schema as :func:`query_region` but restricted to
    exactly the requested positions (the store may contain more positions
    in between that are dropped).
    """
    required = {"chrom", "pos"}
    missing = required - set(sites.columns)
    if missing:
        raise ValueError(f"sites missing required columns: {sorted(missing)}")

    # Group positions by chromosome so we open each chrom partition once.
    by_chrom = sites.group_by("chrom").agg(pl.col("pos").alias("positions"))
    store_p = Path(store)
    samples_resolved = _resolve_samples(store_p, samples)

    parts: list[pl.DataFrame] = []
    for chrom_row in by_chrom.iter_rows(named=True):
        chrom = chrom_row["chrom"]
        positions = chrom_row["positions"]
        if len(positions) == 0:
            continue
        pos_min = int(min(positions))
        pos_max = int(max(positions)) + 1  # half-open upper bound for the prefilter
        pos_set = set(int(p) for p in positions)

        for sample in samples_resolved:
            part_file = store_p / f"sample={sample}" / f"chrom={chrom}" / "part-0.parquet"
            if not part_file.exists():
                continue
            df = (
                pl.scan_parquet(str(part_file))
                .filter((pl.col("pos") >= pos_min) & (pl.col("pos") < pos_max))
                .select(_DEFAULT_COLUMNS)
                .collect()
            )
            if df.height == 0:
                continue
            df = df.filter(pl.col("pos").is_in(pos_set))
            if df.height == 0:
                continue
            df = df.with_columns(
                pl.lit(sample).alias("sample_id"),
                pl.lit(chrom).alias("chrom"),
            )
            parts.append(df)

    if not parts:
        return _empty_frame()
    out = pl.concat(parts, how="vertical_relaxed")
    out = out.with_columns(
        (pl.col("N_meth") / pl.col("coverage").cast(pl.Float64))
        .fill_nan(None).alias("beta")
    )
    return out.select(["sample_id", "chrom", "pos", "strand", "N_meth", "coverage", "beta"]).sort(
        ["sample_id", "chrom", "pos"]
    )


def _empty_frame() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            "sample_id": pl.Utf8,
            "chrom": pl.Utf8,
            "pos": pl.Int32,
            "strand": pl.Utf8,
            "N_meth": pl.Int32,
            "coverage": pl.Int32,
            "beta": pl.Float64,
        }
    )


__all__ = ["query_region", "query_regions", "query_sites"]
