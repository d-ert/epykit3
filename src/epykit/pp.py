from __future__ import annotations

import logging
import shutil
import warnings
from pathlib import Path

import polars as pl

from . import filter as filter_mod
from ._cache import count_store_rows as _count_parquet_rows
from .methyldata import MethylData

logger = logging.getLogger(__name__)


_BED_BASE_COLS = ["chrom", "start", "end"]
_BED_EXTRA_NAMES = ["name", "score", "strand", "thickStart", "thickEnd",
                    "itemRgb", "blockCount", "blockSizes", "blockStarts"]


def _append_store_history(md: MethylData, step: str, path: str, n_sites: int | None) -> None:
    history = md.uns.get("_store_history")
    if not isinstance(history, list):
        history = []
    history.append({"step": step, "path": path, "n_sites": n_sites})
    md.uns["_store_history"] = history


def filter_coverage(
    md: MethylData,
    lo_count: int = 10,
    hi_perc: float = 99.9,
    blacklist_bed: str | None = None,
    output_store: str | None = None,
) -> None:
    """Coverage filtering in-place on a MethylData object."""
    quantile = hi_perc / 100.0 if hi_perc > 1 else hi_perc
    if quantile <= 0 or quantile > 1:
        raise ValueError(f"Invalid hi_perc/quantile value: {hi_perc}")

    # Derive output store path: explicit override, or from analysis_root cache, or legacy behavior
    if output_store:
        out = output_store
    elif md.analysis_root:
        out = str(Path(md.analysis_root) / ".cache" / "filtered")
    else:
        out = f"{md.store}_filtered"

    filter_mod.filter_sites(
        methylstore_path=md.store,
        output_dir=out,
        min_coverage=lo_count,
        max_coverage_quantile=quantile,
        blacklist_bed=blacklist_bed,
    )
    md.store = out
    # _filtered is a derived property; appending to _store_history below is
    # what makes md._filtered evaluate to True.
    md.uns["filter"] = {
        "lo_count": lo_count,
        "hi_perc": hi_perc,
        "blacklist_bed": blacklist_bed,
    }
    n_sites = _count_parquet_rows(out)
    if n_sites is not None:
        md.uns["n_sites_filtered"] = n_sites
        _append_store_history(md, "filtered", out, n_sites)


def normalize_coverage(md: MethylData, method: str = "median") -> None:
    """Per-sample coverage normalisation, in-place on a MethylData.

    Computes a per-sample scaling factor so that each sample's central
    coverage statistic (median by default, or mean) matches a common
    target -- the median (or mean) of the per-sample summaries. Read
    counts are scaled and re-integerised, then ``md.store`` is repointed
    at a new ``.cache/normalized`` (or ``<store>_normalized``) partition.

    This prevents deeper-sequenced samples from dominating pooled-count
    tile / region tests downstream. The per-CpG score test in
    ``ep.tl.dmc`` is much less sensitive to coverage imbalance, but
    ``ep.tl.dmr(method='tile')`` is -- running normalisation between
    coverage filtering and tile aggregation removes that bias.

    .. warning::

        This rescales integer read counts by a per-sample factor and rounds
        back to integers, so it **fabricates / discards reads**: a low-depth
        sample (factor > 1) has its counts scaled *up*, inflating the
        apparent precision the binomial / quasi-binomial count model reads
        off coverage. It is **optional**, not a default step -- reach for it
        when downstream coverage imbalance is a real problem (notably
        ``tl.dmr(method='tile')``), and prefer leaving it off for the
        per-CpG ``tl.dmc`` engines, which already weight by coverage and are
        little affected by imbalance. Same trade-off as ``methylKit``'s
        ``normalizeCoverage``.

    Call order: ``filter_coverage`` -> ``normalize_coverage`` -> ``unite``.

    Parameters
    ----------
    md : MethylData
        Object whose store has been ``filter_coverage``'d.
    method : {"median", "mean"}
        Central statistic to align. ``"median"`` is the robust default
        -- robust to extreme-coverage tails.

    Raises
    ------
    ValueError
        If ``filter_coverage`` has not been called yet.
    """
    if not md._filtered:
        raise ValueError(
            "Run ep.pp.filter_coverage(md) before ep.pp.normalize_coverage(md)."
        )
    if md._united:
        logger.warning(
            "normalize_coverage called after unite(); the recommended order "
            "is filter -> normalize -> unite. Re-running unite() is a no-op "
            "but downstream stats reflect the normalised store."
        )

    if md.analysis_root:
        out = str(Path(md.analysis_root) / ".cache" / "normalized")
    else:
        out = f"{md.store}_normalized"

    factors = filter_mod.normalize_coverage_store(
        methylstore_path=md.store,
        output_dir=out,
        method=method,
    )

    md.store = out
    md.uns["normalize"] = {
        "method": method,
        "factors": factors,
    }
    n_sites = _count_parquet_rows(out)
    if n_sites is not None:
        md.uns["n_sites_normalized"] = n_sites
        _append_store_history(md, "normalized", out, n_sites)


def set_unite_type(md: MethylData, type: str = "union") -> None:
    """Record the site-alignment strategy for downstream DMC processing.

    This is a *state-marker*: it writes ``md.uns["unite"] = {"type": type}``
    and does not materialise the intersection or union. ``ep.tl.dmc`` reads
    this state and passes ``unite=True/False`` to ``process_chromosomes_dmc``,
    which performs the per-chromosome join lazily in O(n_sites) memory --
    identical to the old procedural API.

    Parameters
    ----------
    md : MethylData
    type : {"intersect", "union"}
        Site-set strategy:
          * ``"intersect"`` -- only sites covered in every sample.
          * ``"union"`` -- every site covered in at least one sample.

    Notes
    -----
    Eagerly computing the full intersection (previously stored in
    ``md.uns["site_intersect"]``) caused an OOM on whole-genome data
    because it loaded all 338 M+ rows into RAM at once.
    """
    if type not in {"intersect", "union"}:
        raise ValueError("type must be 'intersect' or 'union'")
    # _united is a derived property -- recording in uns is enough.
    md.uns["unite"] = {"type": type}
    _append_store_history(md, "united", md.store, None)


def unite(md: MethylData, type: str = "union") -> None:
    """Deprecated alias for ``set_unite_type``.

    The original ``unite()`` name suggested a verb that performs the
    union, but this function only writes ``md.uns["unite"]`` and does
    not materialise anything. Use ``set_unite_type()`` for an honest name.

    Scheduled for removal in epykit 2.0.
    """
    warnings.warn(
        "pp.unite() is deprecated; use pp.set_unite_type() for the same "
        "semantics with an honest name. pp.unite() will be removed in 2.0.",
        DeprecationWarning,
        stacklevel=2,
    )
    set_unite_type(md, type=type)


def _read_bed(regions_bed: str, region_id_col: str | None = None) -> pl.DataFrame:
    """Parse a 3/4/6/12-column BED into a DataFrame with columns
    ``chrom, start, end, region_id``.

    BED files are tab-separated, no header. The optional ``region_id_col``
    overrides the default ``name`` (col 4); if absent or the file has fewer
    than 4 columns, region_id falls back to ``"{chrom}:{start}-{end}"``.
    Lines beginning with ``#``, ``track``, or ``browser`` are skipped.
    """
    path = Path(regions_bed)
    if not path.exists():
        raise FileNotFoundError(f"BED file not found: {regions_bed}")

    raw = pl.read_csv(
        str(path),
        separator="\t",
        has_header=False,
        comment_prefix="#",
        infer_schema_length=0,  # everything as Utf8 then cast
        truncate_ragged_lines=True,
    )
    # Drop track/browser lines
    raw = raw.filter(
        ~pl.col("column_1").str.starts_with("track")
        & ~pl.col("column_1").str.starts_with("browser")
    )
    n_cols = raw.width
    if n_cols < 3:
        raise ValueError(
            f"BED file {regions_bed} has only {n_cols} column(s); need at least 3."
        )
    new_names = _BED_BASE_COLS + _BED_EXTRA_NAMES[: max(0, n_cols - 3)]
    raw = raw.rename({f"column_{i+1}": new_names[i] for i in range(min(n_cols, len(new_names)))})

    raw = raw.with_columns([
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
    ])

    if region_id_col and region_id_col in raw.columns:
        raw = raw.with_columns(pl.col(region_id_col).cast(pl.Utf8).alias("region_id"))
    elif "name" in raw.columns:
        raw = raw.with_columns(pl.col("name").cast(pl.Utf8).alias("region_id"))
    else:
        raw = raw.with_columns(
            (pl.col("chrom").cast(pl.Utf8)
             + pl.lit(":")
             + pl.col("start").cast(pl.Utf8)
             + pl.lit("-")
             + pl.col("end").cast(pl.Utf8)
            ).alias("region_id")
        )

    return raw.select(["chrom", "start", "end", "region_id"])


def aggregate_regions(
    md: MethylData,
    regions_bed: str,
    *,
    region_id_col: str | None = None,
    output_store: str | None = None,
    min_cpgs_per_region: int = 1,
) -> None:
    """Aggregate per-CpG methylation counts within user-supplied BED regions.

    After this call, ``md.store`` points at a new partitioned Parquet
    store whose rows are *regions* rather than CpGs, but with a schema
    compatible with the rest of the pipeline so ``ep.tl.dmc(md)`` runs
    unchanged. Each region row carries:

    * ``chrom``, ``pos`` (region midpoint, Int32), ``strand`` ("*")
    * ``context`` (inherited from md.context)
    * ``N_meth`` and ``N_unmeth`` summed across all CpGs in the region
    * ``coverage`` rebuilt as ``N_meth + N_unmeth``
    * ``sample``
    * Extras: ``region_id``, ``start``, ``end``, ``n_cpgs``

    Parameters
    ----------
    md : MethylData
        Object whose store will be re-aggregated.
    regions_bed : str
        Path to a BED file (3/4/6/12 columns). Coordinates 0-based,
        half-open -- same convention as the methylstore.
    region_id_col : str, optional
        Name of the BED column to use as the region identifier. Defaults
        to ``name`` (column 4) if present, otherwise to
        ``"{chrom}:{start}-{end}"``.
    output_store : str, optional
        Override output path. Defaults to ``<analysis_root>/.cache/regions``
        when an analysis root is set, else ``<md.store>_regions``.
    min_cpgs_per_region : int
        Drop region rows with fewer than this many CpGs (per sample) after
        aggregation. Default 1 (keep every overlapped region).

    Notes
    -----
    Regions with zero overlapping CpGs in a sample are not emitted for
    that sample. Use ``ep.pp.unite(md, type="union")`` afterwards if you
    want union semantics across samples.

    The overlap is computed with a sorted-search on per-chromosome CpG /
    region arrays, so no additional dependencies beyond polars + numpy
    are needed (bioframe, used elsewhere in epykit, is *not* required
    here).
    """
    bed = _read_bed(regions_bed, region_id_col=region_id_col)
    n_regions = len(bed)
    if n_regions == 0:
        raise ValueError(f"BED file {regions_bed} contains no regions")

    if output_store:
        out = output_store
    elif md.analysis_root:
        out = str(Path(md.analysis_root) / ".cache" / "regions")
    else:
        out = f"{md.store}_regions"

    out_path = Path(out)
    src = Path(md.store)
    # Clear any stale store from a previous run before writing. filter_sites
    # and normalize_coverage_store both rmtree their output for this reason;
    # without it, re-running aggregate_regions with a different BED (or a
    # stricter min_cpgs_per_region) left prior chrom partitions on disk and
    # md.store was repointed at a store mixing regions from two BEDs (M6).
    # Guard against nuking the source: if out_path IS the current store
    # (a degenerate re-aggregation on an already-aggregated md), refuse
    # rather than delete the data we are about to read.
    if out_path.exists() and out_path.resolve() == src.resolve():
        raise ValueError(
            f"aggregate_regions output {out_path} is the current md.store; "
            f"re-aggregate from the original (pre-regions) store instead."
        )
    if out_path.exists():
        shutil.rmtree(out_path)
    out_path.mkdir(parents=True, exist_ok=True)

    sample_dirs = sorted(src.glob("sample=*"))
    if not sample_dirs:
        raise ValueError(f"No sample=* directories found in {md.store}")

    context_value = md.context or "CpG"

    # Per-chromosome bed dict for fast lookup
    bed_by_chrom: dict[str, pl.DataFrame] = {
        c: bed.filter(pl.col("chrom") == c).sort("start")
        for c in bed["chrom"].unique().to_list()
    }

    for sample_dir in sample_dirs:
        sample = sample_dir.name.removeprefix("sample=")
        for chrom_dir in sorted(sample_dir.glob("chrom=*")):
            chrom = chrom_dir.name.removeprefix("chrom=")
            regions_chrom = bed_by_chrom.get(chrom)
            if regions_chrom is None or len(regions_chrom) == 0:
                continue
            parts = list(chrom_dir.glob("part-*.parquet"))
            if not parts:
                continue
            cpgs = pl.concat([pl.read_parquet(str(p)) for p in parts])
            if len(cpgs) == 0:
                continue

            # Range-join: each CpG is assigned to EVERY region it overlaps
            # (half-open [start, end)). Overlapping / nested BED regions are
            # common (gene bodies, promoters, enhancers) and a CpG inside
            # several regions must contribute to all of them.
            assigned = _assign_cpgs_to_regions(cpgs, regions_chrom)
            if assigned is None or len(assigned) == 0:
                continue

            agg = (
                assigned
                .group_by("region_id", maintain_order=True)
                .agg([
                    pl.col("start").first(),
                    pl.col("end").first(),
                    pl.sum("N_meth").alias("N_meth"),
                    pl.sum("N_unmeth").alias("N_unmeth"),
                    pl.len().alias("n_cpgs"),
                ])
                .filter(pl.col("n_cpgs") >= min_cpgs_per_region)
            )
            if len(agg) == 0:
                continue

            agg = agg.with_columns([
                pl.lit(chrom).alias("chrom"),
                ((pl.col("start") + pl.col("end")) // 2).cast(pl.Int32).alias("pos"),
                pl.lit("*").alias("strand"),
                pl.lit(context_value).alias("context"),
                (pl.col("N_meth") + pl.col("N_unmeth")).cast(pl.Int32).alias("coverage"),
                pl.lit(sample).alias("sample"),
                pl.col("N_meth").cast(pl.Int32),
                pl.col("N_unmeth").cast(pl.Int32),
                pl.col("n_cpgs").cast(pl.Int32),
            ]).select([
                "chrom", "pos", "strand", "context", "N_meth", "N_unmeth",
                "coverage", "sample", "region_id", "start", "end", "n_cpgs",
            ]).sort("pos")

            out_chrom_dir = out_path / f"sample={sample}" / f"chrom={chrom}"
            out_chrom_dir.mkdir(parents=True, exist_ok=True)
            agg.write_parquet(
                str(out_chrom_dir / "part-0.parquet"),
                compression="zstd",
            )

    md.store = str(out_path)
    md.uns["regions"] = {
        "bed": str(Path(regions_bed).resolve()),
        "n_regions": int(n_regions),
        "region_id_col": region_id_col,
        "min_cpgs_per_region": int(min_cpgs_per_region),
    }
    n_rows = _count_parquet_rows(str(out_path))
    if n_rows is not None:
        md.uns["n_sites_regions"] = n_rows
        _append_store_history(md, "regions", str(out_path), n_rows)
    logger.info(
        "aggregate_regions: %d region(s) defined, store at %s", n_regions, out_path
    )


def _assign_cpgs_to_regions(
    cpgs: pl.DataFrame, regions_chrom: pl.DataFrame
) -> pl.DataFrame | None:
    """Assign each CpG to EVERY region it overlaps (half-open intervals).

    Returns a DataFrame with the CpG columns plus ``region_id``, ``start``,
    ``end`` -- one row per (CpG, overlapping region) pair, so a CpG inside
    nested or overlapping regions contributes to all of them. Returns None
    if no CpGs overlap any region.

    The previous implementation kept only the single region at
    ``searchsorted(starts, pos) - 1`` (the last region whose start <= pos),
    which silently dropped CpGs that fell past an inner region's end inside
    an outer region, and gave overlapping regions only one of their shared
    CpGs -- under-counting ``N_meth``/``N_unmeth``/``n_cpgs`` (M7). A range
    join over the half-open intervals is correct for all BED topologies.
    """
    # Inequality (range) join: pos in [start, end). Polars >= 1.7 join_where
    # handles nested/overlapping regions natively (each CpG pairs with every
    # region it falls in). Select only the region columns we need so the
    # shared ``chrom`` column doesn't collide on the join.
    region_cols = regions_chrom.select(["region_id", "start", "end"])
    matched = cpgs.join_where(
        region_cols,
        pl.col("pos") >= pl.col("start"),
        pl.col("pos") < pl.col("end"),
    )
    if len(matched) == 0:
        return None
    return matched


def smooth(
    md: MethylData,
    *,
    method: str = "gaussian",
    bandwidth: int = 1000,
    grid_resolution_bp: int | None = None,
    # BSmooth-only knobs
    ns: int = 70,
    h_bp: int = 1000,
    degree: int = 2,
    min_cpgs_for_smooth: int = 3,
) -> None:
    """Smooth per-sample beta values along the genomic axis.

    Two backends:

    * ``method="gaussian"`` (default) -- coverage-weighted Gaussian kernel
      on a regularised grid. Fast (O(G) per chrom), documented as a
      BSmooth approximation.
    * ``method="bsmooth"`` -- spec-faithful BSmooth (Hansen et al. 2012):
      local weighted-polynomial regression with adaptive bandwidth
      (``max(distance to ns-th CpG, h_bp)``), tricube x coverage weights,
      degree-2 fit by default. Slower than Gaussian but matches the
      Bioconductor ``bsseq``/DSS reference behaviour. Compiled via numba.

    Smoothed values are written to ``md.uns["smooth_path"]`` and are
    accessible from the Parquet store for downstream analyses. Standard
    ``ep.tl.dmc`` / ``ep.tl.dmr`` continue to use raw counts.

    Parameters
    ----------
    md
        MethylData (must have been filtered with ``ep.pp.filter_coverage``).
    method : {"gaussian", "bsmooth"}
        Which smoother to use. Default ``"gaussian"`` preserves prior behaviour.
    bandwidth, grid_resolution_bp
        Gaussian-method options. Ignored when ``method="bsmooth"``.
    ns, h_bp, degree, min_cpgs_for_smooth
        BSmooth-method options. Ignored when ``method="gaussian"``.
        ``ns`` and ``h_bp`` defaults match the Hansen et al. (2012) paper
        (70 CpGs / 1 kb half-window). ``degree=2`` is the canonical
        BSmooth fit; ``degree=1`` is a faster linear fallback.
    """
    if not md._filtered:
        raise ValueError(
            "Run ep.pp.filter_coverage(md) before ep.pp.smooth(md)."
        )

    method = method.lower()
    if method not in ("gaussian", "bsmooth"):
        raise ValueError(
            f"Unknown smoothing method {method!r}. Use 'gaussian' or 'bsmooth'."
        )

    samples = md.obs.get_column("sample_id").to_list()

    if md.analysis_root:
        smooth_path = str(Path(md.analysis_root) / ".cache" / "smoothed")
    else:
        smooth_path = f"{md.store}_smoothed"

    logger.info("Running %s smoothing to %s ...", method, smooth_path)

    if method == "gaussian":
        from .dmr import smooth_methylation_gaussian
        smooth_methylation_gaussian(
            methylstore_path=md.store,
            samples=samples,
            bandwidth=bandwidth,
            grid_resolution_bp=grid_resolution_bp,
            output_path=smooth_path,
        )
        params = {
            "method": "gaussian",
            "bandwidth": bandwidth,
            "grid_resolution_bp": grid_resolution_bp,
        }
    else:  # bsmooth
        from .dmr import smooth_methylation_bsmooth
        smooth_methylation_bsmooth(
            methylstore_path=md.store,
            samples=samples,
            ns=ns,
            h_bp=h_bp,
            degree=degree,
            min_cpgs_for_smooth=min_cpgs_for_smooth,
            output_path=smooth_path,
        )
        params = {
            "method": "bsmooth",
            "ns": ns,
            "h_bp": h_bp,
            "degree": degree,
            "min_cpgs_for_smooth": min_cpgs_for_smooth,
        }

    md.uns["smooth_path"] = smooth_path
    md.uns["smooth_params"] = params
    _append_store_history(md, "smoothed", smooth_path, None)

    logger.info(
        "Smoothing complete (%d samples, method=%s). Results in %s",
        len(samples), method, smooth_path,
    )
