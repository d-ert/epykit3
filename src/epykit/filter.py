"""QC and filtering for Parquet methylation stores.

Notable design points:
  - Blacklist filtering uses one bioframe interval overlap (not a per-region
    .filter() loop).
  - ``filter_sites`` reads each chromosome once: a lightweight coverage-only
    scan computes the genome-wide quantile, then per-chromosome reads do
    filtering + writing.
  - ``intersect_sites`` uses one group_by + count rather than N-1 sequential
    joins.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Optional

import numpy as np
import polars as pl

from . import _cache

logger = logging.getLogger(__name__)

FILTER_MANIFEST_NAME = ".epykit_filter_manifest.json"
NORMALIZE_MANIFEST_NAME = ".epykit_normalize_manifest.json"


def _filter_params_payload(
    min_coverage: int,
    max_coverage_quantile: float,
    blacklist_bed: Optional[str],
) -> dict:
    return {
        "min_coverage": int(min_coverage),
        "max_coverage_quantile": float(max_coverage_quantile),
        "blacklist_bed_sig": (
            _cache.file_signature(Path(blacklist_bed))
            if blacklist_bed
            else None
        ),
    }


def _can_reuse_filtered(
    in_sample_dir: Path, out_sample_dir: Path, params: dict
) -> bool:
    manifest = _cache.load_json(out_sample_dir / FILTER_MANIFEST_NAME)
    if not manifest:
        return False
    if manifest.get("params") != params:
        return False
    if manifest.get("source") != _cache.upstream_sample_signature(in_sample_dir):
        return False
    chroms = manifest.get("chroms")
    if not isinstance(chroms, list):
        return False
    return _cache.sample_is_complete(out_sample_dir, chroms)


# Internal helpers

def _genome_wide_quantile(
    sample_dir: Path,
    quantile: float,
) -> int:
    """Compute a genome-wide coverage quantile by scanning only the coverage
    column across all chromosome Parquet files for one sample.

    Reading a single column is dramatically cheaper than reading the full
    schema, so this first pass is lightweight even on large genomes.
    """
    coverage_series = pl.concat([
        pl.read_parquet(str(part), columns=["coverage"])["coverage"]
        for chrom_dir in sorted(sample_dir.glob("chrom=*"))
        for part in chrom_dir.glob("part-*.parquet")
    ])
    return int(coverage_series.quantile(quantile))


def _apply_blacklist(df: pl.DataFrame, blacklist_bed: str) -> pl.DataFrame:
    """Remove CpG sites that overlap any region in the blacklist BED file.

    Uses bioframe for vectorised interval overlap -- one operation regardless
    of how many regions are in the blacklist.

    PERF-1 fix: replaces the old per-region .filter() loop that built a chain
    of O(n_regions) nested lazy nodes.
    """
    try:
        import bioframe
        import pandas as pd
    except ImportError as exc:
        raise ImportError(
            "bioframe is required for blacklist filtering. "
            "Install it with: pip install bioframe"
        ) from exc

    bl = bioframe.read_table(blacklist_bed, schema="bed3", usecols=[0, 1, 2])

    cpg = pd.DataFrame({
        "chrom": df["chrom"].to_list(),
        "start": df["pos"].to_list(),
        "end":   (df["pos"] + 1).to_list(),  # single-base intervals
    })

    hits = bioframe.overlap(cpg, bl, how="inner", suffixes=("", "_bl"))

    if len(hits) == 0:
        return df

    # Build a set of (chrom, pos) pairs to exclude
    exclude: set[tuple[str, int]] = set(
        zip(hits["chrom"].tolist(), hits["start"].tolist())
    )

    chrom_list = df["chrom"].to_list()
    pos_list = df["pos"].to_list()
    keep = pl.Series(
        [(c, p) not in exclude for c, p in zip(chrom_list, pos_list)]
    )
    return df.filter(keep)


# Public API

def sample_summary(
    methylstore_path: str,
    sample: str,
    output_path: Optional[str] = None,
) -> pl.DataFrame:
    """Compute per-chromosome summary statistics for a single sample.

    Parameters
    ----------
    methylstore_path : str
        Path to the partitioned Parquet methylstore (root containing
        sample=*/chrom=*/)
    sample : str
        Sample identifier (exact match)
    output_path : str, optional
        If provided, write summary as Parquet; default returns to caller

    Returns
    -------
    pl.DataFrame
        Columns: chrom, n_CpGs, mean_coverage, median_coverage,
        global_methylation
    """
    glob_pattern = f"{methylstore_path}/sample={sample}/**/part-*.parquet"
    lf = pl.scan_parquet(glob_pattern)

    stats = (
        lf.group_by("chrom")
        .agg(
            [
                pl.len().alias("n_CpGs"),
                pl.mean("coverage").alias("mean_coverage"),
                pl.median("coverage").alias("median_coverage"),
                (pl.sum("N_meth") / pl.sum("coverage")).alias("global_methylation"),
            ]
        )
        .sort("chrom")
        .collect()
    )

    if output_path:
        stats.write_parquet(output_path)
    return stats


def get_coverage_quantile(
    methylstore_path: str,
    sample: str,
    quantile: float = 0.999,
) -> int:
    """Compute a genome-wide per-sample coverage quantile.

    Reads only the coverage column to minimise I/O.

    Parameters
    ----------
    methylstore_path : str
        Path to the partitioned Parquet methylstore
    sample : str
        Sample identifier
    quantile : float
        Quantile to compute (0.0-1.0); default 0.999

    Returns
    -------
    int
        Coverage value at the specified quantile
    """
    sample_dir = Path(methylstore_path) / f"sample={sample}"
    if not sample_dir.exists():
        raise ValueError(f"Sample directory not found: {sample_dir}")
    return _genome_wide_quantile(sample_dir, quantile)


def filter_sites(
    methylstore_path: str,
    output_dir: str,
    min_coverage: int = 10,
    max_coverage_quantile: float = 0.999,
    blacklist_bed: Optional[str] = None,
    sample: Optional[str] = None,
) -> None:
    """Filter low-quality CpG sites from a Parquet methylstore.

    PERF-2 fix: each sample is now processed chromosome-by-chromosome with a
    lightweight first pass (coverage column only) to determine the genome-wide
    quantile threshold, followed by per-chromosome reads for filtering and
    writing.  The old approach required two full-table scans per sample.

    PERF-1 fix: blacklist filtering uses a single bioframe interval overlap
    instead of a per-region lazy .filter() loop.

    Parameters
    ----------
    methylstore_path : str
        Path to input Parquet methylstore
    output_dir : str
        Path to write filtered Parquet store (mirrors input structure)
    min_coverage : int
        Minimum coverage threshold (default 10)
    max_coverage_quantile : float
        Quantile for the upper coverage cap (default 0.999)
    blacklist_bed : str, optional
        Path to BED file with regions to exclude (chrom, start, end columns)
    sample : str, optional
        If provided, filter only this sample; else all samples in the store

    Returns
    -------
    None
        Writes filtered Parquet store to output_dir
    """
    methylstore_path = Path(methylstore_path)
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    sample_dirs = list(methylstore_path.glob("sample=*"))
    if not sample_dirs:
        raise ValueError(f"No sample=* directories found in {methylstore_path}")

    if sample:
        sample_path = methylstore_path / f"sample={sample}"
        if not sample_path.exists():
            raise ValueError(f"Sample '{sample}' not found in {methylstore_path}")
        samples_to_filter = [sample]
    else:
        samples_to_filter = [d.name.removeprefix("sample=") for d in sorted(sample_dirs)]

    params = _filter_params_payload(min_coverage, max_coverage_quantile, blacklist_bed)

    for samp in samples_to_filter:
        sample_dir = methylstore_path / f"sample={samp}"
        out_sample_dir = output_dir_path / f"sample={samp}"

        if _can_reuse_filtered(sample_dir, out_sample_dir, params):
            logger.info("Filtering sample %s: cached", samp)
            continue

        logger.info("Filtering sample %s...", samp)

        # Wipe stale output before re-running so the manifest's chrom list
        # matches what's on disk.
        if out_sample_dir.exists():
            shutil.rmtree(out_sample_dir)
        out_sample_dir.mkdir(parents=True, exist_ok=True)

        # --- PERF-2: pass 1 -- coverage column only, genome-wide quantile ---
        max_cov = _genome_wide_quantile(sample_dir, max_coverage_quantile)
        logger.info("  Max coverage quantile (%s): %s", max_coverage_quantile, max_cov)

        chrom_dirs = sorted(sample_dir.glob("chrom=*"))
        if not chrom_dirs:
            logger.warning("  no chrom=* dirs found for %s; skipping", samp)
            continue

        # --- PERF-2: pass 2 -- per-chromosome read -> filter -> write ---
        for chrom_dir in chrom_dirs:
            chrom = chrom_dir.name.removeprefix("chrom=")

            parts = list(chrom_dir.glob("part-*.parquet"))
            if not parts:
                continue

            chrom_df = pl.concat([pl.read_parquet(str(p)) for p in parts])

            # Coverage filter
            chrom_df = chrom_df.filter(
                (pl.col("coverage") >= min_coverage)
                & (pl.col("coverage") <= max_cov)
            )

            if len(chrom_df) == 0:
                continue

            # --- PERF-1: blacklist via bioframe interval overlap ---
            if blacklist_bed:
                chrom_df = _apply_blacklist(chrom_df, blacklist_bed)

            if len(chrom_df) == 0:
                continue

            out_chrom_dir = out_sample_dir / f"chrom={chrom}"
            out_chrom_dir.mkdir(parents=True, exist_ok=True)
            chrom_df.write_parquet(
                str(out_chrom_dir / "part-0.parquet"),
                compression="zstd",
            )

        n_out = sum(
            1
            for f in out_sample_dir.rglob("part-*.parquet")
        )
        logger.info("  Written %d chromosome file(s) for %s", n_out, samp)

        _cache.write_json(
            out_sample_dir / FILTER_MANIFEST_NAME,
            {
                "sample_name": samp,
                "source": _cache.upstream_sample_signature(sample_dir),
                "params": params,
                "chroms": _cache.expected_chrom_dirs(out_sample_dir),
            },
        )

    logger.info("Filtered Parquet store written to %s", output_dir_path)


def normalize_coverage_store(
    methylstore_path: str,
    output_dir: str,
    method: str = "median",
) -> dict[str, float]:
    """Per-sample coverage normalisation.

    For each sample, compute a single scalar factor ``s_i`` so that the
    chosen central statistic of coverage matches a common target across
    the cohort:

        method="median": target = median(per-sample medians)
                         s_i    = target / median(cov_i)
        method="mean":   target =   mean(per-sample means)
                         s_i    = target /   mean(cov_i)

    Each row's ``N_meth`` and ``N_unmeth (= coverage - N_meth)`` are then
    scaled by ``s_i`` and rounded to int. ``coverage`` is rebuilt as
    ``N_meth + N_unmeth`` so the equality holds exactly after rounding --
    downstream Polars joins rely on this strict invariant.

    Coverage normalisation should run between ``filter_coverage`` and
    ``tl.dmr(method='tile')`` to prevent deeper-sequenced samples from
    dominating pooled-count tile tests.

    Parameters
    ----------
    methylstore_path : str
        Path to input partitioned Parquet methylstore.
    output_dir : str
        Path to write normalised store (mirrors input partitioning).
    method : {"median", "mean"}
        Central statistic to align.

    Returns
    -------
    dict[str, float]
        Mapping ``{sample_id: scaling_factor}``.
    """
    if method not in ("median", "mean"):
        raise ValueError(f"method must be 'median' or 'mean', got {method!r}")

    src = Path(methylstore_path)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    samples = sorted(
        d.name.removeprefix("sample=") for d in src.glob("sample=*")
    )
    if not samples:
        raise ValueError(f"No sample=* directories found in {src}")

    # --- Cohort cache check ----------------------------------------------
    cohort_sig = {
        samp: _cache.upstream_sample_signature(src / f"sample={samp}")
        for samp in samples
    }
    cohort_params = {"method": method}
    cohort_manifest_path = out / NORMALIZE_MANIFEST_NAME
    cached = _cache.load_json(cohort_manifest_path)
    if (
        cached
        and cached.get("params") == cohort_params
        and cached.get("source") == cohort_sig
        and set(cached.get("factors", {}).keys()) == set(samples)
    ):
        chroms_by_sample = cached.get("chroms_by_sample", {})
        if all(
            _cache.sample_is_complete(out / f"sample={samp}", chroms_by_sample.get(samp, []))
            for samp in samples
        ):
            factors = {s: float(v) for s, v in cached["factors"].items()}
            logger.info("Normalised store cached at %s (method=%s, target=%.2f)",
                        out, method, float(cached.get("target", 0.0)))
            for samp in samples:
                logger.info("  %s: factor=%.4f (cached)", samp, factors[samp])
            return factors

    # --- Pass 1: per-sample central coverage ------------------------------
    summaries: dict[str, float] = {}
    for samp in samples:
        sample_dir = src / f"sample={samp}"
        parts = [
            part
            for chrom_dir in sorted(sample_dir.glob("chrom=*"))
            for part in chrom_dir.glob("part-*.parquet")
        ]
        if not parts:
            raise ValueError(f"Sample {samp!r} has no parquet parts")
        coverage_series = pl.concat([
            pl.read_parquet(str(p), columns=["coverage"])["coverage"]
            for p in parts
        ])
        if len(coverage_series) == 0:
            raise ValueError(f"Sample {samp!r} has zero rows after read")
        summaries[samp] = float(
            coverage_series.median() if method == "median"
            else coverage_series.mean()
        )

    summary_values = np.array(list(summaries.values()), dtype=np.float64)
    target = float(
        np.median(summary_values) if method == "median" else summary_values.mean()
    )

    factors: dict[str, float] = {
        samp: (target / s if s > 0 else 1.0) for samp, s in summaries.items()
    }

    logger.info("Coverage normalisation (method=%s, target=%.2f):", method, target)
    for samp in samples:
        logger.info(
            "  %s: %s_cov=%.2f, factor=%.4f",
            samp, method, summaries[samp], factors[samp],
        )

    # --- Pass 2: scale and write ------------------------------------------
    # Wipe any stale output so the cohort manifest's chrom listing matches
    # what's actually on disk after this run.
    for samp in samples:
        out_samp = out / f"sample={samp}"
        if out_samp.exists():
            shutil.rmtree(out_samp)

    for samp in samples:
        s = factors[samp]
        sample_dir = src / f"sample={samp}"
        for chrom_dir in sorted(sample_dir.glob("chrom=*")):
            chrom = chrom_dir.name.removeprefix("chrom=")
            parts = list(chrom_dir.glob("part-*.parquet"))
            if not parts:
                continue

            chrom_df = pl.concat([pl.read_parquet(str(p)) for p in parts])

            # Capture original N_unmeth before either column is overwritten,
            # so the two scaled counts remain consistent and coverage is
            # rebuilt from them exactly.
            scaled = (
                chrom_df
                .with_columns(
                    (pl.col("coverage") - pl.col("N_meth")).alias("_N_unmeth_orig")
                )
                .with_columns([
                    (pl.col("N_meth").cast(pl.Float64) * s)
                        .round().cast(pl.Int32).alias("N_meth"),
                    (pl.col("_N_unmeth_orig").cast(pl.Float64) * s)
                        .round().cast(pl.Int32).alias("_N_unmeth"),
                ])
                .with_columns(
                    (pl.col("N_meth") + pl.col("_N_unmeth"))
                    .cast(pl.Int32).alias("coverage")
                )
                .drop(["_N_unmeth_orig", "_N_unmeth"])
            )

            out_chrom_dir = out / f"sample={samp}" / f"chrom={chrom}"
            out_chrom_dir.mkdir(parents=True, exist_ok=True)
            scaled.write_parquet(
                str(out_chrom_dir / "part-0.parquet"),
                compression="zstd",
            )

    _cache.write_json(
        cohort_manifest_path,
        {
            "params": cohort_params,
            "source": cohort_sig,
            "target": target,
            "factors": factors,
            "chroms_by_sample": {
                samp: _cache.expected_chrom_dirs(out / f"sample={samp}")
                for samp in samples
            },
        },
    )

    logger.info("Normalised Parquet store written to %s", out)
    return factors


def intersect_sites(
    methylstore_path: str,
    samples: list[str],
    output_path: Optional[str] = None,
) -> pl.DataFrame:
    """Find CpG sites present in all specified samples.

    PERF-3 fix: replaced the N-1 sequential inner-join chain with a single
    lazy scan -> group_by -> count approach.  All samples are processed in one
    pass; a site is retained only when its per-sample count equals len(samples).

    Parameters
    ----------
    methylstore_path : str
        Path to Parquet methylstore
    samples : list[str]
        Sample identifiers to intersect
    output_path : str, optional
        If provided, write (chrom, pos, strand) intersection to Parquet

    Returns
    -------
    pl.DataFrame
        Columns: chrom, pos, strand -- sites present in every listed sample
    """
    if not samples:
        raise ValueError("Must provide at least one sample")

    lf = pl.scan_parquet(
        f"{methylstore_path}/sample=*/chrom=*/part-*.parquet"
    )

    intersect = (
        lf.filter(pl.col("sample").is_in(samples))
        .select(["chrom", "pos", "strand", "sample"])
        .unique()                                   # one row per (site, sample)
        .group_by(["chrom", "pos", "strand"])
        .agg(pl.len().alias("n_samples"))
        .filter(pl.col("n_samples") == len(samples))
        .drop("n_samples")
        .sort(["chrom", "pos"])
        .collect()
    )

    if output_path:
        intersect.write_parquet(output_path)

    return intersect


def load_chromosome_data(
    methylstore_path: str,
    chrom: str,
    samples: list[str],
    site_intersect: Optional[pl.DataFrame] = None,
) -> pl.DataFrame:
    """Load all data for a specific chromosome and set of samples.

    Parameters
    ----------
    methylstore_path : str
        Path to Parquet methylstore
    chrom : str
        Chromosome name (e.g. "chr1")
    samples : list[str]
        Sample identifiers to load
    site_intersect : pl.DataFrame, optional
        DataFrame with columns (chrom, pos, strand). When provided, only sites
        in this set are returned (unite behaviour).

    Returns
    -------
    pl.DataFrame
        Columns: chrom, pos, strand, N_meth, N_unmeth, coverage, sample
    """
    glob_pattern = (
        f"{methylstore_path}/sample=*/chrom={chrom}/part-*.parquet"
    )
    lf = pl.scan_parquet(glob_pattern)
    lf = lf.filter(pl.col("sample").is_in(samples))

    if site_intersect is not None:
        site_intersect_chrom = site_intersect.filter(pl.col("chrom") == chrom)
        if len(site_intersect_chrom) == 0:
            return lf.filter(pl.lit(False)).collect()
        lf = lf.join(
            site_intersect_chrom.lazy(),
            on=["chrom", "pos", "strand"],
            how="inner",
        )

    return lf.collect()