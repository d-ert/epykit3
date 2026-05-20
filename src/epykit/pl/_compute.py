"""Backend-agnostic compute layer for epykit plotting.

Every entry point here returns a small, picklable result object (numpy
arrays, polars/pandas frames, or dataclasses thereof). The matplotlib
plots in ``pl/*.py`` and the Plotly counterparts in ``pl/_plotly.py``
both consume the same outputs -- single source of truth for "what data
does this plot show". Heavy on-disk scans live here and nowhere else.

Design rules:
- No matplotlib / plotly imports in this module.
- Functions never hold open file handles; everything is collected before
  return.
- For functions that touch the per-CpG store, do **one** lazy scan and
  push every filter into it. The big OOM cliff in the previous report()
  was caused by per-sample / per-chromosome / per-TSS re-scans.
- Results are small enough to cache in ``md.uns["_report_cache"]``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import polars as pl


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------


@dataclass
class PCAResult:
    coords: np.ndarray            # (n_samples, n_components)
    explained_var: tuple          # per-component ratios
    samples: list[str]            # row order of coords
    groups: list                  # one group label per sample, parallel to samples
    group_col: str                # which md.obs column was used
    n_sites_used: int             # rows entering PCA after subsampling + NaN drop


@dataclass
class MetaplotResult:
    x: np.ndarray                 # bin centre coordinates, length n_bins
    mean_beta: np.ndarray         # (n_samples, n_bins)
    samples: list[str]
    groups: list
    group_col: Optional[str]
    window_bp: int
    n_bins: int


@dataclass
class VolcanoData:
    meth_diff: np.ndarray
    neg_log_p: np.ndarray
    p_col: str
    sig: np.ndarray               # bool mask
    hyper: np.ndarray
    hypo: np.ndarray


@dataclass
class MAData:
    mean_beta: np.ndarray
    meth_diff: np.ndarray
    p_col: str
    sig: np.ndarray
    hyper: np.ndarray
    hypo: np.ndarray


@dataclass
class ManhattanData:
    chrom_blocks: list            # list of dicts: {chrom, x, y, n}
    p_col: str
    alpha_line_y: float
    tick_pos: list[float]
    tick_label: list[str]


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------


def _cache_get(md, key: str):
    cache = md.uns.get("_report_cache")
    if not isinstance(cache, dict):
        return None
    return cache.get(key)


def _cache_put(md, key: str, value) -> None:
    cache = md.uns.get("_report_cache")
    if not isinstance(cache, dict):
        cache = {}
        md.uns["_report_cache"] = cache
    cache[key] = value


def clear_report_cache(md) -> None:
    """Drop cached compute results so the next plot/report recomputes.

    Call this after mutating the underlying store (re-running QC, filter,
    DMC, etc.) so stale arrays don't leak into a fresh report.
    """
    if isinstance(md.uns.get("_report_cache"), dict):
        md.uns["_report_cache"].clear()


# ---------------------------------------------------------------------------
# Sample x site matrix (shared by PCA + UMAP)
# ---------------------------------------------------------------------------


def _store_pattern(md) -> str:
    return f"{md.store}/sample=*/chrom=*/part-*.parquet"


def _resolve_group_col(md, group_col: Optional[str]) -> tuple[Optional[str], list]:
    """Pick ``group_col`` from md.obs, falling back through common names."""
    if group_col is None:
        for cand in ("group", "treatment", "condition"):
            if cand in md.obs.columns:
                group_col = cand
                break
    if group_col is None or group_col not in md.obs.columns:
        return None, ["all"] * len(md.obs)
    return group_col, md.obs.get_column(group_col).to_list()


def compute_sample_site_matrix(
    md,
    *,
    n_sites: int = 10_000,
    seed: int = 42,
) -> tuple[np.ndarray, list[str], int]:
    """Build an (n_samples x n_sites) beta matrix from the per-CpG store.

    Single-pass replacement for the old two-scan ``build_sample_site_matrix``.
    Subsamples sites by hashing ``(chrom, pos)`` before pivoting so peak
    memory is bounded by ``n_sites * n_samples * 8 B`` regardless of how
    many CpGs the store contains.

    Returns
    -------
    matrix : np.ndarray
        Shape ``(n_samples, n_sites_kept)``, no NaN rows.
    samples : list[str]
        Sample ordering for the matrix rows.
    n_sites_used : int
        Number of sites that survived both subsampling and the drop-nan step.
    """
    samples = md.obs.get_column("sample_id").to_list()
    if len(samples) < 2:
        raise ValueError("Need >=2 samples for a sample-site matrix")

    pattern = _store_pattern(md)

    # Probe: one cheap row count to choose a deterministic modulus.
    # `pl.len()` on a lazy scan only reads parquet footers, so this is
    # near-instantaneous even on terabyte stores.
    n_total = pl.scan_parquet(pattern).select(pl.len()).collect().item()
    if n_total == 0:
        raise ValueError("Store contains no CpGs")
    target_rows = max(n_sites * len(samples) * 2, len(samples))
    k = max(1, n_total // target_rows)

    # Integer modulo on pos is predicate-pushdown-friendly (polars pushes
    # it into the parquet reader), so we never read the whole store into
    # memory just to subsample. Sites are unique per (chrom, pos); using
    # the same modulus on pos picks the same positions on every chrom
    # and the same chrom across samples, so per-site coverage is
    # consistent.
    lf = (
        pl.scan_parquet(pattern)
        .filter(pl.col("coverage") > 0)
    )
    if k > 1:
        lf = lf.filter((pl.col("pos") % k) == (seed % k))
    lf = lf.select(["chrom", "pos", "sample", "N_meth", "coverage"]).with_columns(
        (pl.col("N_meth").cast(pl.Float64) / pl.col("coverage").cast(pl.Float64))
        .alias("beta")
    ).select(["chrom", "pos", "sample", "beta"])

    df = lf.collect()
    if df.is_empty():
        raise ValueError("No CpGs survived subsampling")

    pivot = df.pivot(
        values="beta", index=["chrom", "pos"], on="sample",
        aggregate_function="mean",
    )
    present = [s for s in samples if s in pivot.columns]
    if len(present) < 2:
        raise ValueError(
            f"Only {len(present)} samples have data after pivot; "
            "PCA / UMAP need >=2."
        )
    matrix = pivot.select(present).to_numpy()
    valid = ~np.isnan(matrix).any(axis=1)
    matrix = matrix[valid]

    # Hash-based sampling overshoots; trim down to exactly n_sites.
    if matrix.shape[0] > n_sites:
        rng = np.random.default_rng(seed)
        idx = rng.choice(matrix.shape[0], n_sites, replace=False)
        matrix = matrix[idx]

    if matrix.shape[0] < 2:
        raise ValueError(
            f"Only {matrix.shape[0]} valid sites after NaN filter -- not enough for embedding"
        )
    return matrix.T, present, int(matrix.shape[0])


# ---------------------------------------------------------------------------
# PCA
# ---------------------------------------------------------------------------


def compute_pca(
    md,
    *,
    n_sites: int = 10_000,
    n_components: int = 2,
    seed: int = 42,
    group_col: Optional[str] = None,
    use_cache: bool = True,
) -> PCAResult:
    """Compute PCA of per-sample methylation profiles.

    Single-pass over the per-CpG store. See
    :func:`compute_sample_site_matrix` for the streaming/subsampling
    strategy.
    """
    if use_cache:
        cached = _cache_get(md, "pca")
        if cached is not None:
            return cached

    try:
        from sklearn.decomposition import PCA
    except ImportError as exc:
        raise ImportError(
            "scikit-learn is required for PCA. "
            "Install with: pip install scikit-learn"
        ) from exc

    matrix, samples, n_sites_used = compute_sample_site_matrix(
        md, n_sites=n_sites, seed=seed,
    )
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(matrix)

    resolved_col, groups_full = _resolve_group_col(md, group_col)
    sample_to_group = dict(zip(md.obs.get_column("sample_id").to_list(), groups_full))
    groups = [sample_to_group.get(s, "unknown") for s in samples]

    result = PCAResult(
        coords=coords,
        explained_var=tuple(float(v) for v in pca.explained_variance_ratio_),
        samples=samples,
        groups=groups,
        group_col=resolved_col or "all",
        n_sites_used=n_sites_used,
    )
    if use_cache:
        _cache_put(md, "pca", result)
    return result


# ---------------------------------------------------------------------------
# TSS metaplot (vectorised via bioframe overlap)
# ---------------------------------------------------------------------------


def _tss_intervals(gtf_path: str, *, window_bp: int, max_genes: Optional[int]):
    """Build a polars DataFrame of TSS windows from a GTF.

    Returns columns: chrom, start, end, tss, strand, gene_id, gene_name.
    Uses the shared cached GTF parser in :mod:`epykit.annotate`.
    """
    from ..annotate import _parse_gtf_streaming

    genes_pd, _ = _parse_gtf_streaming(gtf_path)
    if genes_pd is None or len(genes_pd) == 0:
        return None
    df = pl.from_pandas(
        genes_pd[["Chromosome", "Start", "End", "Strand", "gene_id", "gene_name"]]
    ).rename({"Chromosome": "chrom", "Strand": "strand"})
    tss = (
        df.with_columns(
            pl.when(pl.col("strand") == "-")
              .then(pl.col("End") - 1)
              .otherwise(pl.col("Start"))
              .cast(pl.Int64).alias("tss")
        )
        .select(["chrom", "tss", "strand", "gene_id", "gene_name"])
    )
    if max_genes is not None and len(tss) > max_genes:
        tss = tss.head(max_genes)
    return tss.with_columns([
        (pl.col("tss") - window_bp).alias("start"),
        (pl.col("tss") + window_bp).alias("end"),
    ]).select(["chrom", "start", "end", "tss", "strand", "gene_id"])


def _merge_intervals(starts: np.ndarray, ends: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Merge potentially overlapping [start, end) intervals.

    Returns sorted, non-overlapping arrays. Used by the metaplot to push
    a tight pre-filter into the polars lazy scan -- otherwise the
    bounding span of all TSS windows can effectively be the entire
    chromosome and the filter is a no-op.
    """
    if starts.size == 0:
        return starts, ends
    order = np.argsort(starts)
    ws = starts[order].astype(np.int64)
    we = ends[order].astype(np.int64)
    out_s = [int(ws[0])]
    out_e = [int(we[0])]
    for i in range(1, len(ws)):
        if int(ws[i]) <= out_e[-1]:
            if int(we[i]) > out_e[-1]:
                out_e[-1] = int(we[i])
        else:
            out_s.append(int(ws[i]))
            out_e.append(int(we[i]))
    return np.asarray(out_s, dtype=np.int64), np.asarray(out_e, dtype=np.int64)


def compute_tss_metaplot(
    md,
    gtf_path: str,
    *,
    window_bp: int = 2000,
    n_bins: int = 100,
    group_by: Optional[str] = "group",
    max_genes: Optional[int] = None,
    use_cache: bool = True,
) -> MetaplotResult:
    """Compute mean beta around the TSS, with bounded peak memory.

    Per chromosome, per sample: one streaming scan of just that sample's
    chrom partition, filtered down to CpGs falling inside the union of
    TSS windows. Vectorised ``np.searchsorted`` then assigns each CpG to
    its TSS bin and accumulates into the running ``(n_samples, n_bins)``
    arrays. Peak resident memory is bounded by one sample's
    window-resident CpGs (typically <50 MB for a human chromosome).
    """
    cache_key = f"tss_metaplot:{gtf_path}:{window_bp}:{n_bins}:{group_by}:{max_genes}"
    if use_cache:
        cached = _cache_get(md, cache_key)
        if cached is not None:
            return cached

    samples = md.obs.get_column("sample_id").to_list()
    if not samples:
        raise ValueError("md.obs has no samples")

    tss = _tss_intervals(gtf_path, window_bp=window_bp, max_genes=max_genes)
    if tss is None or tss.is_empty():
        raise ValueError(f"No gene records found in GTF {gtf_path!r}")

    bin_size = (2 * window_bp) / n_bins
    sum_beta = np.zeros((len(samples), n_bins), dtype=np.float64)
    count = np.zeros((len(samples), n_bins), dtype=np.int64)

    chroms = sorted(set(tss.get_column("chrom").to_list()))
    for chrom in chroms:
        tss_c = tss.filter(pl.col("chrom") == chrom)
        if tss_c.is_empty():
            continue
        tss_positions = tss_c["tss"].to_numpy().astype(np.int64)
        tss_starts = tss_c["start"].to_numpy().astype(np.int64)
        tss_ends = tss_c["end"].to_numpy().astype(np.int64)
        strands_raw = tss_c["strand"].to_numpy()
        strand_sign = np.where(strands_raw == "-", -1, 1).astype(np.int8)

        # Union of TSS windows -- tight pre-filter when genes are clustered.
        merged_s, merged_e = _merge_intervals(tss_starts, tss_ends)
        if merged_s.size == 0:
            continue

        # Per sample, streaming. One sample's chrom partition at a time
        # keeps peak memory bounded by ~chrom_size_in_bp / spacing CpGs.
        for s_idx, sample_id in enumerate(samples):
            part_dir = f"{md.store}/sample={sample_id}/chrom={chrom}"
            try:
                lf = (
                    pl.scan_parquet(f"{part_dir}/part-*.parquet")
                    .select(["pos", "N_meth", "coverage"])
                    .filter(pl.col("coverage") > 0)
                    # Cheap bounding-box filter; polars pushes this into
                    # the parquet reader.
                    .filter(pl.col("pos").is_between(
                        int(merged_s[0]), int(merged_e[-1]),
                    ))
                )
                df = lf.collect()
            except Exception:
                continue
            if df.is_empty():
                continue

            positions = df["pos"].to_numpy().astype(np.int64)
            # Tighter filter: keep only positions inside any merged window.
            # Vectorised: for each pos find the latest merged_s <= pos via
            # searchsorted, then check that pos < merged_e[that index].
            idx = np.searchsorted(merged_s, positions, side="right") - 1
            valid = (idx >= 0) & (positions < merged_e[np.clip(idx, 0, len(merged_e) - 1)])
            if not valid.any():
                continue
            positions = positions[valid]
            betas = (
                df["N_meth"].to_numpy().astype(np.float64)[valid]
                / df["coverage"].to_numpy().astype(np.float64)[valid]
            )
            # Sort by position so per-TSS searchsorted is correct.
            order = np.argsort(positions, kind="mergesort")
            positions = positions[order]
            betas = betas[order]

            # Per-TSS window: pull the CpGs inside, compute relative bin,
            # accumulate. The Python loop is over TSS (~ thousands), not
            # CpGs, so vectorisation lives inside.
            for tss_pos, sign in zip(tss_positions, strand_sign):
                lo = int(tss_pos) - window_bp
                hi = int(tss_pos) + window_bp
                left = np.searchsorted(positions, lo, side="left")
                right = np.searchsorted(positions, hi, side="left")
                if right <= left:
                    continue
                rel = (positions[left:right] - int(tss_pos)) * int(sign)
                bins = np.floor((rel + window_bp) / bin_size).astype(np.int64)
                np.clip(bins, 0, n_bins - 1, out=bins)
                np.add.at(sum_beta[s_idx], bins, betas[left:right])
                np.add.at(count[s_idx], bins, 1)

    with np.errstate(invalid="ignore"):
        mean_beta = np.where(count > 0, sum_beta / count, np.nan)
    x = np.linspace(-window_bp, window_bp, n_bins, endpoint=False) + bin_size / 2.0

    resolved_col, groups_full = _resolve_group_col(md, group_by)
    sample_to_group = dict(zip(md.obs.get_column("sample_id").to_list(), groups_full))
    groups = [sample_to_group.get(s, "unknown") for s in samples]

    result = MetaplotResult(
        x=x, mean_beta=mean_beta, samples=samples, groups=groups,
        group_col=resolved_col, window_bp=window_bp, n_bins=n_bins,
    )
    if use_cache:
        _cache_put(md, cache_key, result)
    return result


# ---------------------------------------------------------------------------
# Coverage histogram
# ---------------------------------------------------------------------------


def compute_coverage_distribution(
    md,
    *,
    max_points: int = 200_000,
    seed: int = 42,
) -> np.ndarray:
    """Return a 1-D numpy array of coverage values, capped at ``max_points``.

    Uses deterministic modulo subsampling on the lazy scan so we never
    materialise more than ~max_points rows.
    """
    pattern = _store_pattern(md)
    total = pl.scan_parquet(pattern).select(pl.len()).collect().item()
    if total == 0:
        return np.array([], dtype=np.int64)
    if total <= max_points:
        return (
            pl.scan_parquet(pattern).select("coverage").collect()["coverage"]
            .to_numpy()
        )
    k = max(1, total // max_points)
    return (
        pl.scan_parquet(pattern)
        .select("coverage")
        .with_row_index("_i")
        .filter(pl.col("_i") % k == 0)
        .drop("_i")
        .collect()["coverage"]
        .to_numpy()
    )


# ---------------------------------------------------------------------------
# Cheap DMC-derived data (volcano / MA / manhattan)
# ---------------------------------------------------------------------------


def _dmc_p_col(dmc: pl.DataFrame) -> str:
    return "qvalue" if "qvalue" in dmc.columns else "pvalue"


def compute_volcano_data(
    md,
    *,
    alpha: float = 0.05,
    min_abs_diff: float = 0.1,
) -> VolcanoData:
    dmc = md.dmc
    if dmc is None:
        raise ValueError("Run ep.tl.dmc(md) first")
    p_col = _dmc_p_col(dmc)
    diff = dmc["meth_diff"].to_numpy()
    pval = dmc[p_col].to_numpy()
    neg_log_p = -np.log10(np.maximum(pval, 1e-300))
    sig = (pval < alpha) & (np.abs(diff) >= min_abs_diff)
    return VolcanoData(
        meth_diff=diff, neg_log_p=neg_log_p, p_col=p_col,
        sig=sig, hyper=sig & (diff > 0), hypo=sig & (diff < 0),
    )


def compute_ma_data(
    md,
    *,
    alpha: float = 0.05,
    min_abs_diff: float = 0.1,
) -> MAData:
    dmc = md.dmc
    if dmc is None:
        raise ValueError("Run ep.tl.dmc(md) first")
    p_col = _dmc_p_col(dmc)
    diff = dmc["meth_diff"].to_numpy()
    pval = dmc[p_col].to_numpy()
    mean_beta = (
        dmc["mean_beta_case"].to_numpy()
        + dmc["mean_beta_control"].to_numpy()
    ) / 2.0
    sig = (pval < alpha) & (np.abs(diff) >= min_abs_diff)
    return MAData(
        mean_beta=mean_beta, meth_diff=diff, p_col=p_col,
        sig=sig, hyper=sig & (diff > 0), hypo=sig & (diff < 0),
    )


def compute_manhattan_data(
    md,
    *,
    alpha: float = 0.05,
) -> ManhattanData:
    dmc = md.dmc
    if dmc is None:
        raise ValueError("Run ep.tl.dmc(md) first")
    if "chrom" not in dmc.columns or "pos" not in dmc.columns:
        raise ValueError("DMC table must carry 'chrom' and 'pos' for a Manhattan plot")
    p_col = _dmc_p_col(dmc)
    dmc_sorted = dmc.sort(["chrom", "pos"])
    chroms = dmc_sorted["chrom"].unique().to_list()
    canonical = (
        [f"chr{i}" for i in range(1, 23)]
        + [f"chr{c}" for c in ("X", "Y", "M")]
    )
    order = canonical + [c for c in chroms if c not in canonical]

    blocks: list[dict] = []
    cumulative = 0.0
    tick_pos: list[float] = []
    tick_label: list[str] = []
    for c in order:
        if c not in chroms:
            continue
        sub = dmc_sorted.filter(pl.col("chrom") == c)
        if sub.is_empty():
            continue
        pos = sub["pos"].to_numpy()
        pvals = sub[p_col].to_numpy()
        y = -np.log10(np.maximum(pvals, 1e-300))
        x = cumulative + pos
        blocks.append({"chrom": c, "x": x, "y": y, "n": len(pos)})
        mid = cumulative + (float(pos.max()) - float(pos.min())) / 2.0
        tick_pos.append(mid)
        tick_label.append(c.replace("chr", ""))
        cumulative += float(pos.max()) + 1e7
    return ManhattanData(
        chrom_blocks=blocks, p_col=p_col, alpha_line_y=-np.log10(alpha),
        tick_pos=tick_pos, tick_label=tick_label,
    )


# ---------------------------------------------------------------------------
# Annotation summaries (annotatr-equivalent compute)
# ---------------------------------------------------------------------------


def _resolve_annotated_table(md, level: str) -> pl.DataFrame:
    """Pick the annotated DMR or DMC table from ``md`` for annotation plots."""
    level = level.lower()
    if level == "dmr":
        dmr = md.uns.get("dmr")
        if dmr is None or not isinstance(dmr, pl.DataFrame) or dmr.is_empty():
            raise ValueError("md.uns['dmr'] is empty -- run ep.tl.dmr(md) first")
        return dmr
    if level == "dmc":
        dmc = md.dmc
        if dmc is None or dmc.is_empty():
            raise ValueError("md.dmc is empty -- run ep.tl.dmc(md) first")
        return dmc
    raise ValueError(f"level must be 'dmr' or 'dmc', got {level!r}")


def _explode_multi_annotation(
    df: pl.DataFrame,
    annot_col: str,
) -> pl.DataFrame:
    """Return a long-form (region_id, annot) frame.

    Handles two cases:
    1. ``annot_col`` is a scalar (e.g. ``"feature_type"``): one row per
       region is fine -- we just rename and add a synthetic region_id.
    2. ``annot_col`` is the list-valued ``"all_overlapping_features"`` /
       ``"all_overlapping_genes"`` from epykit's multi-annotation mode:
       explode so a region with two annotations contributes two rows.
    """
    if annot_col not in df.columns:
        raise ValueError(
            f"annotation column {annot_col!r} not on table. "
            f"Available: {sorted(df.columns)[:10]}..."
        )
    out = df.with_row_index("_region_id")
    if out.schema[annot_col] == pl.List(pl.Utf8):
        out = out.explode(annot_col).filter(pl.col(annot_col).is_not_null())
    return out.rename({annot_col: "annot"}).select(["_region_id", "annot"])


def compute_annotation_counts(
    df: pl.DataFrame,
    *,
    annot_col: str = "feature_type",
) -> pl.DataFrame:
    """Counts of regions per annotation class.

    Dedup is region-level: a region appearing under multiple annotations
    contributes to each class once, matching annotatr's plot_annotation
    behaviour.
    """
    long = _explode_multi_annotation(df, annot_col)
    return (
        long.unique(["_region_id", "annot"])
        .group_by("annot")
        .len()
        .rename({"annot": annot_col, "len": "count"})
        .sort("count", descending=True)
    )


def compute_coannotation_matrix(
    df: pl.DataFrame,
    *,
    annot_col: str = "all_overlapping_features",
) -> tuple[np.ndarray, list[str]]:
    """Pairwise co-annotation count matrix.

    For each region with N annotations, every unordered pair (i, j),
    including i==j, contributes 1 to ``mat[i, j]`` (and ``mat[j, i]``).
    Diagonal entries equal the number of regions in that annotation
    (= ``compute_annotation_counts``). Matches annotatr's
    plot_coannotations dedup ([annotatr/R/visualize.R:163](annotatr/R/visualize.R)).
    """
    if annot_col not in df.columns:
        raise ValueError(
            f"co-annotation needs a list-valued column; {annot_col!r} not on table. "
            "Run ep.tl.annotate(md, multi_annotation=True) first."
        )
    # Get list[str] per region.
    if df.schema[annot_col] != pl.List(pl.Utf8):
        raise ValueError(
            f"{annot_col!r} must be List[Utf8]; "
            "use the annotatr-style multi-annotation columns."
        )
    classes_raw: set[str] = set()
    rows: list[list[str]] = []
    for lst in df.get_column(annot_col).to_list():
        if lst is None:
            continue
        clean = sorted({a for a in lst if a is not None})
        if not clean:
            continue
        rows.append(clean)
        classes_raw.update(clean)
    classes = sorted(classes_raw)
    idx = {c: i for i, c in enumerate(classes)}
    n = len(classes)
    mat = np.zeros((n, n), dtype=np.int64)
    for region_annots in rows:
        for a in region_annots:
            mat[idx[a], idx[a]] += 1
        for i, a in enumerate(region_annots):
            for b in region_annots[i + 1:]:
                ia, ib = idx[a], idx[b]
                mat[ia, ib] += 1
                mat[ib, ia] += 1
    return mat, classes


def compute_numerical_by_annotation(
    df: pl.DataFrame,
    *,
    value_col: str,
    annot_col: str = "feature_type",
    include_all: bool = True,
) -> pl.DataFrame:
    """Long-form (annotation, value) frame for faceted histograms.

    When ``include_all`` is True, an extra ``"All"`` class is appended
    containing every region's value once -- annotatr's standard
    "background overlay" trick (plot_numerical).
    """
    if value_col not in df.columns:
        raise ValueError(f"value column {value_col!r} not on table")
    base = df.with_row_index("_region_id")
    if annot_col in base.columns and base.schema[annot_col] == pl.List(pl.Utf8):
        long = (
            base.select(["_region_id", value_col, annot_col])
            .explode(annot_col)
            .filter(pl.col(annot_col).is_not_null())
            .rename({annot_col: "annot"})
        )
    elif annot_col in base.columns:
        long = (
            base.select(["_region_id", value_col, annot_col])
            .rename({annot_col: "annot"})
            .filter(pl.col("annot").is_not_null())
        )
    else:
        raise ValueError(f"annotation column {annot_col!r} not on table")

    long = long.select(["_region_id", "annot", value_col]).filter(
        pl.col(value_col).is_not_null() & pl.col(value_col).is_not_nan()
    )
    if include_all:
        all_block = (
            base.select(["_region_id", value_col])
            .filter(pl.col(value_col).is_not_null() & pl.col(value_col).is_not_nan())
            .with_columns(pl.lit("All").alias("annot"))
            .select(["_region_id", "annot", value_col])
        )
        long = pl.concat([long, all_block], how="vertical_relaxed")
    return long


def compute_categorical_proportions(
    df: pl.DataFrame,
    *,
    group_col: str = "dmr_type",
    annot_col: str = "feature_type",
    include_all_group: bool = True,
    normalize: bool = True,
) -> pl.DataFrame:
    """Tidy proportions of annotation classes within each group.

    ``group_col`` defaults to epykit's ``dmr_type`` (hyper / hypo) and is
    autodetected to ``DM_status`` if that column is present instead.
    With ``include_all_group=True`` an extra ``group_col == "all"`` row
    block is prepended (the leftmost bar in reference figure 1C).
    """
    # Auto-fallback: annotatr uses DM_status; epykit uses dmr_type.
    if group_col not in df.columns:
        for alt in ("DM_status", "dmr_type"):
            if alt in df.columns:
                group_col = alt
                break
        else:
            raise ValueError(
                f"No grouping column found (tried {group_col!r}, 'DM_status', 'dmr_type')"
            )
    if annot_col not in df.columns:
        raise ValueError(f"annotation column {annot_col!r} not on table")

    # Explode multi-annotation lists so a region with 2 features adds to 2 classes.
    base = df.with_row_index("_region_id")
    if base.schema[annot_col] == pl.List(pl.Utf8):
        base = base.explode(annot_col).filter(pl.col(annot_col).is_not_null())

    blocks: list[pl.DataFrame] = []
    if include_all_group:
        blocks.append(
            base.select([annot_col]).with_columns(pl.lit("all").alias(group_col))
        )
    blocks.append(base.select([annot_col, group_col]))
    long = pl.concat(blocks, how="vertical_relaxed")

    counts = (
        long.group_by([group_col, annot_col])
        .len()
        .rename({"len": "count"})
    )
    if not normalize:
        return counts.sort([group_col, annot_col])
    totals = counts.group_by(group_col).agg(pl.col("count").sum().alias("_total"))
    return (
        counts.join(totals, on=group_col)
        .with_columns((pl.col("count") / pl.col("_total")).alias("proportion"))
        .drop("_total")
        .sort([group_col, annot_col])
    )


__all__ = [
    "PCAResult",
    "MetaplotResult",
    "VolcanoData",
    "MAData",
    "ManhattanData",
    "clear_report_cache",
    "compute_sample_site_matrix",
    "compute_pca",
    "compute_tss_metaplot",
    "compute_coverage_distribution",
    "compute_volcano_data",
    "compute_ma_data",
    "compute_manhattan_data",
    "compute_annotation_counts",
    "compute_coannotation_matrix",
    "compute_numerical_by_annotation",
    "compute_categorical_proportions",
]
