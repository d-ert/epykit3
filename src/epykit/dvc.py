"""Differentially Variable CpG (DVC) calling.

iEVORA-style: test for variance differences between groups at each CpG.
Sites are flagged DVC when:

    q_variance < alpha   AND   p_mean   > mean_filter_alpha

i.e. the between-group variance differs significantly while the means do
not -- the signature of an outlier-driven shift in variability (common in
cancer / aging methylomes) rather than a simple mean shift.

Memory and I/O follow the same per-chromosome streaming layout as
``dmc.process_chromosomes_dmc``: one sample is loaded at a time, per-site
Welford accumulators give variance, no (n_sites x n_replicates) matrix
is ever built.
"""

from __future__ import annotations

import gc
import logging
import tempfile
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats as sp_stats

from .dmc import (
    _detect_chromosomes,
    _intersect_chrom,
    _load_sample_chrom,
    _union_chrom,
    _welford_init,
    _welford_update,
)

logger = logging.getLogger(__name__)

_DVC_EMPTY_SCHEMA = {
    "chrom": pl.Utf8,
    "pos": pl.Int32,
    "strand": pl.Utf8,
    "n_treatment": pl.Int32,
    "n_control": pl.Int32,
    "var_treatment": pl.Float64,
    "var_control": pl.Float64,
    "var_log_ratio": pl.Float64,
    "p_variance": pl.Float64,
    "q_variance": pl.Float64,
    "p_mean": pl.Float64,
    "q_mean": pl.Float64,
    "is_dvc": pl.Boolean,
}


def _per_site_variance_test(group_a: np.ndarray, group_b: np.ndarray) -> tuple[float, float]:
    """Brown-Forsythe variance equality test (median-centred Levene).

    Equivalent to ``scipy.stats.levene(group_a, group_b, center='median')``.
    Robust to non-normality; correct for bounded U-shaped beta values.

    Parameters
    ----------
    group_a, group_b:
        1-D float arrays of per-replicate beta values (NaN-ignored).

    Returns
    -------
    (f_stat, p_val) -- both NaN if either group has fewer than 2 finite values.
    """
    a = np.asarray(group_a, dtype=np.float64)
    b = np.asarray(group_b, dtype=np.float64)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    # Pass 1: median-centred absolute deviations (Brown-Forsythe).
    z_a = np.abs(a - np.median(a))
    z_b = np.abs(b - np.median(b))
    # Pass 2: one-way F-test on the deviation arrays.
    f_stat, p_val = sp_stats.f_oneway(z_a, z_b)
    return float(f_stat), float(p_val)


def _brown_forsythe_vectorised(
    betas_a: list[np.ndarray],
    betas_b: list[np.ndarray],
    n_sites: int,
) -> np.ndarray:
    """Brown-Forsythe p-values computed vectorially across all sites.

    Parameters
    ----------
    betas_a:
        List of length n_samples_a; each element is a float64 array of
        length n_sites with NaN where coverage == 0.
    betas_b:
        Same for the control group.
    n_sites:
        Number of CpG sites (length of each array in betas_a/betas_b).

    Returns
    -------
    p_var : np.ndarray of shape (n_sites,), dtype float64.
        NaN where either group has fewer than 3 finite observations.

    Notes
    -----
    Requires >=3 per group. At n=2 the two median-centred absolute deviations
    are identical (median is the midpoint), so the within-group sum of squares
    is mathematically 0; floating-point residuals then make it a tiny positive
    number, which inflates the F-statistic and yields spuriously significant
    p-values (anti-conservative). Excluding n<3 sites keeps the test honest.
    """
    # Stack to (n_samples x n_sites) matrices; NaN = missing.
    mat_a = np.stack(betas_a, axis=0)  # shape (na, n_sites)
    mat_b = np.stack(betas_b, axis=0)  # shape (nb, n_sites)

    # Median per site ignoring NaN.
    with np.errstate(invalid="ignore"):
        med_a = np.nanmedian(mat_a, axis=0)  # (n_sites,)
        med_b = np.nanmedian(mat_b, axis=0)

    # Absolute deviations from group median.
    z_a = np.abs(mat_a - med_a[np.newaxis, :])  # (na, n_sites)
    z_b = np.abs(mat_b - med_b[np.newaxis, :])  # (nb, n_sites)

    # Count finite observations per site per group.
    finite_a = np.isfinite(mat_a)  # (na, n_sites)
    finite_b = np.isfinite(mat_b)
    n_a = finite_a.sum(axis=0).astype(np.float64)  # (n_sites,)
    n_b = finite_b.sum(axis=0).astype(np.float64)

    # Brown-Forsythe needs >=3 finite obs per group: at n=2 the median-centred
    # deviations are identical, so within-group SS is ~0 (a floating-point
    # residual) and the F-statistic explodes into spurious significance.
    valid = (n_a >= 3) & (n_b >= 3)
    p_var = np.full(n_sites, np.nan, dtype=np.float64)
    if not np.any(valid):
        return p_var

    # NaN deviations for masked sites (coverage 0).
    z_a = np.where(finite_a, z_a, np.nan)
    z_b = np.where(finite_b, z_b, np.nan)

    # Vectorised one-way F-test over deviation arrays per site.
    # mean_i = nanmean of z_i; grand_mean = weighted combo.
    mean_a = np.nanmean(z_a, axis=0)  # (n_sites,)
    mean_b = np.nanmean(z_b, axis=0)
    N = n_a + n_b
    grand_mean = (n_a * mean_a + n_b * mean_b) / np.maximum(N, 1)

    # Between-group SS (df=1 for two groups).
    ss_between = n_a * (mean_a - grand_mean) ** 2 + n_b * (mean_b - grand_mean) ** 2
    # Within-group SS (df = N - 2 for two groups).
    ss_within_a = np.nansum((z_a - mean_a[np.newaxis, :]) ** 2, axis=0)
    ss_within_b = np.nansum((z_b - mean_b[np.newaxis, :]) ** 2, axis=0)
    ss_within = ss_within_a + ss_within_b
    df_within = N - 2

    with np.errstate(invalid="ignore", divide="ignore"):
        ms_between = ss_between  # df_between = 1
        ms_within = np.where(df_within > 0, ss_within / df_within, np.nan)
        f_stat = np.where(ms_within > 0, ms_between / ms_within, np.nan)

    p_var[valid] = sp_stats.f.sf(f_stat[valid], dfn=1, dfd=df_within[valid])
    return p_var


def _process_one_chromosome_dvc(
    methylstore_path: Path,
    chrom: str,
    canonical_df: pl.DataFrame,
    samples_treatment: list[str],
    samples_control: list[str],
    test: str,
    mean_filter_alpha: float,
    alpha: float,
    min_coverage: int = 1,
) -> pl.DataFrame:
    n_sites = len(canonical_df)
    if n_sites == 0:
        return pl.DataFrame(schema=_DVC_EMPTY_SCHEMA)
    canonical_pos = canonical_df.select("pos")

    # Accumulate per-sample beta arrays for Brown-Forsythe (median-centred
    # Levene). Each list holds one float64 array of length n_sites; NaN
    # where coverage == 0. Memory overhead is O(n_sites * n_replicates)
    # which is tolerable (n_replicates typically 2-20).
    betas_t: list[np.ndarray] = []
    betas_c: list[np.ndarray] = []

    # We still use Welford accumulators for the mean-shift (Welch t) filter —
    # no need to iterate again.
    mean_t, M2_t, n_t = _welford_init(n_sites)
    mean_c, M2_c, n_c = _welford_init(n_sites)

    for s in samples_treatment:
        meth, cov = _load_sample_chrom(methylstore_path, chrom, s, canonical_pos)
        with np.errstate(invalid="ignore", divide="ignore"):
            beta = np.where(cov >= min_coverage, meth.astype(np.float64) / cov, np.nan)
        betas_t.append(beta)
        _welford_update(mean_t, M2_t, n_t, meth, cov)
        del meth, cov
    for s in samples_control:
        meth, cov = _load_sample_chrom(methylstore_path, chrom, s, canonical_pos)
        with np.errstate(invalid="ignore", divide="ignore"):
            beta = np.where(cov >= min_coverage, meth.astype(np.float64) / cov, np.nan)
        betas_c.append(beta)
        _welford_update(mean_c, M2_c, n_c, meth, cov)
        del meth, cov

    var_t = np.where(n_t > 1, M2_t / np.maximum(n_t - 1, 1), np.nan)
    var_c = np.where(n_c > 1, M2_c / np.maximum(n_c - 1, 1), np.nan)

    if test in ("brown_forsythe", "bartlett"):
        # "bartlett" accepted as a legacy alias; always run Brown-Forsythe.
        p_var = _brown_forsythe_vectorised(betas_t, betas_c, n_sites)
    else:
        raise ValueError(
            f"DVC test={test!r} not supported. "
            "Use 'brown_forsythe' (default) or legacy alias 'bartlett'."
        )

    # Welch t on means (mean filter)
    vm_t = np.where(n_t > 1, var_t / np.maximum(n_t, 1), np.nan)
    vm_c = np.where(n_c > 1, var_c / np.maximum(n_c, 1), np.nan)
    se = np.sqrt(np.where((vm_t > 0) | (vm_c > 0), vm_t + vm_c, np.nan))
    with np.errstate(invalid="ignore", divide="ignore"):
        t_stat = np.where(se > 0, (mean_t - mean_c) / se, np.nan)
        dof_num = (vm_t + vm_c) ** 2
        dof_den = np.where(n_t > 1, vm_t**2 / np.maximum(n_t - 1, 1), 0.0) + np.where(
            n_c > 1, vm_c**2 / np.maximum(n_c - 1, 1), 0.0
        )
        dof = np.where(dof_den > 0, dof_num / dof_den, 1.0)
        dof = np.maximum(dof, 1.0)
    p_mean = 2.0 * sp_stats.t.sf(np.abs(t_stat), df=dof)

    with np.errstate(invalid="ignore", divide="ignore"):
        var_log_ratio = np.where(
            (var_t > 0) & (var_c > 0),
            np.log2(var_t / var_c),
            np.nan,
        )

    return pl.DataFrame(
        {
            "chrom": pl.Series([chrom] * n_sites, dtype=pl.Utf8),
            "pos": canonical_df["pos"],
            "strand": canonical_df["strand"],
            "n_treatment": pl.Series(n_t.astype(np.int32)),
            "n_control": pl.Series(n_c.astype(np.int32)),
            "var_treatment": pl.Series(var_t),
            "var_control": pl.Series(var_c),
            "var_log_ratio": pl.Series(var_log_ratio),
            "p_variance": pl.Series(p_var),
            "p_mean": pl.Series(p_mean),
        }
    ).sort("pos")


def process_chromosomes_dvc(
    methylstore_path: str,
    samples_treatment: list[str],
    samples_control: list[str],
    *,
    test: str = "brown_forsythe",
    chromosomes: list[str] | None = None,
    unite: bool = True,
    mean_filter_alpha: float = 0.05,
    alpha: float = 0.05,
    min_coverage: int = 1,
    backend: str = "sequential",
    n_workers: int | None = None,
) -> pl.DataFrame:
    """Run DVC analysis across all chromosomes.

    Returns a DataFrame in ``_DVC_EMPTY_SCHEMA`` with both the variance
    test p/q-values and the mean-test p/q-values (so the caller can apply
    the iEVORA signature filter at any threshold).

    The variance equality test is Brown-Forsythe (median-centred Levene),
    robust to the U-shaped, bounded [0,1] distribution of beta values.
    ``test='bartlett'`` is a deprecated alias (Bartlett's test is not
    implemented) and emits a ``UserWarning``.

    ``min_coverage`` masks per-replicate beta below the threshold before the
    variance test. At very low coverage every beta is forced toward {0, 1}
    with large, coverage-dependent binomial noise; if coverage is imbalanced
    between groups that noise reads as differential biological variance.
    Raise ``min_coverage`` (e.g. 5-10) on cohorts with uneven depth.
    """
    import warnings

    if test not in ("brown_forsythe", "bartlett"):
        raise ValueError(
            f"DVC test={test!r} not supported. "
            "Use 'brown_forsythe' (default) or legacy alias 'bartlett'."
        )
    if test == "bartlett":
        warnings.warn(
            "DVC test='bartlett' is a deprecated alias: Bartlett's test is "
            "not implemented and Brown-Forsythe (median-centred Levene) is "
            "run instead. Pass test='brown_forsythe' to silence this.",
            UserWarning,
            stacklevel=2,
        )
    min_n = min(len(samples_treatment), len(samples_control))
    if min_n < 3:
        warnings.warn(
            f"DVC: Brown-Forsythe needs >=3 replicates per group to estimate "
            f"the within-group spread of deviations, but the smaller group has "
            f"{min_n}. At n=2 the within-group sum of squares is exactly 0, so "
            f"the variance p-values are NaN and no DVCs will be called. Add "
            f"replicates or interpret n_dvc=0 as 'test could not run'.",
            UserWarning,
            stacklevel=2,
        )
    min_coverage = max(1, int(min_coverage))
    store = Path(methylstore_path)
    all_samples = samples_treatment + samples_control
    if chromosomes is None:
        chromosomes = _detect_chromosomes(store)
        logger.info("DVC: auto-detected %d chromosomes", len(chromosomes))

    from ._compute import run_chrom_pipeline

    def _dvc_chrom_handler(chrom: str) -> pl.DataFrame | None:
        canonical_df = (
            _intersect_chrom(store, chrom, all_samples)
            if unite
            else _union_chrom(store, chrom, all_samples)
        )
        if len(canonical_df) == 0:
            return None
        return _process_one_chromosome_dvc(
            store,
            chrom,
            canonical_df,
            samples_treatment,
            samples_control,
            test=test,
            mean_filter_alpha=mean_filter_alpha,
            alpha=alpha,
            min_coverage=min_coverage,
        )

    with tempfile.TemporaryDirectory(prefix="epykit_dvc_") as tmpdir:
        tmp = Path(tmpdir)
        written: list[Path] = []
        for chrom, chrom_result in run_chrom_pipeline(
            chromosomes,
            _dvc_chrom_handler,
            backend=backend,
            n_workers=n_workers,
            label="DVC",
        ):
            tmp_file = tmp / f"{chrom}.parquet"
            chrom_result.write_parquet(str(tmp_file))
            written.append(tmp_file)
            del chrom_result
            gc.collect()

        if not written:
            return pl.DataFrame(schema=_DVC_EMPTY_SCHEMA)
        combined = pl.concat([pl.read_parquet(str(f)) for f in written])

    # BH-correct p_variance and p_mean separately.
    from statsmodels.stats.multitest import multipletests

    def _bh(p: np.ndarray) -> np.ndarray:
        finite = np.isfinite(p)
        q = np.full_like(p, np.nan, dtype=np.float64)
        if finite.any():
            _, q_finite, _, _ = multipletests(p[finite], method="fdr_bh")
            q[finite] = q_finite
        return q

    p_var = combined.get_column("p_variance").to_numpy()
    p_mean_arr = combined.get_column("p_mean").to_numpy()
    q_var = _bh(p_var)
    q_mean = _bh(p_mean_arr)
    is_dvc = (q_var < alpha) & (p_mean_arr > mean_filter_alpha)

    return combined.with_columns(
        [
            pl.Series("q_variance", q_var),
            pl.Series("q_mean", q_mean),
            pl.Series("is_dvc", is_dvc),
        ]
    )


# Differentially Variable Regions (DVR) -- density-based aggregation


def call_dvr_density(
    dvc_df: pl.DataFrame,
    *,
    tile_size_bp: int = 1000,
    min_cpgs_per_tile: int = 5,
    alpha: float = 0.05,
) -> pl.DataFrame:
    """Region-level DVC aggregation via density enrichment.

    Why not Stouffer's combining of per-CpG variance p-values? Variance
    test statistics aren't symmetric around zero and don't combine
    cleanly under independence assumptions (the Fisher / Stouffer
    combiners both assume z-transformable null distributions). Counting
    DVCs per tile and asking "is this tile enriched for DVCs vs the
    genome-wide rate?" sidesteps that and gives a region call with an
    interpretable test statistic:

      H0: DVCs distribute uniformly across CpGs at rate
          p_0 = total_DVCs / total_CpGs.
      Per-tile test: Binomial(n=n_cpgs_in_tile, p=p_0); call is the
          number of DVCs in the tile; one-sided enrichment p-value.

    The tile direction (``dvr_type``) is the majority sign of
    ``var_log_ratio`` over DVC sites in the tile.

    Parameters
    ----------
    dvc_df : pl.DataFrame
        The DVC table produced by :func:`process_chromosomes_dvc`. Must
        carry ``chrom``, ``pos``, ``var_log_ratio``, ``is_dvc``.
    tile_size_bp : int
        Tile width along the genome. Default 1 kb.
    min_cpgs_per_tile : int
        Tiles with fewer than this many DVC-tested CpGs are dropped
        before the binomial test. Default 5.
    alpha : float
        BH-q threshold for the returned ``is_dvr`` flag. Default 0.05.

    Returns
    -------
    pl.DataFrame
        Columns: ``chrom``, ``start``, ``end``, ``n_cpgs``, ``n_dvc``,
        ``frac_dvc``, ``pvalue``, ``qvalue``, ``mean_var_log_ratio``,
        ``dvr_type``  in  {var_up, var_down, mixed}, ``is_dvr``.
    """
    needed = {"chrom", "pos", "var_log_ratio", "is_dvc"}
    missing = needed - set(dvc_df.columns)
    if missing:
        raise ValueError(
            f"call_dvr_density: dvc_df missing columns {sorted(missing)}. Run ep.tl.dvc(md) first."
        )
    if len(dvc_df) == 0:
        return pl.DataFrame(
            schema={
                "chrom": pl.Utf8,
                "start": pl.Int64,
                "end": pl.Int64,
                "n_cpgs": pl.Int64,
                "n_dvc": pl.Int64,
                "frac_dvc": pl.Float64,
                "pvalue": pl.Float64,
                "qvalue": pl.Float64,
                "mean_var_log_ratio": pl.Float64,
                "dvr_type": pl.Utf8,
                "is_dvr": pl.Boolean,
            }
        )

    # Genome-wide DVC rate (the background). Treat null is_dvc / NaN
    # values as not-DVC; they still contribute to the denominator.
    total_cpgs = int(dvc_df.height)
    total_dvc = int(dvc_df.get_column("is_dvc").fill_null(False).cast(pl.Int64).sum())
    p0 = total_dvc / total_cpgs if total_cpgs > 0 else 0.0
    if p0 == 0.0:
        logger.warning("call_dvr_density: no DVCs in the input -- every tile is null.")

    # Per-tile aggregation.
    tiles = (
        dvc_df.lazy()
        .with_columns(
            [
                (pl.col("pos") // tile_size_bp).alias("_tile"),
                pl.col("is_dvc").fill_null(False).cast(pl.Int64).alias("_is_dvc"),
                # Only count signs of DVCs themselves to determine direction.
                pl.when(pl.col("is_dvc").fill_null(False))
                .then(pl.col("var_log_ratio"))
                .otherwise(None)
                .alias("_dvc_vlr"),
            ]
        )
        .group_by(["chrom", "_tile"])
        .agg(
            [
                pl.col("pos").min().alias("_pos_min"),
                pl.col("pos").max().alias("_pos_max"),
                pl.len().alias("n_cpgs"),
                pl.col("_is_dvc").sum().alias("n_dvc"),
                pl.col("_dvc_vlr").mean().alias("mean_var_log_ratio"),
                pl.col("_dvc_vlr").sign().sum().alias("_vlr_sign_sum"),
            ]
        )
        .filter(pl.col("n_cpgs") >= min_cpgs_per_tile)
        .collect()
    )

    if tiles.is_empty():
        return pl.DataFrame(
            schema={
                "chrom": pl.Utf8,
                "start": pl.Int64,
                "end": pl.Int64,
                "n_cpgs": pl.Int64,
                "n_dvc": pl.Int64,
                "frac_dvc": pl.Float64,
                "pvalue": pl.Float64,
                "qvalue": pl.Float64,
                "mean_var_log_ratio": pl.Float64,
                "dvr_type": pl.Utf8,
                "is_dvr": pl.Boolean,
            }
        )

    n_cpgs = tiles.get_column("n_cpgs").to_numpy()
    n_dvc = tiles.get_column("n_dvc").to_numpy()
    sign_sum = tiles.get_column("_vlr_sign_sum").to_numpy()

    # One-sided upper-tail binomial: P(X >= n_dvc | n_cpgs, p0).
    if p0 > 0.0:
        pvals = sp_stats.binom.sf(n_dvc - 1, n_cpgs, p0)
    else:
        # No background DVCs => any tile with >=1 DVC has p~=0;
        # tiles with 0 DVCs have p=1.
        pvals = np.where(n_dvc > 0, 0.0, 1.0)
    pvals = np.clip(pvals, np.finfo(float).tiny, 1.0)

    from statsmodels.stats.multitest import multipletests

    finite = np.isfinite(pvals)
    qvals = np.full_like(pvals, np.nan, dtype=np.float64)
    if finite.any():
        _, q_finite, _, _ = multipletests(pvals[finite], method="fdr_bh")
        qvals[finite] = q_finite

    # Tile direction: majority sign of DVC var_log_ratio.
    direction = np.where(
        sign_sum > 0,
        "var_up",
        np.where(sign_sum < 0, "var_down", "mixed"),
    )
    # Tiles with 0 DVCs have sign_sum = 0; mark them mixed (or null).
    direction = np.where(n_dvc == 0, "none", direction)

    is_dvr = (qvals < alpha) & (n_dvc > 0)
    frac_dvc = n_dvc / np.maximum(n_cpgs, 1)

    out = (
        tiles.with_columns(
            [
                # End is half-open by convention; +1 so a single-CpG tile has end>start.
                (pl.col("_pos_min")).alias("start"),
                (pl.col("_pos_max") + 1).alias("end"),
                pl.Series("frac_dvc", frac_dvc),
                pl.Series("pvalue", pvals),
                pl.Series("qvalue", qvals),
                pl.Series("dvr_type", direction),
                pl.Series("is_dvr", is_dvr),
            ]
        )
        .select(
            [
                "chrom",
                "start",
                "end",
                "n_cpgs",
                "n_dvc",
                "frac_dvc",
                "pvalue",
                "qvalue",
                "mean_var_log_ratio",
                "dvr_type",
                "is_dvr",
            ]
        )
        .sort(["chrom", "start"])
    )

    return out
