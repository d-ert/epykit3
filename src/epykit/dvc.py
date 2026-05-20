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
from typing import Optional

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
    "chrom":          pl.Utf8,
    "pos":            pl.Int32,
    "strand":         pl.Utf8,
    "n_treatment":    pl.Int32,
    "n_control":      pl.Int32,
    "var_treatment":  pl.Float64,
    "var_control":    pl.Float64,
    "var_log_ratio":  pl.Float64,
    "p_variance":     pl.Float64,
    "q_variance":     pl.Float64,
    "p_mean":         pl.Float64,
    "q_mean":         pl.Float64,
    "is_dvc":         pl.Boolean,
}


def _bartlett_per_site(
    var_a: np.ndarray, n_a: np.ndarray,
    var_b: np.ndarray, n_b: np.ndarray,
) -> np.ndarray:
    """Vectorised Bartlett test for equal variances across two groups.

    Returns p-values per site. NaN where either group has fewer than 2
    observations or zero variance.
    """
    n_a_safe = np.maximum(n_a - 1, 0)
    n_b_safe = np.maximum(n_b - 1, 0)
    N = n_a_safe + n_b_safe
    with np.errstate(invalid="ignore", divide="ignore"):
        pooled = (n_a_safe * var_a + n_b_safe * var_b) / np.maximum(N, 1)
        chi2 = (
            N * np.log(np.maximum(pooled, 1e-300))
            - n_a_safe * np.log(np.maximum(var_a, 1e-300))
            - n_b_safe * np.log(np.maximum(var_b, 1e-300))
        )
        # Bartlett correction term
        c = 1.0 + (1.0 / (3.0 * 1.0)) * (
            (1.0 / np.maximum(n_a_safe, 1)) + (1.0 / np.maximum(n_b_safe, 1))
            - (1.0 / np.maximum(N, 1))
        )
        chi2_corrected = chi2 / c
    valid = (n_a >= 2) & (n_b >= 2) & (var_a > 0) & (var_b > 0)
    pvals = np.full_like(chi2_corrected, np.nan, dtype=np.float64)
    pvals[valid] = sp_stats.chi2.sf(chi2_corrected[valid], df=1)
    return pvals


def _process_one_chromosome_dvc(
    methylstore_path: Path,
    chrom: str,
    canonical_df: pl.DataFrame,
    samples_treatment: list[str],
    samples_control: list[str],
    test: str,
    mean_filter_alpha: float,
    alpha: float,
) -> pl.DataFrame:
    n_sites = len(canonical_df)
    if n_sites == 0:
        return pl.DataFrame(schema=_DVC_EMPTY_SCHEMA)
    canonical_pos = canonical_df.select("pos")

    mean_t, M2_t, n_t = _welford_init(n_sites)
    mean_c, M2_c, n_c = _welford_init(n_sites)

    for s in samples_treatment:
        meth, cov = _load_sample_chrom(methylstore_path, chrom, s, canonical_pos)
        _welford_update(mean_t, M2_t, n_t, meth, cov)
        del meth, cov
    for s in samples_control:
        meth, cov = _load_sample_chrom(methylstore_path, chrom, s, canonical_pos)
        _welford_update(mean_c, M2_c, n_c, meth, cov)
        del meth, cov

    var_t = np.where(n_t > 1, M2_t / np.maximum(n_t - 1, 1), np.nan)
    var_c = np.where(n_c > 1, M2_c / np.maximum(n_c - 1, 1), np.nan)

    if test == "bartlett":
        p_var = _bartlett_per_site(var_t, n_t, var_c, n_c)
    elif test in ("levene", "brown_forsythe"):
        raise ValueError(
            f"DVC test={test!r} requires per-replicate centered deviations, "
            "which the Welford streaming budget doesn't keep. Use "
            "test='bartlett' (closed-form on Welford accumulators) or "
            "implement a per-replicate accumulator."
        )
    else:
        raise ValueError(
            f"DVC test={test!r} not supported. Use 'bartlett'."
        )

    # Welch t on means (mean filter)
    vm_t = np.where(n_t > 1, var_t / np.maximum(n_t, 1), np.nan)
    vm_c = np.where(n_c > 1, var_c / np.maximum(n_c, 1), np.nan)
    se = np.sqrt(np.where((vm_t > 0) | (vm_c > 0), vm_t + vm_c, np.nan))
    with np.errstate(invalid="ignore", divide="ignore"):
        t_stat = np.where(se > 0, (mean_t - mean_c) / se, np.nan)
        dof_num = (vm_t + vm_c) ** 2
        dof_den = (
            np.where(n_t > 1, vm_t ** 2 / np.maximum(n_t - 1, 1), 0.0)
            + np.where(n_c > 1, vm_c ** 2 / np.maximum(n_c - 1, 1), 0.0)
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

    return pl.DataFrame({
        "chrom":         pl.Series([chrom] * n_sites, dtype=pl.Utf8),
        "pos":           canonical_df["pos"],
        "strand":        canonical_df["strand"],
        "n_treatment":   pl.Series(n_t.astype(np.int32)),
        "n_control":     pl.Series(n_c.astype(np.int32)),
        "var_treatment": pl.Series(var_t),
        "var_control":   pl.Series(var_c),
        "var_log_ratio": pl.Series(var_log_ratio),
        "p_variance":    pl.Series(p_var),
        "p_mean":        pl.Series(p_mean),
    }).sort("pos")


def process_chromosomes_dvc(
    methylstore_path: str,
    samples_treatment: list[str],
    samples_control: list[str],
    *,
    test: str = "bartlett",
    chromosomes: Optional[list[str]] = None,
    unite: bool = True,
    mean_filter_alpha: float = 0.05,
    alpha: float = 0.05,
    backend: str = "sequential",
    n_workers: Optional[int] = None,
) -> pl.DataFrame:
    """Run DVC analysis across all chromosomes.

    Returns a DataFrame in ``_DVC_EMPTY_SCHEMA`` with both the variance
    test p/q-values and the mean-test p/q-values (so the caller can apply
    the iEVORA signature filter at any threshold).
    """
    if test != "bartlett":
        raise ValueError(
            f"DVC test must be 'bartlett' (the only Welford-compatible "
            f"variance-equality test); got {test!r}. Levene / Brown-Forsythe "
            "need per-replicate centered deviations that aren't kept under "
            "the streaming accumulator."
        )
    store = Path(methylstore_path)
    all_samples = samples_treatment + samples_control
    if chromosomes is None:
        chromosomes = _detect_chromosomes(store)
        logger.info("DVC: auto-detected %d chromosomes", len(chromosomes))

    from ._compute import run_chrom_pipeline

    def _dvc_chrom_handler(chrom: str) -> Optional[pl.DataFrame]:
        canonical_df = (
            _intersect_chrom(store, chrom, all_samples)
            if unite else _union_chrom(store, chrom, all_samples)
        )
        if len(canonical_df) == 0:
            return None
        return _process_one_chromosome_dvc(
            store, chrom, canonical_df,
            samples_treatment, samples_control,
            test=test, mean_filter_alpha=mean_filter_alpha, alpha=alpha,
        )

    with tempfile.TemporaryDirectory(prefix="epykit_dvc_") as tmpdir:
        tmp = Path(tmpdir)
        written: list[Path] = []
        for chrom, chrom_result in run_chrom_pipeline(
            chromosomes, _dvc_chrom_handler,
            backend=backend, n_workers=n_workers, label="DVC",
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

    return combined.with_columns([
        pl.Series("q_variance", q_var),
        pl.Series("q_mean",     q_mean),
        pl.Series("is_dvc",     is_dvc),
    ])


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
            f"call_dvr_density: dvc_df missing columns {sorted(missing)}. "
            "Run ep.tl.dvc(md) first."
        )
    if len(dvc_df) == 0:
        return pl.DataFrame(schema={
            "chrom": pl.Utf8, "start": pl.Int64, "end": pl.Int64,
            "n_cpgs": pl.Int64, "n_dvc": pl.Int64,
            "frac_dvc": pl.Float64, "pvalue": pl.Float64,
            "qvalue": pl.Float64, "mean_var_log_ratio": pl.Float64,
            "dvr_type": pl.Utf8, "is_dvr": pl.Boolean,
        })

    # Genome-wide DVC rate (the background). Treat null is_dvc / NaN
    # values as not-DVC; they still contribute to the denominator.
    total_cpgs = int(dvc_df.height)
    total_dvc = int(
        dvc_df.get_column("is_dvc").fill_null(False).cast(pl.Int64).sum()
    )
    p0 = total_dvc / total_cpgs if total_cpgs > 0 else 0.0
    if p0 == 0.0:
        logger.warning(
            "call_dvr_density: no DVCs in the input -- every tile is null."
        )

    # Per-tile aggregation.
    tiles = (
        dvc_df.lazy()
        .with_columns([
            (pl.col("pos") // tile_size_bp).alias("_tile"),
            pl.col("is_dvc").fill_null(False).cast(pl.Int64).alias("_is_dvc"),
            # Only count signs of DVCs themselves to determine direction.
            pl.when(pl.col("is_dvc").fill_null(False))
            .then(pl.col("var_log_ratio"))
            .otherwise(None)
            .alias("_dvc_vlr"),
        ])
        .group_by(["chrom", "_tile"])
        .agg([
            pl.col("pos").min().alias("_pos_min"),
            pl.col("pos").max().alias("_pos_max"),
            pl.len().alias("n_cpgs"),
            pl.col("_is_dvc").sum().alias("n_dvc"),
            pl.col("_dvc_vlr").mean().alias("mean_var_log_ratio"),
            pl.col("_dvc_vlr").sign().sum().alias("_vlr_sign_sum"),
        ])
        .filter(pl.col("n_cpgs") >= min_cpgs_per_tile)
        .collect()
    )

    if tiles.is_empty():
        return pl.DataFrame(schema={
            "chrom": pl.Utf8, "start": pl.Int64, "end": pl.Int64,
            "n_cpgs": pl.Int64, "n_dvc": pl.Int64,
            "frac_dvc": pl.Float64, "pvalue": pl.Float64,
            "qvalue": pl.Float64, "mean_var_log_ratio": pl.Float64,
            "dvr_type": pl.Utf8, "is_dvr": pl.Boolean,
        })

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
        sign_sum > 0, "var_up",
        np.where(sign_sum < 0, "var_down", "mixed"),
    )
    # Tiles with 0 DVCs have sign_sum = 0; mark them mixed (or null).
    direction = np.where(n_dvc == 0, "none", direction)

    is_dvr = (qvals < alpha) & (n_dvc > 0)
    frac_dvc = n_dvc / np.maximum(n_cpgs, 1)

    out = tiles.with_columns([
        # End is half-open by convention; +1 so a single-CpG tile has end>start.
        (pl.col("_pos_min")).alias("start"),
        (pl.col("_pos_max") + 1).alias("end"),
        pl.Series("frac_dvc", frac_dvc),
        pl.Series("pvalue", pvals),
        pl.Series("qvalue", qvals),
        pl.Series("dvr_type", direction),
        pl.Series("is_dvr", is_dvr),
    ]).select([
        "chrom", "start", "end", "n_cpgs", "n_dvc", "frac_dvc",
        "pvalue", "qvalue", "mean_var_log_ratio", "dvr_type", "is_dvr",
    ]).sort(["chrom", "start"])

    return out
