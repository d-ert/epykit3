"""Epigenetic age clocks and reference-based cell-type deconvolution.

These are *infrastructure* -- generic runners that consume user-supplied
coefficient / reference tables -- not bundled scientific data. The
coefficient tables for the published clocks (Horvath 2013, Hannum 2013,
PhenoAge, DunedinPACE) and the EpiDISH / CIBERSORT reference matrices
have their own licences and citation requirements; we don't redistribute
them. The user points the runner at a CSV / Parquet table and gets a
per-sample age estimate or cell-type composition.

Both runners take a CpG -> coefficient table and a CpG -> (chrom, pos)
manifest. The clocks were trained on Illumina array CpG IDs (cgXXXXXXXX)
which don't carry genomic coordinates directly, so the user-supplied
manifest is what wires probe IDs to the WGBS coordinate system. A
manifest is typically the Illumina HumanMethylation450 / EPIC v1.0
annotation distributed by the array vendor (or the bioconductor packages
``IlluminaHumanMethylation450kanno.ilmn12.hg19`` /
``IlluminaHumanMethylationEPICanno.ilm10b4.hg19``).

The two runners share the same input shape so a user can ship a "clock
table" and a "reference matrix" through the same plumbing:

* Clock table: ``cpg_id``, ``coefficient``; an optional intercept passed
  as a kwarg.
* Reference matrix: ``cpg_id`` + one column per cell type -- the rest of
  the values are the reference beta profile per type.
"""

from __future__ import annotations

import logging

import numpy as np
import polars as pl

logger = logging.getLogger(__name__)


def _build_sample_beta_at_clock_cpgs(
    md,
    coords: pl.DataFrame,
    *,
    coord_cpg_col: str = "cpg_id",
    chrom_col: str = "chrom",
    pos_col: str = "pos",
) -> tuple[np.ndarray, list[str], list[str], np.ndarray]:
    """Build a (n_samples x n_cpgs) beta matrix at the clock / reference
    CpGs.

    Returns
    -------
    beta : np.ndarray, shape (n_samples, n_cpgs)
        beta at clock CpGs in coords' row order. NaN where the sample has
        no coverage at that site (or the site is missing from the
        store).
    sample_ids : list[str]
    cpg_ids : list[str]
        Identifiers of the columns of ``beta``; useful for downstream
        reporting / debugging.
    missing_frac_per_sample : np.ndarray, shape (n_samples,)
        Fraction of clock CpGs that resolved to NaN per sample.
    """
    needed = {coord_cpg_col, chrom_col, pos_col}
    missing_cols = needed - set(coords.columns)
    if missing_cols:
        raise ValueError(
            f"coord table missing columns {sorted(missing_cols)}. "
            f"Expected at minimum: cpg_id, chrom, pos."
        )

    # Materialise the target site set as a polars table for join.
    target = coords.select(
        [
            pl.col(coord_cpg_col).alias("cpg_id"),
            pl.col(chrom_col).cast(pl.Utf8).alias("chrom"),
            pl.col(pos_col).cast(pl.Int64).alias("pos"),
        ]
    )
    cpg_ids = target.get_column("cpg_id").to_list()
    n_cpgs = target.height

    sample_ids = md.obs.get_column("sample_id").to_list()
    n_samples = len(sample_ids)
    beta = np.full((n_samples, n_cpgs), np.nan, dtype=np.float64)

    # Index into target by (chrom, pos) for fast lookup.
    target_keyed = target.with_columns(pl.arange(0, n_cpgs).alias("_row"))
    for s_idx, sample in enumerate(sample_ids):
        sample_df = (
            pl.scan_parquet(f"{md.store}/sample={sample}/chrom=*/part-*.parquet")
            .select(["chrom", "pos", "N_meth", "coverage"])
            .filter(pl.col("coverage") > 0)
            .with_columns(pl.col("pos").cast(pl.Int64))
            .collect()
        )
        joined = target_keyed.join(
            sample_df,
            on=["chrom", "pos"],
            how="inner",
        )
        if joined.is_empty():
            continue
        rows = joined.get_column("_row").to_numpy()
        N_meth = joined.get_column("N_meth").to_numpy().astype(np.float64)
        cov = joined.get_column("coverage").to_numpy().astype(np.float64)
        beta[s_idx, rows] = N_meth / np.maximum(cov, 1.0)

    missing_frac = np.isnan(beta).mean(axis=1)
    return beta, sample_ids, cpg_ids, missing_frac


def _horvath_anti_transform(linear, adult_age: float = 20.0):
    """Horvath (2013) age anti-transformation (``anti.trafo``).

    ``F(x) = (1 + adult_age) * exp(x) - 1``           for ``x < 0``
    ``F(x) = (1 + adult_age) * x   + adult_age``       for ``x >= 0``

    With the canonical ``adult_age = 20`` this is ``21*exp(x) - 1`` and
    ``21*x + 20``. The linear predictor is negative for any predicted age
    below ``adult_age``, so the ``(1 + adult_age)`` factor on the negative
    branch is load-bearing for pediatric / young samples.
    """
    linear = np.asarray(linear, dtype=float)
    return np.where(
        linear < 0,
        (1.0 + adult_age) * np.exp(linear) - 1.0,
        (1.0 + adult_age) * linear + adult_age,
    )


def age_clock(
    md,
    coefficients: pl.DataFrame,
    manifest: pl.DataFrame,
    *,
    intercept: float = 0.0,
    transform: str | None = None,
    coef_cpg_col: str = "cpg_id",
    coef_value_col: str = "coefficient",
    manifest_cpg_col: str = "cpg_id",
    impute_missing: bool = True,
    name: str = "age_clock",
) -> pl.DataFrame:
    """Linear epigenetic-age clock runner.

    Computes one age estimate per sample as

        age_hat = transform( intercept + Sigma_i coefficient_i * beta_i )

    where the sum is over CpGs in ``coefficients`` that are also in
    ``manifest`` (which maps probe IDs to genomic coordinates) and
    actually present in ``md.store``.

    Parameters
    ----------
    md : MethylData
        Methylation store.
    coefficients : pl.DataFrame
        Two-column table from the published clock: a CpG ID column and a
        coefficient column. The intercept, if any, is supplied
        separately via ``intercept``.
    manifest : pl.DataFrame
        Probe -> (chrom, pos) lookup. Must carry the ``manifest_cpg_col``
        column plus ``chrom`` and ``pos``. Typically the array vendor's
        annotation file.
    intercept : float
        Linear-model intercept. Default 0.0 (e.g. Hannum's blood clock
        -- Horvath needs ~0.696).
    transform : {"horvath", None}, optional
        Post-linear transformation. ``"horvath"`` applies the standard
        anti-transform for samples >=20 years old::

            age = (1 + adult_age)*exp(linear) - 1   if linear < 0
            age = (1 + adult_age)*linear + adult_age  otherwise   (adult_age = 20)

        ``None`` (default) returns the raw linear combination.
    coef_cpg_col, coef_value_col : str
        Column names in ``coefficients``. Defaults assume the standard
        ``cpg_id`` / ``coefficient`` layout.
    manifest_cpg_col : str
        Probe-ID column name in ``manifest``. Default ``cpg_id``.
    impute_missing : bool
        If True (default), missing beta values are replaced with the mean
        beta at that CpG across the samples that *do* have coverage. If no
        sample covers a CpG, that CpG drops out (its coefficient is
        ignored) -- its weight is redistributed implicitly by being
        absent from the sum, which biases the result; the returned
        table flags how many CpGs were dropped per sample so the user
        can decide whether the estimate is trustworthy.
    name : str
        Stored under ``md.obs[<name>]`` and used in the returned table
        column. Default ``"age_clock"``.

    Returns
    -------
    pl.DataFrame
        One row per sample with columns ``sample_id``, ``<name>``,
        ``n_cpgs_used``, ``n_cpgs_missing``, ``missing_frac``.
    """
    if coef_cpg_col not in coefficients.columns or coef_value_col not in coefficients.columns:
        raise ValueError(
            f"coefficients table must carry '{coef_cpg_col}' and '{coef_value_col}' columns."
        )
    # Resolve probe -> genomic coordinates.
    coords = (
        manifest.select(
            [
                pl.col(manifest_cpg_col).alias("cpg_id"),
                pl.col("chrom"),
                pl.col("pos"),
            ]
        )
        .join(
            coefficients.select(
                [
                    pl.col(coef_cpg_col).alias("cpg_id"),
                    pl.col(coef_value_col).alias("coefficient"),
                ]
            ),
            on="cpg_id",
            how="inner",
        )
        .drop_nulls(subset=["chrom", "pos"])
    )
    if coords.is_empty():
        raise ValueError(
            "No clock CpGs resolved to genomic coordinates. Check that "
            "the manifest's probe IDs overlap with the coefficient "
            "table's IDs."
        )
    logger.info(
        "age_clock: resolved %d / %d clock CpGs to coordinates.",
        coords.height,
        coefficients.height,
    )

    beta, sample_ids, _cpg_ids, missing_frac = _build_sample_beta_at_clock_cpgs(
        md,
        coords,
    )
    coefs = coords.get_column("coefficient").to_numpy().astype(np.float64)

    # Per-CpG mean across samples (ignoring NaN) for optional imputation.
    if impute_missing:
        with np.errstate(invalid="ignore"):
            cpg_mean = np.nanmean(beta, axis=0)
        beta_imp = np.where(np.isnan(beta), cpg_mean[None, :], beta)
    else:
        beta_imp = beta

    # CpGs with no coverage at *any* sample collapse out of the sum
    # regardless of impute_missing.
    cpg_has_any = np.any(~np.isnan(beta), axis=0)
    n_cpgs_used = int(cpg_has_any.sum())
    n_cpgs_missing = int(coords.height - n_cpgs_used)

    # M-SEC1: with impute_missing=False a CpG missing in one sample (but seen
    # in others) contributes 0 to that sample's linear predictor -- i.e. the
    # age is a sum over observed CpGs only. (Note: zero-filling and excluding
    # the term are mathematically identical; the issue is the *bias* from the
    # missing contributions, not the zero-fill itself.) A clock's intercept is
    # calibrated assuming all CpGs are present, so this biases absolute ages.
    if not impute_missing:
        per_sample_missing = np.isnan(beta) & cpg_has_any[None, :]
        if bool(per_sample_missing.any()):
            import warnings

            warnings.warn(
                "age_clock(impute_missing=False): some clock CpGs are missing "
                "in individual samples and contribute 0 to those samples' "
                "linear predictor (age summed over observed CpGs only). The "
                "clock intercept assumes full coverage, so absolute ages are "
                "biased by the missing CpGs' typical contributions. Use "
                "impute_missing=True (imputes to the per-CpG cohort mean), and "
                "prefer a clock with complete coverage; treat ages as "
                "unreliable when missing_frac is non-trivial.",
                UserWarning,
                stacklevel=2,
            )

    # Replace NaN coefficientxbeta products with 0 for the sum (post-impute
    # any remaining NaN is exactly the cpg_has_any==False columns).
    beta_clean = np.where(np.isnan(beta_imp), 0.0, beta_imp)
    coefs_clean = np.where(cpg_has_any, coefs, 0.0)
    linear = intercept + beta_clean @ coefs_clean

    if transform is None:
        age = linear
    elif transform == "horvath":
        age = _horvath_anti_transform(linear)
    else:
        raise ValueError(f"Unknown transform {transform!r}; use None or 'horvath'.")

    out = pl.DataFrame(
        {
            "sample_id": sample_ids,
            name: age,
            "n_cpgs_used": [n_cpgs_used] * len(sample_ids),
            "n_cpgs_missing": [n_cpgs_missing] * len(sample_ids),
            "missing_frac": missing_frac,
        }
    )
    return out


def deconvolve(
    md,
    reference: pl.DataFrame,
    manifest: pl.DataFrame,
    *,
    method: str = "nnls",
    cell_types: list[str] | None = None,
    manifest_cpg_col: str = "cpg_id",
    ref_cpg_col: str = "cpg_id",
) -> pl.DataFrame:
    """Reference-based cell-type deconvolution.

    Solves, per sample,

        beta_sample ~= R * pi

    where ``R`` is the (n_cpgs x n_cell_types) reference beta matrix and
    ``pi`` is the (n_cell_types,) composition vector. By default the
    solve is non-negative least squares (``method="nnls"``, the
    Houseman / EpiDISH "CP" estimator); the composition is then
    re-normalised to sum to 1.

    The standard published references -- EpiDISH ``centDHSbloodDMC.m``
    for whole blood, the saliva and breast reference panels, etc. -- are
    distributed as Illumina-array beta tables and licensed separately; you
    supply the reference matrix and the array manifest, this runner
    does the math.

    Parameters
    ----------
    md : MethylData
        Sample beta matrix source.
    reference : pl.DataFrame
        ``ref_cpg_col`` column + one column per cell type. Values are
        per-CpG mean beta in that cell type.
    manifest : pl.DataFrame
        Probe -> (chrom, pos) lookup with the column named by
        ``manifest_cpg_col``.
    method : {"nnls"}
        Solver. Only ``"nnls"`` is implemented (the standard Houseman
        estimator). RPC / CIBERSORT-style robust regression are future
        work.
    cell_types : list[str], optional
        Subset / re-order the reference columns. Default uses every
        non-CpG-id column of ``reference``.
    manifest_cpg_col, ref_cpg_col : str
        Probe-ID column names in the manifest and reference tables.

    Returns
    -------
    pl.DataFrame
        One row per sample with columns ``sample_id``, ``cell_type``,
        ``proportion``. Long format so it joins cleanly onto ``md.obs``
        (filter to a specific cell type and pivot if you need a wide
        layout).
    """
    if method != "nnls":
        raise ValueError(
            f"deconvolve: method={method!r} not implemented. Only 'nnls' "
            "(non-negative least squares, the Houseman / EpiDISH CP "
            "estimator) is available."
        )
    if ref_cpg_col not in reference.columns:
        raise ValueError(f"reference table missing CpG column '{ref_cpg_col}'.")
    if cell_types is None:
        cell_types = [c for c in reference.columns if c != ref_cpg_col]
    if not cell_types:
        raise ValueError("No cell-type columns identified in reference.")
    for ct in cell_types:
        if ct not in reference.columns:
            raise ValueError(f"cell type {ct!r} not in reference columns")

    # Resolve probe -> coords. Inner join on CpG ID then on manifest
    # coords; rows where any join misses are dropped.
    coords = (
        manifest.select(
            [
                pl.col(manifest_cpg_col).alias("cpg_id"),
                pl.col("chrom"),
                pl.col("pos"),
            ]
        )
        .join(
            reference.rename({ref_cpg_col: "cpg_id"}),
            on="cpg_id",
            how="inner",
        )
        .drop_nulls(subset=["chrom", "pos"])
    )
    if coords.is_empty():
        raise ValueError("No reference CpGs resolved to genomic coordinates.")
    logger.info(
        "deconvolve: resolved %d reference CpGs across %d cell types.",
        coords.height,
        len(cell_types),
    )

    beta, sample_ids, _cpg_ids, _miss = _build_sample_beta_at_clock_cpgs(
        md,
        coords.select(["cpg_id", "chrom", "pos"]),
    )
    R = coords.select(cell_types).to_numpy().astype(np.float64)  # (n_cpgs, K)

    # Per-sample NNLS. Drop CpGs the sample doesn't cover.
    from scipy.optimize import nnls

    rows = []
    for s_idx, sample in enumerate(sample_ids):
        sample_beta = beta[s_idx]
        mask = ~np.isnan(sample_beta)
        if mask.sum() < R.shape[1] + 5:
            logger.warning(
                "deconvolve: sample %s has only %d usable CpGs for %d cell "
                "types; result will be unreliable.",
                sample,
                int(mask.sum()),
                R.shape[1],
            )
        if mask.sum() == 0:
            pi = np.full(R.shape[1], np.nan, dtype=np.float64)
        else:
            R_sub = R[mask]
            y_sub = sample_beta[mask]
            try:
                pi, _resid = nnls(R_sub, y_sub)
                total = pi.sum()
                if total > 0:
                    pi = pi / total
            except Exception as exc:
                logger.warning(
                    "deconvolve: NNLS failed for sample %s: %s",
                    sample,
                    exc,
                )
                pi = np.full(R.shape[1], np.nan, dtype=np.float64)
        for ct, p in zip(cell_types, pi, strict=True):
            rows.append(
                {
                    "sample_id": sample,
                    "cell_type": ct,
                    "proportion": float(p),
                }
            )
    return pl.DataFrame(rows)


__all__ = ["age_clock", "deconvolve"]
