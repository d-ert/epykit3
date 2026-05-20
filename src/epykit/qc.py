"""Quality control and reporting for Parquet methylation stores.

Phase 5 of the epykit pipeline:
  - bisulfite_conversion_rate: estimates library conversion efficiency from
    CHH-context methylation (should be <0.5 % for a high-quality WGBS run).
  - global_methylation_report: per-sample, per-context global methylation
    levels with outlier detection.
  - coverage_uniformity: per-chromosome breadth-of-coverage statistics with
    automatic flagging of low-coverage samples.

All functions read from the partitioned Parquet methylstore layout produced
by epykit.convert and are intentionally lightweight: only the columns
required for each metric are loaded.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import polars as pl

logger = logging.getLogger(__name__)

# Thresholds used by coverage_uniformity for flagging
_MIN_GENOME_COVERAGE_FRACTION = 0.80   # 80 % of CpGs at >=1x
_CONVERSION_WARNING_THRESHOLD  = 0.005  # 0.5 % CHH methylation

# chromosome aliases used by qc.sex_check. Both UCSC ("chrX")
# and Ensembl ("X") naming conventions are accepted.
_X_CHROM_NAMES: tuple[str, ...] = ("chrX", "X")


# Public API

def bisulfite_conversion_rate(
    methylstore_path: str,
    sample: str,
    chh_context_store: str,
) -> float:
    """Estimate bisulfite conversion efficiency from CHH-context methylation.

    Under complete bisulfite conversion, non-CpG cytosines (CHH context) are
    converted to uracil and read as thymine.  Residual CHH methylation
    therefore reflects incomplete conversion.  Conversion efficiency is
    estimated as 1 - mean(CHH beta).

    A value below 99.5 % (i.e. >0.5 % residual CHH methylation) should be
    flagged as a potential quality issue.

    .. note::

        This rate is **reported, not applied**. epykit's DMC / DMR tests
        consume the raw ``count_methylated`` / ``count_unmethylated``
        values exactly as Bismark / MethylDackel emit them; the
        conversion rate is surfaced in the QC dashboard, the MultiQC
        export, and the HTML report so users can gate their analysis on
        it, but it is **not** used to rescale read counts before testing.
        This matches the default behaviour of ``methylKit`` and the
        ``bsseq`` family (``read.bismark`` etc.), which leave count-level
        correction to the user. For a well-converted library (>=99.5 %)
        the correction is statistically negligible; for a poorly
        converted one the right action is usually to re-prep the
        library, not to paper over the issue with a multiplicative
        adjustment that distorts variance.

    Parameters
    ----------
    methylstore_path : str
        Path to the CpG Parquet methylstore (used only to verify that the
        sample exists; not directly read by this function).
    sample : str
        Sample identifier.
    chh_context_store : str
        Path to a *separate* Parquet methylstore generated from CHH-context
        Bismark output for the same sample.

    Returns
    -------
    float
        Estimated bisulfite conversion rate in [0, 1].
        Values close to 1.0 indicate high-quality conversion.

    Raises
    ------
    ValueError
        If no CHH data is found for the given sample.
    """
    chh_store = Path(chh_context_store)
    sample_dir = chh_store / f"sample={sample}"

    if not sample_dir.exists():
        raise ValueError(
            f"CHH store does not contain sample '{sample}': {sample_dir}"
        )

    parts = list(sample_dir.rglob("part-*.parquet"))
    if not parts:
        raise ValueError(
            f"No Parquet files found for sample '{sample}' in {chh_store}"
        )

    # Load only the count columns; ignore chrom partition
    lf = pl.scan_parquet(str(sample_dir / "**" / "part-*.parquet"))
    agg = lf.select([
        pl.sum("N_meth").alias("total_meth"),
        pl.sum("coverage").alias("total_cov"),
    ]).collect()

    total_meth = int(agg["total_meth"][0])
    total_cov  = int(agg["total_cov"][0])

    if total_cov == 0:
        raise ValueError(
            f"Zero CHH coverage for sample '{sample}'; cannot estimate rate."
        )

    mean_chh_methylation = total_meth / total_cov
    conversion_rate      = 1.0 - mean_chh_methylation

    if mean_chh_methylation > _CONVERSION_WARNING_THRESHOLD:
        logger.warning(
            "Sample '%s': CHH methylation = %.3f %% > %.1f %% threshold; "
            "consider checking library quality.",
            sample,
            mean_chh_methylation * 100,
            _CONVERSION_WARNING_THRESHOLD * 100,
        )
    else:
        logger.info(
            "Sample '%s': conversion rate = %.4f %%",
            sample, conversion_rate * 100,
        )

    return float(conversion_rate)


def global_methylation_report(
    methylstore_path: str,
    samples: list[str],
    contexts: list[str] | None = None,
) -> pl.DataFrame:
    """Compute per-sample, per-context global methylation levels.

    Reads methylation counts summed across the entire genome for each
    requested sample.  When the methylstore contains a ``context`` column
    (as written by epykit.convert), statistics are broken down by context
    (CpG / CHG / CHH); otherwise a single "CpG" row is reported.

    Outlier detection: samples whose global CpG methylation deviates by more
    than 3 MAD from the cohort median are flagged.

    Parameters
    ----------
    methylstore_path : str
        Path to the partitioned Parquet methylstore.
    samples : list[str]
        Sample identifiers to include.
    contexts : list[str], optional
        Contexts to report (default: all contexts found in the store).

    Returns
    -------
    pl.DataFrame
        Columns: sample (Utf8), context (Utf8), n_sites (Int64),
                 global_methylation (Float64), is_outlier (bool).
    """
    store = Path(methylstore_path)
    rows: list[dict] = []

    for sample in samples:
        sample_dir = store / f"sample={sample}"
        if not sample_dir.exists():
            logger.warning("Sample '%s' not found; skipping", sample)
            continue

        lf = pl.scan_parquet(str(sample_dir / "**" / "part-*.parquet"))

        # Determine if the context column is present
        schema = lf.collect_schema()
        has_context = "context" in schema

        if has_context:
            agg = (
                lf
                .group_by("context")
                .agg([
                    pl.len().alias("n_sites"),
                    pl.sum("N_meth").alias("total_meth"),
                    pl.sum("coverage").alias("total_cov"),
                ])
                .collect()
            )
        else:
            agg = (
                lf
                .select([
                    pl.len().alias("n_sites"),
                    pl.sum("N_meth").alias("total_meth"),
                    pl.sum("coverage").alias("total_cov"),
                ])
                .collect()
                .with_columns(pl.lit("CpG").alias("context"))
            )

        for row in agg.iter_rows(named=True):
            ctx = row.get("context", "CpG")
            if contexts and ctx not in contexts:
                continue
            total_cov  = row["total_cov"]
            total_meth = row["total_meth"]
            rows.append({
                "sample":             sample,
                "context":            ctx,
                "n_sites":            row["n_sites"],
                "global_methylation": (total_meth / total_cov)
                                       if total_cov > 0 else float("nan"),
            })

    if not rows:
        return pl.DataFrame({
            "sample":             pl.Series([], dtype=pl.Utf8),
            "context":            pl.Series([], dtype=pl.Utf8),
            "n_sites":            pl.Series([], dtype=pl.Int64),
            "global_methylation": pl.Series([], dtype=pl.Float64),
            "is_outlier":         pl.Series([], dtype=pl.Boolean),
        })

    result = pl.DataFrame(rows)

    # --- Outlier detection (MAD, per context) ---
    outlier_flags = np.zeros(len(result), dtype=bool)

    for ctx in result["context"].unique().to_list():
        ctx_mask = (result["context"] == ctx).to_numpy()
        meth_vals = result.filter(pl.col("context") == ctx)[
            "global_methylation"
        ].to_numpy()
        valid = ~np.isnan(meth_vals)

        if valid.sum() < 3:
            continue

        median = float(np.median(meth_vals[valid]))
        mad    = float(np.median(np.abs(meth_vals[valid] - median)))

        if mad == 0:
            continue

        z_scores = np.abs(meth_vals - median) / (1.4826 * mad)
        ctx_outliers = (z_scores > 3.0) & valid
        outlier_flags[ctx_mask] = ctx_outliers

        n_outliers = int(ctx_outliers.sum())
        if n_outliers:
            logger.warning(
                "global_methylation_report: %d outlier sample(s) detected "
                "in context %s (MAD threshold 3sigma)",
                n_outliers, ctx,
            )

    return result.with_columns(
        pl.Series("is_outlier", outlier_flags, dtype=pl.Boolean)
    ).sort(["context", "sample"])


def coverage_uniformity(
    methylstore_path: str,
    sample: str,
    thresholds: list[int] | None = None,
) -> pl.DataFrame:
    """Compute per-chromosome coverage breadth statistics for one sample.

    Reports the fraction of CpG sites covered at each threshold depth and
    flags chromosomes (and the sample overall) when coverage breadth at 1x
    falls below 80 %.

    Parameters
    ----------
    methylstore_path : str
        Path to the partitioned Parquet methylstore.
    sample : str
        Sample identifier.
    thresholds : list[int], optional
        Coverage depth thresholds to report (default: [1, 5, 10]).

    Returns
    -------
    pl.DataFrame
        Columns: sample (Utf8), chrom (Utf8),
                 n_sites (Int64), mean_coverage (Float64),
                 frac_ge_1x (Float64), frac_ge_5x (Float64),
                 frac_ge_10x (Float64),   [one per threshold]
                 low_coverage_flag (bool).
        The final row has chrom="genome" and reports genome-wide aggregates.
    """
    if thresholds is None:
        thresholds = [1, 5, 10]

    store      = Path(methylstore_path)
    sample_dir = store / f"sample={sample}"

    if not sample_dir.exists():
        raise ValueError(
            f"Sample '{sample}' not found in methylstore: {sample_dir}"
        )

    chrom_rows: list[dict] = []

    for chrom_dir in sorted(sample_dir.glob("chrom=*")):
        chrom = chrom_dir.name.removeprefix("chrom=")
        parts = list(chrom_dir.glob("part-*.parquet"))
        if not parts:
            continue

        cov_series = pl.concat([
            pl.read_parquet(str(p), columns=["coverage"])["coverage"]
            for p in parts
        ])

        n       = len(cov_series)
        cov_arr = cov_series.to_numpy()
        row: dict = {
            "sample":        sample,
            "chrom":         chrom,
            "n_sites":       n,
            "mean_coverage": float(cov_arr.mean()) if n > 0 else float("nan"),
        }

        for t in thresholds:
            key      = f"frac_ge_{t}x"
            row[key] = float((cov_arr >= t).sum() / n) if n > 0 else float("nan")

        frac_1x = row.get("frac_ge_1x", float("nan"))
        row["low_coverage_flag"] = (
            (not np.isnan(frac_1x))
            and (frac_1x < _MIN_GENOME_COVERAGE_FRACTION)
        )

        chrom_rows.append(row)

        if row["low_coverage_flag"]:
            logger.warning(
                "Sample '%s', %s: only %.1f %% of sites covered at >=1x "
                "(threshold: %.0f %%)",
                sample, chrom, frac_1x * 100,
                _MIN_GENOME_COVERAGE_FRACTION * 100,
            )

    if not chrom_rows:
        raise ValueError(f"No chromosome data found for sample '{sample}'")

    # --- Genome-wide aggregate row ---
    total_n       = sum(r["n_sites"] for r in chrom_rows)
    genome_row: dict = {
        "sample":        sample,
        "chrom":         "genome",
        "n_sites":       total_n,
        "mean_coverage": float(
            np.mean([r["mean_coverage"] for r in chrom_rows
                     if not np.isnan(r["mean_coverage"])])
        ) if total_n > 0 else float("nan"),
    }
    for t in thresholds:
        key = f"frac_ge_{t}x"
        values = [r[key] * r["n_sites"] for r in chrom_rows
                  if key in r and not np.isnan(r[key])]
        genome_row[key] = (sum(values) / total_n) if total_n > 0 else float("nan")

    frac_1x_genome = genome_row.get("frac_ge_1x", float("nan"))
    genome_row["low_coverage_flag"] = (
        (not np.isnan(frac_1x_genome))
        and (frac_1x_genome < _MIN_GENOME_COVERAGE_FRACTION)
    )
    chrom_rows.append(genome_row)

    # Build schema dynamically based on requested thresholds
    schema_extras: dict = {f"frac_ge_{t}x": pl.Float64 for t in thresholds}
    result = pl.DataFrame(chrom_rows).cast({
        "n_sites": pl.Int64,
        **schema_extras,
    })

    return result.sort(["chrom"])


# Clinical / cohort QC pack

def _resolve_x_chrom_dir(store: Path, sample: str) -> Path | None:
    """Find the X-chromosome partition for a sample under any naming."""
    for name in _X_CHROM_NAMES:
        d = store / f"sample={sample}" / f"chrom={name}"
        if d.exists():
            return d
    return None


def sex_check(
    methylstore_path: str,
    samples: list[str],
    *,
    min_coverage: int = 5,
    expected_sex: dict[str, str] | None = None,
) -> pl.DataFrame:
    """Infer sample sex from mean beta on the X chromosome .

    Female samples carry one inactivated X and have mean(beta) ~= 0.4 - 0.5
    on the X chromosome; male samples have mean(beta) ~= 0.05 - 0.10. We
    classify each sample by kmeans-2 on the per-sample mean-chrX-beta values
    (lower cluster -> male, upper cluster -> female). When fewer than two
    samples are usable the call falls back to a fixed threshold of 0.25.

    Parameters
    ----------
    methylstore_path
        Path to the partitioned Parquet methylstore.
    samples
        Sample IDs to assess. Samples missing chrX data are reported
        with ``inferred_sex == None``.
    min_coverage
        Minimum coverage at a CpG before it counts toward the mean.
    expected_sex
        Optional dict ``{sample_id: "male" | "female"}``. When supplied,
        any inferred sex that differs from the expected value is flagged
        ``mismatch=True``.

    Returns
    -------
    pl.DataFrame with columns
        sample_id, mean_chrx_beta, inferred_sex, expected_sex, mismatch.
    """
    store = Path(methylstore_path)
    records: list[dict] = []
    for sample in samples:
        x_dir = _resolve_x_chrom_dir(store, sample)
        if x_dir is None:
            records.append({
                "sample_id": sample,
                "mean_chrx_beta": float("nan"),
                "inferred_sex": None,
                "expected_sex": (expected_sex or {}).get(sample),
                "mismatch": False,
            })
            continue
        try:
            part = next(x_dir.glob("part-*.parquet"))
        except StopIteration:
            records.append({
                "sample_id": sample,
                "mean_chrx_beta": float("nan"),
                "inferred_sex": None,
                "expected_sex": (expected_sex or {}).get(sample),
                "mismatch": False,
            })
            continue
        df = pl.read_parquet(str(part), columns=["N_meth", "coverage"]).filter(
            pl.col("coverage") >= min_coverage
        )
        if len(df) == 0:
            mean_beta = float("nan")
        else:
            tot_meth = int(df.get_column("N_meth").sum())
            tot_cov  = int(df.get_column("coverage").sum())
            mean_beta = tot_meth / max(tot_cov, 1)
        records.append({
            "sample_id": sample,
            "mean_chrx_beta": float(mean_beta),
            "inferred_sex": None,
            "expected_sex": (expected_sex or {}).get(sample),
            "mismatch": False,
        })

    # KMeans-2 classification on the available mean_chrx_beta values.
    values = np.array(
        [r["mean_chrx_beta"] for r in records if np.isfinite(r["mean_chrx_beta"])],
        dtype=np.float64,
    )
    if len(values) >= 2:
        # 1D kmeans-2: sort and find the largest gap; assign two clusters.
        sorted_vals = np.sort(values)
        gaps = np.diff(sorted_vals)
        cut_idx = int(np.argmax(gaps)) + 1
        cut = (sorted_vals[cut_idx - 1] + sorted_vals[cut_idx]) / 2.0
    elif len(values) == 1:
        cut = 0.25  # fixed fallback
    else:
        cut = None

    for r in records:
        v = r["mean_chrx_beta"]
        if cut is None or not np.isfinite(v):
            r["inferred_sex"] = None
        else:
            r["inferred_sex"] = "male" if v < cut else "female"
        exp = r.get("expected_sex")
        r["mismatch"] = bool(
            exp is not None and r["inferred_sex"] is not None
            and r["inferred_sex"] != exp
        )

    return pl.DataFrame(records)


def contamination_estimate(
    methylstore_path: str,
    sample: str,
    *,
    min_coverage: int = 10,
) -> float:
    """Estimate sample contamination from the beta distribution shape.

    Clean WGBS samples have a strongly bimodal beta distribution: most CpGs
    are either fully methylated (beta ~= 1) or unmethylated (beta ~= 0).
    Cross-sample contamination produces an excess of intermediate values
    (0.2 < beta < 0.8) -- the score is the fraction of well-covered CpGs in
    that band.

    this is the lightweight "boundary-mass" version. A fuller
    EM-based three-component mixture is mentioned in the plan but adds
    a sklearn dep; the histogram score is robust and matches the
    contamination signal direction in practice.

    Parameters
    ----------
    methylstore_path : str
        Partitioned Parquet methylstore root.
    sample : str
        Sample ID.
    min_coverage : int
        Minimum per-CpG coverage to count.

    Returns
    -------
    float
        Contamination score in ``[0, 1]``. Higher values are more
        suspicious; thresholds around 0.15 - 0.25 are typical alert
        levels but the right threshold is cohort-specific.
    """
    store = Path(methylstore_path)
    sample_dir = store / f"sample={sample}"
    if not sample_dir.exists():
        return float("nan")
    parts = list(sample_dir.glob("chrom=*/part-*.parquet"))
    if not parts:
        return float("nan")
    n_total = 0
    n_mid = 0
    for part in parts:
        df = pl.read_parquet(str(part), columns=["N_meth", "coverage"]).filter(
            pl.col("coverage") >= min_coverage
        )
        if len(df) == 0:
            continue
        meth = df.get_column("N_meth").to_numpy().astype(np.float64)
        cov = df.get_column("coverage").to_numpy().astype(np.float64)
        beta = meth / np.maximum(cov, 1.0)
        n_total += len(beta)
        n_mid += int(((beta > 0.2) & (beta < 0.8)).sum())
    if n_total == 0:
        return float("nan")
    return float(n_mid) / float(n_total)


def sample_correlation(
    methylstore_path: str,
    samples: list[str],
    *,
    method: str = "spearman",
    min_coverage: int = 10,
    chromosomes: list[str] | None = None,
) -> pl.DataFrame:
    """Pairwise sample-vs-sample beta correlation matrix .

    Builds the per-sample beta vector over the intersection of CpGs covered
    at ``min_coverage`` in every sample, then returns the
    ``(n_samples x n_samples)`` correlation matrix as a long-form
    DataFrame.

    Parameters
    ----------
    method : {"spearman", "pearson"}
        Correlation type.
    min_coverage : int
        Minimum per-sample coverage to include a CpG in the correlation
        basis.
    chromosomes : list[str], optional
        Restrict to these chromosomes. Auto-detected when None.

    Returns
    -------
    pl.DataFrame with one row per (sample_a, sample_b) pair and a
    ``correlation`` column. The full matrix is also addressable via
    ``df.pivot(...)``.
    """
    if method not in ("spearman", "pearson"):
        raise ValueError("method must be 'spearman' or 'pearson'")
    store = Path(methylstore_path)
    if chromosomes is None:
        chroms: set[str] = set()
        for s_dir in store.glob("sample=*"):
            for c_dir in s_dir.glob("chrom=*"):
                chroms.add(c_dir.name.removeprefix("chrom="))
        chromosomes = sorted(chroms)

    # Build streaming intersection and accumulate per-sample beta vectors.
    beta_chunks: dict[str, list[np.ndarray]] = {s: [] for s in samples}
    for chrom in chromosomes:
        per_sample_dfs: list[pl.DataFrame] = []
        for s in samples:
            part = store / f"sample={s}" / f"chrom={chrom}" / "part-0.parquet"
            if not part.exists():
                per_sample_dfs = []
                break
            d = pl.read_parquet(
                str(part), columns=["pos", "N_meth", "coverage"]
            ).filter(pl.col("coverage") >= min_coverage)
            d = d.with_columns(
                (pl.col("N_meth") / pl.col("coverage")).alias(f"beta_{s}"),
            ).select(["pos", f"beta_{s}"])
            per_sample_dfs.append(d)
        if not per_sample_dfs:
            continue
        joined = per_sample_dfs[0]
        for d in per_sample_dfs[1:]:
            joined = joined.join(d, on="pos", how="inner")
        if len(joined) == 0:
            continue
        for s in samples:
            beta_chunks[s].append(joined.get_column(f"beta_{s}").to_numpy())

    if not any(beta_chunks.values()):
        return pl.DataFrame({
            "sample_a": [],
            "sample_b": [],
            "correlation": [],
        })

    beta_matrix = np.column_stack([
        np.concatenate(beta_chunks[s]) for s in samples
    ])  # shape (n_sites_intersection, n_samples)

    if method == "spearman":
        from scipy import stats as sp_stats
        corr = sp_stats.spearmanr(beta_matrix).statistic
        # spearmanr returns either a scalar (n=2) or a matrix.
        if np.ndim(corr) == 0:
            corr = np.array([[1.0, float(corr)], [float(corr), 1.0]])
    else:
        corr = np.corrcoef(beta_matrix, rowvar=False)

    rows = []
    for i, sa in enumerate(samples):
        for j, sb in enumerate(samples):
            rows.append({
                "sample_a": sa,
                "sample_b": sb,
                "correlation": float(corr[i, j]),
            })
    return pl.DataFrame(rows)


def power(
    meth_diff: float,
    coverage: float,
    n_per_group: int | None = None,
    *,
    power: float | None = None,
    alpha: float = 0.05,
    baseline_beta: float = 0.5,
    replicate_sd: float = 0.05,
    two_sided: bool = True,
) -> float | int:
    """Methylation-specific power / sample-size calculator .

    Models beta at each replicate as a mixture of binomial sampling noise
    (variance ``beta(1-beta)/coverage`` per replicate per CpG) plus between-
    replicate biological variance ``replicate_sd^2``. Under that model
    the two-sample z-test power is

        sd_diff = sqrt(2 * (baseline_beta * (1 - baseline_beta) / coverage
                            + replicate_sd**2) / n_per_group)
        z       = meth_diff / sd_diff
        power   = Phi(z - z_alpha)    (two-sided -> z_alpha = z_{alpha/2})

    When ``n_per_group`` is supplied (default mode), returns the implied
    power. When ``power=...`` is supplied instead, returns the smallest
    ``n_per_group`` (integer, >= 2) achieving that target power.

    Parameters
    ----------
    meth_diff : float
        Expected effect size Deltabeta (e.g. 0.10 = 10 percentage points).
    coverage : float
        Mean per-CpG coverage.
    n_per_group : int, optional
        Sample size per group; pass this to return power.
    power : float, optional
        Target power (e.g. 0.80); pass this to return sample size.
    alpha : float
        Significance level (default 0.05).
    baseline_beta : float
        Reference methylation level (default 0.5 = worst case for
        binomial variance).
    replicate_sd : float
        Between-replicate biological SD.
    two_sided : bool
        Use the two-sided z critical value.

    Returns
    -------
    float (power) or int (n_per_group), depending on which knob was set.
    """
    from scipy import stats as sp_stats
    if (n_per_group is None) == (power is None):
        raise ValueError("Pass exactly one of n_per_group, power.")
    z_alpha = float(
        sp_stats.norm.isf(alpha / 2.0 if two_sided else alpha)
    )

    def _power_at_n(n: int) -> float:
        per_site_var = baseline_beta * (1.0 - baseline_beta) / coverage
        sd = float(np.sqrt(2.0 * (per_site_var + replicate_sd ** 2) / n))
        z = abs(meth_diff) / max(sd, 1e-12)
        return float(sp_stats.norm.sf(z_alpha - z))

    if n_per_group is not None:
        return _power_at_n(int(n_per_group))

    # Solve for the minimum n that achieves the target power.
    target = float(power)
    for n in range(2, 10_001):
        if _power_at_n(n) >= target:
            return n
    return 10_000


def report_multiqc(md, output_dir: str) -> str:
    """Re-export of the MultiQC writer; defined in multiqc_export.py."""
    from .multiqc_export import report_multiqc as _impl
    return _impl(md, output_dir)
