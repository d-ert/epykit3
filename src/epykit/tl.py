"""High-level orchestrators for the standard WGBS analysis flow.

Public entry points: ``qc``, ``dmc``, ``dmr``, ``annotate``,
``diagnose_dmr_calling``. Each mutates a ``MethylData`` in place (except
``diagnose_dmr_calling`` which returns a dict) -- results land in
``md.obs`` / ``md.varm`` / ``md.uns``. See the module docstrings of
``dmc.py`` and ``dmr.py`` for the underlying engines.
"""

from __future__ import annotations

import gc
import logging
from typing import Any, Literal
import polars as pl

logger = logging.getLogger(__name__)

from .annotate import annotate_cpg_islands, annotate_features, _GTF_CACHE
from .dmc import (
    apply_multiple_testing_correction,
    empirical_fdr_for_dmc,
    process_chromosomes_dmc,
)
from .dmr import (
    DMR_PRESETS,  # noqa: F401 -- public re-export: `from epykit.tl import DMR_PRESETS`
    call_dmr_chain_merge,
    call_dmr_sliding_window,
    call_dmr_tile_based,
    empirical_fdr_for_dmr,
)
from .methyldata import MethylData
from .qc import bisulfite_conversion_rate, coverage_uniformity, global_methylation_report


def _auto_test(
    md: MethylData,
    design: str | None = None,
    covariates: list[str] | None = None,
    allow_n1: bool = False,
) -> str:
    """Auto-dispatcher: post-0.7.5 surface is closed to {fisher, lr}.
    fisher at n=1 (only engine that works), lr at n>=2.

    When a covariate design is supplied, we MUST use the binomial GLM path
    (``"glm"``) because the closed-form ``lr`` path doesn't admit
    covariates. The choice is therefore unconditional whenever the user
    asks for adjustment.

    ``allow_n1`` is forwarded to :func:`_auto_test_simple` and only takes
    effect when there are fewer than 2 replicates per group.
    """
    if design is not None or (covariates is not None and len(covariates) > 0):
        return "glm"
    return _auto_test_simple(md, allow_n1=allow_n1)


# One-shot warning gate: ``tl.dmc`` emits a UserWarning the first time a
# user explicitly selects ``test="fisher"`` in a session. We don't want to
# spam them across thousands of chromosomes/sites -- once is enough.
_FISHER_WARNED = False
_DMR_FDR_NOTED = False


def _note_dmr_fdr_calibration_once() -> None:
    """One-shot note that asymptotic DMR q-values are not a calibrated FDR.

    The region ``combined_qvalue`` / tile ``qvalue`` are BH corrections of
    Stouffer-combined per-CpG p-values; under CpG spatial correlation they
    are anti-conservative (a ranking signal, not a calibrated region-level
    FDR). Emitted at INFO (not a warning) so it informs interactive users
    without polluting stderr. Fires only when ``empirical_fdr=False``.
    """
    global _DMR_FDR_NOTED
    if _DMR_FDR_NOTED:
        return
    _DMR_FDR_NOTED = True
    logger.info(
        "DMR q-values (combined_qvalue / tile qvalue) are asymptotic BH on "
        "Stouffer-combined p-values: a well-calibrated *ranking* signal, but "
        "anti-conservative as a region-level FDR under CpG spatial "
        "correlation. For calibrated region inference pass empirical_fdr=True "
        "(tile / sliding_window) and threshold empirical_qvalue."
    )


def _warn_fisher_once() -> None:
    """Emit a one-shot UserWarning when the user explicitly picks fisher."""
    global _FISHER_WARNED
    if _FISHER_WARNED:
        return
    _FISHER_WARNED = True
    import warnings
    warnings.warn(
        "test='fisher' ignores between-replicate variance; p-values are "
        "anti-conservative. Prefer test='lr' at n >= 2.",
        UserWarning, stacklevel=3,
    )


def _check_n1_and_union_footgun(
    md: MethylData,
    allow_n1: bool,
    min_samples_treatment: int,
    min_samples_control: int,
    unit: str = "sites",
) -> None:
    """Enforce n>=2 per group (unless allow_n1) and warn on union+0/0."""
    if min(len(md.treatment_ids), len(md.control_ids)) < 2 and not allow_n1:
        _auto_test_simple(md, allow_n1=False)  # raises ValueError
    unite_info = md.uns.get("unite")
    if (
        unite_info is not None
        and unite_info.get("type") == "union"
        and min_samples_treatment == 0
        and min_samples_control == 0
    ):
        import warnings
        warnings.warn(
            f"unite='union' with min_samples_treatment=min_samples_control=0 "
            f"will test {unit} covered in only one sample per group. "
            f"Recommended: both >= 2 (or unite='intersect').",
            UserWarning, stacklevel=3,
        )


def _auto_test_simple(md: MethylData, allow_n1: bool = False) -> str:
    """Pick a sensible test based on group size.

    Current default at n>=2: ``"lr"`` -- the quasi-binomial likelihood-ratio
    chi-square with per-site McCullagh-Nelder dispersion. Closed-form on
    the streaming (S0_g, S1_g, Sigmam^2/n_g) accumulators we already keep for
    the score test. LR is closer to nominal type-I error than the score
    test at the small samples (n=6) and boundary proportions typical in
    WGBS.

    The default returned here MUST match the CLI ``--test`` default (lr) and
    the ``--test`` default for ``dmr`` (lr). See cli.py for the single source
    of truth.

    At n=1 (single replicate per group) there is no between-replicate
    variability for phi to estimate. By default this is treated as a hard error
    (statistical inference is not credible). Pass ``allow_n1=True`` to opt
    into the Fisher exact fallback (anti-conservative; warns at runtime).
    """
    min_group = min(len(md.treatment_ids), len(md.control_ids))
    if min_group < 2:
        if not allow_n1:
            raise ValueError(
                f"At least 2 replicates per group are required for valid "
                f"statistical inference (got treatment={len(md.treatment_ids)}, "
                f"control={len(md.control_ids)}). To proceed anyway with "
                f"Fisher exact on pooled reads (no between-replicate variance), "
                f"pass allow_n1=True to ep.tl.dmc(). Be aware p-values from "
                f"this path are anti-conservative and should not be reported "
                f"as evidence of differential methylation."
            )
        import warnings
        warnings.warn(
            "n<2 per group: falling back to Fisher exact on pooled reads. "
            "Between-replicate variance is ignored and p-values are anti-conservative.",
            UserWarning,
            stacklevel=3,
        )
        return "fisher"
    return "lr"


def _resolve_tsv_output(
    tsv: str | None,
    csv: str | None,
    *,
    tsv_full: bool = False,
    csv_full: bool = False,
    tsv_alpha: float = 0.05,
    csv_alpha: float = 0.05,
) -> tuple[str | None, bool, float]:
    """Resolve the ``tsv*`` output kwargs, honouring the deprecated ``csv*`` aliases.

    ``csv`` / ``csv_full`` / ``csv_alpha`` are deprecated aliases for ``tsv`` /
    ``tsv_full`` / ``tsv_alpha``. epykit writes tab-delimited TSV by default
    (pass a path ending in ``.csv`` for comma output), so the ``csv*`` names
    were misleading. Passing any ``csv*`` value emits a ``DeprecationWarning``
    and is treated as the matching ``tsv*`` value -- the new ``tsv*`` name wins
    if both are given. Returns ``(path, full, alpha)``.
    """
    if csv is not None or csv_full or csv_alpha != 0.05:
        import warnings
        warnings.warn(
            "The `csv` / `csv_full` / `csv_alpha` arguments are deprecated "
            "aliases for `tsv` / `tsv_full` / `tsv_alpha` and will be removed in "
            "a future release. epykit writes tab-delimited TSV by default; pass "
            "a path ending in `.csv` for comma-separated output.",
            DeprecationWarning,
            stacklevel=3,
        )
    path = tsv if tsv is not None else csv
    full = tsv_full or csv_full
    alpha = tsv_alpha if tsv_alpha != 0.05 else csv_alpha
    return path, full, alpha


def _resolve_auto_tsv(
    md: MethylData,
    tsv,
    csv,
    *,
    default_name: str,
    tsv_full: bool = False,
    csv_full: bool = False,
    tsv_alpha: float = 0.05,
    csv_alpha: float = 0.05,
) -> tuple[str | None, bool, float, bool]:
    """Resolve TSV output for the **auto-emit** analyses (dmc / dmr / annotate).

    Unlike :func:`_resolve_tsv_output` (opt-in), these write a human-readable
    table by default: with no explicit path the table is auto-emitted to
    ``<md.analysis_root>/results/<default_name>`` -- the same place ``md.save``
    writes. ``tsv`` semantics:

    * ``None`` / ``True`` -> auto-emit (the default)
    * ``False``           -> don't write
    * ``"path.tsv"``      -> write there (``.csv`` suffix -> comma-delimited)

    ``csv`` / ``csv_full`` / ``csv_alpha`` are deprecated aliases (warn).
    Auto-emit is globally suppressed by ``EPYKIT_NO_AUTO_TSV`` and skipped
    silently when there is no ``analysis_root`` to anchor on (e.g. an in-memory
    MethylData built without ``store_dir``).

    Returns ``(path_or_None, full, alpha, is_auto)``. ``is_auto`` is True only
    when the path was auto-derived, so the caller can make the auto write
    best-effort while letting an explicit path surface errors.
    """
    if csv is not None or csv_full or csv_alpha != 0.05:
        import warnings
        warnings.warn(
            "The `csv` / `csv_full` / `csv_alpha` arguments are deprecated "
            "aliases for `tsv` / `tsv_full` / `tsv_alpha` and will be removed "
            "in a future release. epykit writes tab-delimited TSV by default; "
            "pass a path ending in `.csv` for comma-separated output.",
            DeprecationWarning,
            stacklevel=3,
        )
    full = tsv_full or csv_full
    alpha = tsv_alpha if tsv_alpha != 0.05 else csv_alpha

    explicit = tsv if isinstance(tsv, str) else (csv if isinstance(csv, str) else None)
    if explicit is not None:
        return explicit, full, alpha, False
    if tsv is False:
        return None, full, alpha, False

    import os
    if os.environ.get("EPYKIT_NO_AUTO_TSV") in ("1", "true", "True"):
        return None, full, alpha, False
    root = getattr(md, "analysis_root", None)
    if not root:
        return None, full, alpha, False
    from pathlib import Path
    return str(Path(root) / "results" / default_name), full, alpha, True


def _emit_result_tsv(write_fn, path: str, *, is_auto: bool) -> None:
    """Run an export writer, swallowing failures only for the auto-emit path.

    An explicit ``tsv=`` path surfaces errors (the user asked for the file); an
    auto-emitted default is best-effort, so a write hiccup never aborts the
    analysis -- the results are already on ``md`` regardless.
    """
    try:
        write_fn()
    except Exception as exc:  # noqa: BLE001 -- auto-emit must not break the run
        if not is_auto:
            raise
        logger.warning("Auto TSV emit to %s skipped: %s", path, exc)


def qc(
    md: MethylData,
    chh_context_store: str | None = None,
    *,
    run_sex_check: bool = False,
    run_contamination: bool = False,
    run_sample_correlation: bool = False,
    correlation_method: str = "spearman",
    expected_sex_col: str | None = None,
    tsv: str | None = None,
    csv: str | None = None,  # deprecated alias for tsv (epykit writes TSV)
) -> None:
    """Populate md.obs with per-sample QC metrics and cache QC tables in md.uns.

    additions are opt-in via the ``run_*`` flags so the default
    ``tl.qc(md)`` keeps the existing fast subset.
    """
    samples = md.obs.get_column("sample_id").to_list()

    global_report = global_methylation_report(md.store, samples)
    cov_reports = [coverage_uniformity(md.store, sample) for sample in samples]
    cov_report = pl.concat(cov_reports) if cov_reports else pl.DataFrame()

    obs = md.obs
    cpG_report = global_report.filter(pl.col("context") == "CpG")
    if len(cpG_report) > 0:
        obs = obs.join(
            cpG_report.select([pl.col("sample").alias("sample_id"), "global_methylation"]),
            on="sample_id",
            how="left",
        )

    if len(cov_report) > 0:
        cov_genome = cov_report.filter(pl.col("chrom") == "genome")
        obs = obs.join(
            cov_genome.select([
                "sample",
                "mean_coverage",
                "frac_ge_1x",
                "frac_ge_5x",
                "frac_ge_10x",
                "low_coverage_flag",
            ]).rename({"sample": "sample_id"}),
            on="sample_id",
            how="left",
        )

    if chh_context_store:
        conv = []
        for sample in samples:
            try:
                rate = bisulfite_conversion_rate(md.store, sample, chh_context_store)
            except Exception as exc:
                logger.warning(
                    "Bisulfite conversion rate unavailable for sample %s "
                    "(CHH store %s): %s",
                    sample, chh_context_store, exc,
                )
                rate = None
            conv.append({"sample_id": sample, "bisulfite_conversion_rate": rate})
        obs = obs.join(pl.DataFrame(conv), on="sample_id", how="left")

    # --- clinical QC additions ----------------------------------
    if run_sex_check:
        from .qc import sex_check as _sex_check
        expected = None
        if expected_sex_col and expected_sex_col in obs.columns:
            expected = {
                row["sample_id"]: row[expected_sex_col]
                for row in obs.iter_rows(named=True)
                if row.get(expected_sex_col) is not None
            }
        sex_df = _sex_check(md.store, samples, expected_sex=expected)
        md.uns["qc_sex_check"] = sex_df
        obs = obs.join(
            sex_df.select(["sample_id", "inferred_sex", "mismatch"]).rename(
                {"mismatch": "sex_mismatch"}
            ),
            on="sample_id", how="left",
        )

    if run_contamination:
        from .qc import contamination_estimate as _contam
        scores = [
            {"sample_id": s, "contamination_score": float(_contam(md.store, s))}
            for s in samples
        ]
        obs = obs.join(pl.DataFrame(scores), on="sample_id", how="left")

    if run_sample_correlation:
        from .qc import sample_correlation as _samp_corr
        corr_df = _samp_corr(md.store, samples, method=correlation_method)
        md.uns["qc_sample_correlation"] = corr_df
        if len(corr_df) > 0:
            # Per-sample min off-diagonal correlation (low -> likely swap).
            off_diag = corr_df.filter(pl.col("sample_a") != pl.col("sample_b"))
            min_corr = (
                off_diag.group_by("sample_a")
                .agg(pl.min("correlation").alias("min_pairwise_corr"))
                .rename({"sample_a": "sample_id"})
            )
            obs = obs.join(min_corr, on="sample_id", how="left")

    md.obs = obs
    md.uns["qc_global_methylation"] = global_report
    md.uns["qc_coverage_uniformity"] = cov_report

    tsv, _, _ = _resolve_tsv_output(tsv, csv)
    if tsv is not None:
        from .export import qc_to_tsv
        qc_to_tsv(md, tsv)


def dmc(
    md: MethylData,
    test: str = "auto",
    chromosomes: list[str] | None = None,
    min_samples_treatment: int | None = None,
    min_samples_control: int = 0,
    dispersion: str = "eb",
    reference: str = "adaptive",
    allow_n1: bool = False,
    # Section 1 of multi-group / continuous-covariate contrasts
    formula: str | None = None,
    contrast=None,
    covariates: list[str] | None = None,
    treatment_col: str = "treatment",
    # permutation-based empirical FDR (binary path only) ----------
    empirical_fdr: bool = False,
    n_perm: int = 100,
    perm_seed: int = 42,
    perm_n_jobs: int = 1,
    *,
    backend: str = "sequential",
    n_workers: int | None = None,
    glm_backend: str = "cpu",
    resumable: bool = False,
    # Memory contract (since 1.0.1): materialize=True assembles the full
    # per-CpG result onto md.varm (the historical behaviour, ~700 MB at
    # 22M CpGs); materialize=False keeps only the streaming DMCStore handle
    # so peak memory stays O(largest chromosome) end-to-end.
    materialize: bool = True,
    use_smoothed: bool = False,
    smoothing: bool = False,
    smoothing_span_bp: int = 500,
    # FDR procedure (since 0.7.1) ----------------------------------------
    fdr_method: str = "fdr_bh",
    # Neighbour-aware p-value combining (since 0.7.1) --------------------
    neighbour_combine: bool = False,
    neighbour_bp: int = 500,
    # Separation-aware Fisher fallback (since 0.7.1) ---------------------
    sep_fallback: bool = False,
    sep_threshold: float = 0.9,
    # Power stack convenience (since 0.7.2, updated 1.0) ------------------
    power_stack: Literal["auto", "lr+", "conservative", "off"] | bool = "off",
    # Explicit patsy reference level for categorical factors (since 0.7.5) -
    reference_level: str | None = None,
    # TSV output (since 1.0). Auto-emits significant DMCs to
    # <analysis_root>/results/dmc.significant.tsv by default; tsv=False disables,
    # tsv="path" overrides. csv* are deprecated aliases (epykit writes TSV).
    tsv: str | bool | None = None,
    tsv_full: bool = False,
    tsv_alpha: float = 0.05,
    csv: str | None = None,
    csv_full: bool = False,
    csv_alpha: float = 0.05,
) -> None:
    """Run DMC calling and store result in md.varm['dmc_<test>'].

    Smoothed-input mode (``use_smoothed=True``)
    -------------------------------------------
    Requires that ``ep.pp.smooth(md)`` has been run first (either method,
    ``"gaussian"`` or ``"bsmooth"``). The DMC test then runs on
    pseudo-counts derived from the smoothed beta values:

      ``N_meth' = round(beta_smooth * coverage)``,
      ``N_unmeth' = coverage - N_meth'``

    rather than the raw per-CpG counts. All 8 test backends work
    unchanged. Results land in ``md.varm["dmc_<test>_smoothed"]`` so the
    smoothed and raw runs don't collide in the same session.

    .. warning::
       This is **not** equivalent to DSS's ``DMLfit.multiFactor(smoothing=TRUE)``.
       DSS uses BSmooth-smoothed estimates only in the per-CpG **dispersion**
       step; the per-CpG GLM still fits raw counts. The pseudo-count
       approach here is more aggressive -- it replaces the count signal
       entirely with the locally-averaged version, which can wash out
       genuine per-CpG signal at default BSmooth parameters (ns=70,
       h_bp=1000). For replicating DSS-style analyses, prefer
       ``use_smoothed=False`` (the default) with ``test="lr"`` -- the
       quasi-binomial LR with McCullagh-Nelder dispersion is the closest
       per-CpG match to DSS's count model in epykit. Reach for
       ``use_smoothed=True`` only when you genuinely want a strongly
       regularised local-mean test, with loosened smoother parameters
       (e.g. ``ns=20``).

    Backend selection
    -----------------
    ``backend="sequential"`` (default) walks chromosomes one at a time on
    the calling process -- bit-identical to the pre-0.4 behaviour.
    ``backend="dask"`` and ``backend="ray"`` submit one task per
    chromosome to a worker pool; both require the corresponding optional
    extra (``pip install 'epykit[distributed]'`` for Dask, ``epykit[ray]``
    for Ray). ``n_workers`` controls the pool size; ``None`` lets the
    backend pick a default.

    ``glm_backend="cpu"`` (default) runs the batched IRLS used by
    ``test="glm"`` / formula contrasts on the CPU via numpy. Setting
    ``glm_backend="gpu"`` routes the IRLS through CuPy (requires
    ``pip install 'epykit[gpu]'``) -- only affects the GLM hot path; the
    closed-form ``lr`` / ``score`` tests stay CPU-only.

    ``resumable=True`` (default False) participates in the 0.4.0
    pipeline manifest: when an earlier ``tl.dmc`` call with the same
    inputs + params is recorded in ``<store>/.epykit_manifest.json``,
    the prior result is loaded from its sidecar parquet and the
    computation is skipped. Inputs include the methylstore fingerprint,
    treatment / control sample lists, ``test``, ``formula``,
    ``contrast``, ``covariates``, ``min_samples_*``, ``dispersion``,
    ``reference``. Set False (the default) to preserve pre-0.4
    behaviour -- no manifest read, no skip, no sidecar write.

    Memory contract (``materialize``)
    ---------------------------------
    DMC *computation* is O(largest chromosome): the engine streams
    chromosome-by-chromosome and BH correction rewrites per-chrom parquet
    in place. ``materialize=True`` (the default) then assembles the full
    per-CpG table onto ``md.varm`` for plotting / report / export (a
    whole-genome operation; ~700 MB at 22M CpGs). Pass
    ``materialize=False`` to keep peak memory at O(largest chromosome)
    end-to-end: the result lives only as the on-disk ``DMCStore`` (reachable
    via ``md.dmc_store`` and streamed by ``tl.dmr``); ``md.dmc`` then
    materialises on demand when accessed. ``materialize=False`` is
    incompatible with the eager-only post-processors ``neighbour_combine``,
    ``empirical_fdr`` and ``use_smoothed`` (they raise ``ValueError``), and
    suppresses TSV auto-export.

    Parameters
    ----------
    md : MethylData
        Analysis object containing the methylstore path and the
        treatment/control sample lists.
    test : str
        One of ``"auto"``, ``"lr"``, ``"welch_t"``, ``"fisher"``, or
        ``"glm"``. ``"auto"`` resolves to ``"fisher"`` at n<2 and ``"lr"``
        (the recommended default) at n>=2.

        Engines removed in 0.7.5 (raise ``ValueError`` with a migration
        hint): ``"logit_t"`` (use ``"welch_t"``), ``"bb_lr"`` (use
        ``"lr"``), ``"score"`` (use ``"lr"``), ``"cmh"`` (use
        ``formula='~ group + batch'``).

        When ``formula`` and/or ``contrast`` are supplied, the test is
        forced to a GLM-based path regardless of ``test=``.
    formula : str, optional
        patsy formula on ``md.obs`` columns, e.g. ``"~ group"`` for a
        multi-group test or ``"~ age + sex"`` for a continuous-covariate
        primary effect. When supplied with ``contrast``, the engine fits
        the GLM once per site and runs a Wald / joint-F test against the
        contrast.
    contrast : str or np.ndarray, optional
        Contrast specification. Accepts:
        - a column name in the resolved design (``"age"`` for a continuous
          covariate primary effect; produces a single-coef Wald-z^2 test
          with meth-scale CIs);
        - a factor name (``"group"``); every dummy of that factor is
          included -> joint F-test (multi-group);
        - a patsy linear-combination string
          (``"group[T.KO] - group[T.WT]"``); produces a single-row contrast;
        - a raw ``(k, p)`` matrix.
    covariates : list[str], optional
        Convenience list of column names to include as nuisance terms.
        Combined with ``formula`` and the resolved ``treatment_col``.
    treatment_col : str, default ``"treatment"``
        Name of the binary 0/1 column in ``md.obs`` used by the legacy
        binary path. Ignored when ``contrast`` is supplied and resolves
        without it.
    dispersion : {"site", "chrom", "shrink", "eb"}
        McCullagh-Nelder dispersion strategy used by the ``"lr"`` test.
        Default ``"eb"`` shrinks the per-site Pearson
        residual estimate toward a chromosome-wide pool via empirical-Bayes
        weights (stable at low n / low coverage). Alternatives:
        ``"site"`` uses the noisy per-site estimate only; ``"chrom"`` uses
        the chromosome-pooled phi for every site; ``"shrink"`` is a James-Stein-style weighted average of per-site and
        chromosome estimates with a fixed pseudo-df weight. See
        :func:`_score_finalize` in ``dmc.py`` for the math.
    chromosomes : list[str], optional
        Restrict to a subset of chromosomes. Auto-detected when None.
    min_samples_treatment, min_samples_control : int
        per-site minimum number of samples with non-zero coverage in
        each group. Sites that fail are NaN'd out before FDR correction.
        Primarily useful when ``ep.pp.unite(..., type="union")`` was used so
        that union-introduced zero-coverage rows aren't treated as real
        observations.

    """
    if min_samples_treatment is None:
        min_samples_treatment = 0

    tsv, tsv_full, tsv_alpha, _tsv_is_auto = _resolve_auto_tsv(
        md, tsv, csv, default_name="dmc.significant.tsv",
        tsv_full=tsv_full, csv_full=csv_full,
        tsv_alpha=tsv_alpha, csv_alpha=csv_alpha,
    )

    if test == "logit_t":
        raise ValueError(
            "test='logit_t' was removed in 0.7.5 (miscalibrated near β=0/1). "
            "Use test='welch_t' for the replicate-aware β-mean test or "
            "test='lr' for the recommended default."
        )

    if test == "bb_lr":
        raise ValueError(
            "test='bb_lr' was removed in 0.7.5 (TPR < 8% at n ≤ 4 + a "
            "dispersion-df bug). Use test='lr' (recommended) which uses "
            "the same quasi-binomial dispersion but pools counts per group "
            "for higher power at small n."
        )

    if test == "score":
        raise ValueError(
            "test='score' was removed in 0.7.5 (strictly dominated by "
            "test='lr' in finite samples; asymptotically equivalent under "
            "H0). Switch test='score' -> test='lr'; output schema is "
            "identical."
        )

    if test == "cmh":
        raise ValueError(
            "test='cmh' was removed in 0.7.5 (stratification semantics "
            "confusing; dominated by GLM with batch covariate). For "
            "stratified analysis use tl.dmc(formula='~ group + batch'), "
            "which gives proper dispersion correction and handles "
            "continuous covariates."
        )

    # --- New contrast / multi-group path -------------------------------------
    if formula is not None or contrast is not None:
        if empirical_fdr:
            # Same refusal as the DMR path: label shuffling invalidates
            # the stratified design that formula= encodes.
            raise ValueError(
                "empirical_fdr=True is not supported with the contrast / "
                "multi-group DMC path (label shuffling invalidates the "
                "stratified design). Use the binary treatment / control "
                "path or implement a custom stratified permutation."
            )
        _run_dmc_contrast(
            md, test=test, formula=formula, contrast=contrast,
            covariates=covariates, treatment_col=treatment_col,
            chromosomes=chromosomes,
            min_samples_treatment=min_samples_treatment,
            min_samples_control=min_samples_control,
            dispersion=dispersion, reference=reference,
            reference_level=reference_level,
        )
        # P1-11 deprecation notice for GLM / contrast path.
        import warnings as _warnings
        _warnings.warn(
            "The 'log2_odds_ratio' column is deprecated and is slated for "
            "removal in 1.1. Use 'log2_odds_ratio_pooled' for pooled-count "
            "tests (lr, fisher) or 'coef_treatment_log2' for the glm backend. "
            "The transitional column is NaN-filled.",
            FutureWarning, stacklevel=2,
        )
        if tsv is not None:
            from .export import dmc_to_tsv
            _emit_result_tsv(
                lambda: dmc_to_tsv(md, tsv, alpha=tsv_alpha, full=tsv_full),
                tsv, is_auto=_tsv_is_auto,
            )
        return

    # Unconditional n=1 guard: applies whether test is "auto" or explicit.
    # _auto_test_simple raises ValueError when allow_n1=False; trigger that
    # check up front so explicit test="lr"/"fisher" with n<2 also gets
    # caught instead of silently running on degenerate data.
    _check_n1_and_union_footgun(
        md, allow_n1, min_samples_treatment, min_samples_control,
    )
    selected_test = _auto_test(md, allow_n1=allow_n1) if test == "auto" else test

    # --- lr+ power-stack dispatch (1.0) ---
    # The four knobs the stack controls:
    #   neighbour_combine, fdr_method, sep_fallback, dispersion ("eb"
    #   is already the default in this function).
    # power_stack values:
    #   "lr+" / True / "auto"  -> engage all four at any n
    #   "conservative"          -> engage only at n <= 2 (pre-1.0 behavior)
    #   "off" / False           -> leave knobs at user-passed values
    if isinstance(power_stack, bool):
        power_stack = "lr+" if power_stack else "off"
    if power_stack not in {"auto", "lr+", "conservative", "off"}:
        raise ValueError(
            f"power_stack must be one of {{'auto','lr+','conservative','off'}} "
            f"or a bool; got {power_stack!r}"
        )

    if selected_test == "lr" and power_stack in {"auto", "lr+", "conservative"}:
        _min_n = min(len(md.treatment_ids), len(md.control_ids))
        engage = (power_stack != "conservative") or (_min_n <= 2)
        if engage:
            if not neighbour_combine:
                neighbour_combine = True
                logger.info(
                    "Auto-enabling neighbour_combine (lr+ stack, "
                    "power_stack=%s, n=%d). Pass power_stack='off' to "
                    "disable.", power_stack, _min_n,
                )
            if fdr_method == "fdr_bh":
                fdr_method = "fdr_tsbh"
                logger.info(
                    "Auto-switching fdr_method 'fdr_bh' -> 'fdr_tsbh' "
                    "(lr+ stack, power_stack=%s).", power_stack,
                )
            if not sep_fallback:
                sep_fallback = True
                logger.info(
                    "Auto-enabling sep_fallback (lr+ stack, "
                    "power_stack=%s).", power_stack,
                )

    if selected_test == "fisher":
        _warn_fisher_once()
    unite_info = md.uns.get("unite")
    unite = (unite_info is not None) and (unite_info.get("type") == "intersect")

    # materialize=False keeps only the streaming DMCStore handle, so it
    # cannot run the eager-only post-processors that need the full in-memory
    # result table. Reject the combination explicitly rather than silently
    # producing different output. (neighbour_combine may have been
    # auto-enabled just above by power_stack='lr+'; use_smoothed routes
    # through a temp store that is cleaned up on return, so its DMCStore
    # would not survive.)
    if not materialize:
        _incompat = [
            name for name, on in (
                ("neighbour_combine", neighbour_combine),
                ("empirical_fdr", empirical_fdr),
                ("use_smoothed", use_smoothed),
            ) if on
        ]
        if _incompat:
            raise ValueError(
                "materialize=False keeps only the streaming DMCStore handle "
                "and cannot run features that post-process the full in-memory "
                f"result table: {', '.join(_incompat)}. Re-run with "
                "materialize=True (the default) for these, or disable them "
                "(neighbour_combine may have been auto-enabled by "
                "power_stack='lr+')."
            )

    # 0.4.0 checkpoint/resume: when resumable=True, fingerprint the input
    # + params and look up a prior matching run in the pipeline manifest.
    resume_root = None
    resume_sig = None
    resume_stage_name = None
    if resumable:
        from .dmc import _canonicalise_test_name as _canon
        from ._cache import input_signature, manifest_find, manifest_append
        from pathlib import Path
        canonical_for_key = _canon(selected_test)
        resume_stage_name = f"dmc_{canonical_for_key}"
        resume_root = md.analysis_root or md.store
        resume_sig = input_signature(
            md.store,
            sorted(md.treatment_ids),
            sorted(md.control_ids),
            {
                "test": selected_test,
                "chromosomes": chromosomes,
                "unite": unite,
                "min_samples_treatment": min_samples_treatment,
                "min_samples_control": min_samples_control,
                "dispersion": dispersion,
                "reference": reference,
                "empirical_fdr": empirical_fdr,
                "n_perm": n_perm if empirical_fdr else None,
                "perm_seed": perm_seed if empirical_fdr else None,
                # lr+ stack knobs: changing any of these changes engine
                # output, so a cached result computed at different values
                # must NOT be reused. Omitting any of these from the
                # signature causes silent cache hits in parameter sweeps.
                "power_stack": power_stack,
                "sep_fallback": sep_fallback,
                "sep_threshold": sep_threshold,
                "neighbour_combine": neighbour_combine,
                "neighbour_bp": neighbour_bp,
                "fdr_method": fdr_method,
            },
        )
        if resume_root:
            prior = manifest_find(resume_root, resume_stage_name)
            if prior is not None and prior.get("input_sig") == resume_sig:
                sidecar = Path(prior["output_path"])
                if not sidecar.is_absolute():
                    sidecar = Path(resume_root) / sidecar
                if sidecar.exists():
                    import logging as _lg
                    _lg.getLogger(__name__).info(
                        "[resume] %s: loading cached result from %s",
                        resume_stage_name, sidecar,
                    )
                    md.varm[resume_stage_name] = pl.read_parquet(str(sidecar))
                    md.uns["dmc"] = {
                        "test_requested": test,
                        "test_used": canonical_for_key,
                        "n_sites": len(md.varm[resume_stage_name]),
                        "unite": unite,
                        "min_samples_treatment": min_samples_treatment,
                        "min_samples_control": min_samples_control,
                        "dispersion": dispersion,
                        "reference": reference,
                        "empirical_fdr": empirical_fdr,
                        "n_perm": n_perm if empirical_fdr else None,
                        "perm_seed": perm_seed if empirical_fdr else None,
                        "power_stack": power_stack,
                        "sep_fallback": sep_fallback,
                        "sep_threshold": sep_threshold,
                        "neighbour_combine": neighbour_combine,
                        "neighbour_bp": neighbour_bp,
                        "fdr_method": fdr_method,
                        "last_key": resume_stage_name,
                        "resumed": True,
                    }
                    if tsv is not None:
                        from .export import dmc_to_tsv
                        _emit_result_tsv(
                            lambda: dmc_to_tsv(md, tsv, alpha=tsv_alpha, full=tsv_full),
                            tsv, is_auto=_tsv_is_auto,
                        )
                    return

    # Smoothed-input mode: materialise a temp methylstore of smoothed
    # pseudo-counts and route the DMC engine at it. Cleans up automatically
    # when the with-block exits.
    import tempfile as _tempfile
    from pathlib import Path as _Path
    if use_smoothed:
        import warnings as _warnings
        _warnings.warn(
            "use_smoothed=True (pseudo-count transform of raw reads via "
            "BSmooth) is NOT equivalent to DSS's smoothing=TRUE -- it's "
            "too aggressive (washes out per-CpG signal at default BSmooth "
            "parameters). For DSS-style behavior, use smoothing=True "
            "(applies DSS's uniform-box +/-smoothing_span_bp//2 moving "
            "average to each sample's raw counts before they hit the "
            "test, matching DMLfit.multiFactor(smoothing=TRUE)). The "
            "use_smoothed pseudo-count path will be removed in a future "
            "minor release.",
            DeprecationWarning, stacklevel=2,
        )
        if "smooth_path" not in md.uns:
            raise ValueError(
                "use_smoothed=True requires ep.pp.smooth(md) first "
                "(either method='gaussian' or 'bsmooth'). The smoothed "
                "sidecar at md.uns['smooth_path'] is the input to the "
                "pseudo-count transform that feeds the DMC test."
            )
        smooth_path = md.uns["smooth_path"]
        _smoothed_tmp = _tempfile.TemporaryDirectory(prefix="epykit_dmc_smoothed_")
        _dmc_store = _smoothed_tmp.name
        from ._smoothed_store import build_smoothed_pseudo_count_store
        build_smoothed_pseudo_count_store(
            raw_store=_Path(md.store),
            smooth_store=_Path(smooth_path),
            samples=md.obs.get_column("sample_id").to_list(),
            out_dir=_Path(_dmc_store),
        )
    else:
        _smoothed_tmp = None
        _dmc_store = md.store

    dmc_store = None
    try:
        # Use return_store=True so the per-chrom parquet directory
        # becomes the source of truth. Both BH correction and
        # downstream DMR can then stream chromosomes from disk and
        # avoid materialising the full 22M-row table in memory.
        dmc_store = process_chromosomes_dmc(
            methylstore_path=_dmc_store,
            samples_treatment=md.treatment_ids,
            samples_control=md.control_ids,
            test=selected_test,
            chromosomes=chromosomes,
            unite=unite,
            min_samples_treatment=min_samples_treatment,
            min_samples_control=min_samples_control,
            dispersion=dispersion,
            reference=reference,
            backend=backend,
            n_workers=n_workers,
            glm_backend=glm_backend,
            return_store=True,
            smoothing=smoothing,
            smoothing_span_bp=smoothing_span_bp,
            sep_fallback=sep_fallback,
            sep_threshold=sep_threshold,
        )
        dmc_store = apply_multiple_testing_correction(dmc_store, method=fdr_method)

        if materialize:
            # Materialise the full DataFrame for md.varm back-compat
            # (plot.py / report.py / pl modules consume md.dmc as a
            # DataFrame). With chrom/strand stored as pl.Enum this is
            # roughly 700 MB at 22M rows vs. ~2 GB before -- manageable
            # alongside the per-chrom DMR working set.
            result = dmc_store.to_dataframe()
        else:
            # Keep only the streaming DMCStore handle (O(largest chromosome)
            # end-to-end). neighbour_combine / empirical_fdr were rejected
            # above, so the post-processing blocks below are inert
            # (their guards short-circuit on the False flag before touching
            # `result`). md.dmc materialises on demand from store_path.
            result = None

        # Neighbour-aware p-value combining (RADMeth-style, since 0.7.1).
        # When enabled, run the signed-Stouffer combiner over the per-CpG
        # raw p-values and emit the combined p-value as a sibling column
        # `pvalue_combined` (with `qvalue_combined` from BH/Storey on the
        # combined values, and `pvalue_combined_n_neighbours` /
        # `qvalue_combined_reject` as audit columns). The canonical
        # `pvalue` / `qvalue` columns remain the raw per-CpG values --
        # downstream consumers that want the combined values must opt in
        # by reading the `_combined` columns explicitly. Sites without
        # enough neighbours fall back to their raw p-value identity.
        if neighbour_combine and len(result) > 0:
            from .dmc import combine_neighbour_pvalues
            result = combine_neighbour_pvalues(result, neighbour_bp=neighbour_bp)
            # Keep `pvalue` / `qvalue` as the raw per-CpG values.
            # `combine_neighbour_pvalues` already added a `pvalue_combined`
            # column; produce `qvalue_combined` next to it via BH on the
            # combined p-values. Downstream consumers that want the
            # combined values must opt in by reading `pvalue_combined` /
            # `qvalue_combined`.
            result = apply_multiple_testing_correction(
                result, method=fdr_method,
                pvalue_col="pvalue_combined", qvalue_col="qvalue_combined",
            )

        if empirical_fdr and len(result) > 0:
            result = empirical_fdr_for_dmc(
                methylstore_path=_dmc_store,
                samples_treatment=md.treatment_ids,
                samples_control=md.control_ids,
                observed_dmc=result,
                n_perm=n_perm,
                seed=perm_seed,
                n_jobs=perm_n_jobs,
                test=selected_test,
                chromosomes=chromosomes,
                unite=unite,
                min_samples_treatment=min_samples_treatment,
                min_samples_control=min_samples_control,
                dispersion=dispersion,
                reference=reference,
                # M3: engine knobs that can overwrite the per-site p-value
                # MUST be applied identically in observed and null runs,
                # otherwise the Westfall-Young statistic compares deflated
                # observed p-values against an un-deflated null pool.
                sep_fallback=sep_fallback,
                sep_threshold=sep_threshold,
                smoothing=smoothing,
                smoothing_span_bp=smoothing_span_bp,
            )
    finally:
        if _smoothed_tmp is not None:
            _smoothed_tmp.cleanup()

    # Canonicalise key name (test_used reflects the canonical name post-rename)
    from .dmc import _canonicalise_test_name
    canonical_used = _canonicalise_test_name(selected_test)
    key = f"dmc_{canonical_used}_smoothed" if use_smoothed else f"dmc_{canonical_used}"
    if materialize:
        md.varm[key] = result
        n_sites = len(result)
    else:
        # No eager table on varm; the DMCStore is the source of truth.
        # md.get_dmc() / md.dmc resolve it on demand from store_path.
        n_sites = dmc_store.total_sites if dmc_store is not None else 0
    md.uns["dmc"] = {
        "test_requested": test,
        "test_used": canonical_used,
        "n_sites": n_sites,
        "materialized": bool(materialize),
        "unite": unite,
        "min_samples_treatment": min_samples_treatment,
        "min_samples_control": min_samples_control,
        "dispersion": dispersion,
        "reference": reference,
        "empirical_fdr": empirical_fdr,
        "n_perm": n_perm if empirical_fdr else None,
        "perm_seed": perm_seed if empirical_fdr else None,
        "power_stack": power_stack,
        "sep_fallback": sep_fallback,
        "sep_threshold": sep_threshold,
        "neighbour_combine": neighbour_combine,
        "neighbour_bp": neighbour_bp,
        "fdr_method": fdr_method,
        # explicit pointer so MethylData.get_dmc() / .dmc resolve to the
        # table the user just wrote, regardless of which other tests have
        # been run in the same session.
        "last_key": key,
        "use_smoothed": use_smoothed,
        "smooth_method": (
            md.uns.get("smooth_params", {}).get("method") if use_smoothed else None
        ),
        # DSS-style count smoothing (DMLfit.multiFactor(smoothing=TRUE)
        # analogue). Surface params for both modes so the metadata
        # round-trips through save / report consistently; the span is
        # only meaningful when smoothing == True.
        "smoothing": bool(smoothing),
        "smoothing_span_bp": int(smoothing_span_bp) if smoothing else None,
        # Path to the persistent per-chrom DMC store. Lets downstream
        # tools (esp. tl.dmr(method='sliding_window')) stream
        # chromosomes from disk instead of holding the materialised
        # DataFrame plus a per-chrom working set in memory at the
        # same time.
        "store_path": str(dmc_store.path) if dmc_store is not None else None,
    }

    # 0.4.0 checkpoint/resume: persist a sidecar parquet + manifest entry
    # so a subsequent resumable=True call with the same fingerprint can
    # skip the computation entirely. Best-effort: a failed write logs but
    # does not propagate (the in-memory result is still valid).
    if resumable and resume_root and resume_sig and resume_stage_name and result is not None:
        try:
            from ._cache import manifest_append
            from pathlib import Path
            sidecar_dir = Path(resume_root) / ".epykit_results"
            sidecar_dir.mkdir(parents=True, exist_ok=True)
            sidecar = sidecar_dir / f"{resume_stage_name}.parquet"
            result.write_parquet(str(sidecar))
            manifest_append(
                resume_root, resume_stage_name,
                params={
                    "test": canonical_used,
                    "unite": unite,
                    "min_samples_treatment": min_samples_treatment,
                    "min_samples_control": min_samples_control,
                    "dispersion": dispersion,
                    "reference": reference,
                    "empirical_fdr": empirical_fdr,
                },
                input_sig=resume_sig,
                output_path=str(sidecar),
                extra={"n_sites": len(result)},
            )
        except OSError as exc:
            import logging as _lg
            _lg.getLogger(__name__).warning(
                "[resume] failed to persist %s sidecar: %s",
                resume_stage_name, exc,
            )

    # P1-11 deprecation notice – emitted once per tl.dmc call (not per-row,
    # not per-chromosome).  The transitional 'log2_odds_ratio' column is
    # NaN-filled and slated for removal in 1.1.
    import warnings as _warnings
    _warnings.warn(
        "The 'log2_odds_ratio' column is deprecated and is slated for "
        "removal in 1.1. Use 'log2_odds_ratio_pooled' for pooled-count "
        "tests (lr, fisher) or 'coef_treatment_log2' for the glm backend. "
        "The transitional column is NaN-filled.",
        FutureWarning, stacklevel=2,
    )

    if tsv is not None:
        if materialize:
            from .export import dmc_to_tsv
            _emit_result_tsv(
                lambda: dmc_to_tsv(md, tsv, alpha=tsv_alpha, full=tsv_full),
                tsv, is_auto=_tsv_is_auto,
            )
        else:
            logger.info(
                "materialize=False: skipping DMC TSV auto-export (no "
                "in-memory result table). Export later via "
                "ep.export.dmc_to_tsv(md) or re-run with materialize=True."
            )


def _run_dmc_contrast(
    md: MethylData,
    *,
    test: str,
    formula: str | None,
    contrast,
    covariates: list[str] | None,
    treatment_col: str,
    chromosomes: list[str] | None,
    min_samples_treatment: int,
    min_samples_control: int,
    dispersion: str,
    reference: str,
    reference_level: str | None = None,
) -> None:
    """Internal: multi-group / continuous-covariate primary-effect DMC.

    Always uses test='glm_contrast' internally. Uses ALL samples in
    md.obs order (not the binary case/control split), so the design
    matrix matches md.obs row-for-row.
    """
    from .dmc import process_chromosomes_dmc
    from ._glm import build_design, resolve_contrast

    if not md.obs.height:
        raise ValueError("md.obs is empty; cannot build a design matrix.")
    samples_all = md.obs.get_column("sample_id").to_list()

    # Build design -- without requiring a treatment column if we have a
    # formula that doesn't reference one. The user's `treatment_col`
    # default ("treatment") is *only* required when the existing binary
    # path would have used it; here we let the formula speak.
    need_treatment = (treatment_col in md.obs.columns) and (
        formula is None or treatment_col in formula
    )
    design_full, _design_reduced, coef_idx, term_names, formula_used, design_info = (
        build_design(
            md.obs,
            samples_ordered=samples_all,
            formula=formula,
            covariates=covariates,
            treatment_col=treatment_col,
            require_treatment_col=need_treatment,
            return_design_info=True,
            reference_level=reference_level,
        )
    )

    # Resolve the contrast against the design.
    if contrast is None:
        # Default: a single-coef contrast on `treatment_col` (this happens
        # when the user supplies a `formula=` for covariate adjustment but
        # no explicit contrast).
        contrast = treatment_col
    C, contrast_label = resolve_contrast(contrast, term_names, design_info=design_info)

    # Build per-sample level labels for the multi-group output schema:
    # take the FIRST term that is a factor of `contrast` (when contrast is
    # a factor name), otherwise no per-level breakdown.
    group_labels: list[str] | None = None
    if isinstance(contrast, str) and contrast in md.obs.columns:
        # Either a continuous column (single coef) or a categorical column
        # (joint test). For both, emit per-level labels for downstream
        # mean_beta_<level> columns when the column is categorical.
        col = md.obs.get_column(contrast)
        if col.dtype == pl.Utf8 or col.dtype == pl.Categorical:
            group_labels = col.cast(pl.Utf8).to_list()

    # Determine which samples are "case" vs "control" for the
    # backwards-compatible binary columns. If treatment_col is on obs and
    # carries a numeric 0/1 signal, use it; otherwise leave both empty so
    # mean_beta_case/control remain NaN (uninterpretable for multi-group).
    samples_case_local: list[str] = []
    samples_control_local: list[str] = []
    if treatment_col in md.obs.columns:
        try:
            mask_treat = (
                md.obs.get_column(treatment_col).cast(pl.Float64, strict=False) == 1
            ).to_list()
            samples_case_local = [s for s, m in zip(samples_all, mask_treat) if m]
            samples_control_local = [s for s, m in zip(samples_all, mask_treat) if not m]
        except Exception as exc:
            logger.warning(
                "Could not derive case/control split from treatment column "
                "%r; mean_beta_case/control will be NaN: %s",
                treatment_col, exc,
            )

    unite_info = md.uns.get("unite")
    unite = (unite_info is not None) and (unite_info.get("type") == "intersect")

    dmc_store_contrast = process_chromosomes_dmc(
        methylstore_path=md.store,
        samples_treatment=samples_case_local,
        samples_control=samples_control_local,
        test="glm_contrast",
        chromosomes=chromosomes,
        unite=unite,
        min_samples_treatment=min_samples_treatment,
        min_samples_control=min_samples_control,
        dispersion=dispersion,
        reference=reference,
        design_full=design_full,
        contrast_matrix=C,
        contrast_label=contrast_label,
        samples_all_ordered=samples_all,
        group_labels_per_sample=group_labels,
        return_store=True,
    )
    dmc_store_contrast = apply_multiple_testing_correction(
        dmc_store_contrast, method="fdr_bh"
    )
    result = dmc_store_contrast.to_dataframe()

    key = "dmc_glm_contrast"
    md.varm[key] = result
    md.uns["dmc"] = {
        "test_requested": test,
        "test_used": "glm_contrast",
        "n_sites": len(result),
        "unite": unite,
        "formula": formula_used,
        "contrast": contrast_label,
        "design_terms": term_names,
        "covariates": list(covariates) if covariates else None,
        "treatment_col": treatment_col,
        "min_samples_treatment": min_samples_treatment,
        "min_samples_control": min_samples_control,
        "dispersion": dispersion,
        "reference": reference,
        "last_key": key,
        "store_path": str(dmc_store_contrast.path),
    }


def dmr(
    md: MethylData,
    method: str = "chain_merge",
    # Parameter preset bundle (chain_merge only; see DMR_PRESETS) ----------
    preset: str | None = None,
    # Tile-method options ---------------------------------------------------
    tile_size_bp: int = 1000,
    min_cpgs_per_tile: int = 5,
    test: str = "auto",
    chromosomes: list[str] | None = None,
    min_samples_treatment: int | None = None,
    min_samples_control: int = 0,
    dispersion: str = "site",
    reference: str = "adaptive",
    # Covariate design (tile-method only) ----------------------------------
    design: str | None = None,
    covariates: list[str] | None = None,
    treatment_col: str = "treatment",
    # Sliding-window options ------------------------------------------------
    window_bp: int = 500,
    step_bp: int = 250,
    min_cpgs: int = 5,
    min_sites_significant: int = 3,
    # Chain-merge options (DSS callDMR semantics) --------------------------
    dis_merge_bp: int = 500,
    pct_sig: float = 0.5,
    minlen_bp: int = 50,
    use_q_for_sig: bool = False,
    # Shared filters --------------------------------------------------------
    alpha: float = 0.05,
    min_abs_meth_diff: float = 0.1,
    min_mean_qvalue: float | None = 0.05,
    # Replicate-count guard --------------------------------------------------
    allow_n1: bool = False,
    # permutation-based empirical FDR (tile method only) --------
    empirical_fdr: bool = False,
    n_perm: int = 100,
    perm_seed: int = 42,
    perm_n_jobs: int = 1,
    empirical_strata: str | None = None,
    *,
    backend: str = "sequential",
    n_workers: int | None = None,
    merge_adjacent: bool = True,
    # Auto-emits the DMR table to <analysis_root>/results/dmr.tsv by default;
    # tsv=False disables, tsv="path" overrides. csv= is a deprecated alias.
    tsv: str | bool | None = None,
    csv: str | None = None,
) -> None:
    """Run DMR calling and store result in ``md.uns['dmr']``.

    Backend selection (``method="tile"`` only)
    ------------------------------------------
    ``backend="sequential"`` (default) walks chromosomes one at a time.
    ``backend="dask"`` / ``"ray"`` parallelise the per-chromosome DMC
    pass run on the tiled methylstore; requires ``[distributed]`` or
    ``[ray]``. ``n_workers`` controls pool size. The
    ``method="sliding_window"`` path operates on the in-memory DMC table
    and is unaffected by this setting.

    Four methods are supported (see ``method`` below for the full list):

    * ``method="chain_merge"`` (default, recommended) -- DSS-style
      chaining of significant CpGs into regions; tune via ``preset``.
      The recommended caller for most WGBS analyses.
    * ``method="tile"`` -- aggregates read counts within fixed tiles and
      runs a single test per tile. Requires direct access to ``md.store``
      and the per-sample methylstore; does not need a prior DMC table.
    * ``method="segment"`` -- rule-based 3-state segmentation over the
      ``meth_diff`` signal with Stouffer-combined per-segment p-values.
    * ``method="sliding_window"`` -- the legacy in-tree method: takes the
      DMC result on ``md`` and combines per-CpG p-values within overlapping
      windows with signed Stouffer's Z. Faster (no extra I/O) but
      substantially lower-power than tile-based aggregation at typical
      WGBS coverage.

    Region-level FDR calibration (important)
    ----------------------------------------
    The DMR-level ``combined_qvalue`` (chain_merge / sliding_window /
    segment) and tile ``qvalue`` are BH corrections of **asymptotic** region
    p-values built by combining per-CpG p-values with signed Stouffer's Z.
    Stouffer's combination assumes the per-CpG statistics are independent;
    adjacent WGBS CpGs are positively correlated, so the combined null
    variance exceeds 1 and these region q-values are **anti-conservative**
    (too small) in CpG-dense regions. Treat them as a well-calibrated
    *ranking* signal, **not** a calibrated region-level FDR. For trustworthy
    region-level inference pass ``empirical_fdr=True`` (``method="tile"`` or
    ``"sliding_window"``), which re-runs the engine on shuffled labels and
    adds permutation ``empirical_pvalue`` / ``empirical_qvalue`` columns;
    threshold ``empirical_qvalue`` instead of ``combined_qvalue`` for FDR
    control. (Running without ``empirical_fdr`` logs a one-time INFO note to
    this effect.)

    Parameters
    ----------
    method : {"tile", "sliding_window", "segment", "chain_merge"}
        Which DMR algorithm to run.
    preset : {"strict", "default", "permissive"}, optional
        Parameter preset bundle for ``method="chain_merge"``. Applies
        ``(alpha, min_abs_meth_diff, dis_merge_bp, min_cpgs, pct_sig,
        minlen_bp)`` from :data:`epykit.tl.DMR_PRESETS` for the chosen
        bundle. Any explicit kwarg passed alongside ``preset`` overrides
        the bundled value. Ignored for other methods.

        Preset summary:

        * ``"strict"`` -- validation-ready DMRs only (alpha=1e-6,
          min_cpgs=5, min_abs_meth_diff=0.20). DSS-strict alpha.
        * ``"default"`` -- balanced (alpha=1e-4, min_abs_meth_diff=0.10).
          One order looser on alpha than DSS to capture real-but-moderate
          signal; keeps the 10% per-CpG effect-size floor. Recommended
          starting point for general WGBS analyses.
        * ``"permissive"`` -- recall-oriented (alpha=1e-4,
          dis_merge_bp=200, min_abs_meth_diff=0.05). Expect lower PPV.
    tile_size_bp, min_cpgs_per_tile : int
        Tile-method options. Default ``tile_size_bp=1000``.
    test : str
        Statistical test for tile-method (ignored when ``method="sliding_window"``).
        ``"auto"`` resolves the same way as in :func:`dmc`.
    chromosomes : list[str], optional
        Restrict tile-method processing to these chromosomes.
    min_samples_treatment, min_samples_control : int
        Per-tile sample-count guard for tile-method.
    window_bp, step_bp, min_cpgs, min_sites_significant : int
        Sliding-window method options.
    alpha : float
        q-value threshold for "significant" at the DMC / tile level
        (used by both methods, with different downstream meanings).
    min_abs_meth_diff : float
        Minimum |meth_diff| for a DMC / tile to count.
    min_mean_qvalue : float or None
        post-hoc filter on the DMR-level **q-value**
        (``combined_qvalue`` for sliding-window, ``qvalue`` for tile).
        DMRs with q >= ``min_mean_qvalue`` are dropped. Set to None to keep
        all candidate DMRs. Default 0.05.

        This parameter was previously named ``min_mean_pvalue`` and applied
        to the uncorrected p-value, which was not FDR-controlled across the
        DMR set.
    """
    if min_samples_treatment is None:
        min_samples_treatment = 0
    tsv, _, _, _tsv_is_auto = _resolve_auto_tsv(md, tsv, csv, default_name="dmr.tsv")
    if not empirical_fdr:
        # Asymptotic DMR q-values are a ranking signal, not a calibrated
        # region-level FDR under CpG spatial correlation (M5). Point users at
        # empirical_fdr=True for calibrated inference.
        _note_dmr_fdr_calibration_once()
    if method == "tile":
        _check_n1_and_union_footgun(
            md, allow_n1, min_samples_treatment, min_samples_control, unit="tiles",
        )
        selected_test = (
            _auto_test(md, design=design, covariates=covariates, allow_n1=allow_n1)
            if test == "auto"
            else test
        )
        if selected_test == "fisher":
            _warn_fisher_once()
        unite_info = md.uns.get("unite")
        unite = (unite_info is not None) and (unite_info.get("type") == "intersect")

        # ---- Build covariate design matrix when requested -----------------
        design_full = None
        design_reduced = None
        coef_idx = None
        term_names: list[str] = []
        formula_used: str | None = None
        if selected_test == "glm" or design is not None or covariates is not None:
            from ._glm import build_design
            samples_ordered = md.treatment_ids + md.control_ids
            design_full, design_reduced, coef_idx, term_names, formula_used = build_design(
                md.obs,
                samples_ordered=samples_ordered,
                formula=design,
                covariates=covariates,
                treatment_col=treatment_col,
            )
            # Force GLM regardless of what 'auto' resolved to: covariates only
            # work with the GLM path.
            selected_test = "glm"

        dmr_df = call_dmr_tile_based(
            methylstore_path=md.store,
            samples_treatment=md.treatment_ids,
            samples_control=md.control_ids,
            tile_size_bp=tile_size_bp,
            test=selected_test,
            chromosomes=chromosomes,
            min_cpgs_per_tile=min_cpgs_per_tile,
            alpha=alpha,
            min_abs_meth_diff=min_abs_meth_diff,
            unite=unite,
            min_samples_treatment=min_samples_treatment,
            min_samples_control=min_samples_control,
            dispersion=dispersion,
            reference=reference,
            design_full=design_full,
            design_reduced=design_reduced,
            coef_idx=coef_idx,
            backend=backend,
            n_workers=n_workers,
            merge_adjacent=merge_adjacent,
        )

        # Optional q-value post-filter (the tile path already filtered at
        # `alpha`, but a stricter user threshold is allowed here).
        if len(dmr_df) > 0 and min_mean_qvalue is not None and "qvalue" in dmr_df.columns:
            dmr_df = dmr_df.filter(pl.col("qvalue") < min_mean_qvalue)

        # permutation FDR. Refuses to run when a covariate design
        # is in play (shuffling treatment labels invalidates the assumed
        # covariate structure).
        if empirical_fdr:
            if design is not None or (covariates is not None and len(covariates) > 0):
                raise ValueError(
                    "empirical_fdr=True is not supported with covariate "
                    "designs (label-shuffling invalidates stratification). "
                    "Use a stratified-permutation scheme manually if needed."
                )
            # Build strata map from obs column when empirical_strata= supplied.
            strata_map: dict[str, list[str]] | None = None
            if empirical_strata is not None and empirical_strata in md.obs.columns:
                all_samples = list(md.treatment_ids) + list(md.control_ids)
                obs_indexed = md.obs.filter(
                    pl.col("sample_id").is_in(all_samples)
                )
                strata_map = {}
                for row in obs_indexed.iter_rows(named=True):
                    strata_map.setdefault(row[empirical_strata], []).append(
                        row["sample_id"]
                    )
            if len(dmr_df) > 0:
                dmr_df = empirical_fdr_for_dmr(
                    methylstore_path=md.store,
                    samples_treatment=md.treatment_ids,
                    samples_control=md.control_ids,
                    observed_dmr=dmr_df,
                    n_perm=n_perm,
                    seed=perm_seed,
                    n_jobs=perm_n_jobs,
                    empirical_strata=strata_map,
                    tile_size_bp=tile_size_bp,
                    test=selected_test,
                    chromosomes=chromosomes,
                    min_cpgs_per_tile=min_cpgs_per_tile,
                    alpha=alpha,
                    min_abs_meth_diff=min_abs_meth_diff,
                    unite=unite,
                    min_samples_treatment=min_samples_treatment,
                    min_samples_control=min_samples_control,
                    dispersion=dispersion,
                    reference=reference,
                )

        md.uns["dmr"] = dmr_df
        md.uns["dmr_params"] = {
            "method": "tile",
            "tile_size_bp": tile_size_bp,
            "min_cpgs_per_tile": min_cpgs_per_tile,
            "test": selected_test,
            "alpha": alpha,
            "min_abs_meth_diff": min_abs_meth_diff,
            "min_mean_qvalue": min_mean_qvalue,
            "min_samples_treatment": min_samples_treatment,
            "min_samples_control": min_samples_control,
            "unite": unite,
            "dispersion": dispersion,
            "reference": reference,
            "design": design,
            "covariates": list(covariates) if covariates else None,
            "treatment_col": treatment_col,
            "formula_used": formula_used,
            "design_terms": term_names if term_names else None,
            "empirical_fdr": empirical_fdr,
            "n_perm": n_perm if empirical_fdr else None,
            "perm_seed": perm_seed if empirical_fdr else None,
        }
        if tsv is not None:
            from .export import dmr_to_tsv
            _emit_result_tsv(lambda: dmr_to_tsv(md, tsv), tsv, is_auto=_tsv_is_auto)
        return

    if method == "sliding_window":
        # Prefer streaming from the persistent DMC store when ep.tl.dmc
        # has staged one -- keeps DMR peak memory at O(largest chrom)
        # instead of holding the full 22M-row DataFrame plus a per-chrom
        # working set in memory at the same time.
        dmc_store_path = md.uns.get("dmc", {}).get("store_path")
        dmc_input: object
        if dmc_store_path:
            from ._dmc_store import DMCStore
            from pathlib import Path as _Path
            store_path = _Path(dmc_store_path)
            if (store_path / ".epykit_dmc_manifest.json").exists():
                dmc_input = DMCStore.open(store_path)
            else:
                dmc_input = md.dmc
        else:
            dmc_input = md.dmc

        if dmc_input is None:
            raise ValueError(
                "No DMC results available. Run ep.tl.dmc(md) first, or use "
                "method='tile' which goes directly to the methylstore."
            )

        dmr_df = call_dmr_sliding_window(
            dmc_results=dmc_input,
            window_bp=window_bp,
            step_bp=step_bp,
            min_cpgs=min_cpgs,
            min_sites_significant=min_sites_significant,
            alpha=alpha,
            min_abs_meth_diff=min_abs_meth_diff,
        )
        # filter on the BH-corrected DMR-level q-value, not the raw
        # combined p-value. ``call_dmr_sliding_window`` now adds
        # ``combined_qvalue`` itself.
        if len(dmr_df) > 0 and min_mean_qvalue is not None:
            q_col = "combined_qvalue" if "combined_qvalue" in dmr_df.columns else "combined_pvalue"
            dmr_df = dmr_df.filter(pl.col(q_col) < min_mean_qvalue)

        md.uns["dmr"] = dmr_df
        md.uns["dmr_params"] = {
            "method": "sliding_window",
            "window_bp": window_bp,
            "step_bp": step_bp,
            "min_cpgs": min_cpgs,
            "min_sites_significant": min_sites_significant,
            "alpha": alpha,
            "min_abs_meth_diff": min_abs_meth_diff,
            "min_mean_qvalue": min_mean_qvalue,
        }
        if tsv is not None:
            from .export import dmr_to_tsv
            _emit_result_tsv(lambda: dmr_to_tsv(md, tsv), tsv, is_auto=_tsv_is_auto)
        return

    if method == "segment":
        from .dmr_segment import call_dmr_rule_segment
        dmc_df = md.dmc
        if dmc_df is None:
            raise ValueError(
                "method='segment' needs a DMC table on md. Run ep.tl.dmc(md) first."
            )
        dmr_df = call_dmr_rule_segment(
            dmc_df,
            min_cpgs=min_cpgs,
            min_abs_meth_diff=min_abs_meth_diff,
            alpha=alpha,
        )
        md.uns["dmr"] = dmr_df
        md.uns["dmr_params"] = {
            "method": "segment",
            "min_cpgs": min_cpgs,
            "min_abs_meth_diff": min_abs_meth_diff,
            "alpha": alpha,
        }
        if tsv is not None:
            from .export import dmr_to_tsv
            _emit_result_tsv(lambda: dmr_to_tsv(md, tsv), tsv, is_auto=_tsv_is_auto)
        return

    if method == "chain_merge":
        # DSS::callDMR semantics -- chain contiguous sig CpGs whose gap is
        # <= dis_merge_bp, then filter by minlen_bp / min_cpgs / pct_sig.
        # Reuses the same DMC store as sliding_window when available so a
        # 22M-CpG run stays streaming-friendly.
        dmc_store_path = md.uns.get("dmc", {}).get("store_path")
        if dmc_store_path:
            from ._dmc_store import DMCStore
            from pathlib import Path as _Path
            store_path = _Path(dmc_store_path)
            if (store_path / ".epykit_dmc_manifest.json").exists():
                dmc_input = DMCStore.open(store_path)
            else:
                dmc_input = md.dmc
        else:
            dmc_input = md.dmc

        if dmc_input is None:
            raise ValueError(
                "method='chain_merge' needs a DMC table on md. "
                "Run ep.tl.dmc(md) first."
            )

        dmr_df = call_dmr_chain_merge(
            dmc_input,
            preset=preset,
            alpha=alpha,
            min_abs_meth_diff=min_abs_meth_diff,
            dis_merge_bp=dis_merge_bp,
            min_cpgs=min_cpgs,
            pct_sig=pct_sig,
            minlen_bp=minlen_bp,
            use_q_for_sig=use_q_for_sig,
        )

        # Same post-hoc q-value filter as the sliding-window path: drop
        # candidate DMRs whose BH-corrected combined q-value isn't sig.
        if len(dmr_df) > 0 and min_mean_qvalue is not None:
            q_col = "combined_qvalue" if "combined_qvalue" in dmr_df.columns else "combined_pvalue"
            dmr_df = dmr_df.filter(pl.col(q_col) < min_mean_qvalue)

        md.uns["dmr"] = dmr_df
        md.uns["dmr_params"] = {
            "method": "chain_merge",
            "alpha": alpha,
            "min_abs_meth_diff": min_abs_meth_diff,
            "dis_merge_bp": dis_merge_bp,
            "min_cpgs": min_cpgs,
            "pct_sig": pct_sig,
            "minlen_bp": minlen_bp,
            "use_q_for_sig": use_q_for_sig,
            "min_mean_qvalue": min_mean_qvalue,
        }
        if tsv is not None:
            from .export import dmr_to_tsv
            _emit_result_tsv(lambda: dmr_to_tsv(md, tsv), tsv, is_auto=_tsv_is_auto)
        return

    raise ValueError(
        f"Unknown DMR method '{method}'. Expected 'tile', 'sliding_window', "
        f"'segment', or 'chain_merge'."
    )


def diagnose_dmr_calling(
    md: MethylData,
    reference_dmrs: pl.DataFrame,
    *,
    dmc_key: str | None = None,
    chromosomes: list[str] | None = None,
    alpha_threshold: float = 1e-5,
) -> dict:
    """Classify reference DMRs by recovery status to debug missing-DMR causes.

    Given a reference DMR set (e.g. from a published paper or another
    pipeline), bucket each reference DMR by *why* epykit's current DMC +
    DMR-calling produced (or didn't produce) an overlapping DMR. Lets
    you diagnose a low recall number into actionable categories instead
    of guessing which parameter to tune.

    Five buckets:

    * ``SUCCESS_OVERLAP`` -- our DMR set already contains a region that
      overlaps this reference DMR. Nothing to fix.
    * ``H1_NO_CPGS`` -- 0 of our united CpGs fall inside the reference
      coordinates. The coverage filter or unite step dropped them;
      relax ``min_coverage`` or use ``type="union"`` for ``pp.unite``.
    * ``H2_NO_SIG_CPGS`` -- at least one CpG present but none reach
      ``q < 0.05``. The DMC test statistic is too conservative for this
      region. The only fix is a more powerful test (e.g. a spatial-
      covariance Wald test); no DMR-caller tuning can recover this.
    * ``H3a_WEAK_ALPHA`` -- has sig CpGs at ``q < 0.05`` but none reach
      ``alpha_threshold`` (default 1e-5, matching DSS callDMR). Recover
      by loosening ``alpha`` in :func:`ep.tl.dmr` (e.g. to 1e-4 or
      1e-3, ideally via ``preset="permissive"``).
    * ``H3b_STRUCTURE`` -- CpGs at ``q < alpha_threshold`` exist in the
      region but no DMR was called. Chain-merge dropped the candidate
      on a structural filter: ``min_cpgs``, ``pct_sig``, ``minlen_bp``,
      or ``dis_merge_bp``. Loosen ``dis_merge_bp`` first (highest
      Pareto leverage), then ``min_cpgs``.

    Parameters
    ----------
    md : MethylData
        Must have a DMC table populated (``ep.tl.dmc`` already run) and
        ``md.uns['dmr']`` populated (``ep.tl.dmr`` already run).
    reference_dmrs : polars DataFrame
        Reference DMR set with at least ``chrom``, ``start``, ``end``
        columns. Coordinates assumed 0-based half-open (BED convention).
    dmc_key : str, optional
        Specific key in ``md.varm`` to use as the DMC table. Defaults to
        ``md.uns['dmc']['last_key']``.
    chromosomes : list of str, optional
        Restrict analysis to these chromosomes (e.g. main chroms only,
        skipping ``_random`` / alt contigs). Defaults to all chromosomes
        present in both the reference and the DMC table.
    alpha_threshold : float, default 1e-5
        The per-CpG significance cutoff that was used in chain-merge.
        Determines the H3a vs H3b boundary. Pass the same value you used
        in ``ep.tl.dmr(alpha=...)`` so the diagnosis matches your run.

    Returns
    -------
    dict with keys:

    * ``"counts"`` (dict[str, int]) -- count per bucket
    * ``"bucket_indices"`` (dict[str, list[int]]) -- 0-based row indices
      into ``reference_dmrs`` for each bucket
    * ``"n_reference"`` (int) -- total reference DMRs analyzed
    * ``"summary"`` (str) -- multi-line human-readable summary
    * ``"alpha_threshold"`` (float) -- the threshold used
    """
    import numpy as np
    from collections import defaultdict

    # ---- Resolve inputs ----
    if dmc_key is None:
        dmc_key = md.uns.get("dmc", {}).get("last_key")
    if dmc_key is None or dmc_key not in md.varm:
        raise ValueError(
            "No DMC table found on md. Run ep.tl.dmc(md, ...) first, or pass "
            "dmc_key= explicitly."
        )
    dmc = md.varm[dmc_key]

    if "dmr" not in md.uns:
        raise ValueError(
            "No DMR table found at md.uns['dmr']. Run ep.tl.dmr(md, ...) "
            "first so the diagnostic knows which reference DMRs we recovered."
        )
    ours_dmr = md.uns["dmr"]

    # Normalize reference column types
    ref = reference_dmrs.with_columns([
        pl.col("chrom").cast(pl.Utf8),
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
    ])

    # Optional chromosome filter
    if chromosomes is not None:
        chrom_set = set(chromosomes)
        ref = ref.filter(pl.col("chrom").is_in(chrom_set))

    # Pick the q-value column on the DMC table (prefer qvalue over pvalue)
    qcol = "qvalue" if "qvalue" in dmc.columns else (
        "pvalue" if "pvalue" in dmc.columns else None
    )
    if qcol is None:
        raise ValueError(
            f"DMC table at varm[{dmc_key!r}] has neither 'qvalue' nor "
            f"'pvalue' columns; cannot diagnose."
        )
    pos_col = "pos" if "pos" in dmc.columns else "start"

    # ---- Build indexed lookups ----
    # Our DMRs: chrom -> sorted [(start, end), ...]
    ours_by_chr: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for r in ours_dmr.iter_rows(named=True):
        ours_by_chr[str(r["chrom"])].append((int(r["start"]), int(r["end"])))
    for c in ours_by_chr:
        ours_by_chr[c].sort()

    def _has_overlap(intervals, s, e):
        # Linear scan is fine for typical DMR-set sizes; bisect would help
        # only with 100k+ DMRs which is unusual.
        for ps, pe in intervals:
            if pe < s: continue
            if ps > e: break
            return True
        return False

    # DMC table: per-chrom (sorted positions array, q-values array). We
    # convert to numpy for searchsorted speed; the cost is O(n_cpgs) per
    # chromosome but only paid once across all reference DMRs.
    dmc_pd = dmc.to_pandas()
    dmc_pd["chrom"] = dmc_pd["chrom"].astype(str)
    dmc_by_chr: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for chrom, g in dmc_pd.groupby("chrom"):
        if chromosomes is not None and chrom not in chrom_set:
            continue
        g = g.sort_values(pos_col)
        dmc_by_chr[chrom] = (
            g[pos_col].to_numpy(dtype=np.int64),
            g[qcol].to_numpy(dtype=np.float64),
        )

    # ---- Classify ----
    buckets: dict[str, list[int]] = {
        "SUCCESS_OVERLAP": [], "H1_NO_CPGS": [], "H2_NO_SIG_CPGS": [],
        "H3a_WEAK_ALPHA": [], "H3b_STRUCTURE": [],
    }

    for idx, r in enumerate(ref.iter_rows(named=True)):
        chrom = str(r["chrom"])
        rs, re_ = int(r["start"]), int(r["end"])

        # Bucket 0: already recovered?
        if _has_overlap(ours_by_chr.get(chrom, []), rs, re_):
            buckets["SUCCESS_OVERLAP"].append(idx)
            continue

        # Bucket H1: no CpGs in the region
        entry = dmc_by_chr.get(chrom)
        if entry is None:
            buckets["H1_NO_CPGS"].append(idx)
            continue
        pos_arr, q_arr = entry
        lo = int(np.searchsorted(pos_arr, rs, side="left"))
        hi = int(np.searchsorted(pos_arr, re_, side="right"))
        if lo == hi:
            buckets["H1_NO_CPGS"].append(idx)
            continue

        region_q = q_arr[lo:hi]
        min_q = float(np.nanmin(region_q))

        if min_q >= 0.05:
            buckets["H2_NO_SIG_CPGS"].append(idx)
        elif min_q >= alpha_threshold:
            buckets["H3a_WEAK_ALPHA"].append(idx)
        else:
            buckets["H3b_STRUCTURE"].append(idx)

    # ---- Assemble result ----
    counts = {k: len(v) for k, v in buckets.items()}
    n_ref = sum(counts.values())

    # Build a human-readable summary
    lines = [
        f"DMR-calling diagnostic on {n_ref} reference DMRs:",
        f"  alpha_threshold = {alpha_threshold:.0e}  (matches chain-merge alpha)",
        "",
    ]
    bucket_help = {
        "SUCCESS_OVERLAP":  "already recovered (no fix needed)",
        "H1_NO_CPGS":       "no CpGs present  -> coverage/unite issue",
        "H2_NO_SIG_CPGS":   "no sig CpGs at q<0.05  -> need better test stat (e.g. Wald-smoothed)",
        "H3a_WEAK_ALPHA":   f"sig CpGs at q<0.05 but not q<{alpha_threshold:.0e}  -> loosen alpha (preset='permissive')",
        "H3b_STRUCTURE":    "sig CpGs exist but chain-merge dropped  -> loosen dis_merge_bp first",
    }
    for name in ("SUCCESS_OVERLAP", "H1_NO_CPGS", "H2_NO_SIG_CPGS",
                 "H3a_WEAK_ALPHA", "H3b_STRUCTURE"):
        n = counts[name]
        pct = n / max(n_ref, 1)
        lines.append(f"  {name:<18} {n:>5} ({pct:>5.1%})  -- {bucket_help[name]}")

    return {
        "counts": counts,
        "bucket_indices": buckets,
        "n_reference": n_ref,
        "summary": "\n".join(lines),
        "alpha_threshold": alpha_threshold,
    }


def dvc(
    md: MethylData,
    test: str = "brown_forsythe",
    chromosomes: list[str] | None = None,
    alpha: float = 0.05,
    mean_filter_alpha: float = 0.05,
    min_coverage: int = 1,
    *,
    backend: str = "sequential",
    n_workers: int | None = None,
    tsv: str | None = None,
    tsv_full: bool = False,
    tsv_alpha: float = 0.05,
    # Deprecated aliases for the tsv* args above (epykit writes TSV by default).
    csv: str | None = None,
    csv_full: bool = False,
    csv_alpha: float = 0.05,
) -> None:
    """Differential-Variability CpG calling (iEVORA-style).

    Identifies CpGs whose between-replicate variance differs significantly
    between the treatment and control groups *while* the means do not --
    the signature of an outlier-driven shift in variability that purely
    mean-based DMC analysis misses (cancer / aging methylomes).

    Result is stored at ``md.varm["dvc"]`` with columns:
        chrom, pos, strand, n_treatment, n_control,
        var_treatment, var_control, var_log_ratio,
        p_variance, q_variance, p_mean, q_mean, is_dvc

    Parameters
    ----------
    test : {"brown_forsythe", "bartlett"}
        Variance-equality test. ``"brown_forsythe"`` (default; median-centred
        Levene) is robust to the bounded, U-shaped beta distribution and is
        the test actually run. ``"bartlett"`` is a deprecated alias that runs
        Brown-Forsythe and emits a ``UserWarning`` (Bartlett's test is not
        implemented). Brown-Forsythe needs >=3 replicates per group; at n=2
        the within-group sum of squares is 0 so the variance p-values are NaN
        and no DVCs are called (a warning is emitted).
    min_coverage : int
        Mask per-replicate beta below this coverage before the variance test
        (default 1, i.e. any covered site). Raise it on cohorts with
        imbalanced sequencing depth, where low-coverage binomial noise can
        masquerade as differential biological variance.
    alpha : float
        q-value cutoff on the variance test for the ``is_dvc`` flag.
    mean_filter_alpha : float
        Sites are flagged DVC only when ``p_mean > mean_filter_alpha`` --
        i.e. variance changes that aren't accompanied by mean changes.
    """
    from .dvc import process_chromosomes_dvc
    unite_info = md.uns.get("unite")
    unite = (unite_info is not None) and (unite_info.get("type") == "intersect")
    result = process_chromosomes_dvc(
        methylstore_path=md.store,
        samples_treatment=md.treatment_ids,
        samples_control=md.control_ids,
        test=test,
        chromosomes=chromosomes,
        unite=unite,
        mean_filter_alpha=mean_filter_alpha,
        alpha=alpha,
        min_coverage=min_coverage,
        backend=backend,
        n_workers=n_workers,
    )
    md.varm["dvc"] = result
    md.uns["dvc"] = {
        "test": test,
        "alpha": alpha,
        "mean_filter_alpha": mean_filter_alpha,
        "n_sites": len(result),
        "n_dvc": int(result.get_column("is_dvc").sum()) if len(result) else 0,
        "unite": unite,
    }

    tsv, tsv_full, tsv_alpha = _resolve_tsv_output(
        tsv, csv, tsv_full=tsv_full, csv_full=csv_full,
        tsv_alpha=tsv_alpha, csv_alpha=csv_alpha,
    )
    if tsv is not None:
        from .export import dvc_to_tsv
        dvc_to_tsv(md, tsv, alpha=tsv_alpha, full=tsv_full)


def dvr(
    md: MethylData,
    *,
    tile_size_bp: int = 1000,
    min_cpgs_per_tile: int = 5,
    alpha: float = 0.05,
) -> None:
    """Differentially Variable Regions -- density-based aggregation of DVC.

    Requires ``ep.tl.dvc(md)`` to have been run first; reads
    ``md.varm['dvc']`` and writes the region call to ``md.uns['dvr']``.
    See :func:`epykit.dvc.call_dvr_density` for the statistical model
    (per-tile binomial enrichment vs the genome-wide DVC rate).

    Parameters
    ----------
    tile_size_bp : int
        Tile width in bp. Default 1 kb.
    min_cpgs_per_tile : int
        Tiles below this size are dropped. Default 5.
    alpha : float
        BH q-value threshold for the ``is_dvr`` flag. Default 0.05.
    """
    if "dvc" not in md.varm or md.varm["dvc"] is None:
        raise ValueError(
            "md.varm['dvc'] is missing. Run ep.tl.dvc(md) before ep.tl.dvr(md)."
        )
    from .dvc import call_dvr_density
    dvr_df = call_dvr_density(
        md.varm["dvc"],
        tile_size_bp=tile_size_bp,
        min_cpgs_per_tile=min_cpgs_per_tile,
        alpha=alpha,
    )
    md.uns["dvr"] = dvr_df
    md.uns["dvr_params"] = {
        "method": "density",
        "tile_size_bp": tile_size_bp,
        "min_cpgs_per_tile": min_cpgs_per_tile,
        "alpha": alpha,
        "n_regions": int(dvr_df.height),
        "n_dvr": int(dvr_df.get_column("is_dvr").sum()) if dvr_df.height else 0,
    }


def age_clock(
    md: MethylData,
    coefficients,
    manifest,
    *,
    intercept: float = 0.0,
    transform: str | None = None,
    impute_missing: bool = True,
    name: str = "age_clock",
) -> None:
    """Run a linear epigenetic-age clock and write per-sample ages to ``md.obs``.

    Thin orchestrator over :func:`epykit.clocks.age_clock` that joins the
    resulting per-sample age onto ``md.obs`` so downstream code can use it
    as a regular obs column (e.g. as a covariate via ``tl.dmc(formula=...)``).
    The full per-CpG diagnostic table is parked at
    ``md.uns[f'{name}_diagnostics']`` for QC.
    """
    from .clocks import age_clock as _age_clock
    result = _age_clock(
        md, coefficients, manifest,
        intercept=intercept, transform=transform,
        impute_missing=impute_missing, name=name,
    )
    md.obs = md.obs.join(
        result.select(["sample_id", name]),
        on="sample_id", how="left",
    )
    md.uns[f"{name}_diagnostics"] = result


def deconvolve(
    md: MethylData,
    reference,
    manifest,
    *,
    method: str = "nnls",
    cell_types: list[str] | None = None,
    uns_key: str = "deconvolution",
) -> None:
    """Reference-based cell-type deconvolution; results to ``md.uns[uns_key]``.

    Thin orchestrator over :func:`epykit.clocks.deconvolve`. The long-form
    proportions table is stored at ``md.uns[uns_key]``; the wide pivot
    (one column per cell type) is left-joined onto ``md.obs`` so cell-type
    proportions can be used as covariates.
    """
    from .clocks import deconvolve as _deconvolve
    long = _deconvolve(
        md, reference, manifest, method=method, cell_types=cell_types,
    )
    md.uns[uns_key] = long
    if long.is_empty():
        return
    wide = long.pivot(
        values="proportion", index="sample_id", on="cell_type",
        aggregate_function="first",
    )
    # Prefix columns so cell-type names like 'CD4T' don't collide with
    # existing obs columns.
    rename_map = {c: f"frac_{c}" for c in wide.columns if c != "sample_id"}
    wide = wide.rename(rename_map)
    md.obs = md.obs.join(wide, on="sample_id", how="left")


def annotate(
    md: MethylData,
    gtf: str | None = None,
    cpg_islands: str | None = None,
    significant_only: bool = True,
    alpha: float = 0.05,
    promoter_upstream_bp: int = 2000,
    promoter_downstream_bp: int = 200,
    clear_gtf_cache: bool = True,
    multi_annotation: bool = True,
    *,
    refgene: str | None = None,
    gene_type_filter: str | list[str] | tuple[str, ...] | None = None,
    features: list[str] | tuple[str, ...] | None = None,
    # Auto-emits the annotated DMC table to
    # <analysis_root>/results/dmc_annotated.tsv by default; tsv=False disables,
    # tsv="path" overrides. csv= is a deprecated alias.
    tsv: str | bool | None = None,
    csv: str | None = None,
) -> None:
    """Annotate DMC/DMR outputs.

    By default only significant DMCs are annotated to avoid OOM. Set
    `significant_only=False` to annotate all sites (not recommended for
    whole-genome datasets).

    Parameters
    ----------
    gtf, refgene : str or None
        Annotation source. Provide exactly one of ``gtf`` (GENCODE / Ensembl
        GTF) or ``refgene`` (UCSC ``refGene.txt(.gz)`` -- HOMER's default
        catalog, gives the highest paper-gene recall for methylation work
        because it's curated and protein-coding-biased).
    gene_type_filter : str or list of str or None, keyword-only
        Restrict the gene catalog before building overlap intervals and the
        nearest-TSS index. Typical: ``"protein_coding"`` to drop lincRNAs /
        pseudogenes / novel predictions. Works on both sources.
    clear_gtf_cache : bool, optional
        If True (default), clear the GTF cache and run garbage collection
        after annotation. Set to False if you plan to call annotate()
        multiple times to reuse the cached GTF.
    multi_annotation : bool, optional
        If True (default), populate annotatr-style columns on every
        annotated table: ``nearest_tss_gene`` / ``nearest_tss_distance``
        (HOMER's nearest-TSS rule), plus ``all_overlapping_genes`` /
        ``all_overlapping_features`` (one-to-many). Set False to skip them
        and keep only the legacy single-best gene-name columns. See
        :func:`epykit.annotate.annotate_features` for details.
    features : sequence of str or None, keyword-only
        Override the feature classes built by
        :func:`epykit.annotate.annotate_features`. Default ``None`` lets
        the lower-level function pick its own default (the full HOMER set:
        promoter / 5UTR / exon / intron / 3UTR / TTS / noncoding). Pass a
        narrower tuple to skip categories -- e.g.
        ``("promoter", "exon", "intron")`` for the pre-0.x coarse
        breakdown when downstream code expects only those buckets.
    """
    if gtf and refgene:
        raise ValueError("Provide only one of gtf or refgene, not both")
    feature_source_present = bool(gtf or refgene)
    if not feature_source_present and not cpg_islands:
        raise ValueError("Provide at least one of gtf / refgene / cpg_islands")

    for key, df in list(md.varm.items()):
        if not key.startswith("dmc"):
            continue

        # Match old behavior: annotate only significant sites to avoid OOM
        if significant_only:
            p_col = "qvalue" if "qvalue" in df.columns else "pvalue"
            ann = df.filter(pl.col(p_col) < alpha)
        else:
            ann = df

        if len(ann) == 0:
            continue

        if feature_source_present:
            # Pass through to the new annotate_features API. `gtf` and
            # `refgene` are the wrapper-level convenience kwargs; under
            # the hood there's one ``annotation`` argument with explicit
            # source. We forward whichever the caller set.
            annotation_path = gtf if gtf is not None else refgene
            forwarded_source = "gtf" if gtf is not None else "refgene"
            af_kwargs: dict[str, Any] = dict(
                source=forwarded_source,
                promoter_upstream_bp=promoter_upstream_bp,
                promoter_downstream_bp=promoter_downstream_bp,
                multi_annotation=multi_annotation,
                gene_type_filter=gene_type_filter,
            )
            if features is not None:
                af_kwargs["features"] = features
            ann = annotate_features(ann, annotation_path, **af_kwargs)
        if cpg_islands:
            ann = annotate_cpg_islands(ann, cpg_island_bed=cpg_islands)

        # Store as separate key so full DMC results are preserved
        md.varm[f"{key}_annotated"] = ann

    if "dmr" in md.uns and isinstance(md.uns["dmr"], pl.DataFrame) and feature_source_present:
        annotation_path = gtf if gtf is not None else refgene
        forwarded_source = "gtf" if gtf is not None else "refgene"
        af_kwargs = dict(
            source=forwarded_source,
            promoter_upstream_bp=promoter_upstream_bp,
            promoter_downstream_bp=promoter_downstream_bp,
            multi_annotation=multi_annotation,
            gene_type_filter=gene_type_filter,
        )
        if features is not None:
            af_kwargs["features"] = features
        md.uns["dmr"] = annotate_features(md.uns["dmr"], annotation_path, **af_kwargs)

    md.uns["annotation"] = {
        "gtf": gtf,
        "refgene": refgene,
        "cpg_islands": cpg_islands,
        "significant_only": significant_only,
        "alpha": alpha,
        "promoter_upstream_bp": promoter_upstream_bp,
        "promoter_downstream_bp": promoter_downstream_bp,
        "multi_annotation": multi_annotation,
        "gene_type_filter": gene_type_filter,
        "features": list(features) if features is not None else None,
    }

    # Clear GTF cache if requested (default: True)
    if clear_gtf_cache and gtf:
        _GTF_CACHE.clear()
        gc.collect()

    # Auto-emit the annotated DMC table as a human-readable TSV (default on).
    tsv_out, _, _, _tsv_is_auto = _resolve_auto_tsv(
        md, tsv, csv, default_name="dmc_annotated.tsv"
    )
    if tsv_out is not None and md.get_dmc(annotated=True) is not None:
        from .export import dmc_to_tsv
        _emit_result_tsv(
            lambda: dmc_to_tsv(md, tsv_out, full=True), tsv_out, is_auto=_tsv_is_auto,
        )


def asm(
    md: MethylData,
    *,
    bam,
    vcf,
    min_reads_per_haplotype: int = 10,
    min_phased_snvs: int = 1,
    chromosomes: list[str] | None = None,
    caller: str = "bismark",
) -> None:
    """Allele-specific methylation (ASM) caller -- 0.5.0.

    See :func:`epykit.asm.call_asm` for the algorithm. Per-CpG ASM tests
    are stored at ``md.varm["asm"]`` with the same column names as the
    ``dmc_*`` family (``pvalue``, ``qvalue``, ``meth_diff``) so volcano
    / Manhattan plots work without modification.

    Parameters
    ----------
    bam : mapping
        ``{sample_id -> bam_path}``. BAMs must be coordinate-sorted and
        indexed; per-base methylation calls come from Bismark ``XM``
        tags or MethylDackel ``MM/ML`` tags.
    vcf : str | Path
        Per-individual VCF (bgzipped + tabix preferred). Heterozygous
        biallelic SNVs are used as phasing anchors.
    """
    from .asm import asm as _asm
    _asm(md, bam=bam, vcf=vcf,
         min_reads_per_haplotype=min_reads_per_haplotype,
         min_phased_snvs=min_phased_snvs,
         chromosomes=chromosomes,
         caller=caller)


def entropy(
    md: MethylData,
    *,
    bam,
    window_cpgs: int = 4,
    min_reads: int = 10,
    chromosomes: list[str] | None = None,
    caller: str = "bismark",
) -> None:
    """Methylation entropy caller -- 0.5.0.

    See :func:`epykit.entropy.call_entropy` for the algorithm. Per-CpG-
    window Shannon entropy is stored at ``md.varm["entropy"]``.
    """
    from .entropy import entropy as _entropy
    _entropy(md, bam=bam, window_cpgs=window_cpgs, min_reads=min_reads,
             chromosomes=chromosomes, caller=caller)


def pmd(
    md: MethylData,
    *,
    samples: list[str] | None = None,
    bandwidth_bp: float = 10_000,
    beta_threshold: float = 0.55,
    min_pmd_bp: int = 100_000,
    chromosomes: list[str] | None = None,
    backend: str = "sequential",
    n_workers: int | None = None,
) -> None:
    """Partially methylated domain (PMD) caller -- 0.6.0.

    See :func:`epykit.pmd.call_pmd_one_sample` for the algorithm.
    Per-sample, megabase-scale 2-state HMM segmentation; results land
    in ``md.uns["pmd"]``.
    """
    from .pmd import pmd as _pmd
    _pmd(md, samples=samples, bandwidth_bp=bandwidth_bp,
         beta_threshold=beta_threshold, min_pmd_bp=min_pmd_bp,
         chromosomes=chromosomes, backend=backend, n_workers=n_workers)


def hmr(
    md: MethylData,
    *,
    samples: list[str] | None = None,
    hmr_threshold: float = 0.30,
    lmr_max_density: float = 0.020,
    min_cpgs: int = 4,
    chromosomes: list[str] | None = None,
    backend: str = "sequential",
    n_workers: int | None = None,
) -> None:
    """HMR / LMR caller (MethylSeekR-style) -- 0.6.0.

    See :func:`epykit.hmr.call_hmr_one_sample` for the algorithm.
    Two-state HMM per sample over raw per-CpG beta; results land in
    ``md.uns["hmr"]`` and ``md.uns["lmr"]``.
    """
    from .hmr import hmr as _hmr
    _hmr(md, samples=samples, hmr_threshold=hmr_threshold,
         lmr_max_density=lmr_max_density, min_cpgs=min_cpgs,
         chromosomes=chromosomes, backend=backend, n_workers=n_workers)