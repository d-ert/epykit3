"""The stages of ``tl.dmc``.

``tl.dmc`` builds one :class:`DMCConfig` and then runs the stages in this
module in a fixed order. Every stage takes the resolved config through a
:class:`DMCPlan` rather than loose arguments, so no stage re-derives a
knob. The three ways a result can arrive (the formula / contrast run, the
``resumable=True`` cache hit and the ordinary binary run) all produce a
:class:`DMCOutcome` of the same shape, so :func:`publish` and
:func:`finish` do not branch on how the result was produced.

Rules the split keeps:

- ``md.uns["dmc"]`` is written in exactly one place, :func:`publish`,
  through :meth:`DMCConfig.to_uns`;
- nothing new is materialised: :func:`run_engine` returns the streaming
  :class:`DMCStore` and :func:`post_process` is the only binary-path stage
  that may hold the full table;
- warnings raised inside a stage carry the ``stacklevel`` that points at
  the caller of ``tl.dmc``. A stage is one frame below ``tl.dmc``, the
  ``open_input_store`` generator body sits one frame below the context
  manager's ``__enter__``, and the shared ``tl`` helpers take the depth
  they need as an argument.

The shared TSV, warning and sample-count helpers stay in ``tl`` because
``tl.dmr`` uses them too. Stages import them inside the function body,
after ``tl`` has finished importing this module, so there is no import
cycle.
"""

from __future__ import annotations

import logging
import tempfile
import warnings
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import polars as pl

from ._dmc_config import DMCConfig
from .dmc import (
    DMCStore,
    apply_multiple_testing_correction,
    combine_neighbour_pvalues,
    empirical_fdr_for_dmc,
    process_chromosomes_dmc,
)
from .methyldata import MethylData

logger = logging.getLogger(__name__)

Mode = Literal["binary", "contrast"]

_LOG2_ODDS_RATIO_NOTICE = (
    "The 'log2_odds_ratio' column is deprecated and is slated for removal "
    "in a future major release. Use 'log2_odds_ratio_pooled' for "
    "pooled-count tests (lr, fisher) or 'coef_treatment_log2' for the glm "
    "backend. The transitional column is NaN-filled."
)


@dataclass(frozen=True)
class TsvPlan:
    """Where the significant-DMC TSV goes, resolved once.

    ``None`` in :attr:`DMCPlan.tsv` means no export. ``is_auto`` keeps the
    "auto-emitted, so a failure only logs" behaviour of ``_emit_result_tsv``.
    """

    path: str
    full: bool
    alpha: float
    is_auto: bool


@dataclass(frozen=True)
class DMCPlan:
    """Everything the run needs, decided before any chromosome is read.

    Built once by :func:`plan_run`; read-only afterwards. ``cfg`` is the
    config after :meth:`DMCConfig.apply_power_stack`, so the stages never
    re-derive a knob. ``key`` is the ``md.varm`` key the result lands under
    (``dmc_<test>``, ``dmc_<test>_smoothed`` or ``dmc_glm_contrast``).
    """

    cfg: DMCConfig
    mode: Mode
    selected_test: str
    """The engine after ``"auto"`` resolution, ``"glm_contrast"`` in contrast
    mode; the ``test_used`` recorded in uns and the resume stage name."""
    unite: bool
    smooth_method: str | None
    key: str
    tsv: TsvPlan | None


@dataclass(frozen=True)
class ContrastDesign:
    """What the contrast run resolved from ``md.obs``.

    Carried on the outcome so :func:`publish` can record ``formula``,
    ``contrast`` and ``design_terms`` without re-resolving them.
    """

    formula_used: str
    contrast_label: str
    design_terms: list[str]


@dataclass(frozen=True)
class ResumeTicket:
    """The ``resumable=True`` lookup, hit or miss.

    ``hit`` is the cached table when the manifest has a matching entry and
    its sidecar exists; ``None`` otherwise. The other three fields are what
    :func:`persist_resume` needs to write the sidecar after a miss.
    """

    stage_name: str
    root: str
    signature: str
    hit: pl.DataFrame | None


@dataclass(frozen=True)
class DMCOutcome:
    """One shape for the three ways a result can arrive.

    ``result`` is the full table when the run materialised (binary with
    ``materialize=True``, every contrast run, every resume hit) and ``None``
    when only the streaming store exists. ``store`` is ``None`` on a resume
    hit, which is why ``store_path`` is ``None`` in uns on that path.
    """

    result: pl.DataFrame | None
    store: DMCStore | None
    resumed: bool = False
    design: ContrastDesign | None = None

    @property
    def n_sites(self) -> int:
        """``len(result)`` when materialised, else the store's site count."""
        if self.result is not None:
            return len(self.result)
        return self.store.total_sites if self.store is not None else 0


def plan_run(md: MethylData, cfg: DMCConfig) -> DMCPlan:
    """Validate the request and fix every run-time choice.

    In order: TSV resolution, ``cfg.validate``, the formula / contrast
    dispatch (with the refusals that path cannot honour), the n=1 and
    union guards, ``"auto"`` test selection, ``apply_power_stack``,
    ``validate_resolved``, the one-time Fisher warning, ``unite`` from
    ``md.uns["unite"]`` and ``smooth_method`` from
    ``md.uns["smooth_params"]``.

    Raises ``ValueError`` in the same order as the pre-split ``tl.dmc``.
    Emits no result; everything after this stage is mechanical.
    """
    # Local imports: these helpers are shared with tl.dmr and live in tl,
    # which imports this module. Warnings from them are three frames below
    # the caller of tl.dmc (helper -> plan_run -> tl.dmc -> caller).
    from .tl import (
        _auto_test_simple,
        _check_n1_and_union_footgun,
        _resolve_auto_tsv,
        _warn_fisher_once,
    )

    tsv_path, tsv_full, tsv_alpha, tsv_is_auto = _resolve_auto_tsv(
        md,
        cfg.tsv,
        cfg.csv,
        default_name="dmc.significant.tsv",
        tsv_full=cfg.tsv_full,
        csv_full=cfg.csv_full,
        tsv_alpha=cfg.tsv_alpha,
        csv_alpha=cfg.csv_alpha,
        stacklevel=4,
    )
    tsv = (
        TsvPlan(path=tsv_path, full=tsv_full, alpha=tsv_alpha, is_auto=tsv_is_auto)
        if tsv_path is not None
        else None
    )

    cfg.validate()

    unite_info = md.uns.get("unite")
    unite = (unite_info is not None) and (unite_info.get("type") == "intersect")

    if cfg.formula is not None or cfg.contrast is not None:
        _refuse_contrast_knobs(cfg)
        return DMCPlan(
            cfg=cfg,
            mode="contrast",
            selected_test="glm_contrast",
            unite=unite,
            smooth_method=None,
            key="dmc_glm_contrast",
            tsv=tsv,
        )

    # Unconditional n=1 guard: applies whether test is "auto" or explicit.
    # _auto_test_simple raises ValueError when allow_n1=False; trigger that
    # check up front so explicit test="lr"/"fisher" with n<2 also gets
    # caught instead of silently running on degenerate data.
    _check_n1_and_union_footgun(
        md,
        cfg.allow_n1,
        cfg.min_samples_treatment,
        cfg.min_samples_control,
        stacklevel=4,
    )
    selected_test = (
        _auto_test_simple(md, allow_n1=cfg.allow_n1, stacklevel=4)
        if cfg.test == "auto"
        else cfg.test
    )

    # lr+ power-stack dispatch (1.0): resolves power_stack into
    # neighbour_combine / fdr_method / sep_fallback (dispersion="eb" is
    # already the default), then checks the knobs that only make sense
    # once the stack is applied.
    cfg = cfg.apply_power_stack(selected_test, min(len(md.treatment_ids), len(md.control_ids)))
    cfg.validate_resolved()

    if selected_test == "fisher":
        _warn_fisher_once(stacklevel=4)
    smooth_method = md.uns.get("smooth_params", {}).get("method") if cfg.use_smoothed else None

    key = f"dmc_{selected_test}_smoothed" if cfg.use_smoothed else f"dmc_{selected_test}"
    return DMCPlan(
        cfg=cfg,
        mode="binary",
        selected_test=selected_test,
        unite=unite,
        smooth_method=smooth_method,
        key=key,
        tsv=tsv,
    )


def _refuse_contrast_knobs(cfg: DMCConfig) -> None:
    """Refuse the knobs the formula / contrast path cannot honour.

    ``empirical_fdr`` and ``materialize=False`` raise up front rather than
    being silently ignored; an unknown ``power_stack`` raises as on the
    binary path, and a valid one is ignored with a notice because the GLM
    has none of the ``lr+`` knobs.
    """
    if cfg.empirical_fdr:
        # Same refusal as the DMR path: label shuffling invalidates
        # the stratified design that formula= encodes.
        raise ValueError(
            "empirical_fdr=True is not supported with the contrast / "
            "multi-group DMC path (label shuffling invalidates the "
            "stratified design). Use the binary treatment / control "
            "path or implement a custom stratified permutation."
        )
    if not cfg.materialize:
        # This path always assembles the full result onto md.varm;
        # refuse the argument rather than silently ignore it.
        raise ValueError(
            "materialize=False is not supported on the formula / contrast "
            "path yet: it always assembles the full per-CpG result onto "
            "md.varm. Re-run with materialize=True (the default)."
        )
    cfg.validate_resolved()
    if cfg.power_stack != "off":
        logger.info(
            "power_stack=%r is ignored on the formula / contrast path: the GLM has no lr+ knobs.",
            cfg.power_stack,
        )


def lookup_resume(md: MethylData, plan: DMCPlan) -> ResumeTicket | None:
    """Fingerprint the inputs and look for a prior run.

    Returns ``None`` when ``cfg.resumable`` is off or there is no analysis
    root to anchor the manifest on. On a hit the cached sidecar is already
    loaded into ``ticket.hit``; the orchestrator turns it into a
    :class:`DMCOutcome` with ``resumed=True`` and skips the engine. Only
    defined for binary mode; the contrast path has no resume support.
    """
    if not plan.cfg.resumable:
        return None
    from ._cache import input_signature, manifest_find

    stage_name = f"dmc_{plan.selected_test}"
    root = md.analysis_root or md.store
    if not root:
        return None
    signature = input_signature(
        md.store,
        sorted(md.treatment_ids),
        sorted(md.control_ids),
        plan.cfg.resume_signature_params(selected_test=plan.selected_test, unite=plan.unite),
    )
    hit: pl.DataFrame | None = None
    prior = manifest_find(root, stage_name)
    if prior is not None and prior.get("input_sig") == signature:
        sidecar = Path(prior["output_path"])
        if not sidecar.is_absolute():
            sidecar = Path(root) / sidecar
        if sidecar.exists():
            logger.info("[resume] %s: loading cached result from %s", stage_name, sidecar)
            hit = pl.read_parquet(str(sidecar))
    return ResumeTicket(stage_name=stage_name, root=root, signature=signature, hit=hit)


@contextmanager
def open_input_store(md: MethylData, plan: DMCPlan) -> Iterator[str]:
    """Yield the methylstore path the engine reads.

    ``md.store`` normally. With ``use_smoothed=True`` this builds the
    pseudo-count store in a temporary directory, after the existing
    ``DeprecationWarning`` and the ``smooth_path`` check, and removes it
    when the ``with`` block exits, on both normal exit and exceptions.
    """
    if not plan.cfg.use_smoothed:
        yield md.store
        return

    # Frames above this one: contextmanager.__enter__, tl.dmc, its caller.
    warnings.warn(
        "use_smoothed=True (pseudo-count transform of raw reads via "
        "BSmooth) is NOT equivalent to DSS's smoothing=TRUE -- it's "
        "too aggressive (washes out per-CpG signal at default BSmooth "
        "parameters). For DSS-style behavior, use smoothing=True "
        "(applies DSS's uniform-box +/-smoothing_span_bp//2 moving "
        "average to each sample's raw counts before they hit the "
        "test, matching DMLfit.multiFactor(smoothing=TRUE)). The "
        "use_smoothed pseudo-count path will be removed in a future "
        "minor release.",
        DeprecationWarning,
        stacklevel=4,
    )
    if "smooth_path" not in md.uns:
        raise ValueError(
            "use_smoothed=True requires ep.pp.smooth(md) first "
            "(either method='gaussian' or 'bsmooth'). The smoothed "
            "sidecar at md.uns['smooth_path'] is the input to the "
            "pseudo-count transform that feeds the DMC test."
        )
    from ._smoothed_store import build_smoothed_pseudo_count_store

    with tempfile.TemporaryDirectory(prefix="epykit_dmc_smoothed_") as tmp_dir:
        build_smoothed_pseudo_count_store(
            raw_store=Path(md.store),
            smooth_store=Path(md.uns["smooth_path"]),
            samples=md.obs.get_column("sample_id").to_list(),
            out_dir=Path(tmp_dir),
        )
        yield tmp_dir


def run_engine(md: MethylData, plan: DMCPlan, store_path: str) -> DMCStore:
    """Stream the per-CpG test and apply the FDR correction in place.

    One call to ``process_chromosomes_dmc(..., return_store=True)`` with the
    knobs read from ``plan.cfg``, then ``apply_multiple_testing_correction``
    on the store. Returns the streaming store; nothing is materialised here.
    The per-chrom parquet directory is the source of truth, so both the
    correction and downstream DMR stream chromosomes from disk.
    """
    cfg = plan.cfg
    store = process_chromosomes_dmc(
        methylstore_path=store_path,
        samples_treatment=md.treatment_ids,
        samples_control=md.control_ids,
        test=plan.selected_test,
        chromosomes=cfg.chromosomes,
        unite=plan.unite,
        min_samples_treatment=cfg.min_samples_treatment,
        min_samples_control=cfg.min_samples_control,
        dispersion=cfg.dispersion,
        reference=cfg.reference,
        backend=cfg.backend,
        n_workers=cfg.n_workers,
        glm_backend=cfg.glm_backend,
        return_store=True,
        smoothing=cfg.smoothing,
        smoothing_span_bp=cfg.smoothing_span_bp,
        sep_fallback=cfg.sep_fallback,
        sep_threshold=cfg.sep_threshold,
    )
    return apply_multiple_testing_correction(store, method=cfg.fdr_method)


def post_process(
    md: MethylData, plan: DMCPlan, store: DMCStore, store_path: str
) -> pl.DataFrame | None:
    """Materialise and run the eager post-processors.

    Returns ``store.to_dataframe()`` with ``neighbour_combine`` and
    ``empirical_fdr`` applied when they are on, or ``None`` when
    ``materialize=False`` (``validate_resolved`` already refused the
    combination of ``materialize=False`` with either post-processor).
    ``store_path`` is needed because the permutation null re-reads the
    input store, so this stage runs while a smoothed temp store is alive.
    """
    cfg = plan.cfg
    if not cfg.materialize:
        # Keep only the streaming DMCStore handle (O(largest chromosome)
        # end-to-end). md.dmc materialises on demand from store_path.
        return None

    # Materialise the full DataFrame for md.varm back-compat (plot.py /
    # report.py / pl modules consume md.dmc as a DataFrame). With
    # chrom/strand stored as pl.Enum this is roughly 700 MB at 22M rows.
    result = store.to_dataframe()

    # Neighbour-aware p-value combining (RADMeth-style, since 0.7.1). The
    # canonical `pvalue` / `qvalue` columns remain the raw per-CpG values;
    # the combined values land in `pvalue_combined` / `qvalue_combined`
    # (plus `pvalue_combined_n_neighbours` / `qvalue_combined_reject` as
    # audit columns). Sites without enough neighbours keep their raw
    # p-value identity.
    if cfg.neighbour_combine and len(result) > 0:
        result = combine_neighbour_pvalues(result, neighbour_bp=cfg.neighbour_bp)
        result = apply_multiple_testing_correction(
            result,
            method=cfg.fdr_method,
            pvalue_col="pvalue_combined",
            qvalue_col="qvalue_combined",
        )

    if cfg.empirical_fdr and len(result) > 0:
        result = empirical_fdr_for_dmc(
            methylstore_path=store_path,
            samples_treatment=md.treatment_ids,
            samples_control=md.control_ids,
            observed_dmc=result,
            n_perm=cfg.n_perm,
            seed=cfg.perm_seed,
            n_jobs=cfg.perm_n_jobs,
            test=plan.selected_test,
            chromosomes=cfg.chromosomes,
            unite=plan.unite,
            min_samples_treatment=cfg.min_samples_treatment,
            min_samples_control=cfg.min_samples_control,
            dispersion=cfg.dispersion,
            reference=cfg.reference,
            # Engine knobs that can overwrite the per-site p-value MUST be
            # applied identically in observed and null runs, otherwise the
            # Westfall-Young statistic compares deflated observed p-values
            # against an un-deflated null pool.
            sep_fallback=cfg.sep_fallback,
            sep_threshold=cfg.sep_threshold,
            smoothing=cfg.smoothing,
            smoothing_span_bp=cfg.smoothing_span_bp,
        )
    return result


def run_contrast(md: MethylData, plan: DMCPlan) -> DMCOutcome:
    """The formula / contrast run end to end.

    Always uses ``test="glm_contrast"`` and ALL samples in ``md.obs`` order
    (not the binary case / control split), so the design matrix matches
    ``md.obs`` row-for-row. Builds the design, resolves the contrast,
    derives the per-level labels and the back-compat case / control split,
    runs the engine, applies the FDR correction and materialises. Always
    materialises. Returns the outcome with :class:`ContrastDesign` attached
    so :func:`publish` can record it.
    """
    from ._glm import build_design, resolve_contrast

    cfg = plan.cfg
    if not md.obs.height:
        raise ValueError("md.obs is empty; cannot build a design matrix.")
    samples_all = md.obs.get_column("sample_id").to_list()
    treatment_col = cfg.treatment_col

    # Build design without requiring a treatment column if we have a
    # formula that doesn't reference one. The `treatment_col` default
    # ("treatment") is only required when the binary path would have used
    # it; here we let the formula speak.
    need_treatment = (treatment_col in md.obs.columns) and (
        cfg.formula is None or treatment_col in cfg.formula
    )
    design_full, _design_reduced, _coef_idx, term_names, formula_used, design_info = build_design(
        md.obs,
        samples_ordered=samples_all,
        formula=cfg.formula,
        covariates=cfg.covariates,
        treatment_col=treatment_col,
        require_treatment_col=need_treatment,
        return_design_info=True,
        reference_level=cfg.reference_level,
    )

    # Resolve the contrast against the design. The default is a single-coef
    # contrast on `treatment_col` (a `formula=` for covariate adjustment
    # with no explicit contrast).
    contrast = cfg.contrast if cfg.contrast is not None else treatment_col
    contrast_matrix, contrast_label = resolve_contrast(
        contrast, term_names, design_info=design_info
    )

    # Per-sample level labels for the multi-group output schema: when the
    # contrast names a categorical obs column, emit per-level labels for
    # the downstream mean_beta_<level> columns.
    group_labels: list[str] | None = None
    if isinstance(contrast, str) and contrast in md.obs.columns:
        col = md.obs.get_column(contrast)
        if col.dtype == pl.Utf8 or col.dtype == pl.Categorical:
            group_labels = col.cast(pl.Utf8).to_list()

    # Case / control split for the backwards-compatible binary columns. If
    # treatment_col is on obs and carries a numeric 0/1 signal, use it;
    # otherwise leave both empty so mean_beta_case/control remain NaN
    # (uninterpretable for multi-group).
    samples_case: list[str] = []
    samples_control: list[str] = []
    if treatment_col in md.obs.columns:
        try:
            mask_treat = (
                md.obs.get_column(treatment_col).cast(pl.Float64, strict=False) == 1
            ).to_list()
            # Both lists are columns of md.obs, so they have the same length.
            samples_case = [s for s, m in zip(samples_all, mask_treat, strict=True) if m]
            samples_control = [s for s, m in zip(samples_all, mask_treat, strict=True) if not m]
        except Exception as exc:
            logger.warning(
                "Could not derive case/control split from treatment column "
                "%r; mean_beta_case/control will be NaN: %s",
                treatment_col,
                exc,
            )

    store = process_chromosomes_dmc(
        methylstore_path=md.store,
        samples_treatment=samples_case,
        samples_control=samples_control,
        test="glm_contrast",
        chromosomes=cfg.chromosomes,
        unite=plan.unite,
        min_samples_treatment=cfg.min_samples_treatment,
        min_samples_control=cfg.min_samples_control,
        dispersion=cfg.dispersion,
        reference=cfg.reference,
        design_full=design_full,
        contrast_matrix=contrast_matrix,
        contrast_label=contrast_label,
        samples_all_ordered=samples_all,
        group_labels_per_sample=group_labels,
        return_store=True,
    )
    store = apply_multiple_testing_correction(store, method=cfg.fdr_method)
    return DMCOutcome(
        result=store.to_dataframe(),
        store=store,
        design=ContrastDesign(
            formula_used=formula_used,
            contrast_label=contrast_label,
            design_terms=term_names,
        ),
    )


def publish(md: MethylData, plan: DMCPlan, outcome: DMCOutcome) -> None:
    """Write ``md.varm[plan.key]`` and ``md.uns["dmc"]``.

    The only writer of the uns record. ``md.varm`` gets the table only when
    the run materialised; otherwise the DMCStore is the source of truth and
    ``md.get_dmc()`` / ``md.dmc`` resolve it on demand from ``store_path``.
    """
    if outcome.result is not None:
        md.varm[plan.key] = outcome.result
    design = outcome.design
    md.uns["dmc"] = plan.cfg.to_uns(
        test_used=plan.selected_test,
        n_sites=outcome.n_sites,
        materialized=outcome.result is not None,
        unite=plan.unite,
        last_key=plan.key,
        store_path=str(outcome.store.path) if outcome.store is not None else None,
        resumed=outcome.resumed,
        smooth_method=plan.smooth_method,
        formula=design.formula_used if design is not None else None,
        contrast=design.contrast_label if design is not None else None,
        design_terms=design.design_terms if design is not None else None,
    )


def persist_resume(
    md: MethylData, plan: DMCPlan, ticket: ResumeTicket, outcome: DMCOutcome
) -> None:
    """Write the sidecar parquet and the manifest entry after a miss.

    Only when the run materialised. Best effort: an ``OSError`` logs a
    warning and does not propagate (the in-memory result is still valid).
    Runs after :func:`publish` so the order of side effects is unchanged.
    """
    result = outcome.result
    if result is None:
        return
    from ._cache import manifest_append

    cfg = plan.cfg
    try:
        sidecar_dir = Path(ticket.root) / ".epykit_results"
        sidecar_dir.mkdir(parents=True, exist_ok=True)
        sidecar = sidecar_dir / f"{ticket.stage_name}.parquet"
        result.write_parquet(str(sidecar))
        manifest_append(
            ticket.root,
            ticket.stage_name,
            params={
                "test": plan.selected_test,
                "unite": plan.unite,
                "min_samples_treatment": cfg.min_samples_treatment,
                "min_samples_control": cfg.min_samples_control,
                "dispersion": cfg.dispersion,
                "reference": cfg.reference,
                "empirical_fdr": cfg.empirical_fdr,
            },
            input_sig=ticket.signature,
            output_path=str(sidecar),
            extra={"n_sites": len(result)},
        )
    except OSError as exc:
        logger.warning("[resume] failed to persist %s sidecar: %s", ticket.stage_name, exc)


def finish(md: MethylData, plan: DMCPlan, outcome: DMCOutcome) -> None:
    """The ``log2_odds_ratio`` FutureWarning and the TSV export.

    Emits the deprecation notice once per ``tl.dmc`` call on every path,
    including the resume cache hit, then writes the TSV when ``plan.tsv``
    is set and a table exists, or logs the ``materialize=False`` skip.
    """
    # Frames above this one: tl.dmc, its caller.
    warnings.warn(_LOG2_ODDS_RATIO_NOTICE, FutureWarning, stacklevel=3)

    tsv = plan.tsv
    if tsv is None:
        return
    if outcome.result is None:
        logger.info(
            "materialize=False: skipping DMC TSV auto-export (no "
            "in-memory result table). Export later via "
            "ep.export.dmc_to_tsv(md) or re-run with materialize=True."
        )
        return
    from .export import dmc_to_tsv
    from .tl import _emit_result_tsv

    _emit_result_tsv(
        lambda: dmc_to_tsv(md, tsv.path, alpha=tsv.alpha, full=tsv.full),
        tsv.path,
        is_auto=tsv.is_auto,
    )
