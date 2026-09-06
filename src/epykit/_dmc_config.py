"""The keyword knobs of ``tl.dmc`` carried as one frozen record.

``DMCConfig`` holds every keyword argument ``ep.tl.dmc`` accepts, with the
same names and defaults. The orchestrator builds one per call and uses it
to validate the request, resolve the ``lr+`` power stack, fingerprint the
``resumable=True`` lookup and write the ``md.uns["dmc"]`` metadata record,
so those four pieces of logic have a single owner. The public
``ep.tl.dmc(...)`` signature is unchanged.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, replace
from typing import Any, Literal

logger = logging.getLogger(__name__)

PowerStack = Literal["auto", "lr+", "conservative", "off"]
_POWER_STACK_MODES = frozenset({"auto", "lr+", "conservative", "off"})

# Engines removed in 0.7.5. Each raises ValueError with the migration hint.
_REMOVED_ENGINES = {
    "logit_t": (
        "test='logit_t' was removed in 0.7.5 (miscalibrated near β=0/1). "
        "Use test='welch_t' for the replicate-aware β-mean test or "
        "test='lr' for the recommended default."
    ),
    "bb_lr": (
        "test='bb_lr' was removed in 0.7.5 (TPR < 8% at n ≤ 4 + a "
        "dispersion-df bug). Use test='lr' (recommended) which uses "
        "the same quasi-binomial dispersion but pools counts per group "
        "for higher power at small n."
    ),
    "score": (
        "test='score' was removed in 0.7.5 (strictly dominated by "
        "test='lr' in finite samples; asymptotically equivalent under "
        "H0). Switch test='score' -> test='lr'; output schema is "
        "identical."
    ),
    "cmh": (
        "test='cmh' was removed in 0.7.5 (stratification semantics "
        "confusing; dominated by GLM with batch covariate). For "
        "stratified analysis use tl.dmc(formula='~ group + batch'), "
        "which gives proper dispersion correction and handles "
        "continuous covariates."
    ),
}


@dataclass(frozen=True)
class DMCConfig:
    """Every keyword knob of ``tl.dmc``; see that docstring for meanings.

    Two inputs are normalised on construction, as ``tl.dmc`` has always
    done: ``min_samples_treatment=None`` reads as 0, and a bool
    ``power_stack`` is the ``"lr+"`` / ``"off"`` alias.
    """

    test: str = "auto"
    chromosomes: list[str] | None = None
    min_samples_treatment: int | None = None
    min_samples_control: int = 0
    dispersion: str = "eb"
    reference: str = "adaptive"
    allow_n1: bool = False
    formula: str | None = None
    contrast: Any = None
    covariates: list[str] | None = None
    treatment_col: str = "treatment"
    empirical_fdr: bool = False
    n_perm: int = 100
    perm_seed: int = 42
    perm_n_jobs: int = 1
    backend: str = "sequential"
    n_workers: int | None = None
    glm_backend: str = "cpu"
    resumable: bool = False
    materialize: bool = True
    use_smoothed: bool = False
    smoothing: bool = False
    smoothing_span_bp: int = 500
    fdr_method: str = "fdr_bh"
    neighbour_combine: bool = False
    neighbour_bp: int = 500
    sep_fallback: bool = False
    sep_threshold: float = 0.9
    power_stack: PowerStack | bool = "off"
    reference_level: str | None = None
    tsv: str | bool | None = None
    tsv_full: bool = False
    tsv_alpha: float = 0.05
    csv: str | None = None
    csv_full: bool = False
    csv_alpha: float = 0.05

    def __post_init__(self) -> None:
        if self.min_samples_treatment is None:
            object.__setattr__(self, "min_samples_treatment", 0)
        if isinstance(self.power_stack, bool):
            object.__setattr__(self, "power_stack", "lr+" if self.power_stack else "off")

    def validate(self) -> None:
        """Reject the engines removed in 0.7.5 with their migration hints.

        Runs before the formula / contrast dispatch, so a removed engine
        name is refused on every path.
        """
        message = _REMOVED_ENGINES.get(self.test)
        if message is not None:
            raise ValueError(message)

    def validate_resolved(self) -> None:
        """Checks that need the power stack resolved first (binary path only).

        ``power_stack`` must be a known mode, and ``materialize=False`` cannot
        run the eager-only post-processors that need the full in-memory
        result table. ``neighbour_combine`` may have been switched on by the
        stack itself, which is why this runs after :meth:`apply_power_stack`.
        """
        if self.power_stack not in _POWER_STACK_MODES:
            raise ValueError(
                f"power_stack must be one of {{'auto','lr+','conservative','off'}} "
                f"or a bool; got {self.power_stack!r}"
            )
        if not self.materialize:
            incompatible = [
                name
                for name, on in (
                    ("neighbour_combine", self.neighbour_combine),
                    ("empirical_fdr", self.empirical_fdr),
                    ("use_smoothed", self.use_smoothed),
                )
                if on
            ]
            if incompatible:
                raise ValueError(
                    "materialize=False keeps only the streaming DMCStore handle "
                    "and cannot run features that post-process the full in-memory "
                    f"result table: {', '.join(incompatible)}. Re-run with "
                    "materialize=True (the default) for these, or disable them "
                    "(neighbour_combine may have been auto-enabled by "
                    "power_stack='lr+')."
                )

    def apply_power_stack(self, selected_test: str, min_n: int) -> DMCConfig:
        """Resolve the ``lr+`` power stack into its component knobs.

        Returns a copy with ``neighbour_combine``, ``fdr_method`` and
        ``sep_fallback`` switched on for the ``lr`` engine: ``"lr+"`` and
        ``"auto"`` engage at any n, ``"conservative"`` only at n <= 2,
        ``"off"`` leaves the user's values alone. The fourth component,
        ``dispersion="eb"``, is already the default. Knobs the user set
        themselves are kept.
        """
        if selected_test != "lr" or self.power_stack not in {"auto", "lr+", "conservative"}:
            return self
        if self.power_stack == "conservative" and min_n > 2:
            return self
        changes: dict[str, Any] = {}
        if not self.neighbour_combine:
            changes["neighbour_combine"] = True
            logger.info(
                "Auto-enabling neighbour_combine (lr+ stack, "
                "power_stack=%s, n=%d). Pass power_stack='off' to "
                "disable.",
                self.power_stack,
                min_n,
            )
        if self.fdr_method == "fdr_bh":
            changes["fdr_method"] = "fdr_tsbh"
            logger.info(
                "Auto-switching fdr_method 'fdr_bh' -> 'fdr_tsbh' (lr+ stack, power_stack=%s).",
                self.power_stack,
            )
        if not self.sep_fallback:
            changes["sep_fallback"] = True
            logger.info(
                "Auto-enabling sep_fallback (lr+ stack, power_stack=%s).",
                self.power_stack,
            )
        return replace(self, **changes) if changes else self

    def resume_signature_params(self, *, selected_test: str, unite: bool) -> dict[str, Any]:
        """The params half of the ``resumable=True`` fingerprint.

        Every knob that changes engine output is listed, including the
        ``lr+`` stack knobs: leaving one out would let a parameter sweep
        silently reuse a cached result computed at different values.
        """
        return {
            "test": selected_test,
            "chromosomes": self.chromosomes,
            "unite": unite,
            "min_samples_treatment": self.min_samples_treatment,
            "min_samples_control": self.min_samples_control,
            "dispersion": self.dispersion,
            "reference": self.reference,
            "empirical_fdr": self.empirical_fdr,
            "n_perm": self.n_perm if self.empirical_fdr else None,
            "perm_seed": self.perm_seed if self.empirical_fdr else None,
            "power_stack": self.power_stack,
            "sep_fallback": self.sep_fallback,
            "sep_threshold": self.sep_threshold,
            "neighbour_combine": self.neighbour_combine,
            "neighbour_bp": self.neighbour_bp,
            "fdr_method": self.fdr_method,
        }

    def to_uns(
        self,
        *,
        test_used: str,
        n_sites: int,
        materialized: bool,
        unite: bool,
        last_key: str,
        store_path: str | None,
        resumed: bool = False,
        smooth_method: str | None = None,
        formula: str | None = None,
        contrast: str | None = None,
        design_terms: list[str] | None = None,
    ) -> dict[str, Any]:
        """The ``md.uns["dmc"]`` record, with the same keys on every path.

        Readers: ``MethylData.get_dmc`` / ``.dmc_store`` / ``.save``
        (``last_key``, ``store_path``), ``report.py`` (``test_used``,
        ``test_requested``, ``fdr_method``), ``multiqc_export`` and
        ``pl.dashboard`` (``test_used``, ``n_sites``, ``unite``),
        ``export.export_tables`` and the CLI (``last_key``).

        The formula / contrast path passes ``design_terms`` together with
        the resolved ``formula`` and ``contrast`` label. That path does not
        consume the power stack, permutation FDR or smoothing knobs, so
        they are recorded as ``None`` there; the binary path records the
        contrast fields as ``None``. ``store_path`` is ``None`` when no
        DMCStore was opened (the resume cache hit).
        """
        binary = design_terms is None

        def binary_only(value: Any) -> Any:
            return value if binary else None

        return {
            "test_requested": self.test,
            "test_used": test_used,
            "n_sites": n_sites,
            "materialized": materialized,
            "unite": unite,
            "min_samples_treatment": self.min_samples_treatment,
            "min_samples_control": self.min_samples_control,
            "dispersion": self.dispersion,
            "reference": self.reference,
            "empirical_fdr": binary_only(self.empirical_fdr),
            "n_perm": binary_only(self.n_perm if self.empirical_fdr else None),
            "perm_seed": binary_only(self.perm_seed if self.empirical_fdr else None),
            "power_stack": binary_only(self.power_stack),
            "sep_fallback": binary_only(self.sep_fallback),
            "sep_threshold": binary_only(self.sep_threshold),
            "neighbour_combine": binary_only(self.neighbour_combine),
            "neighbour_bp": binary_only(self.neighbour_bp),
            "fdr_method": self.fdr_method,
            # Explicit pointer so MethylData.get_dmc() / .dmc resolve to the
            # table just written, whatever other tests ran in the session.
            "last_key": last_key,
            "use_smoothed": binary_only(self.use_smoothed),
            "smooth_method": binary_only(smooth_method),
            # DSS-style count smoothing; the span only means something when
            # smoothing is on.
            "smoothing": binary_only(bool(self.smoothing)),
            "smoothing_span_bp": binary_only(
                int(self.smoothing_span_bp) if self.smoothing else None
            ),
            # Persistent per-chromosome DMCStore, so tl.dmr can stream
            # chromosomes from disk instead of holding the full table.
            "store_path": store_path,
            "resumed": resumed,
            "formula": formula,
            "contrast": contrast,
            "design_terms": design_terms,
            "covariates": None if binary else (list(self.covariates) if self.covariates else None),
            "treatment_col": None if binary else self.treatment_col,
        }
