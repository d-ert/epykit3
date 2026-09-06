"""The DMC engine registry: which engines exist and the facts about them.

One frozen :class:`EngineSpec` per engine, holding only the facts other
modules consume: whether the engine is a public ``test=`` choice, whether
the ``lr+`` power stack resolves onto it, and which effect-size column it
emits. The numerical code stays in ``dmc.py``; the sample-size warnings,
the ``"auto"`` selection rule and the design checks stay with the code
that runs them.

This module imports only the standard library, so ``_dmc_config``, the
stages, ``tl`` and ``cli`` can import it without pulling in ``dmc.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType


@dataclass(frozen=True)
class EngineSpec:
    """The facts about one DMC engine that other modules read."""

    name: str
    public: bool
    """Offered as a ``test=`` / ``--test`` choice. ``glm_contrast`` is internal."""
    power_stack_applies: bool
    """Whether :meth:`DMCConfig.apply_power_stack` resolves the ``lr+`` stack onto it."""
    effect_column: str
    """The effect-size column the engine emits: a pooled log2 odds ratio for
    the count engines, the logit coefficient in log2 units for the GLMs."""


# Public engines first, in CLI choice order.
ENGINES = MappingProxyType(
    {
        spec.name: spec
        for spec in (
            EngineSpec(
                name="lr",
                public=True,
                power_stack_applies=True,
                effect_column="log2_odds_ratio_pooled",
            ),
            EngineSpec(
                name="glm",
                public=True,
                power_stack_applies=False,
                effect_column="coef_treatment_log2",
            ),
            EngineSpec(
                name="welch_t",
                public=True,
                power_stack_applies=False,
                effect_column="log2_odds_ratio_pooled",
            ),
            EngineSpec(
                name="fisher",
                public=True,
                power_stack_applies=False,
                effect_column="log2_odds_ratio_pooled",
            ),
            EngineSpec(
                name="glm_contrast",
                public=False,
                power_stack_applies=False,
                effect_column="coef_treatment_log2",
            ),
        )
    }
)

PUBLIC_ENGINES: tuple[str, ...] = tuple(spec.name for spec in ENGINES.values() if spec.public)
"""The ``--test`` choice list of both CLI commands, in the documented order."""

# Engines removed in 0.7.5. Each raises ValueError with the migration hint.
REMOVED_ENGINES = MappingProxyType(
    {
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
)


def engine_spec(test: str) -> EngineSpec:
    """Look up a resolved engine name.

    Raises ``ValueError`` with the migration hint for an engine removed in
    0.7.5, and with the public choice list for any other unknown name.
    ``"auto"`` is not an engine: the high-level path resolves it before
    calling this.
    """
    spec = ENGINES.get(test)
    if spec is not None:
        return spec
    hint = REMOVED_ENGINES.get(test)
    if hint is not None:
        raise ValueError(hint)
    raise ValueError(f"Unknown DMC test {test!r}. Choose one of: {', '.join(PUBLIC_ENGINES)}.")
