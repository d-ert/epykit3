"""PROTOTYPE, not wired in: the stage split of ``tl.dmc``.

This module is a stub to react to for the map ticket "Stage boundaries for
the tl.dmc split". Every function carries its final signature and docstring
and a pointer to the ``tl.py`` lines it absorbs (line numbers as of
``900ea6c``), but no body has moved yet. ``dmc_sketch`` at the bottom is the
orchestrator as it would read after the split; it replaces ``tl.py`` lines
640 to 991 and ``_run_dmc_contrast`` (994 to 1142).

Design rules the split keeps:

- the public ``ep.tl.dmc(...)`` signature and docstring do not change;
- ``DMCConfig`` stays the single carrier of the knobs, and every stage takes
  the resolved config through ``DMCPlan`` rather than loose arguments;
- ``md.uns["dmc"]`` is written in exactly one place, ``publish``, through
  ``DMCConfig.to_uns``;
- the resume cache hit and the formula / contrast run are outcomes of the
  same shape as the binary run, so ``publish`` and ``finish`` do not branch
  on how the result was produced;
- nothing new is materialised: ``run_engine`` returns the streaming
  ``DMCStore`` and ``post_process`` is the only stage that may hold the
  full table, exactly where ``tl.dmc`` holds it today.
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Literal

import polars as pl

from ._dmc_config import DMCConfig
from .dmc import DMCStore
from .methyldata import MethylData

Mode = Literal["binary", "contrast"]


@dataclass(frozen=True)
class TsvPlan:
    """Where the significant-DMC TSV goes, resolved once (tl.py 640 to 649).

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
    """The engine after ``"auto"`` resolution; ``"glm_contrast"`` in contrast mode."""
    canonical_test: str
    """``_canonicalise_test_name(selected_test)``; the ``test_used`` recorded in uns."""
    unite: bool
    smooth_method: str | None
    key: str
    tsv: TsvPlan | None


@dataclass(frozen=True)
class ContrastDesign:
    """What the contrast run resolved from ``md.obs`` (tl.py 1036 to 1102).

    Carried on the outcome so ``publish`` can record ``formula``,
    ``contrast`` and ``design_terms`` without re-resolving them.
    """

    formula_used: str
    contrast_label: str
    design_terms: list[str]


@dataclass(frozen=True)
class ResumeTicket:
    """The ``resumable=True`` lookup, hit or miss (tl.py 702 to 757).

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
    hit, which is why ``store_path`` is ``None`` in uns on that path today.
    """

    result: pl.DataFrame | None
    store: DMCStore | None
    resumed: bool = False
    design: ContrastDesign | None = None

    @property
    def n_sites(self) -> int:
        """``len(result)`` when materialised, else ``store.total_sites`` (tl.py 907 to 913)."""
        raise NotImplementedError


def plan_run(md: MethylData, cfg: DMCConfig) -> DMCPlan:
    """Validate the request and fix every run-time choice.

    Absorbs tl.py 640 to 700 in today's order: TSV resolution, ``cfg.validate``,
    the formula / contrast dispatch, the n=1 and union guards, ``"auto"``
    test selection, ``apply_power_stack``, ``validate_resolved``, the
    one-time Fisher warning, ``unite`` from ``md.uns["unite"]`` and
    ``smooth_method`` from ``md.uns["smooth_params"]``. Contrast mode also
    absorbs the refusals at the top of ``_run_dmc_contrast`` (1010 to 1034).

    Raises the same ``ValueError`` s as today, in the same order. Emits no
    result; everything after this stage is mechanical.
    """
    raise NotImplementedError


def lookup_resume(md: MethylData, plan: DMCPlan) -> ResumeTicket | None:
    """Fingerprint the inputs and look for a prior run (tl.py 702 to 736).

    Returns ``None`` when ``cfg.resumable`` is off. On a hit the cached
    sidecar is already loaded into ``ticket.hit``; the orchestrator turns it
    into a :class:`DMCOutcome` with ``resumed=True`` and skips the engine.
    Only defined for binary mode; the contrast path has no resume support
    today and keeps not having it.
    """
    raise NotImplementedError


@contextmanager
def open_input_store(md: MethylData, plan: DMCPlan) -> Iterator[str]:
    """Yield the methylstore path the engine reads (tl.py 759 to 801, 898 to 900).

    ``md.store`` normally. With ``use_smoothed=True`` this builds the
    pseudo-count store in a temporary directory, after the existing
    ``DeprecationWarning`` and the ``smooth_path`` check, and removes it on
    exit; the ``try / finally`` around the engine call becomes this
    ``with`` block. The warning's ``stacklevel`` moves from 2 to 3 so it
    still points at the caller of ``tl.dmc``.
    """
    raise NotImplementedError
    yield  # pragma: no cover


def run_engine(md: MethylData, plan: DMCPlan, store_path: str) -> DMCStore:
    """Stream the per-CpG test and apply BH in place (tl.py 809 to 829).

    One call to ``process_chromosomes_dmc(..., return_store=True)`` with the
    knobs read from ``plan.cfg``, then ``apply_multiple_testing_correction``
    on the store. Returns the streaming store; nothing is materialised here.
    This is the seam the engine registry ticket will later own.
    """
    raise NotImplementedError


def post_process(
    md: MethylData, plan: DMCPlan, store: DMCStore, store_path: str
) -> pl.DataFrame | None:
    """Materialise and run the eager post-processors (tl.py 831 to 897).

    Returns ``store.to_dataframe()`` with ``neighbour_combine`` and
    ``empirical_fdr`` applied when they are on, or ``None`` when
    ``materialize=False`` (``validate_resolved`` already refused the
    combination of ``materialize=False`` with either post-processor, so the
    guards inside are the same short-circuits as today). ``store_path`` is
    needed because the permutation null re-reads the input store.
    """
    raise NotImplementedError


def run_contrast(md: MethylData, plan: DMCPlan) -> DMCOutcome:
    """The formula / contrast run end to end (tl.py 1036 to 1128).

    Builds the design, resolves the contrast, derives the per-level labels
    and the back-compat case / control split, runs
    ``process_chromosomes_dmc(test="glm_contrast", ...)``, applies BH and
    materialises. Always materialises, as today. Returns the outcome with
    :class:`ContrastDesign` attached so ``publish`` can record it.
    """
    raise NotImplementedError


def publish(md: MethylData, plan: DMCPlan, outcome: DMCOutcome) -> None:
    """Write ``md.varm[plan.key]`` and ``md.uns["dmc"]`` (tl.py 736 to 748, 905 to 922, 1130 to 1142).

    The only writer of the uns record. Calls ``plan.cfg.to_uns`` with
    ``test_used=plan.canonical_test``, ``n_sites=outcome.n_sites``,
    ``materialized=outcome.result is not None``, ``last_key=plan.key``,
    ``store_path`` from ``outcome.store`` (``None`` on a resume hit),
    ``resumed=outcome.resumed``, ``smooth_method=plan.smooth_method`` and the
    three contrast fields from ``outcome.design`` when present.
    """
    raise NotImplementedError


def persist_resume(
    md: MethylData, plan: DMCPlan, ticket: ResumeTicket, outcome: DMCOutcome
) -> None:
    """Write the sidecar parquet and the manifest entry (tl.py 924 to 961).

    Only after a miss with a materialised result. Best effort: an ``OSError``
    logs a warning and does not propagate, as today. Runs after
    :func:`publish` so the order of side effects is unchanged.
    """
    raise NotImplementedError


def finish(md: MethylData, plan: DMCPlan, outcome: DMCOutcome) -> None:
    """The ``log2_odds_ratio`` FutureWarning and the TSV export (tl.py 963 to 991, 656 to 674).

    Emits the deprecation warning once per call (``stacklevel`` 3 from
    here), then writes the TSV when ``plan.tsv`` is set and a table exists,
    or logs the ``materialize=False`` skip. Today the resume hit returns
    before the warning; see the open question on whether it should.
    """
    raise NotImplementedError


def dmc_sketch(md: MethylData, cfg: DMCConfig) -> None:
    """The orchestrator after the split. Replaces tl.py 640 to 991 inside ``tl.dmc``.

    ``tl.dmc`` keeps its 36 explicit parameters and its docstring, builds
    ``cfg`` as it does now (601 to 638), then runs this body.
    """
    plan = plan_run(md, cfg)

    if plan.mode == "contrast":
        outcome = run_contrast(md, plan)
        publish(md, plan, outcome)
        finish(md, plan, outcome)
        return

    ticket = lookup_resume(md, plan)
    if ticket is not None and ticket.hit is not None:
        outcome = DMCOutcome(result=ticket.hit, store=None, resumed=True)
        publish(md, plan, outcome)
        finish(md, plan, outcome)
        return

    with open_input_store(md, plan) as store_path:
        store = run_engine(md, plan, store_path)
        result = post_process(md, plan, store, store_path)
    outcome = DMCOutcome(result=result, store=store)

    publish(md, plan, outcome)
    if ticket is not None:
        persist_resume(md, plan, ticket, outcome)
    finish(md, plan, outcome)
