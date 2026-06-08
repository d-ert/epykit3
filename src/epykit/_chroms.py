"""Canonical chromosome definitions shared across analysis, ingestion, plotting.

Single source of truth for "which contigs are main-assembly chromosomes" so the
DMC/DMR default filter (``canonical_only``), the ingestion opt-in, and the
plotters all agree on one definition. Kept dependency-free on purpose: it is
imported by low-level I/O (``io``/``convert``) as well as the analysis and
plotting layers, so it must not pull in heavy modules.

Naming conventions
------------------
The predicate is naming-convention-agnostic: a leading ``chr`` prefix is
stripped and the remainder compared case-insensitively, so it accepts both the
UCSC convention (``chr1``, ``chrM``) and the Ensembl convention (``1``, ``MT``).

Scope
-----
The canonical set is mammalian numeric naming (human ``1``-``22``; mouse
``1``-``19`` is a subset), the sex chromosomes ``X``/``Y``, and the
mitochondrion ``M``/``MT`` (kept -- see
``docs/superpowers/specs/2026-06-08-canonical-chrom-filter-design.md``).
Non-mammalian assemblies with roman-numeral or named chromosomes are *not*
recognised; callers analysing those must pass ``canonical_only=False`` (or an
explicit ``chromosomes=`` list). The per-run audit log makes an over-eager drop
visible.
"""

from __future__ import annotations

import logging
from typing import Iterable

logger = logging.getLogger(__name__)

# Core chromosome identifiers (no ``chr`` prefix), compared case-insensitively.
# Covers UCSC and Ensembl naming once the optional prefix is stripped.
CANONICAL_CHROM_CORES: frozenset = frozenset(
    {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}
)

# Ordered canonical names in UCSC convention -- used by the plotters for a
# stable genome-wide axis order (chr1..chr22, then X, Y, M).
CANONICAL_CHROMS_UCSC: tuple = tuple(
    [f"chr{i}" for i in range(1, 23)] + [f"chr{c}" for c in ("X", "Y", "M")]
)


def _core(name: str) -> str:
    """Strip an optional leading ``chr`` prefix (any case); return the rest."""
    return name[3:] if name[:3].lower() == "chr" else name


def is_canonical_chrom(name: str) -> bool:
    """Return True for a main-assembly chromosome under UCSC or Ensembl naming.

    Accepts ``chr1``/``1`` ... ``chr22``/``22``, ``chrX``/``X``, ``chrY``/``Y``
    and the mitochondrion ``chrM``/``chrMT``/``M``/``MT``. Rejects unplaced /
    unlocalized / alt contigs (``chr14_KI270722v1_random``, ``chrUn_*``,
    ``GL000216v2``, ``KI270722.1``) and out-of-range names (``chr23``, ``chr0``).
    """
    return _core(str(name)).upper() in CANONICAL_CHROM_CORES


def filter_canonical(chroms: Iterable[str]) -> list:
    """Order-preserving filter of *chroms* down to the canonical chromosomes."""
    return [c for c in chroms if is_canonical_chrom(c)]


def filter_canonical_logged(chroms: Iterable[str], *, context: str = "") -> list:
    """``filter_canonical`` plus one INFO line naming the dropped contigs.

    The shared audit message behind the DMC/DMR ``canonical_only`` default, so a
    dropped scaffold is never silent (and a non-mammalian assembly being
    over-dropped stays visible -- see the module docstring). ``context`` tags
    the log line with the call site (e.g. ``"dmc"``, ``"dmr/tile"``).
    """
    chroms = list(chroms)
    kept = filter_canonical(chroms)
    if len(kept) != len(chroms):
        dropped = [c for c in chroms if not is_canonical_chrom(c)]
        logger.info(
            "canonical_only%s: keeping %d canonical chromosome(s); dropping "
            "%d non-canonical contig(s): %s%s. Pass canonical_only=False "
            "(CLI: --all-contigs) to include them.",
            f" [{context}]" if context else "",
            len(kept),
            len(dropped),
            ", ".join(dropped[:5]),
            ", ..." if len(dropped) > 5 else "",
        )
    return kept
