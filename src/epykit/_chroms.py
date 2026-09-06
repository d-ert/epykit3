"""Canonical chromosome names shared by ingestion, analysis and plotting.

One definition of "main-assembly chromosome" so every ``canonical_only``
option in the pipeline drops the same contigs. The module is dependency-free
on purpose: low-level I/O imports it as well as the analysis and plotting
layers, so it must not pull in polars or numpy.

Naming
------
A leading ``chr`` prefix is optional and compared case-insensitively, so the
UCSC convention (``chr1``, ``chrM``) and the Ensembl convention (``1``,
``MT``) both match.

Scope
-----
The list is a fixed, human-style set: autosomes ``1`` to ``22``, the sex
chromosomes ``X`` and ``Y``, and the mitochondrion as ``M`` or ``MT``. It is
not a species-aware assembly validator. Mouse ``1`` to ``19`` happens to be a
subset; assemblies with roman-numeral or named chromosomes are not
recognised and must not use ``canonical_only=True``. Filtering is opt-in
everywhere, and :func:`filter_canonical_logged` names the dropped contigs so
an over-eager drop is visible in the log.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable

logger = logging.getLogger(__name__)

# Core identifiers without the ``chr`` prefix, compared case-insensitively.
CANONICAL_CHROM_CORES: frozenset[str] = frozenset(
    {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}
)

# The canonical names in UCSC convention and genome order (chr1..chr22, chrX,
# chrY, chrM). Plotters use it for a stable genome-wide axis.
CANONICAL_CHROMS_UCSC: tuple[str, ...] = tuple(
    [f"chr{i}" for i in range(1, 23)] + [f"chr{c}" for c in ("X", "Y", "M")]
)

# How many dropped contigs the audit line lists before summarising the rest.
_MAX_LISTED_DROPPED = 5


def _core(name: str) -> str:
    """Strip an optional leading ``chr`` prefix (any case) and return the rest."""
    return name[3:] if name[:3].lower() == "chr" else name


def is_canonical_chrom(name: str) -> bool:
    """Return True for a main-assembly chromosome under UCSC or Ensembl naming.

    Accepts ``chr1``/``1`` to ``chr22``/``22``, ``chrX``/``X``, ``chrY``/``Y``
    and the mitochondrion as ``chrM``/``chrMT``/``M``/``MT``. Rejects unplaced,
    unlocalised and alt contigs (``chr14_KI270722v1_random``, ``chrUn_*``,
    ``GL000216v2``, ``KI270722.1``) and out-of-range names (``chr23``,
    ``chr0``, ``chr01``).
    """
    return _core(str(name)).upper() in CANONICAL_CHROM_CORES


def filter_canonical(chroms: Iterable[str]) -> list[str]:
    """Return the canonical chromosomes of *chroms* in their original order."""
    return [c for c in chroms if is_canonical_chrom(c)]


def filter_canonical_logged(chroms: Iterable[str], *, context: str = "") -> list[str]:
    """Filter like :func:`filter_canonical` and log one INFO line per call.

    The line is the audit trail behind an opt-in ``canonical_only=True``: it
    names the dropped contigs and says how to keep them, so a dropped scaffold
    is never silent. Nothing is logged when every chromosome is canonical.
    ``context`` tags the line with the call site, for example ``"dmc"`` or
    ``"dmr/tile"``.
    """
    chroms = list(chroms)
    kept = filter_canonical(chroms)
    if len(kept) == len(chroms):
        return kept
    dropped = [c for c in chroms if not is_canonical_chrom(c)]
    listed = ", ".join(dropped[:_MAX_LISTED_DROPPED])
    rest = len(dropped) - _MAX_LISTED_DROPPED
    logger.info(
        "canonical_only%s: keeping %d canonical chromosome(s); dropping %d contig(s): %s%s. "
        "Omit canonical_only=True to retain them.",
        f" [{context}]" if context else "",
        len(kept),
        len(dropped),
        listed,
        f" and {rest} more" if rest > 0 else "",
    )
    return kept
