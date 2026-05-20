"""BAM ingestion for read-level methylation analyses.

This module is the shared input layer for analyses that need read-level
methylation information -- currently :mod:`epykit.asm` and
:mod:`epykit.entropy`. It does **not** replace the Bismark / MethylDackel
``.cov`` -> Parquet flow that drives ordinary DMC / DMR analyses;
``read_bismark`` and ``read_methyldackel`` remain unchanged.

Two BAM dialects are supported:

  * ``"bismark"`` -- methylation calls live in the per-base ``XM`` tag.
    The character codes are: ``Z`` methylated CpG, ``z`` unmethylated
    CpG, ``X``/``x`` methylated/unmethylated CHG, ``H``/``h`` CHH,
    ``.`` no call. See Bismark's docs.
  * ``"methyldackel"`` -- the SAM standard ``MM`` (modified-base
    positions) and ``ML`` (likelihood) tags. ``MM:Z:C+m,...`` flags
    methylated C positions on the forward strand.

Both produce the same long-form output (one row per read x covered
CpG): ``(read_id, chrom, pos, methylation_status, base_qual, mapq,
mate_pair_id, strand, allele_base)``. Downstream code is dialect-agnostic.

pysam is an optional dependency. Install with
``pip install 'epykit[bam]'`` (Linux/macOS only -- pysam has no Windows
wheel).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Iterator, Optional

import polars as pl

logger = logging.getLogger(__name__)


# Output schema: every BAM reader emits rows in this shape so downstream
# code is dialect-agnostic.
_BAM_READ_SCHEMA = {
    "read_id": pl.Utf8,
    "chrom": pl.Utf8,
    "pos": pl.Int32,
    "methylation_status": pl.Int8,   # 1 = methylated, 0 = unmethylated, -1 = no call
    "base_qual": pl.Int8,
    "mapq": pl.Int8,
    "mate_pair_id": pl.Int8,         # 1 / 2 for paired-end, 0 for single-end
    "strand": pl.Utf8,               # "+" / "-"
    "allele_base": pl.Utf8,          # base call at pos, used by ASM for haplotyping
}


def _require_pysam():
    """Lazy import of pysam with a clear install message."""
    try:
        import pysam  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "pysam is required for BAM ingestion. "
            "Install with: pip install 'epykit[bam]' "
            "(Linux/macOS only; pysam has no Windows wheel -- "
            "use 'epykit[methylkit]' for the tabix-only feature set on Windows)."
        ) from exc
    return pysam


# Bismark XM tag dispatch
_BISMARK_METHYLATED = set("ZXH")     # methylated CpG / CHG / CHH
_BISMARK_UNMETHYLATED = set("zxh")   # unmethylated CpG / CHG / CHH


def read_methylation_calls(
    bam_path: str | Path,
    *,
    regions: Optional[Iterable[tuple[str, int, int]]] = None,
    min_mapq: int = 10,
    min_baseq: int = 20,
    caller: str = "bismark",
    context: str = "CpG",
) -> pl.DataFrame:
    """Extract per-base methylation calls from a BAM file.

    Parameters
    ----------
    bam_path
        Path to a coordinate-sorted, indexed BAM. ``.bai`` must be
        alongside (or a ``.csi``).
    regions
        Iterable of ``(chrom, start, end)`` tuples to restrict the
        scan. ``None`` (default) reads the entire BAM.
    min_mapq, min_baseq
        Quality filters. Reads with ``mapq < min_mapq`` and individual
        positions with ``base_qual < min_baseq`` are dropped.
    caller
        BAM dialect: ``"bismark"`` (XM tag, default) or
        ``"methyldackel"`` (MM/ML tags).
    context
        Methylation context to retain: ``"CpG"`` (default), ``"CHG"``,
        ``"CHH"``, or ``"any"``. Bismark XM letters encode context.

    Returns
    -------
    pl.DataFrame
        Long-form frame, one row per (read, covered position). See
        ``_BAM_READ_SCHEMA`` for the exact columns.
    """
    pysam = _require_pysam()

    caller = caller.lower()
    if caller not in ("bismark", "methyldackel"):
        raise ValueError(f"Unknown caller {caller!r}. Use 'bismark' or 'methyldackel'.")

    bam_p = Path(bam_path)
    if not bam_p.exists():
        raise FileNotFoundError(f"BAM not found: {bam_p}")

    rows: list[dict[str, object]] = []
    with pysam.AlignmentFile(str(bam_p), "rb") as bam:
        if regions is None:
            iterator = bam.fetch(until_eof=True)
        else:
            iterator = _multi_region_iter(bam, regions)

        for read in iterator:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_duplicate or read.is_qcfail:
                continue
            if read.mapping_quality < min_mapq:
                continue

            if caller == "bismark":
                rows.extend(_extract_bismark(read, min_baseq=min_baseq, context=context))
            else:
                rows.extend(_extract_methyldackel(read, min_baseq=min_baseq))

    if not rows:
        return pl.DataFrame(schema=_BAM_READ_SCHEMA)

    return pl.DataFrame(rows, schema=_BAM_READ_SCHEMA)


def _multi_region_iter(bam, regions: Iterable[tuple[str, int, int]]) -> Iterator:
    """Yield reads from each (chrom, start, end) in turn."""
    for chrom, start, end in regions:
        try:
            yield from bam.fetch(chrom, start, end)
        except ValueError:
            logger.warning("Region %s:%d-%d skipped (chrom not in BAM)", chrom, start, end)


def _read_meta(read) -> tuple[str, str, int, int]:
    """Return (read_id, strand, mate_pair_id, mapq) shared across dialects."""
    strand = "-" if read.is_reverse else "+"
    if read.is_paired:
        mate = 1 if read.is_read1 else 2
    else:
        mate = 0
    return read.query_name, strand, mate, read.mapping_quality


def _extract_bismark(read, *, min_baseq: int, context: str) -> list[dict[str, object]]:
    """Pull per-base methylation calls from a Bismark XM tag."""
    try:
        xm = read.get_tag("XM")
    except KeyError:
        return []
    if not xm:
        return []

    read_id, strand, mate, mapq = _read_meta(read)
    chrom = read.reference_name
    if chrom is None:
        return []

    aligned_pairs = read.get_aligned_pairs(matches_only=True)
    if not aligned_pairs:
        return []
    quals = read.query_qualities  # may be None for older BAMs
    seq = read.query_sequence

    keep_letters: set[str]
    if context == "CpG":
        keep_letters = {"Z", "z"}
    elif context == "CHG":
        keep_letters = {"X", "x"}
    elif context == "CHH":
        keep_letters = {"H", "h"}
    else:
        keep_letters = _BISMARK_METHYLATED | _BISMARK_UNMETHYLATED

    out: list[dict[str, object]] = []
    for query_idx, ref_pos in aligned_pairs:
        if query_idx is None or query_idx >= len(xm):
            continue
        letter = xm[query_idx]
        if letter not in keep_letters:
            continue
        baseq = int(quals[query_idx]) if quals is not None else 0
        if baseq < min_baseq:
            continue
        meth_status = 1 if letter in _BISMARK_METHYLATED else 0
        allele = seq[query_idx] if seq else ""
        out.append({
            "read_id": read_id,
            "chrom": chrom,
            "pos": int(ref_pos),
            "methylation_status": meth_status,
            "base_qual": min(baseq, 127),
            "mapq": min(mapq, 127),
            "mate_pair_id": mate,
            "strand": strand,
            "allele_base": allele,
        })
    return out


def _extract_methyldackel(read, *, min_baseq: int) -> list[dict[str, object]]:
    """Pull per-base methylation calls from SAM MM/ML tags.

    The MM tag encodes positions of modified Cs (counted in C-residues
    along the read); ML encodes per-modification likelihood (0-255). A
    call is "methylated" when ML >= 128 (>= 0.5 probability), otherwise
    "unmethylated".
    """
    try:
        mm = read.get_tag("MM")
        ml = read.get_tag("ML")
    except KeyError:
        return []
    if not mm:
        return []

    # MM format: "C+m,N1,N2,...;" -- N_i = number of C residues to skip.
    # We support the simple "C+m" prefix (5mC on the forward strand). For
    # robustness, scan the first specifier and bail out on anything else.
    spec, _, body = mm.partition(",")
    if not spec.startswith("C+m"):
        return []
    if not body:
        return []
    try:
        skips = [int(s) for s in body.rstrip(";").split(",") if s]
    except ValueError:
        return []

    read_id, strand, mate, mapq = _read_meta(read)
    chrom = read.reference_name
    if chrom is None:
        return []
    seq = read.query_sequence
    if not seq:
        return []
    quals = read.query_qualities

    # Build a list of C-residue indices along the read.
    c_indices = [i for i, base in enumerate(seq) if base in ("C", "c")]

    # Each skip value points to the (skip+1)-th remaining C residue.
    target_query_idx: list[int] = []
    cursor = 0
    for skip in skips:
        cursor += skip
        if cursor >= len(c_indices):
            break
        target_query_idx.append(c_indices[cursor])
        cursor += 1

    if not target_query_idx:
        return []

    if not isinstance(ml, (list, tuple, bytes, bytearray)):
        return []
    if len(ml) < len(target_query_idx):
        return []

    # Map query indices -> ref positions via the aligned-pairs table.
    aligned = dict(read.get_aligned_pairs(matches_only=True))

    out: list[dict[str, object]] = []
    for q_idx, prob in zip(target_query_idx, ml[: len(target_query_idx)]):
        ref_pos = aligned.get(q_idx)
        if ref_pos is None:
            continue
        baseq = int(quals[q_idx]) if quals is not None else 0
        if baseq < min_baseq:
            continue
        meth_status = 1 if int(prob) >= 128 else 0
        out.append({
            "read_id": read_id,
            "chrom": chrom,
            "pos": int(ref_pos),
            "methylation_status": meth_status,
            "base_qual": min(baseq, 127),
            "mapq": min(mapq, 127),
            "mate_pair_id": mate,
            "strand": strand,
            "allele_base": seq[q_idx] if seq else "",
        })
    return out


__all__ = ["read_methylation_calls"]
