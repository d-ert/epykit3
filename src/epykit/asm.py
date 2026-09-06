"""Allele-specific methylation (ASM) caller.

For each heterozygous SNV in a per-individual VCF, partition the
overlapping reads by their base at the SNV position into haplotype 1 /
haplotype 2, then project each read's methylation calls onto every CpG
it covers. Per CpG, build a 2 x 2 (h1_meth, h1_unmeth, h2_meth,
h2_unmeth) contingency table and test allele-specific methylation
with Fisher's exact test.

Schema compatibility
--------------------
Results land in ``md.varm["asm"]`` with column names matching the
``dmc_*`` family (``pvalue``, ``qvalue``, ``meth_diff``, ``chrom``,
``pos``, ``strand``). This lets ``pl.volcano(md, key="asm")``, the
Manhattan plotter, and other DMC-aware visualisations work without
modification.

Inputs
------
* ``bam``: ``{sample_id -> bam_path}``. BAMs must be coordinate-sorted
  and indexed; reads need either Bismark ``XM`` tags or
  MethylDackel ``MM/ML`` tags.
* ``vcf``: path to a per-individual VCF (bgzipped + tabix-indexed
  preferred). Heterozygous biallelic SNVs are used as phasing anchors.

Bisulfite-safe anchors
----------------------
Bisulfite conversion reads an unmethylated C as T on reads aligned to the
converted top strand (Bismark ``XG:Z:CT``) and an unmethylated G as A on
reads from the converted bottom strand (``XG:Z:GA``). At a heterozygous
SNV whose alleles include that base, the base a read shows at the anchor
depends on its methylation state as well as on its allele, so phasing by
the raw base would fabricate ASM. A read is therefore assigned to an
allele only when its ``XG`` strand cannot convert either allele of the
unordered REF/ALT class:

======  ======  ======  ======================
class   XG=CT   XG=GA   missing or unknown XG
======  ======  ======  ======================
A/T     yes     yes     yes
A/G     yes     no      no
G/T     yes     no      no
C/T     no      yes     no
A/C     no      yes     no
C/G     no      no      no
======  ======  ======  ======================

C/G anchors (and records whose REF or ALT is not A, C, G or T) are
skipped before any read is fetched. Reads without a recognised ``XG``
tag, such as MethylDackel / bwa-meth BAMs or Bismark BAMs with the tag
stripped, can only phase A/T anchors. The rule is deliberately
conservative: it also drops reads that a methylation-aware genotyper
could have rescued, so fewer anchors and reads phase than in releases
before 1.2 and some previously reported sites disappear. The caller has
not been validated on real biological data.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import polars as pl

from .bam_io import _require_pysam, read_methylation_calls
from .dmc import apply_multiple_testing_correction, fisher_exact_vectorized

logger = logging.getLogger(__name__)


_ASM_SCHEMA = {
    "sample_id": pl.Utf8,
    "chrom": pl.Utf8,
    "pos": pl.Int32,
    "strand": pl.Utf8,
    "h1_meth": pl.Int32,
    "h1_unmeth": pl.Int32,
    "h2_meth": pl.Int32,
    "h2_unmeth": pl.Int32,
    "n_h1": pl.Int32,
    "n_h2": pl.Int32,
    "meth_diff": pl.Float32,
    "pvalue": pl.Float64,
}

_VALID_BASES = frozenset("ACGT")
_AT = frozenset("AT")
_CG = frozenset("CG")
# The base that bisulfite conversion can change on reads of each Bismark
# XG genome-conversion strand: an unmethylated C reads as T on CT-strand
# reads, an unmethylated G reads as A on GA-strand reads.
_CONVERTIBLE_BASE = {"CT": "C", "GA": "G"}


def _snv_class(ref: str, alt: str) -> frozenset[str] | None:
    """Unordered class of an upper-cased biallelic SNV; None unless two distinct ACGT bases."""
    pair = frozenset((ref, alt))
    if len(pair) != 2 or not pair <= _VALID_BASES:
        return None
    return pair


def _xg_tag(read) -> str:
    """Bismark genome-conversion strand of a read (``"CT"`` or ``"GA"``); ``""`` when untagged."""
    try:
        return str(read.get_tag("XG"))
    except KeyError:
        return ""


def _strand_can_phase(snv_class: frozenset[str], xg: str) -> bool:
    """True when a read on conversion strand ``xg`` shows both alleles literally.

    A class that contains the strand's convertible base is rejected, because
    the read's base at the anchor then depends on methylation state. Reads
    with a missing or unrecognised tag can only phase A/T anchors.
    """
    convertible = _CONVERTIBLE_BASE.get(xg)
    if convertible is None:
        return snv_class == _AT
    return convertible not in snv_class


def call_asm(
    bam: Mapping[str, str | Path],
    vcf: str | Path,
    *,
    min_reads_per_haplotype: int = 10,
    min_baseq: int = 20,
    min_mapq: int = 10,
    min_phased_snvs: int = 1,
    caller: str = "bismark",
    chromosomes: list[str] | None = None,
) -> pl.DataFrame:
    """Run the ASM caller across one or more samples.

    Returns a long-form polars DataFrame with one row per
    (sample, CpG with enough h1/h2 reads). The DataFrame is BH-corrected
    on ``pvalue`` to produce ``qvalue``.

    Reads are assigned to an allele only at bisulfite-safe anchors for
    their Bismark ``XG`` conversion strand (see the module docstring for
    the table). C/G anchors are never used; reads without a recognised
    ``XG`` tag phase A/T anchors only. One INFO line per sample reports
    how many anchors phased at least one read, how many anchors were
    rejected by class before any read was fetched, and how many
    read-anchor observations the strand rule rejected. This
    strand-aware exclusion is a conservative substitute for
    methylation-aware genotyping, and the caller has not been validated
    on real biological data.
    """
    pysam = _require_pysam()
    vcf_p = Path(vcf)
    if not vcf_p.exists():
        raise FileNotFoundError(f"VCF not found: {vcf_p}")

    sample_frames: list[pl.DataFrame] = []
    for sample_id, bam_path in bam.items():
        logger.info("[ASM] %s: extracting methylation calls", sample_id)
        meth_df = read_methylation_calls(
            bam_path,
            caller=caller,
            min_baseq=min_baseq,
            min_mapq=min_mapq,
        )
        if meth_df.height == 0:
            logger.warning("[ASM] %s: no usable reads; skipping", sample_id)
            continue
        if chromosomes is not None:
            meth_df = meth_df.filter(pl.col("chrom").is_in(chromosomes))
            if meth_df.height == 0:
                continue

        logger.info(
            "[ASM] %s: %d reads x CpG rows; phasing via VCF",
            sample_id,
            meth_df.height,
        )
        df = _call_asm_one_sample(
            pysam,
            meth_df,
            str(vcf_p),
            bam_path,
            sample_id=sample_id,
            min_reads_per_haplotype=min_reads_per_haplotype,
            min_phased_snvs=min_phased_snvs,
            min_mapq=min_mapq,
        )
        if df.height > 0:
            df = df.with_columns(pl.lit(sample_id).alias("sample_id"))
            sample_frames.append(df)

    if not sample_frames:
        return pl.DataFrame(schema=_ASM_SCHEMA).with_columns(
            pl.lit(None).cast(pl.Float64).alias("qvalue"),
        )

    combined = pl.concat(sample_frames, how="vertical_relaxed")
    combined = apply_multiple_testing_correction(combined, method="fdr_bh")
    return combined.sort(["sample_id", "chrom", "pos"])


def _call_asm_one_sample(
    pysam,
    meth_df: pl.DataFrame,
    vcf_path: str,
    bam_path: str | Path,
    *,
    sample_id: str,
    min_reads_per_haplotype: int,
    min_phased_snvs: int,
    min_mapq: int,
) -> pl.DataFrame:
    """Per-sample ASM: phase reads via het SNVs, build 2x2 tables per CpG."""
    # ---- 1. Map each read_id -> haplotype via het SNVs in the VCF ----
    read_haplotype: dict[str, int] = {}  # read_id -> 1 or 2
    read_phased_snv_count: dict[str, int] = {}  # read_id -> # of confirming SNVs
    n_anchors_phased = 0  # het SNVs that assigned at least one read
    n_anchors_rejected_class = 0  # het SNVs skipped before fetch (C/G or non-ACGT)
    n_reads_rejected_xg = 0  # read-anchor observations the strand rule dropped

    with pysam.VariantFile(vcf_path) as vcf, pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for rec in vcf:
            if not rec.alts or len(rec.alts) != 1:
                continue  # skip multi-allelic / non-bi-allelic for simplicity
            ref = rec.ref.upper()
            alt = rec.alts[0].upper()
            if len(ref) != 1 or len(alt) != 1:
                continue  # SNVs only (no indels)

            # Is this SNV heterozygous in the sample?
            if not _is_het(rec):
                continue

            # C/G anchors convert on both strands; non-ACGT records have no class.
            snv_class = _snv_class(ref, alt)
            if snv_class is None or snv_class == _CG:
                n_anchors_rejected_class += 1
                continue

            chrom = rec.chrom
            pos = rec.pos - 1  # VCF is 1-based, BAM is 0-based
            try:
                pile = bam.fetch(chrom, pos, pos + 1)
            except ValueError:
                continue
            anchor_phased = False
            for read in pile:
                if read.is_unmapped or read.mapping_quality < min_mapq:
                    continue
                if read.is_secondary or read.is_supplementary:
                    continue
                # Only phase from a strand on which neither allele converts.
                if not _strand_can_phase(snv_class, _xg_tag(read)):
                    n_reads_rejected_xg += 1
                    continue
                # Find the base in this read at the SNV position.
                aligned = dict(read.get_aligned_pairs(matches_only=True))
                # Reverse lookup: ref_pos -> query_idx.
                qry = None
                for q, r in aligned.items():
                    if r == pos:
                        qry = q
                        break
                if qry is None:
                    continue
                seq = read.query_sequence
                if not seq or qry >= len(seq):
                    continue
                base = seq[qry].upper()
                if base == ref:
                    hap = 1
                elif base == alt:
                    hap = 2
                else:
                    continue
                anchor_phased = True
                prior = read_haplotype.get(read.query_name)
                if prior is None:
                    read_haplotype[read.query_name] = hap
                    read_phased_snv_count[read.query_name] = 1
                elif prior == hap:
                    read_phased_snv_count[read.query_name] += 1
                else:
                    # SNV-disagreeing read: drop from phasing entirely.
                    read_haplotype.pop(read.query_name, None)
                    read_phased_snv_count.pop(read.query_name, None)
            if anchor_phased:
                n_anchors_phased += 1

    logger.info(
        "[ASM] %s: phasing summary: anchors_phased=%d anchors_rejected_class=%d reads_rejected_xg=%d",
        sample_id,
        n_anchors_phased,
        n_anchors_rejected_class,
        n_reads_rejected_xg,
    )
    if not read_haplotype:
        logger.warning(
            "[ASM] %s: no phaseable reads (no bisulfite-safe het SNV covered by any read)",
            sample_id,
        )
        return pl.DataFrame(schema={k: v for k, v in _ASM_SCHEMA.items() if k != "sample_id"})

    # Filter to reads that hit >= min_phased_snvs.
    kept_reads = {rid for rid, n in read_phased_snv_count.items() if n >= min_phased_snvs}
    if not kept_reads:
        return pl.DataFrame(schema={k: v for k, v in _ASM_SCHEMA.items() if k != "sample_id"})

    # ---- 2. Project methylation calls through the haplotype map ----
    # Annotate each (read, pos) row with its haplotype; drop rows whose
    # read isn't phased.
    hap_df = pl.DataFrame(
        {
            "read_id": list(kept_reads),
            "haplotype": [read_haplotype[r] for r in kept_reads],
        },
        schema={"read_id": pl.Utf8, "haplotype": pl.Int32},
    )

    joined = meth_df.join(hap_df, on="read_id", how="inner")
    if joined.height == 0:
        return pl.DataFrame(schema={k: v for k, v in _ASM_SCHEMA.items() if k != "sample_id"})

    # ---- 3. Per CpG: build (h1_meth, h1_unmeth, h2_meth, h2_unmeth) ----
    per_cpg = (
        joined.group_by(["chrom", "pos", "strand"])
        .agg(
            [
                ((pl.col("haplotype") == 1) & (pl.col("methylation_status") == 1))
                .sum()
                .cast(pl.Int32)
                .alias("h1_meth"),
                ((pl.col("haplotype") == 1) & (pl.col("methylation_status") == 0))
                .sum()
                .cast(pl.Int32)
                .alias("h1_unmeth"),
                ((pl.col("haplotype") == 2) & (pl.col("methylation_status") == 1))
                .sum()
                .cast(pl.Int32)
                .alias("h2_meth"),
                ((pl.col("haplotype") == 2) & (pl.col("methylation_status") == 0))
                .sum()
                .cast(pl.Int32)
                .alias("h2_unmeth"),
            ]
        )
        .with_columns(
            [
                (pl.col("h1_meth") + pl.col("h1_unmeth")).alias("n_h1"),
                (pl.col("h2_meth") + pl.col("h2_unmeth")).alias("n_h2"),
            ]
        )
        .filter(
            (pl.col("n_h1") >= min_reads_per_haplotype)
            & (pl.col("n_h2") >= min_reads_per_haplotype)
        )
        .sort(["chrom", "pos"])
    )
    if per_cpg.height == 0:
        logger.warning(
            "[ASM] no CpGs reached %d reads per haplotype.",
            min_reads_per_haplotype,
        )
        return pl.DataFrame(schema={k: v for k, v in _ASM_SCHEMA.items() if k != "sample_id"})

    # ---- 4. Fisher exact + meth_diff per CpG ----
    h1m = per_cpg["h1_meth"].to_numpy()
    h1u = per_cpg["h1_unmeth"].to_numpy()
    h2m = per_cpg["h2_meth"].to_numpy()
    h2u = per_cpg["h2_unmeth"].to_numpy()
    pvals, _log2_or = fisher_exact_vectorized(h1m, h1u, h2m, h2u)

    n_h1 = h1m + h1u
    n_h2 = h2m + h2u
    with np.errstate(invalid="ignore", divide="ignore"):
        beta_h1 = h1m / np.where(n_h1 > 0, n_h1, np.nan)
        beta_h2 = h2m / np.where(n_h2 > 0, n_h2, np.nan)
    meth_diff = (beta_h1 - beta_h2).astype(np.float32)

    return per_cpg.with_columns(
        [
            pl.Series("meth_diff", meth_diff),
            pl.Series("pvalue", pvals),
        ]
    )


def _is_het(record) -> bool:
    """True if the (first / only) sample in the VCF record is heterozygous."""
    if not record.samples:
        return False
    sample = next(iter(record.samples.values()))
    gt = sample.get("GT")
    if gt is None or len(gt) < 2:
        return False
    # Heterozygous = two different non-None allele indices.
    a, b = gt[0], gt[1]
    return a is not None and b is not None and a != b


def asm(
    md,
    *,
    bam: Mapping[str, str | Path],
    vcf: str | Path,
    min_reads_per_haplotype: int = 10,
    min_phased_snvs: int = 1,
    chromosomes: list[str] | None = None,
    caller: str = "bismark",
) -> None:
    """Run ASM and store results in ``md.varm["asm"]``.

    See :func:`call_asm` for the full algorithm. This is the
    ``MethylData``-aware wrapper that mirrors the ``tl.dmc`` / ``tl.dmr``
    convention.
    """
    md_samples = set(md.obs.get_column("sample_id").to_list())
    missing = [s for s in bam if s not in md_samples]
    if missing:
        raise ValueError(
            f"bam keys not in md.obs.sample_id: {missing[:5]}{' ...' if len(missing) > 5 else ''}"
        )
    result = call_asm(
        bam=bam,
        vcf=vcf,
        min_reads_per_haplotype=min_reads_per_haplotype,
        min_phased_snvs=min_phased_snvs,
        chromosomes=chromosomes,
        caller=caller,
    )
    md.varm["asm"] = result
    md.uns["asm"] = {
        "n_sites": int(result.height),
        "n_samples": len(bam),
        "min_reads_per_haplotype": min_reads_per_haplotype,
        "min_phased_snvs": min_phased_snvs,
        "vcf": str(vcf),
    }


__all__ = ["asm", "call_asm"]
