"""ASM (allele-specific methylation) caller tests.

Every test builds a tiny synthetic BAM + VCF at run time: 30 bp reads on
``chr_asm`` that cover one phasing SNV at 0-based position 50 and one
measured CpG at position 60. Each read is described by the base it shows
at the SNV, its Bismark XM letter at the CpG (``Z`` methylated, ``z``
unmethylated) and its Bismark ``XG`` genome-conversion tag (``None``
leaves the read untagged).

The tests cover the planted A/G signal, the bisulfite confound at a C/T
anchor (an unmethylated C converts to T on CT-strand reads, so raw base
matching would fabricate ASM), and the per-class, per-strand anchor rule.

Skips cleanly when pysam isn't available (Windows).
"""

from __future__ import annotations

import logging
from pathlib import Path

import polars as pl
import pytest

pysam = pytest.importorskip("pysam", reason="pysam required for ASM tests.")

from epykit.asm import call_asm  # noqa: E402 -- must follow importorskip

ANCHOR_POS = 50  # 0-based reference position of the phasing SNV (read offset 5)
CPG_POS = 60  # 0-based reference position of the measured CpG (read offset 15)
READ_START = 45
READ_LEN = 30

# One synthetic read: (base at the anchor, XM letter at the CpG, XG tag).
Read = tuple[str, str, str | None]


def _write_bam(bam_path: Path, reads: list[Read]) -> Path:
    """Write a coordinate-sorted, indexed BAM with one 30 bp read per spec."""
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr_asm", "LN": 200}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for i, (anchor_base, cpg_call, xg) in enumerate(reads):
            r = pysam.AlignedSegment(out.header)
            r.query_name = f"read_{i}"
            seq = ["T"] * READ_LEN  # filler
            seq[ANCHOR_POS - READ_START] = anchor_base
            seq[CPG_POS - READ_START] = "C"
            r.query_sequence = "".join(seq)
            r.flag = 0
            r.reference_id = 0
            r.reference_start = READ_START
            r.mapping_quality = 60
            r.cigar = [(0, READ_LEN)]
            r.query_qualities = pysam.qualitystring_to_array("E" * READ_LEN)
            xm = ["."] * READ_LEN
            xm[CPG_POS - READ_START] = cpg_call
            r.set_tag("XM", "".join(xm))
            if xg is not None:
                r.set_tag("XG", xg)
            out.write(r)

    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))
    return bam_path


def _write_vcf(raw_path: Path, records: list[tuple[int, str, str]]) -> Path:
    """Write a bgzipped, tabix-indexed VCF; ``records`` are ``(pos_1based, ref, alt)``, all het."""
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr_asm,length=200>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
    ]
    for pos, ref, alt in records:
        lines.append(f"chr_asm\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT\t0/1")
    raw_path.write_text("\n".join(lines) + "\n")
    vcf_path = raw_path.with_suffix(".vcf.gz")
    pysam.tabix_compress(str(raw_path), str(vcf_path), force=True)
    pysam.tabix_index(str(vcf_path), preset="vcf", force=True)
    return vcf_path


def _run(tmp_path: Path, reads: list[Read], records: list[tuple[int, str, str]]) -> pl.DataFrame:
    bam_path = _write_bam(tmp_path / "asm.bam", reads)
    vcf_path = _write_vcf(tmp_path / "asm.vcf", records)
    return call_asm(
        bam={"sample1": bam_path},
        vcf=vcf_path,
        min_reads_per_haplotype=5,
        min_phased_snvs=1,
    )


def _counts(df: pl.DataFrame) -> tuple[int, int, int, int]:
    """(h1_meth, h1_unmeth, h2_meth, h2_unmeth) at the planted CpG."""
    planted = df.filter(pl.col("pos") == CPG_POS)
    assert planted.height == 1, f"planted CpG missing; got {df}"
    row = planted.row(0, named=True)
    return row["h1_meth"], row["h1_unmeth"], row["h2_meth"], row["h2_unmeth"]


def _summary(caplog) -> str:
    lines = [r.getMessage() for r in caplog.records if "phasing summary" in r.getMessage()]
    assert len(lines) == 1, f"expected one phasing summary; got {lines}"
    return lines[0]


CT_ANCHOR = [(ANCHOR_POS + 1, "C", "T")]


def test_asm_recovers_planted_signal(tmp_path):
    """The A/G anchor stays usable on CT-strand reads (neither allele converts)."""
    reads: list[Read] = [("A", "Z", "CT")] * 10 + [("G", "z", "CT")] * 10
    df = _run(tmp_path, reads, [(ANCHOR_POS + 1, "A", "G")])

    assert _counts(df) == (10, 0, 0, 10)
    row = df.filter(pl.col("pos") == CPG_POS).row(0, named=True)
    # Meth diff should be ~1.0 (fully methylated h1 vs unmethylated h2).
    assert abs(row["meth_diff"] - 1.0) < 1e-6
    # Fisher's exact on a (10,0) vs (0,10) table is highly significant.
    assert row["pvalue"] < 0.01
    assert "qvalue" in df.columns


def test_ct_anchor_on_ct_reads_is_not_phased(tmp_path, caplog):
    """A null C/T anchor seen only on CT-strand reads contributes no ASM result.

    Both alleles are 50% methylated. On CT reads an unmethylated C converts
    to T, so raw base matching would file the ten converted C-allele reads
    under the T allele and report 10/0 versus 10/20 (Fisher p ~ 4e-4).
    """
    caplog.set_level(logging.INFO, logger="epykit.asm")
    reads: list[Read] = (
        [("C", "Z", "CT")] * 10  # C allele, methylated: the anchor C is protected
        + [("T", "z", "CT")] * 10  # C allele, unmethylated: the anchor C reads as T
        + [("T", "Z", "CT")] * 10  # T allele, methylated
        + [("T", "z", "CT")] * 10  # T allele, unmethylated
    )
    df = _run(tmp_path, reads, CT_ANCHOR)

    assert df.height == 0, f"unsafe CT-only anchor was phased: {df}"
    assert _summary(caplog).endswith(
        "anchors_phased=0 anchors_rejected_class=0 reads_rejected_xg=40"
    )


def test_ct_anchor_on_ga_reads_balanced_is_null(tmp_path):
    """On GA-strand reads C is literal, so the balanced anchor gives 10/10 versus 10/10."""
    reads: list[Read] = (
        [("C", "Z", "GA")] * 10
        + [("C", "z", "GA")] * 10
        + [("T", "Z", "GA")] * 10
        + [("T", "z", "GA")] * 10
    )
    df = _run(tmp_path, reads, CT_ANCHOR)

    assert _counts(df) == (10, 10, 10, 10)
    row = df.filter(pl.col("pos") == CPG_POS).row(0, named=True)
    assert row["meth_diff"] == 0.0
    assert row["pvalue"] == pytest.approx(1.0)


def test_ct_anchor_on_ga_reads_detects_signal(tmp_path):
    """A genuine C/T ASM signal is still called from GA-strand reads."""
    reads: list[Read] = [("C", "Z", "GA")] * 10 + [("T", "z", "GA")] * 10
    df = _run(tmp_path, reads, CT_ANCHOR)

    assert _counts(df) == (10, 0, 0, 10)
    assert df.filter(pl.col("pos") == CPG_POS)["pvalue"][0] < 0.01


def test_unsafe_reads_never_contribute(tmp_path, caplog):
    """At a C/T anchor only GA-strand reads are assigned; CT and untagged reads are dropped."""
    caplog.set_level(logging.INFO, logger="epykit.asm")
    reads: list[Read] = (
        [("C", "Z", "GA")] * 10
        + [("T", "z", "GA")] * 10
        + [("T", "Z", "CT")] * 10  # would inflate h2_meth if phased
        + [("C", "z", None)] * 5  # would inflate h1_unmeth if phased
    )
    df = _run(tmp_path, reads, CT_ANCHOR)

    assert _counts(df) == (10, 0, 0, 10)
    assert _summary(caplog).endswith(
        "anchors_phased=1 anchors_rejected_class=0 reads_rejected_xg=15"
    )


def test_cg_and_invalid_anchors_rejected_before_fetch(tmp_path, caplog):
    """C/G and non-ACGT anchors count as class rejections; an untagged A/T anchor still phases."""
    caplog.set_level(logging.INFO, logger="epykit.asm")
    reads: list[Read] = [("A", "Z", None)] * 10 + [("T", "z", None)] * 10
    records = [
        (ANCHOR_POS + 1, "A", "T"),
        (ANCHOR_POS + 6, "C", "G"),
        (ANCHOR_POS + 8, "N", "A"),
    ]
    df = _run(tmp_path, reads, records)

    assert _counts(df) == (10, 0, 0, 10)
    assert _summary(caplog).endswith(
        "anchors_phased=1 anchors_rejected_class=2 reads_rejected_xg=0"
    )


def test_lowercase_vcf_alleles_are_normalised(tmp_path):
    reads: list[Read] = [("A", "Z", None)] * 10 + [("T", "z", None)] * 10
    df = _run(tmp_path, reads, [(ANCHOR_POS + 1, "a", "t")])

    assert _counts(df) == (10, 0, 0, 10)


# The acceptance table: unordered SNV class -> XG strands that may phase it.
# Reads with a missing or unrecognised XG tag can only phase A/T anchors.
_SAFE_ON = {
    "AT": {"CT", "GA"},
    "AG": {"CT"},
    "GT": {"CT"},
    "CT": {"GA"},
    "AC": {"GA"},
    "CG": set(),
}


@pytest.mark.parametrize(
    "xg", ["CT", "GA", None, "XX"], ids=["xg_ct", "xg_ga", "untagged", "xg_bad"]
)
@pytest.mark.parametrize("swap", [False, True], ids=["ref_alt", "alt_ref"])
@pytest.mark.parametrize("snv_class", sorted(_SAFE_ON))
def test_anchor_class_rule(tmp_path, snv_class, swap, xg):
    ref, alt = snv_class
    if swap:
        ref, alt = alt, ref
    reads: list[Read] = [(ref, "Z", xg)] * 10 + [(alt, "z", xg)] * 10
    df = _run(tmp_path, reads, [(ANCHOR_POS + 1, ref, alt)])

    accepted = xg in _SAFE_ON[snv_class] if xg in ("CT", "GA") else snv_class == "AT"
    if accepted:
        assert _counts(df) == (10, 0, 0, 10)
    else:
        assert df.height == 0, f"{ref}/{alt} on XG={xg} must not phase; got {df}"
