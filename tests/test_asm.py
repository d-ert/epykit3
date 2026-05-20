"""ASM (allele-specific methylation) caller tests.

We build a tiny synthetic BAM + VCF where one CpG is fully methylated
on h1 and fully unmethylated on h2 (a planted ASM signal). The test
then verifies the caller flags that CpG as significant.

Skips cleanly when pysam isn't available (Windows).
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

pysam = pytest.importorskip("pysam", reason="pysam required for ASM tests.")

from epykit.asm import call_asm


def _write_synth_bam_and_vcf(tmp_path: Path) -> tuple[Path, Path]:
    """Build a tiny BAM + VCF with one ASM-positive CpG planted."""
    bam_path = tmp_path / "asm.bam"
    vcf_path = tmp_path / "asm.vcf.gz"

    # Reference: 100 bp on chr_asm. Het SNV at pos 50 (REF=A, ALT=G).
    # 10 reads carry A (h1, all methylated at the CpG at pos 60).
    # 10 reads carry G (h2, all UNMETHYLATED at the CpG at pos 60).
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr_asm", "LN": 200}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        # Reads start at pos 45, cover SNV at pos 50 and CpG at pos 60.
        # query bases for the SNV are at offset 5; CpG at offset 15.
        for i in range(20):
            r = pysam.AlignedSegment(out.header)
            r.query_name = f"read_{i}"
            # Build a 30bp sequence with the SNV base at offset 5 and
            # a 'C' at offset 15 (the CpG position).
            snv_base = "A" if i < 10 else "G"
            seq_list = ["N"] * 30
            for j in range(30):
                seq_list[j] = "T"  # filler
            seq_list[5] = snv_base
            seq_list[15] = "C"
            r.query_sequence = "".join(seq_list)
            r.flag = 0
            r.reference_id = 0
            r.reference_start = 45
            r.mapping_quality = 60
            r.cigar = [(0, 30)]
            r.query_qualities = pysam.qualitystring_to_array("E" * 30)
            # XM: 30 chars. Position 15 carries the CpG call.
            xm = ["."] * 30
            xm[15] = "Z" if i < 10 else "z"   # h1 -> methylated, h2 -> unmethylated
            r.set_tag("XM", "".join(xm))
            out.write(r)

    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))

    # VCF with one heterozygous SNV at chr_asm:51 (1-based).
    vcf_text = """##fileformat=VCFv4.2
##contig=<ID=chr_asm,length=200>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr_asm\t51\t.\tA\tG\t60\tPASS\t.\tGT\t0/1
"""
    raw_vcf = tmp_path / "asm.vcf"
    raw_vcf.write_text(vcf_text)
    pysam.tabix_compress(str(raw_vcf), str(vcf_path), force=True)
    pysam.tabix_index(str(vcf_path), preset="vcf", force=True)
    return bam_path, vcf_path


def test_asm_recovers_planted_signal(tmp_path):
    bam_path, vcf_path = _write_synth_bam_and_vcf(tmp_path)
    df = call_asm(
        bam={"sample1": bam_path},
        vcf=vcf_path,
        min_reads_per_haplotype=5,
        min_phased_snvs=1,
    )

    # At least one CpG should pass the filter (the one at pos 60).
    assert df.height >= 1
    # The planted CpG must be there: h1 all-methylated, h2 all-unmethylated.
    planted = df.filter(pl.col("pos") == 60)
    assert planted.height == 1, f"planted CpG missing; got {df}"
    row = planted.row(0, named=True)
    assert row["h1_meth"] == 10 and row["h1_unmeth"] == 0
    assert row["h2_meth"] == 0 and row["h2_unmeth"] == 10
    # Meth diff should be ~1.0 (fully methylated h1 vs unmethylated h2).
    assert abs(row["meth_diff"] - 1.0) < 1e-6
    # Fisher's exact on a (10,0) vs (0,10) table is highly significant.
    assert row["pvalue"] < 0.01
    # BH q-value column present.
    assert "qvalue" in df.columns
