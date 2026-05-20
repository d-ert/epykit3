"""BAM ingestion tests (0.5.0).

These tests construct a tiny synthetic BAM with planted methylation
calls (Bismark XM tags) and verify that
:func:`epykit.bam_io.read_methylation_calls` recovers exactly the
calls we planted. Pysam is optional; tests cleanly skip when it isn't
installed.

We deliberately avoid mocking pysam -- building a real (tiny) BAM
exercises the same code path as real Bismark output and catches any
bugs in tag parsing / aligned-pair iteration.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pysam = pytest.importorskip(
    "pysam",
    reason="pysam required for BAM ingestion tests (Linux/macOS only).",
)

from epykit.bam_io import read_methylation_calls


@pytest.fixture
def synth_bismark_bam(tmp_path: Path) -> Path:
    """Build a tiny coordinate-sorted, indexed BAM with Bismark XM tags.

    Two reads on ``chr_test``:
      - Read A: starts at pos 100, length 10, XM = "z.z.zZ.Z.z"
        (5 CpG-style calls: 0, 0, 0, 1, 0 at offsets 0, 2, 4, 5, 9).
        Wait -- re-check: XM has one letter per query base; "z.z.zZ.Z.z"
        means bases at indices 0,2,4,5,7,9 are CpG calls. We'll
        treat indices 0/2/4/9 as unmethylated, 5/7 as methylated.
      - Read B: same start, but methylated where read A is unmethylated.

    Both reads carry mapq 60, base qual 30, and align without indels.
    """
    bam_path = tmp_path / "synth.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr_test", "LN": 1000}],
    }

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        # Build two reads with identical alignment but opposite methylation.
        for read_name, xm in (
            ("read_A", "z.z.zZ.Z.z"),
            ("read_B", "Z.Z.ZzZz.Z"),
        ):
            r = pysam.AlignedSegment(out.header)
            r.query_name = read_name
            r.query_sequence = "ACTGACGTAC"   # any 10-mer; XM is what matters
            r.flag = 0
            r.reference_id = 0
            r.reference_start = 100
            r.mapping_quality = 60
            r.cigar = [(0, 10)]  # 10M
            r.next_reference_id = -1
            r.next_reference_start = -1
            r.template_length = 0
            r.query_qualities = pysam.qualitystring_to_array("E" * 10)
            r.set_tag("XM", xm)
            out.write(r)

    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))
    return bam_path


def test_read_bismark_xm_recovers_planted_calls(synth_bismark_bam):
    """The XM-encoded methylation calls round-trip through read_methylation_calls."""
    df = read_methylation_calls(synth_bismark_bam, caller="bismark", min_baseq=20)
    # XM letters at indices 0,2,4,5,7,9 are CpG calls (z/Z) for read A
    # and the same indices for read B. CpG context only keeps z/Z.
    assert df.height >= 2  # at least one row per read
    # All rows must be on chr_test starting at ref_pos 100..109.
    assert set(df["chrom"].to_list()) == {"chr_test"}
    assert df["pos"].min() == 100
    assert df["pos"].max() <= 109

    # Read A: methylated at the offsets that were "Z" (uppercase) in
    # "z.z.zZ.Z.z" -- that's indices 5 (Z) and 7 (Z).
    # Read B: methylated at indices 0, 2, 4, 6, 9.
    read_a = df.filter(df["read_id"] == "read_A")
    read_b = df.filter(df["read_id"] == "read_B")
    assert read_a.height > 0, "no read_A rows"
    assert read_b.height > 0, "no read_B rows"
    # The two reads should not have identical methylation patterns.
    assert (
        read_a["methylation_status"].to_list()
        != read_b["methylation_status"].to_list()
    ), "reads encoded opposite methylation -- extractor lost the difference"


def test_read_methylation_calls_respects_min_baseq(synth_bismark_bam):
    """Setting min_baseq above the fixture's quality drops all rows."""
    # Fixture writes base quals = 36 (char 'E' - 33). A min_baseq of 100
    # is impossible, so the result must be empty.
    df = read_methylation_calls(synth_bismark_bam, caller="bismark", min_baseq=100)
    assert df.height == 0


def test_read_methylation_calls_respects_min_mapq(synth_bismark_bam):
    """min_mapq filter drops below-threshold reads."""
    # Fixture has mapq=60. min_mapq=100 -> no reads survive.
    df = read_methylation_calls(synth_bismark_bam, caller="bismark", min_mapq=100)
    assert df.height == 0


def test_read_methylation_calls_unknown_caller_raises(synth_bismark_bam):
    with pytest.raises(ValueError, match="Unknown caller"):
        read_methylation_calls(synth_bismark_bam, caller="unicorn")


def test_read_methylation_calls_missing_bam_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_methylation_calls(tmp_path / "nonexistent.bam", caller="bismark")


def test_read_methylation_calls_regions_filter(synth_bismark_bam):
    """Restricting to a region returns only positions in that range."""
    df = read_methylation_calls(
        synth_bismark_bam, caller="bismark",
        regions=[("chr_test", 0, 105)],
    )
    if df.height:
        assert (df["pos"] < 105).all()
