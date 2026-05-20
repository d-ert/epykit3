"""Methylation entropy caller tests.

Two synthetic regions are planted:
  - An "ordered" region: all reads carry the same methylation pattern
    over the window -> entropy ~= 0.
  - A "disordered" region: reads carry every possible pattern with
    roughly equal frequency -> normalised entropy ~= 1.

Skips when pysam isn't available.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

pysam = pytest.importorskip("pysam", reason="pysam required for entropy tests.")

from epykit.entropy import call_entropy


def _write_synth_bam(tmp_path: Path, ordered: bool) -> Path:
    """Write a tiny BAM with planted methylation patterns.

    Four CpGs placed at offsets 0, 4, 8, 12 in a 16-bp read. 16 reads.
    - ordered=True  -> every read has methylation pattern 1100 ->
                      one pattern dominates -> low entropy.
    - ordered=False -> reads enumerate all 16 patterns -> max entropy.
    """
    bam_path = tmp_path / ("ordered.bam" if ordered else "disordered.bam")
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr_ent", "LN": 200}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for i in range(16):
            r = pysam.AlignedSegment(out.header)
            r.query_name = f"r{i}"
            r.query_sequence = "C" * 16
            r.flag = 0
            r.reference_id = 0
            r.reference_start = 0
            r.mapping_quality = 60
            r.cigar = [(0, 16)]
            r.query_qualities = pysam.qualitystring_to_array("E" * 16)
            # Methylation pattern: a 4-bit number -> bits placed at offsets 0,4,8,12.
            pattern = 0b1100 if ordered else i  # 0..15 enumerates all patterns
            xm = ["."] * 16
            for cpg_idx, offset in enumerate((0, 4, 8, 12)):
                meth = (pattern >> cpg_idx) & 1
                xm[offset] = "Z" if meth else "z"
            r.set_tag("XM", "".join(xm))
            out.write(r)

    pysam.sort("-o", str(bam_path), str(bam_path))
    pysam.index(str(bam_path))
    return bam_path


def test_entropy_low_for_ordered_region(tmp_path):
    bam_path = _write_synth_bam(tmp_path, ordered=True)
    df = call_entropy(
        bam={"s1": bam_path}, window_cpgs=4, min_reads=8,
    )
    assert df.height >= 1
    win = df.filter((pl.col("start") == 0) & (pl.col("end") == 13))
    assert win.height == 1
    # All 16 reads carry identical pattern -> entropy = 0.
    assert win["entropy"][0] == pytest.approx(0.0, abs=1e-9)
    assert win["normalised_entropy"][0] == pytest.approx(0.0, abs=1e-9)


def test_entropy_high_for_disordered_region(tmp_path):
    bam_path = _write_synth_bam(tmp_path, ordered=False)
    df = call_entropy(
        bam={"s1": bam_path}, window_cpgs=4, min_reads=8,
    )
    assert df.height >= 1
    win = df.filter((pl.col("start") == 0) & (pl.col("end") == 13))
    assert win.height == 1
    # 16 reads, 16 distinct patterns, uniform -> entropy = log2(16) = 4.
    assert win["entropy"][0] == pytest.approx(4.0, abs=1e-9)
    assert win["normalised_entropy"][0] == pytest.approx(1.0, abs=1e-9)


def test_entropy_rejects_invalid_window_cpgs(tmp_path):
    bam_path = _write_synth_bam(tmp_path, ordered=True)
    with pytest.raises(ValueError, match="window_cpgs"):
        call_entropy(bam={"s1": bam_path}, window_cpgs=1)
    with pytest.raises(ValueError, match="window_cpgs"):
        call_entropy(bam={"s1": bam_path}, window_cpgs=9)


def test_entropy_skips_windows_below_min_reads(tmp_path):
    bam_path = _write_synth_bam(tmp_path, ordered=True)
    df = call_entropy(bam={"s1": bam_path}, window_cpgs=4, min_reads=100)
    # 16 reads, threshold 100 -> no windows pass.
    assert df.height == 0
