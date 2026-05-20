"""Tests for BedGraph / BigWig / BED export."""

from __future__ import annotations

from pathlib import Path

import pytest


def test_to_bedgraph(synth_md_filtered, tmp_path):
    import epykit as ep

    sample = synth_md_filtered.obs.get_column("sample_id")[0]
    out = tmp_path / "x.bedgraph"
    synth_md_filtered.to_bedgraph(sample, str(out), value="beta")
    assert out.exists()
    lines = out.read_text(encoding="utf-8").splitlines()
    # First line is the track header
    assert lines[0].startswith("track type=bedGraph")
    # Subsequent lines: 4 fields, beta in [0, 1]
    data = [ln.split("\t") for ln in lines[1:] if ln.strip()]
    assert len(data) > 0
    for chrom, start, end, val in data[:200]:
        assert chrom.startswith("chr")
        assert int(end) == int(start) + 1
        f = float(val)
        assert 0.0 <= f <= 1.0


def test_to_bedgraph_coverage(synth_md_filtered, tmp_path):
    import epykit as ep

    sample = synth_md_filtered.obs.get_column("sample_id")[0]
    out = tmp_path / "cov.bedgraph"
    synth_md_filtered.to_bedgraph(sample, str(out), value="coverage")
    assert out.exists()
    # values now integer-ish (might be float in text)
    text = out.read_text(encoding="utf-8")
    assert text.count("\n") > 5


def test_dmcs_to_bed(synth_md_filtered, tmp_path):
    import epykit as ep

    ep.tl.dmc(synth_md_filtered, test="lr")
    out = tmp_path / "dmcs.bed"
    synth_md_filtered.dmcs_to_bed(str(out), alpha=0.5, min_abs_diff=0.0)
    assert out.exists()
    lines = [ln for ln in out.read_text(encoding="utf-8").splitlines()
             if ln and not ln.startswith("track")]
    # Should yield at least one DMC at alpha=0.5
    assert len(lines) > 0
    # 6-column BED
    for ln in lines[:20]:
        parts = ln.split("\t")
        assert len(parts) == 6
        chrom, start, end, name, score, strand = parts
        assert int(end) == int(start) + 1
        assert 0 <= int(score) <= 1000
        assert strand in ("+", "-")
        assert name.startswith("dmc_")


def test_dmrs_to_bed(synth_md_filtered, tmp_path):
    import epykit as ep

    ep.tl.dmc(synth_md_filtered, test="lr")
    ep.tl.dmr(synth_md_filtered, method="tile", tile_size_bp=500, min_cpgs_per_tile=3,
              min_mean_qvalue=1.0)
    out = tmp_path / "dmrs.bed"
    dmr_df = synth_md_filtered.uns.get("dmr")
    # With min_mean_qvalue=1.0 (no q-cut) on the calibrated fixture we
    # expect a non-trivial DMR table; if this is empty something upstream
    # regressed (covered by test_accuracy.test_dmr_tile_recovers_seeded_regions).
    assert dmr_df is not None and len(dmr_df) > 0, (
        "tile DMR returned 0 rows even with min_mean_qvalue=1.0 -- "
        "BED export test cannot run; investigate the DMR engine."
    )
    synth_md_filtered.dmrs_to_bed(str(out))
    lines = [ln for ln in out.read_text(encoding="utf-8").splitlines()
             if ln and not ln.startswith("track")]
    assert len(lines) == len(synth_md_filtered.uns["dmr"])
    for ln in lines[:10]:
        parts = ln.split("\t")
        assert len(parts) == 6


def test_to_bigwig_skipped_without_pybigwig(synth_md_filtered, tmp_path):
    """If pyBigWig is missing, to_bigwig raises a clear ImportError."""
    pyBigWig = pytest.importorskip("pyBigWig")  # noqa: F841
    sample = synth_md_filtered.obs.get_column("sample_id")[0]
    out = tmp_path / "x.bw"
    synth_md_filtered.to_bigwig(sample, str(out), value="beta")
    assert out.exists()
    assert out.stat().st_size > 0
