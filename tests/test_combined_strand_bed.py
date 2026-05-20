"""Tests for the 12-column strand-collapsed methylation BED adapter.

The format (one row per CpG, tab-separated, no header)::

    chrom  start  end  fwd_M  fwd_T  fwd_pct  rev_M  rev_T  rev_pct  M  T  pct

The adapter must:
  1. Parse the 12 columns.
  2. Project to the canonical Bismark-layout methylstore using the
     combined-strand triplet (cols 10-12).
  3. Round-trip downstream: a converted file feeds ep.tl.dmc cleanly.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

import epykit as ep
from epykit.convert import convert_sample


def _write_synth_bed(path: Path, rows: list[tuple]) -> None:
    """Write a 12-col TSV. Each row: (chrom, start, end, fM, fT, fpct, rM, rT, rpct, M, T, pct)."""
    with path.open("w") as f:
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def test_convert_sample_combined_strand_bed_basic(tmp_path):
    """Six rows of the GSE263850 schema convert cleanly to the methylstore."""
    bed = tmp_path / "sample1.bed"
    _write_synth_bed(bed, [
        # chrom start end | fwd_M fwd_T fwd%   | rev_M rev_T rev%  | M T pct
        ("chr1", 10542, 10543,  4,  7, 57.14,  1,  1, 100.00,  5,  8, 62.50),
        ("chr1", 10563, 10564,  8,  8, 100.00, 1,  1, 100.00,  9,  9, 100.00),
        ("chr1", 10571, 10572,  8,  8, 100.00, 1,  1, 100.00,  9,  9, 100.00),
        ("chr1", 10577, 10578,  4,  8,  50.00, 1,  1, 100.00,  5,  9,  55.56),
        ("chr1", 10579, 10580,  9,  9, 100.00, 1,  1, 100.00, 10, 10, 100.00),
        ("chr2",  1000,  1001,  3, 10,  30.00, 0,  0,   0.00,  3, 10,  30.00),
    ])
    out_dir = tmp_path / "store"
    convert_sample(
        str(bed), "sample1", str(out_dir),
        format="combined_strand_bed",
    )
    chr1_parquet = out_dir / "sample=sample1" / "chrom=chr1" / "part-0.parquet"
    chr2_parquet = out_dir / "sample=sample1" / "chrom=chr2" / "part-0.parquet"
    assert chr1_parquet.exists()
    assert chr2_parquet.exists()

    chr1 = pl.read_parquet(str(chr1_parquet))
    # Row 1: pos=10542, N_meth=5, coverage=8, N_unmeth=3
    row0 = chr1.filter(pl.col("pos") == 10542).row(0, named=True)
    assert row0["N_meth"] == 5
    assert row0["coverage"] == 8
    assert row0["N_unmeth"] == 3

    # Row 4 (chr1): the asymmetric case 4/8 vs 1/1 -- combined is 5/9.
    row3 = chr1.filter(pl.col("pos") == 10577).row(0, named=True)
    assert row3["N_meth"] == 5
    assert row3["coverage"] == 9
    assert row3["N_unmeth"] == 4


def test_combined_strand_bed_zero_coverage_strand(tmp_path):
    """A row with one strand at coverage 0 must still parse and produce the right combined counts."""
    bed = tmp_path / "s.bed"
    _write_synth_bed(bed, [
        ("chr1", 100, 101,  3, 10, 30.00,  0, 0, 0.00,  3, 10, 30.00),
    ])
    convert_sample(str(bed), "s", str(tmp_path / "store"),
                   format="combined_strand_bed")
    df = pl.read_parquet(str(tmp_path / "store" / "sample=s" / "chrom=chr1" / "part-0.parquet"))
    row = df.row(0, named=True)
    assert row["N_meth"] == 3 and row["coverage"] == 10 and row["N_unmeth"] == 7


def test_combined_strand_bed_end_to_end_dmc(tmp_path):
    """A two-sample, two-condition synth BED runs through read -> filter -> dmc."""
    # Two samples per group, very simple methylation: treatment is hyper at 5 of 10 sites.
    chroms = ["chr1"] * 10
    positions = list(range(1000, 1100, 10))
    base_M = [3, 3, 3, 3, 3, 8, 8, 8, 8, 8]  # half hypo, half hyper (treatment)
    rev_pad = [(1, 1, 100.0)] * 10

    def write_sample(name: str, M_modifier: int):
        bed = tmp_path / f"{name}.bed"
        rows = []
        for c, p, mb, rp in zip(chroms, positions, base_M, rev_pad):
            M = max(0, min(10, mb + M_modifier))
            T = 10
            pct = 100.0 * M / T
            rows.append(
                (c, p, p + 1, M, T, pct, rp[0], rp[1], rp[2],
                 M + rp[0], T + rp[1], 100.0 * (M + rp[0]) / (T + rp[1]))
            )
        _write_synth_bed(bed, rows)
        return bed

    ctrl_a = write_sample("ctrl_a", 0)
    ctrl_b = write_sample("ctrl_b", 0)
    tr_a   = write_sample("tr_a", +2)   # treatment slightly more methylated
    tr_b   = write_sample("tr_b", +2)

    samplesheet = tmp_path / "sheet.csv"
    samplesheet.write_text(
        "sample_id,group,path\n"
        f"ctrl_a,control,{ctrl_a}\n"
        f"ctrl_b,control,{ctrl_b}\n"
        f"tr_a,treatment,{tr_a}\n"
        f"tr_b,treatment,{tr_b}\n"
    )
    md = ep.read_combined_strand_bed(
        str(samplesheet),
        treatment_group="treatment",
        control_group="control",
        store_dir=str(tmp_path / "ms"),
    )
    ep.pp.filter_coverage(md, lo_count=5, hi_perc=99.9)
    ep.pp.unite(md, type="intersect")
    ep.tl.dmc(md, test="lr")
    dmc = md.varm["dmc_lr"]
    assert dmc.height > 0
    assert "pvalue" in dmc.columns
    assert "meth_diff" in dmc.columns


def test_read_combined_strand_bed_exposed_on_package():
    assert hasattr(ep, "read_combined_strand_bed")
    assert callable(ep.read_combined_strand_bed)


def test_convert_sample_unknown_format_still_rejects(tmp_path):
    """Sanity check: unknown formats raise (back-compat with existing behaviour)."""
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    with pytest.raises(ValueError, match="Unknown format"):
        convert_sample(str(bed), "s", str(tmp_path / "out"), format="excel")
