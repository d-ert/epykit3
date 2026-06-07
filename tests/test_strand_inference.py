"""Strand inference from a reference FASTA must use the 0-based ``pos``.

Regression test for the C1-followup bug: Bismark ``.cov`` is 1-based, and
``_infer_strand`` previously indexed the 0-based reference array by the raw
1-based ``start``, reading the base one position 3' of the cytosine. That
landed on the forward G of the CpG dinucleotide and mislabelled every
+ strand CpG as ``-`` (and pushed the - strand C onto whatever followed the
dyad, often ``*``). The coordinate (C1) fix corrected the ``pos`` column but
not this lookup.

Reference chr1 = A A C G T T   (1-based positions 1..6)
                 1 2 3 4 5 6
A CpG dyad sits at 1-based 3(C)/4(G):
  * the + strand cytosine is the C reported at 1-based start = 3,
  * the - strand cytosine pairs with the forward G reported at start = 4.
Correct inference is therefore ['+', '-'].
"""
from __future__ import annotations

import polars as pl
import pytest

pytest.importorskip("pyfaidx")

from epykit.convert import _infer_strand, convert_sample


@pytest.fixture
def ref_fasta(tmp_path):
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\nAACGTT\n")
    return str(fa)


def test_infer_strand_uses_zero_based_pos(ref_fasta):
    # df carries the internal 0-based `pos` (1-based 3/4 -> 0-based 2/3).
    df = pl.DataFrame({"chrom": ["chr1", "chr1"], "pos": [2, 3]})
    assert _infer_strand(df, ref_fasta).to_list() == ["+", "-"]


def test_convert_sample_strand_end_to_end(tmp_path, ref_fasta):
    cov = tmp_path / "s1.cov"
    # 1-based Bismark: start == end. + strand C at 3, - strand C at 4.
    cov.write_text("chr1\t3\t3\t100\t5\t0\n" "chr1\t4\t4\t100\t3\t0\n")

    out = tmp_path / "store"
    resolved = convert_sample(
        str(cov),
        "s1",
        str(out),
        coordinate_base="one_based",
        reference_fasta=ref_fasta,
        merge_strands=False,  # isolate raw per-site strand (no pair merge)
    )
    assert resolved == "one_based"

    part = out / "sample=s1" / "chrom=chr1" / "part-0.parquet"
    res = pl.read_parquet(part).select(["pos", "strand"]).sort("pos")

    # pos is 0-based (1-based 3/4 -> 2/3); strand follows the cytosine base.
    assert res["pos"].to_list() == [2, 3]
    assert res["strand"].to_list() == ["+", "-"]
