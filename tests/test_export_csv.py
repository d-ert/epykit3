"""Tests for CSV/TSV export of epykit result tables."""
from __future__ import annotations

import csv
from pathlib import Path

import polars as pl
from epykit import export


def _make_dmr_frame() -> pl.DataFrame:
    return pl.DataFrame({
        "chrom": ["chr1", "chr1", "chr2"],
        "start": [200, 100, 50],
        "end":   [300, 200, 150],
        "meth_diff": [0.3, -0.4, 0.2],
        "qvalue":  [0.01, 0.02, 0.5],
        "dmr_type": ["hyper", "hypo", "hyper"],
    })


def _stub_md_with_dmr(tmp_path: Path):
    """Build a minimal MethylData carrying a DMR table in md.uns['dmr']."""
    from epykit.methyldata import MethylData
    obs = pl.DataFrame({
        "sample_id": ["s1", "s2"],
        "group": ["case", "ctrl"],
    })
    store = tmp_path / "store"
    store.mkdir()
    md = MethylData(obs=obs, store=str(store))
    md.uns["dmr"] = _make_dmr_frame()
    return md


def test_dmr_to_tsv_writes_full_table_chrom_start_sorted(tmp_path):
    md = _stub_md_with_dmr(tmp_path)
    out = tmp_path / "dmr.tsv"
    export.dmr_to_tsv(md, str(out))

    text = out.read_text(encoding="utf-8")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    assert len(rows) == 3
    assert [(r["chrom"], int(r["start"])) for r in rows] == [
        ("chr1", 100), ("chr1", 200), ("chr2", 50),
    ]


def test_dmr_to_csv_uses_comma_for_csv_suffix(tmp_path):
    md = _stub_md_with_dmr(tmp_path)
    out = tmp_path / "dmr.csv"
    export.dmr_to_tsv(md, str(out))

    text = out.read_text(encoding="utf-8")
    # Header line must contain commas, no tabs.
    header = text.splitlines()[0]
    assert "," in header
    assert "\t" not in header
