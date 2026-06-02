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


def _make_dmc_frame(with_combined: bool = False) -> pl.DataFrame:
    rows = {
        "chrom": ["chr1", "chr1", "chr1", "chr2", "chr2"],
        "pos":   [100, 200, 300, 50, 150],
        "meth_diff": [0.3, -0.4, 0.05, 0.2, -0.1],
        "pvalue": [1e-6, 1e-5, 0.5, 1e-3, 0.4],
        "qvalue": [1e-5, 1e-4, 0.6, 1e-2, 0.5],
    }
    if with_combined:
        # Flip combined values so significance differs from raw qvalue.
        rows["pvalue_combined"] = [0.5, 1e-7, 0.5, 0.5, 1e-9]
        rows["qvalue_combined"] = [0.6, 1e-6, 0.6, 0.6, 1e-8]
    return pl.DataFrame(rows)


def _stub_md_with_dmc(tmp_path: Path, *, key: str = "dmc_lr",
                       with_combined: bool = False):
    from epykit.methyldata import MethylData
    obs = pl.DataFrame({
        "sample_id": ["s1", "s2"], "group": ["case", "ctrl"],
    })
    store = tmp_path / "store"
    store.mkdir()
    md = MethylData(obs=obs, store=str(store))
    md.varm[key] = _make_dmc_frame(with_combined=with_combined)
    md.uns["dmc"] = {"last_key": key}
    return md


def test_dmc_to_tsv_significant_only_qvalue_asc(tmp_path):
    md = _stub_md_with_dmc(tmp_path)
    out = tmp_path / "dmc.significant.tsv"
    export.dmc_to_tsv(md, str(out))  # default alpha=0.05, full=False

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    # 3 rows have qvalue < 0.05 (1e-5, 1e-4, 1e-2); the other 2 (q=0.5, 0.6) are dropped.
    assert len(rows) == 3
    # qvalue ascending
    qvalues = [float(r["qvalue"]) for r in rows]
    assert qvalues == sorted(qvalues)


def test_dmc_to_tsv_full_writes_all_rows_genomic_order(tmp_path):
    md = _stub_md_with_dmc(tmp_path)
    out = tmp_path / "dmc.tsv"
    export.dmc_to_tsv(md, str(out), full=True)

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    assert len(rows) == 5
    assert [(r["chrom"], int(r["pos"])) for r in rows] == [
        ("chr1", 100), ("chr1", 200), ("chr1", 300),
        ("chr2", 50),  ("chr2", 150),
    ]


def test_dmc_to_tsv_alpha_override(tmp_path):
    md = _stub_md_with_dmc(tmp_path)
    out = tmp_path / "dmc.strict.tsv"
    export.dmc_to_tsv(md, str(out), alpha=1e-3)

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    # Only rows with qvalue < 1e-3 -> qvalue 1e-5 and 1e-4
    assert len(rows) == 2
    assert all(float(r["qvalue"]) < 1e-3 for r in rows)


def test_dmc_to_tsv_uses_qvalue_combined_when_present(tmp_path):
    md = _stub_md_with_dmc(tmp_path, with_combined=True)
    out = tmp_path / "dmc.combined.tsv"
    export.dmc_to_tsv(md, str(out))

    rows = list(csv.DictReader(
        out.read_text(encoding="utf-8").splitlines(), delimiter="\t",
    ))
    # qvalue_combined < 0.05 for rows with combined values 1e-6 and 1e-8.
    assert len(rows) == 2
    qc = sorted(float(r["qvalue_combined"]) for r in rows)
    assert qc == [1e-8, 1e-6]
