"""Coordinate-convention regression tests (C1).

Pin the 1-based Bismark .cov -> 0-based store ``pos`` contract end-to-end so a
future regression (treating 1-based input as 0-based, or vice versa) fails
loudly instead of silently shifting every CpG by 1 bp -- a plausible-but-wrong
error that survives per-CpG statistics (a uniform shift cancels) but corrupts
annotation, cross-format unite, and exported BED coordinates.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl

import epykit as ep


def _write_cov(path: Path, rows) -> None:
    """rows: iterable of (chrom, start, end, pct, m, u)."""
    with gzip.open(path, "wt", newline="") as fh:
        for chrom, start, end, pct, m, u in rows:
            fh.write(f"{chrom}\t{start}\t{end}\t{pct:.2f}\t{m}\t{u}\n")


def _write_bedgraph(path: Path, rows) -> None:
    with gzip.open(path, "wt", newline="") as fh:
        fh.write('track type="bedGraph" description="CpG methylation levels"\n')
        for chrom, start, end, pct, m, u in rows:
            fh.write(f"{chrom}\t{start}\t{end}\t{pct:.2f}\t{m}\t{u}\n")


def _sheet(path: Path, sample_paths: dict) -> None:
    path.write_text(
        "sample_id,group,path\n"
        + "\n".join(f"{sid},{grp},{p}" for sid, (grp, p) in sample_paths.items())
        + "\n"
    )


def _store_pos(md, sample: str, chrom: str) -> list[int]:
    part = Path(md.store) / f"sample={sample}" / f"chrom={chrom}" / "part-0.parquet"
    return pl.read_parquet(part).sort("pos")["pos"].to_list()


def test_bismark_one_based_cov_maps_to_zero_based_pos(tmp_path):
    """Real Bismark .cov is 1-based (start == end). A cytosine at 1-based 1000
    must land at 0-based pos 999."""
    _write_cov(tmp_path / "t1.bismark.cov.gz",
               [("chr1", 1000, 1000, 80.0, 8, 2), ("chr1", 2000, 2000, 50.0, 5, 5)])
    _write_cov(tmp_path / "c1.bismark.cov.gz",
               [("chr1", 1000, 1000, 80.0, 8, 2), ("chr1", 2000, 2000, 50.0, 5, 5)])
    sheet = tmp_path / "sheet.csv"
    _sheet(sheet, {"t1": ("treatment", tmp_path / "t1.bismark.cov.gz"),
                   "c1": ("control", tmp_path / "c1.bismark.cov.gz")})
    md = ep.read_bismark(str(sheet), treatment_group="treatment",
                         control_group="control", store_dir=str(tmp_path / "store"))
    assert _store_pos(md, "t1", "chr1") == [999, 1999]


def test_methyldackel_zero_based_unchanged(tmp_path):
    """MethylDackel bedGraph is 0-based half-open; pos == start."""
    _write_bedgraph(tmp_path / "t1.bedGraph.gz", [("chr1", 1000, 1001, 80.0, 8, 2)])
    _write_bedgraph(tmp_path / "c1.bedGraph.gz", [("chr1", 1000, 1001, 80.0, 8, 2)])
    sheet = tmp_path / "sheet.csv"
    _sheet(sheet, {"t1": ("treatment", tmp_path / "t1.bedGraph.gz"),
                   "c1": ("control", tmp_path / "c1.bedGraph.gz")})
    md = ep.read_methyldackel(str(sheet), treatment_group="treatment",
                              control_group="control", store_dir=str(tmp_path / "store"))
    assert _store_pos(md, "t1", "chr1") == [1000]


def test_same_cytosine_aligns_across_formats(tmp_path):
    """Bismark (1-based start=1001) and MethylDackel (0-based start=1000) for
    the SAME physical cytosine must produce the same store pos (1000)."""
    _write_cov(tmp_path / "bi_t.bismark.cov.gz", [("chr1", 1001, 1001, 80.0, 8, 2)])
    _write_cov(tmp_path / "bi_c.bismark.cov.gz", [("chr1", 1001, 1001, 80.0, 8, 2)])
    bi_sheet = tmp_path / "bi.csv"
    _sheet(bi_sheet, {"t1": ("treatment", tmp_path / "bi_t.bismark.cov.gz"),
                      "c1": ("control", tmp_path / "bi_c.bismark.cov.gz")})
    md_bi = ep.read_bismark(str(bi_sheet), treatment_group="treatment",
                            control_group="control", store_dir=str(tmp_path / "bi_store"))

    _write_bedgraph(tmp_path / "md_t.bedGraph.gz", [("chr1", 1000, 1001, 80.0, 8, 2)])
    _write_bedgraph(tmp_path / "md_c.bedGraph.gz", [("chr1", 1000, 1001, 80.0, 8, 2)])
    md_sheet = tmp_path / "md.csv"
    _sheet(md_sheet, {"t1": ("treatment", tmp_path / "md_t.bedGraph.gz"),
                      "c1": ("control", tmp_path / "md_c.bedGraph.gz")})
    md_md = ep.read_methyldackel(str(md_sheet), treatment_group="treatment",
                                 control_group="control", store_dir=str(tmp_path / "md_store"))

    assert _store_pos(md_bi, "t1", "chr1") == _store_pos(md_md, "t1", "chr1") == [1000]


def test_coordinate_base_override_forces_zero_based(tmp_path):
    """coordinate_base='zero_based' bypasses auto-detect on a start==end file."""
    _write_cov(tmp_path / "t1.bismark.cov.gz", [("chr1", 1000, 1000, 80.0, 8, 2)])
    _write_cov(tmp_path / "c1.bismark.cov.gz", [("chr1", 1000, 1000, 80.0, 8, 2)])
    sheet = tmp_path / "sheet.csv"
    _sheet(sheet, {"t1": ("treatment", tmp_path / "t1.bismark.cov.gz"),
                   "c1": ("control", tmp_path / "c1.bismark.cov.gz")})
    md = ep.read_bismark(str(sheet), treatment_group="treatment",
                         control_group="control", store_dir=str(tmp_path / "store"),
                         coordinate_base="zero_based")
    assert _store_pos(md, "t1", "chr1") == [1000]  # forced: no shift


def test_one_based_cov_annotates_to_correct_gene(tmp_path):
    """End-to-end: a 1-based Bismark CpG -> 0-based store pos -> annotate maps
    it to the gene it physically sits in (pins the ingestion->annotation chain)."""
    from epykit.annotate import annotate_features

    # Cytosine at 1-based 5001 -> 0-based store pos 5000.
    _write_cov(tmp_path / "t1.bismark.cov.gz", [("chr1", 5001, 5001, 80.0, 8, 2)])
    _write_cov(tmp_path / "c1.bismark.cov.gz", [("chr1", 5001, 5001, 80.0, 8, 2)])
    sheet = tmp_path / "sheet.csv"
    _sheet(sheet, {"t1": ("treatment", tmp_path / "t1.bismark.cov.gz"),
                   "c1": ("control", tmp_path / "c1.bismark.cov.gz")})
    md = ep.read_bismark(str(sheet), treatment_group="treatment",
                         control_group="control", store_dir=str(tmp_path / "store"))
    assert _store_pos(md, "t1", "chr1") == [5000]

    gtf = tmp_path / "g.gtf"
    gtf.write_text(
        'chr1\tt\tgene\t4001\t6000\t.\t+\t.\tgene_id "g1"; gene_name "G1";\n'
        'chr1\tt\texon\t4001\t6000\t.\t+\t.\tgene_id "g1"; gene_name "G1";\n'
    )
    sites = pl.DataFrame({"chrom": ["chr1"], "pos": [5000]})
    out = annotate_features(sites, str(gtf))
    assert out["gene_name"][0] == "G1"
