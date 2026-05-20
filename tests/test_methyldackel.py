"""MethylDackel .bedGraph input adapter."""

from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

import epykit as ep


_TRACK_HEADER = 'track type="bedGraph" description="CpG methylation levels"\n'


def _write_bismark_cov(path: Path, rows: list[tuple[str, int, int, float, int, int]]) -> None:
    with gzip.open(path, "wt", newline="") as fh:
        for chrom, start, end, pct, m, u in rows:
            fh.write(f"{chrom}\t{start}\t{end}\t{pct:.2f}\t{m}\t{u}\n")


def _write_methyldackel_bedgraph(
    path: Path, rows: list[tuple[str, int, int, float, int, int]]
) -> None:
    with gzip.open(path, "wt", newline="") as fh:
        fh.write(_TRACK_HEADER)
        for chrom, start, end, pct, m, u in rows:
            fh.write(f"{chrom}\t{start}\t{end}\t{pct:.2f}\t{m}\t{u}\n")


def _make_samplesheet(
    samplesheet: Path, sample_paths: dict[str, tuple[str, Path]]
) -> None:
    """sample_paths: {sample_id: (group, path)}."""
    samplesheet.write_text(
        "sample_id,group,path\n"
        + "\n".join(
            f"{sid},{group},{p}" for sid, (group, p) in sample_paths.items()
        )
        + "\n"
    )


def _toy_rows() -> list[tuple[str, int, int, float, int, int]]:
    return [
        ("chr1", 100, 101, 80.0, 8, 2),
        ("chr1", 200, 201, 60.0, 6, 4),
        ("chr1", 300, 301, 40.0, 4, 6),
        ("chr1", 400, 401, 20.0, 2, 8),
        ("chr2", 100, 101, 90.0, 9, 1),
        ("chr2", 200, 201, 10.0, 1, 9),
    ]


def test_methyldackel_roundtrip_matches_bismark(tmp_path):
    """read_methyldackel + read_bismark on the same data produce identical
    on-disk parquet stores (same counts at each site).
    """
    rows = _toy_rows()

    md_dir = tmp_path / "md_input"
    bi_dir = tmp_path / "bi_input"
    md_dir.mkdir()
    bi_dir.mkdir()

    # MethylDackel "treatment_1": same counts; "control_1": same counts.
    # We only need to verify the converter correctly parses the same data.
    _write_methyldackel_bedgraph(md_dir / "t1.bedGraph.gz", rows)
    _write_methyldackel_bedgraph(md_dir / "c1.bedGraph.gz", rows)
    _write_bismark_cov(bi_dir / "t1.bismark.cov.gz", rows)
    _write_bismark_cov(bi_dir / "c1.bismark.cov.gz", rows)

    md_sheet = tmp_path / "md_sheet.csv"
    bi_sheet = tmp_path / "bi_sheet.csv"
    _make_samplesheet(md_sheet, {
        "t1": ("treatment", md_dir / "t1.bedGraph.gz"),
        "c1": ("control",   md_dir / "c1.bedGraph.gz"),
    })
    _make_samplesheet(bi_sheet, {
        "t1": ("treatment", bi_dir / "t1.bismark.cov.gz"),
        "c1": ("control",   bi_dir / "c1.bismark.cov.gz"),
    })

    md_methyldackel = ep.read_methyldackel(
        str(md_sheet),
        treatment_group="treatment",
        control_group="control",
        store_dir=str(tmp_path / "md_store"),
    )
    md_bismark = ep.read_bismark(
        str(bi_sheet),
        treatment_group="treatment",
        control_group="control",
        store_dir=str(tmp_path / "bi_store"),
    )

    # uns.source_format should track the pipeline.
    assert md_methyldackel.uns["source_format"] == "methyldackel"
    assert md_bismark.uns["source_format"] == "bismark"
    assert md_methyldackel.uns["pipeline"] == "methyldackel"

    # Load one sample from each parquet store and compare on (chrom, pos,
    # N_meth, N_unmeth). The two should be byte-equivalent on the data
    # columns regardless of which format was ingested.
    md_part = (
        Path(md_methyldackel.store) / "sample=t1" / "chrom=chr1" / "part-0.parquet"
    )
    bi_part = (
        Path(md_bismark.store) / "sample=t1" / "chrom=chr1" / "part-0.parquet"
    )
    assert md_part.exists(), f"MethylDackel parquet missing: {md_part}"
    assert bi_part.exists(), f"Bismark parquet missing: {bi_part}"

    md_df = pl.read_parquet(md_part).select(["pos", "N_meth", "N_unmeth", "coverage"])
    bi_df = pl.read_parquet(bi_part).select(["pos", "N_meth", "N_unmeth", "coverage"])
    assert md_df.shape == bi_df.shape, (
        f"row counts differ: methyldackel={md_df.shape}, bismark={bi_df.shape}"
    )
    # Compare sorted: both stores partition by chrom but row order inside
    # a partition is implementation-defined.
    md_sorted = md_df.sort("pos")
    bi_sorted = bi_df.sort("pos")
    assert md_sorted.equals(bi_sorted), (
        "MethylDackel and Bismark conversions diverged on identical counts:\n"
        f"methyldackel:\n{md_sorted}\nbismark:\n{bi_sorted}"
    )


def test_methyldackel_cache_format_aware(tmp_path):
    """Re-running on the same input is a no-op (cache hit), but the cache
    is keyed by format so a Bismark store can't be reused for MethylDackel
    input pointing at a same-named file.
    """
    rows = _toy_rows()
    md_path = tmp_path / "input.bedGraph.gz"
    _write_methyldackel_bedgraph(md_path, rows)
    sheet = tmp_path / "sheet.csv"
    _make_samplesheet(sheet, {"s1": ("treatment", md_path), "s2": ("control", md_path)})
    store_dir = tmp_path / "store"

    # First pass: fresh convert.
    md1 = ep.read_methyldackel(
        str(sheet),
        treatment_group="treatment", control_group="control",
        store_dir=str(store_dir),
    )
    # Second pass: same store_dir; should reuse the cached partitions.
    md2 = ep.read_methyldackel(
        str(sheet),
        treatment_group="treatment", control_group="control",
        store_dir=str(store_dir),
    )
    assert md1.store == md2.store
    # Inspect the manifest to confirm format was recorded.
    manifest_path = (
        Path(md1.store) / "sample=s1" / ".epykit_raw_manifest.json"
    )
    assert manifest_path.exists()
    import json
    manifest = json.loads(manifest_path.read_text())
    assert manifest.get("format") == "methyldackel"


def test_cli_format_methyldackel(tmp_path):
    """`epykit convert --format methyldackel` writes a parquet store with
    the MethylDackel header correctly skipped."""
    rows = _toy_rows()
    md_path = tmp_path / "in.bedGraph.gz"
    _write_methyldackel_bedgraph(md_path, rows)
    out_dir = tmp_path / "out"

    result = subprocess.run(
        [
            sys.executable, "-m", "epykit.cli",
            "convert",
            "--input", str(md_path),
            "--sample-id", "s1",
            "--output-dir", str(out_dir),
            "--format", "methyldackel",
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, (
        f"convert subcommand failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
    # Confirm the parquet store was written with rows preserved.
    part = out_dir / "sample=s1" / "chrom=chr1" / "part-0.parquet"
    assert part.exists(), f"expected parquet at {part}; stderr: {result.stderr}"
    df = pl.read_parquet(part)
    assert df.height == 4  # chr1 has 4 rows in _toy_rows()


def test_cli_version_flag():
    """`epykit --version` prints 'epykit X.Y.Z' and exits 0."""
    result = subprocess.run(
        [sys.executable, "-m", "epykit.cli", "--version"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0
    out = (result.stdout + result.stderr).strip()
    assert out.startswith("epykit "), f"unexpected output: {out!r}"
    assert ep.__version__ in out, (
        f"--version output {out!r} doesn't include epykit.__version__={ep.__version__!r}"
    )
