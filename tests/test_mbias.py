"""Bismark M-bias parser + pl.qc.mbias_plot smoke test."""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from epykit.nfcore_qc import parse_bismark_mbias
import epykit as ep


_TINY_PE_REPORT = """CpG context (R1)
================
position\tcount methylated\tcount unmethylated\t% methylation\tcoverage
1\t100\t50\t66.67\t150
2\t102\t48\t68.00\t150
3\t98\t52\t65.33\t150
4\t101\t49\t67.33\t150
5\t99\t51\t66.00\t150

CpG context (R2)
================
position\tcount methylated\tcount unmethylated\t% methylation\tcoverage
1\t90\t60\t60.00\t150
2\t92\t58\t61.33\t150
3\t94\t56\t62.67\t150
4\t95\t55\t63.33\t150
5\t96\t54\t64.00\t150

CHG context (R1)
================
position\tcount methylated\tcount unmethylated\t% methylation\tcoverage
1\t1\t199\t0.50\t200
2\t0\t200\t0.00\t200

CHH context (R1)
================
position\tcount methylated\tcount unmethylated\t% methylation\tcoverage
1\t1\t499\t0.20\t500
2\t2\t498\t0.40\t500
"""


def _write_report(path: Path) -> None:
    path.write_text(_TINY_PE_REPORT, encoding="utf-8")


def test_parse_bismark_mbias_round_trip(tmp_path):
    """Parser pulls every context/read panel into a long table."""
    p = tmp_path / "sample1.M-bias.txt"
    _write_report(p)

    df = parse_bismark_mbias(str(p))
    assert df.height == 5 + 5 + 2 + 2  # CpG R1, CpG R2, CHG R1, CHH R1
    assert set(df.get_column("context").unique().to_list()) == {"CpG", "CHG", "CHH"}
    assert set(df.get_column("read").unique().to_list()) == {"R1", "R2"}

    # Spot-check a few values.
    cpg_r1 = df.filter((pl.col("context") == "CpG") & (pl.col("read") == "R1"))
    assert cpg_r1.height == 5
    row1 = cpg_r1.row(0, named=True)
    assert row1["position"] == 1
    assert row1["n_meth"] == 100
    assert row1["n_unmeth"] == 50
    assert abs(row1["percent"] - 66.67) < 0.01
    assert row1["coverage"] == 150


def test_parse_bismark_mbias_empty_file_returns_empty_schema(tmp_path):
    p = tmp_path / "empty.M-bias.txt"
    p.write_text("not an mbias report\n", encoding="utf-8")
    df = parse_bismark_mbias(str(p))
    assert df.is_empty()
    # Schema columns are still defined so consumers can concat safely.
    for col in (
        "position", "context", "read", "n_meth", "n_unmeth", "percent", "coverage",
    ):
        assert col in df.columns


def test_mbias_plot_renders_from_paths(tmp_path):
    """pl.qc.mbias_plot accepts a dict of {sample: Path} and parses inline."""
    p1 = tmp_path / "s1.M-bias.txt"
    p2 = tmp_path / "s2.M-bias.txt"
    _write_report(p1)
    _write_report(p2)
    fig, ax = ep.pl.mbias_plot(
        {"s1": p1, "s2": str(p2)}, context="CpG",
    )
    # Two samples x R1 + R2 = 4 lines.
    assert len(ax.get_lines()) == 4
    assert ax.get_xlabel() == "Read position (bp)"
    # y-axis should be locked to 0..100 (percent methylation).
    ymin, ymax = ax.get_ylim()
    assert ymin == 0 and ymax == 100


def test_mbias_plot_renders_from_parsed_dataframes(tmp_path):
    """Pre-parsed DataFrames are accepted directly."""
    p1 = tmp_path / "s1.M-bias.txt"
    _write_report(p1)
    df = parse_bismark_mbias(str(p1))
    fig, ax = ep.pl.mbias_plot({"s1": df}, context="CHG")
    # CHG has only R1, 2 rows -> 1 line.
    assert len(ax.get_lines()) == 1


def test_mbias_plot_raises_on_missing_context(tmp_path):
    """Asking for a context not present in any sample is a hard error."""
    p1 = tmp_path / "s1.M-bias.txt"
    p1.write_text(
        "CpG context (R1)\n"
        "================\n"
        "position\tcount methylated\tcount unmethylated\t% methylation\tcoverage\n"
        "1\t1\t1\t50.00\t2\n",
        encoding="utf-8",
    )
    df = parse_bismark_mbias(str(p1))
    with pytest.raises(ValueError, match="No rows for context"):
        ep.pl.mbias_plot({"s1": df}, context="CHG")
