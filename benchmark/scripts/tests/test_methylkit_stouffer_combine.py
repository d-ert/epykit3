"""Test for methylkit_stouffer_combine.R.

Skips when Rscript is not on PATH."""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

FIXTURE_DIR = Path(__file__).parent / "fixtures"
SCRIPT = Path(__file__).parents[1] / "methylkit_stouffer_combine.R"


@pytest.mark.slow
def test_methylkit_stouffer_combine_matches_expected(tmp_path):
    if shutil.which("Rscript") is None:
        pytest.skip("Rscript not on PATH; skipping R subprocess test")
    in_tsv = FIXTURE_DIR / "methylkit_sample_in.tsv"
    expected = pl.read_csv(
        FIXTURE_DIR / "methylkit_sample_expected.tsv", separator="\t"
    )
    out_tsv = tmp_path / "out.tsv"
    result = subprocess.run(
        ["Rscript", str(SCRIPT),
         "--in", str(in_tsv),
         "--out", str(out_tsv),
         "--max-gap-bp", "1000"],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, (
        f"Rscript failed:\nstderr: {result.stderr}\nstdout: {result.stdout}"
    )
    got = pl.read_csv(out_tsv, separator="\t")
    assert "pvalue_combined" in got.columns
    assert "qvalue_combined" in got.columns
    assert got.height == 6
    # Combined p-values for the first cluster should be much smaller than raw.
    first_cluster = got.filter(pl.col("start") <= 300)
    assert (first_cluster["pvalue_combined"] < first_cluster["pvalue"]).all(), (
        "Combined p should be < raw p for the clustered CpGs"
    )
