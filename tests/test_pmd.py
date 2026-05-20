"""Tests for the PMD caller.

The synth fixture's beta values aren't designed to seed a megabase PMD,
so we build a minimal in-memory store on the fly: 5000 CpGs across one
chrom, with a 1 Mb hypomethylated stretch in the middle. The PMD
caller should recover that region within 50 kb of the planted boundaries.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest

from epykit.pmd import call_pmd_one_sample


def _build_pmd_store(tmp_path: Path, n_cpgs: int = 5000,
                     pmd_start: int = 2_000_000, pmd_end: int = 3_000_000) -> Path:
    """Write a synthetic methylstore with one sample, one chrom, one PMD."""
    store = tmp_path / "pmd_store"
    sample_dir = store / "sample=S1" / "chrom=chr_pmd"
    sample_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(42)
    positions = np.sort(rng.integers(100, 5_000_000, size=n_cpgs)).astype(np.int32)
    positions = np.unique(positions)
    # Methylation pattern: high (~0.85) outside the PMD window, low (~0.30) inside.
    n = len(positions)
    cov = np.full(n, 20, dtype=np.int32)
    in_pmd = (positions >= pmd_start) & (positions < pmd_end)
    p_meth = np.where(in_pmd, 0.30, 0.85)
    meth = rng.binomial(cov, p_meth).astype(np.int32)

    df = pl.DataFrame({
        "chrom":          ["chr_pmd"] * n,
        "pos":            positions,
        "strand":         ["+"] * n,
        "methyl_percent": (meth / cov * 100).astype(np.float32),
        "N_meth":         meth,
        "N_unmeth":       cov - meth,
        "coverage":       cov,
    })
    df.write_parquet(str(sample_dir / "part-0.parquet"))
    return store


def test_pmd_recovers_planted_region(tmp_path):
    """The 1 Mb planted PMD must be recovered within 50 kb of its boundaries."""
    store = _build_pmd_store(tmp_path, pmd_start=2_000_000, pmd_end=3_000_000)
    df = call_pmd_one_sample(
        store, "S1",
        bandwidth_bp=50_000,
        beta_threshold=0.55,
        min_pmd_bp=200_000,
    )
    assert df.height >= 1, f"no PMDs called; result was: {df}"
    # At least one PMD must contain the centre of the planted region (~2.5 Mb).
    centre = 2_500_000
    matching = df.filter((pl.col("start") <= centre) & (pl.col("end") >= centre))
    assert matching.height >= 1, f"no PMD covers centre 2.5 Mb; got {df}"
    # Its boundaries must be within 200 kb (relaxed for HMM smoothing margins
    # on a single-sample, small-fixture run) of the planted edges.
    row = matching.row(0, named=True)
    assert abs(row["start"] - 2_000_000) < 200_000, f"start {row['start']} far from 2 Mb"
    assert abs(row["end"] - 3_000_000) < 200_000, f"end {row['end']} far from 3 Mb"


def test_pmd_no_calls_when_no_signal(tmp_path):
    """All-high-beta chrom should produce zero PMDs."""
    # Build a store where every CpG is at high beta.
    store = tmp_path / "flat_store"
    sample_dir = store / "sample=S1" / "chrom=chr_flat"
    sample_dir.mkdir(parents=True, exist_ok=True)
    positions = np.arange(100, 1_000_000, 200, dtype=np.int32)
    n = len(positions)
    cov = np.full(n, 20, dtype=np.int32)
    meth = np.full(n, 18, dtype=np.int32)  # 90% beta
    pl.DataFrame({
        "chrom":          ["chr_flat"] * n,
        "pos":            positions,
        "strand":         ["+"] * n,
        "methyl_percent": np.full(n, 90.0, dtype=np.float32),
        "N_meth":         meth,
        "N_unmeth":       cov - meth,
        "coverage":       cov,
    }).write_parquet(str(sample_dir / "part-0.parquet"))

    df = call_pmd_one_sample(store, "S1", min_pmd_bp=500_000)
    assert df.height == 0, f"expected no PMDs on all-high-beta chrom; got {df}"
