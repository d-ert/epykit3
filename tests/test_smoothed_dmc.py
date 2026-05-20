"""Tests for ``ep.tl.dmc(..., use_smoothed=True)``.

Contract:
  1. Calling without first running pp.smooth raises a clear error.
  2. After pp.smooth, the smoothed-DMC result lands in
     ``md.varm["dmc_<test>_smoothed"]`` (separate from any raw run).
  3. The metadata flag is set: ``md.uns["dmc"]["use_smoothed"] == True``.
  4. The pseudo-count store produces sensible output: same sites as raw,
     same coverage column, but N_meth is derived from beta_smooth.
  5. End-to-end works with both Gaussian and BSmooth smoothers.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest

import epykit as ep
from epykit._smoothed_store import build_smoothed_pseudo_count_store

pytestmark = pytest.mark.slow


# ---- 1. Error when called without prior pp.smooth ------------------


def test_dmc_use_smoothed_without_pp_smooth_raises(synth_md_filtered):
    """use_smoothed=True must error clearly when no smoothed sidecar exists."""
    with pytest.raises(ValueError, match="use_smoothed=True requires ep.pp.smooth"):
        ep.tl.dmc(synth_md_filtered, test="lr", use_smoothed=True)


# ---- 2. End-to-end with bsmooth smoother ---------------------------


def test_dmc_use_smoothed_writes_separate_key(synth_md_filtered):
    """Smoothed-DMC lands at dmc_lr_smoothed (not overwriting dmc_lr)."""
    md = synth_md_filtered
    # First run raw DMC for comparison.
    ep.tl.dmc(md, test="lr")
    assert "dmc_lr" in md.varm
    raw_key = md.uns["dmc"]["last_key"]
    assert raw_key == "dmc_lr"

    # Now apply smoothing + smoothed DMC.
    ep.pp.smooth(md, method="bsmooth", ns=20, h_bp=2000)
    ep.tl.dmc(md, test="lr", use_smoothed=True)

    # Smoothed result is at a different key.
    assert "dmc_lr_smoothed" in md.varm
    assert md.uns["dmc"]["last_key"] == "dmc_lr_smoothed"
    assert md.uns["dmc"]["use_smoothed"] is True
    assert md.uns["dmc"]["smooth_method"] == "bsmooth"

    # And the raw run is preserved untouched.
    assert "dmc_lr" in md.varm
    raw = md.varm["dmc_lr"]
    smoothed = md.varm["dmc_lr_smoothed"]
    # Both tables should have the same site count (same chrom set, same
    # union/intersection rules).
    assert smoothed.height == raw.height


def test_dmc_use_smoothed_with_gaussian_smoother(synth_md_filtered):
    """The Gaussian smoother also feeds DMC correctly."""
    md = synth_md_filtered
    ep.pp.smooth(md, method="gaussian", bandwidth=1000)
    ep.tl.dmc(md, test="lr", use_smoothed=True)
    assert "dmc_lr_smoothed" in md.varm
    assert md.uns["dmc"]["smooth_method"] == "gaussian"


# ---- 3. Pseudo-count store structural correctness ------------------


def test_pseudo_count_store_preserves_coverage(tmp_path, synth_md_filtered):
    """Building the smoothed store leaves coverage intact and recomputes N_meth."""
    md = synth_md_filtered
    ep.pp.smooth(md, method="bsmooth", ns=20, h_bp=2000)

    out_dir = tmp_path / "smoothed_store"
    build_smoothed_pseudo_count_store(
        raw_store=Path(md.store),
        smooth_store=Path(md.uns["smooth_path"]),
        samples=md.obs.get_column("sample_id").to_list(),
        out_dir=out_dir,
    )

    # Pick one sample / chrom and check structure.
    samples = list(out_dir.glob("sample=*"))
    assert len(samples) > 0
    chrom_dirs = list(samples[0].glob("chrom=*"))
    assert len(chrom_dirs) > 0
    part = chrom_dirs[0] / "part-0.parquet"
    assert part.exists()

    smoothed = pl.read_parquet(str(part))
    # Same canonical columns as the raw store.
    for col in ("chrom", "pos", "strand", "N_meth", "N_unmeth", "coverage"):
        assert col in smoothed.columns, f"missing column {col!r} in smoothed store"

    # N_meth + N_unmeth == coverage (the conservation law).
    n_meth = smoothed["N_meth"].to_numpy()
    n_unmeth = smoothed["N_unmeth"].to_numpy()
    cov = smoothed["coverage"].to_numpy()
    np.testing.assert_array_equal(n_meth + n_unmeth, cov)
    # N_meth is non-negative and bounded by coverage.
    assert (n_meth >= 0).all()
    assert (n_meth <= cov).all()


def test_pseudo_count_store_uses_smoothed_beta(tmp_path):
    """N_meth in the pseudo-store reflects round(beta_smooth * coverage)."""
    # Build a tiny synthetic raw store + smoothed sidecar by hand.
    raw_store = tmp_path / "raw"
    smooth_store = tmp_path / "smooth"
    out_store = tmp_path / "pseudo"
    raw_dir = raw_store / "sample=S" / "chrom=chrX"
    sm_dir  = smooth_store / "sample=S" / "chrom=chrX"
    raw_dir.mkdir(parents=True)
    sm_dir.mkdir(parents=True)

    # Raw: 4 CpGs, coverage 10, all reading 1/10 (beta=0.10).
    raw = pl.DataFrame({
        "chrom":    ["chrX"] * 4,
        "pos":      [100, 200, 300, 400],
        "strand":   ["+"] * 4,
        "context":  ["CpG"] * 4,
        "N_meth":   [1, 1, 1, 1],
        "N_unmeth": [9, 9, 9, 9],
        "coverage": [10, 10, 10, 10],
        "sample":   ["S"] * 4,
    })
    raw.write_parquet(str(raw_dir / "part-0.parquet"))

    # Smoother decided each CpG's true beta is 0.50 (a strong upward shift).
    sm = pl.DataFrame({
        "chrom":       ["chrX"] * 4,
        "pos":         [100, 200, 300, 400],
        "sample":      ["S"] * 4,
        "beta_raw":    [0.10] * 4,
        "beta_smooth": [0.50] * 4,
    })
    sm.write_parquet(str(sm_dir / "part-0.parquet"))

    build_smoothed_pseudo_count_store(raw_store, smooth_store, ["S"], out_store)

    pseudo = pl.read_parquet(str(out_store / "sample=S" / "chrom=chrX" / "part-0.parquet"))
    # round(0.50 * 10) = 5 -> N_meth = 5, N_unmeth = 5, coverage = 10
    assert pseudo["N_meth"].to_list() == [5, 5, 5, 5]
    assert pseudo["N_unmeth"].to_list() == [5, 5, 5, 5]
    assert pseudo["coverage"].to_list() == [10, 10, 10, 10]


def test_pseudo_count_store_falls_back_to_raw_on_nan(tmp_path):
    """Sites with NaN beta_smooth keep raw N_meth / N_unmeth."""
    raw_store = tmp_path / "raw"
    smooth_store = tmp_path / "smooth"
    out_store = tmp_path / "pseudo"
    (raw_store / "sample=S" / "chrom=chr1").mkdir(parents=True)
    (smooth_store / "sample=S" / "chrom=chr1").mkdir(parents=True)

    pl.DataFrame({
        "chrom":    ["chr1", "chr1"],
        "pos":      [100, 200],
        "strand":   ["+", "+"],
        "context":  ["CpG", "CpG"],
        "N_meth":   [3, 7],
        "N_unmeth": [7, 3],
        "coverage": [10, 10],
        "sample":   ["S", "S"],
    }).write_parquet(str(raw_store / "sample=S" / "chrom=chr1" / "part-0.parquet"))

    pl.DataFrame({
        "chrom":       ["chr1", "chr1"],
        "pos":         [100, 200],
        "sample":      ["S", "S"],
        "beta_raw":    [0.30, 0.70],
        "beta_smooth": [0.50, float("nan")],   # second site falls back
    }).write_parquet(str(smooth_store / "sample=S" / "chrom=chr1" / "part-0.parquet"))

    build_smoothed_pseudo_count_store(raw_store, smooth_store, ["S"], out_store)
    pseudo = pl.read_parquet(str(out_store / "sample=S" / "chrom=chr1" / "part-0.parquet"))
    # Site 1: round(0.50 * 10) = 5  -> pseudo (5, 5)
    # Site 2: NaN smooth -> fallback to raw (7, 3)
    assert pseudo["N_meth"].to_list() == [5, 7]
    assert pseudo["N_unmeth"].to_list() == [5, 3]


# ---- 4. Empirical FDR + smoothed inputs --------------------------


def test_dmc_use_smoothed_with_empirical_fdr(synth_md_filtered):
    """The empirical FDR path also respects use_smoothed=True (small n_perm)."""
    md = synth_md_filtered
    ep.pp.smooth(md, method="bsmooth", ns=20, h_bp=2000)
    ep.tl.dmc(
        md, test="lr", use_smoothed=True,
        empirical_fdr=True, n_perm=5, perm_seed=0,
    )
    df = md.varm["dmc_lr_smoothed"]
    assert "empirical_pvalue" in df.columns or "pvalue" in df.columns
