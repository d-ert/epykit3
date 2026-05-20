"""HMR / LMR caller tests.

A synthetic store with two planted hypo regions of different CpG
density: one dense (HMR / CpG-island-like) and one sparse (LMR).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

from epykit.hmr import call_hmr_one_sample


def _build_hmr_store(tmp_path: Path) -> Path:
    """Build a synthetic store with:
       - a dense HMR (50 CpGs in 500 bp) at pos 5000-5500
       - a sparse LMR (4 CpGs spread across 2000 bp) at pos 50000-52000
       - high-beta filler everywhere else.
    """
    store = tmp_path / "hmr_store"
    sample_dir = store / "sample=S1" / "chrom=chr_hmr"
    sample_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(7)
    # Filler: 200 evenly-spaced high-beta CpGs across the chrom, EXCLUDING
    # the planted hypo regions so they aren't contaminated. Leave a
    # clear gap on either side of the LMR so the HMM can transition.
    filler_lo = np.linspace(100, 4000, 50).astype(np.int32)
    filler_mid = np.linspace(7000, 45_000, 100).astype(np.int32)
    filler_hi = np.linspace(56_000, 60_000, 50).astype(np.int32)
    # Dense HMR region: 50 CpGs in 500 bp.
    hmr = np.linspace(5000, 5500, 50).astype(np.int32)
    # Sparse LMR region: 8 CpGs across 2800 bp (density ~0.003, below 0.020).
    lmr = np.linspace(49_000, 51_800, 8).astype(np.int32)

    positions = np.unique(np.sort(np.concatenate([filler_lo, filler_mid, filler_hi, hmr, lmr])))

    in_hmr = (positions >= 5000) & (positions <= 5500)
    in_lmr = (positions >= 49_000) & (positions <= 51_800)
    cov = np.full(len(positions), 25, dtype=np.int32)
    p_meth = np.where(in_hmr | in_lmr, 0.15, 0.85)
    meth = rng.binomial(cov, p_meth).astype(np.int32)

    pl.DataFrame({
        "chrom":          ["chr_hmr"] * len(positions),
        "pos":            positions,
        "strand":         ["+"] * len(positions),
        "methyl_percent": (meth / cov * 100).astype(np.float32),
        "N_meth":         meth,
        "N_unmeth":       cov - meth,
        "coverage":       cov,
    }).write_parquet(str(sample_dir / "part-0.parquet"))
    return store


def test_hmr_recovers_dense_region(tmp_path):
    store = _build_hmr_store(tmp_path)
    df = call_hmr_one_sample(
        store, "S1",
        hmr_threshold=0.30,
        lmr_max_density=0.020,  # CpGs per bp
        min_cpgs=3,
    )
    # The dense HMR should be in the output, tagged kind="HMR".
    hmr_rows = df.filter(
        (df["start"] <= 5300) & (df["end"] >= 5300) & (df["kind"] == "HMR")
    )
    assert hmr_rows.height >= 1, f"dense HMR missing; got {df}"


def test_hmr_distinguishes_lmr_by_density(tmp_path):
    store = _build_hmr_store(tmp_path)
    df = call_hmr_one_sample(
        store, "S1",
        hmr_threshold=0.30,
        lmr_max_density=0.020,
        min_cpgs=3,
    )
    # The sparse LMR region (8 CpGs in 2800 bp = density ~0.003, < 0.020)
    # must be tagged kind="LMR".
    lmr_rows = df.filter(
        (df["start"] <= 50_000) & (df["end"] >= 50_000) & (df["kind"] == "LMR")
    )
    assert lmr_rows.height >= 1, (
        f"sparse region must be tagged LMR; got types {df['kind'].to_list()}"
    )


def test_hmr_no_calls_on_high_beta_only(tmp_path):
    """All-high-beta store yields zero HMRs/LMRs."""
    store = tmp_path / "flat"
    sd = store / "sample=S1" / "chrom=chr_flat"
    sd.mkdir(parents=True, exist_ok=True)
    positions = np.linspace(100, 100_000, 500).astype(np.int32)
    cov = np.full(len(positions), 20, dtype=np.int32)
    pl.DataFrame({
        "chrom": ["chr_flat"] * len(positions),
        "pos": positions,
        "strand": ["+"] * len(positions),
        "methyl_percent": np.full(len(positions), 90.0, dtype=np.float32),
        "N_meth": np.full(len(positions), 18, dtype=np.int32),
        "N_unmeth": np.full(len(positions), 2, dtype=np.int32),
        "coverage": cov,
    }).write_parquet(str(sd / "part-0.parquet"))

    df = call_hmr_one_sample(store, "S1", hmr_threshold=0.30, min_cpgs=4)
    assert df.height == 0, f"expected no calls on flat chrom; got {df}"
