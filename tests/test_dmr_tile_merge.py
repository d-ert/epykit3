"""Tests for _merge_adjacent_tiles in dmr.py."""
import math

import numpy as np
import polars as pl
import scipy.stats as sp_stats

from epykit.dmr import _merge_adjacent_tiles

_SCHEMA = {
    "chrom": pl.Utf8,
    "start": pl.Int32,
    "end": pl.Int32,
    "n_cpgs": pl.Int32,
    "n_case": pl.Int32,
    "n_control": pl.Int32,
    "mean_beta_case": pl.Float64,
    "mean_beta_control": pl.Float64,
    "meth_diff": pl.Float64,
    "log2_odds_ratio": pl.Float64,
    "pvalue": pl.Float64,
    "qvalue": pl.Float64,
    "dmr_type": pl.Utf8,
}


def _mk(rows: list[dict]) -> pl.DataFrame:
    defaults = {
        "n_case": 3, "n_control": 3,
        "mean_beta_case": 0.8, "mean_beta_control": 0.4,
        "log2_odds_ratio": 1.0, "qvalue": 0.01,
    }
    full = [{**defaults, **r} for r in rows]
    return pl.DataFrame(full, schema=_SCHEMA)


def test_tile_merges_adjacent_significant():
    """Two adjacent hyper tiles should merge into one."""
    df = _mk([
        {"chrom": "chr1", "start": 0,    "end": 1000, "n_cpgs": 10,
         "meth_diff": 0.3, "pvalue": 1e-5, "dmr_type": "hyper"},
        {"chrom": "chr1", "start": 1000, "end": 2000, "n_cpgs": 8,
         "meth_diff": 0.4, "pvalue": 1e-6, "dmr_type": "hyper"},
    ])
    out = _merge_adjacent_tiles(df)
    assert len(out) == 1
    assert out["start"][0] == 0
    assert out["end"][0] == 2000
    assert out["n_cpgs"][0] == 18
    assert out["dmr_type"][0] == "hyper"
    assert out["pvalue"][0] < 1e-5


def test_tile_no_merge_different_direction():
    """Adjacent hyper + hypo tiles should NOT merge."""
    df = _mk([
        {"chrom": "chr1", "start": 0,    "end": 1000, "n_cpgs": 10,
         "meth_diff": 0.3, "pvalue": 1e-5, "dmr_type": "hyper"},
        {"chrom": "chr1", "start": 1000, "end": 2000, "n_cpgs": 8,
         "meth_diff": -0.3, "pvalue": 1e-5, "dmr_type": "hypo"},
    ])
    out = _merge_adjacent_tiles(df)
    assert len(out) == 2


def test_tile_no_merge_gap():
    """Tiles separated by a gap should NOT merge."""
    df = _mk([
        {"chrom": "chr1", "start": 0,    "end": 1000, "n_cpgs": 10,
         "meth_diff": 0.3, "pvalue": 1e-5, "dmr_type": "hyper"},
        {"chrom": "chr1", "start": 2000, "end": 3000, "n_cpgs": 8,
         "meth_diff": 0.3, "pvalue": 1e-5, "dmr_type": "hyper"},
    ])
    out = _merge_adjacent_tiles(df)
    assert len(out) == 2


def test_tile_no_merge_different_chrom():
    """Tiles on different chromosomes should NOT merge."""
    df = _mk([
        {"chrom": "chr1", "start": 0,    "end": 1000, "n_cpgs": 10,
         "meth_diff": 0.3, "pvalue": 1e-5, "dmr_type": "hyper"},
        {"chrom": "chr2", "start": 0,    "end": 1000, "n_cpgs": 8,
         "meth_diff": 0.3, "pvalue": 1e-5, "dmr_type": "hyper"},
    ])
    out = _merge_adjacent_tiles(df)
    assert len(out) == 2


def test_tile_three_way_merge():
    """Three consecutive tiles should merge into one."""
    df = _mk([
        {"chrom": "chr1", "start": 0,    "end": 1000, "n_cpgs": 5,
         "meth_diff": 0.2, "pvalue": 1e-4, "dmr_type": "hyper"},
        {"chrom": "chr1", "start": 1000, "end": 2000, "n_cpgs": 7,
         "meth_diff": 0.3, "pvalue": 1e-5, "dmr_type": "hyper"},
        {"chrom": "chr1", "start": 2000, "end": 3000, "n_cpgs": 6,
         "meth_diff": 0.4, "pvalue": 1e-6, "dmr_type": "hyper"},
    ])
    out = _merge_adjacent_tiles(df)
    assert len(out) == 1
    assert out["start"][0] == 0
    assert out["end"][0] == 3000
    assert out["n_cpgs"][0] == 18


def test_tile_empty_input():
    """Empty input returns empty frame with same schema."""
    df = _mk([])
    out = _merge_adjacent_tiles(df)
    assert len(out) == 0
    assert set(out.columns) == set(df.columns)


def test_three_tile_chain_uses_correct_stouffer_denominator():
    """Three adjacent hyper tiles: combined p must be Stouffer-correct,
    not iterative-pairwise with a stale sqrt(2)."""
    # Three tiles with the same one-sided p (chosen so the math is easy).
    # Each two-sided p = 0.01 -> abs_z = norm.isf(0.005) ~ 2.5758
    # Sum z = 3 * 2.5758 = 7.7274
    # Correct Stouffer-combined z (n=3) = 7.7274 / sqrt(3) ~ 4.4615
    # Two-sided p = 2 * norm.sf(4.4615) ~ 8.13e-6
    p_each = 0.01
    abs_z = sp_stats.norm.isf(p_each / 2.0)
    expected_z = 3 * abs_z / math.sqrt(3)
    expected_p = 2.0 * sp_stats.norm.sf(expected_z)

    df = _mk([
        {"chrom": "chr1", "start":    0, "end": 1000, "n_cpgs": 10,
         "meth_diff": 0.3, "pvalue": p_each, "dmr_type": "hyper"},
        {"chrom": "chr1", "start": 1000, "end": 2000, "n_cpgs":  9,
         "meth_diff": 0.3, "pvalue": p_each, "dmr_type": "hyper"},
        {"chrom": "chr1", "start": 2000, "end": 3000, "n_cpgs": 11,
         "meth_diff": 0.3, "pvalue": p_each, "dmr_type": "hyper"},
    ])
    out = _merge_adjacent_tiles(df)
    assert len(out) == 1
    got_p = out["pvalue"][0]
    assert np.isclose(got_p, expected_p, rtol=1e-6), (
        f"Combined p={got_p:.3e}, expected Stouffer-3 {expected_p:.3e}. "
        f"Current code uses iterative pairwise /sqrt(2) and one-sided "
        f"isf on a two-sided p; both are wrong."
    )


def test_two_tile_merge_uses_two_sided_isf():
    """Even a 2-tile merge needs isf(p/2), not isf(p), because the
    input pvalue is two-sided."""
    p_each = 0.001
    abs_z = sp_stats.norm.isf(p_each / 2.0)
    expected_z = (abs_z + abs_z) / math.sqrt(2)
    expected_p = 2.0 * sp_stats.norm.sf(expected_z)

    df = _mk([
        {"chrom": "chr1", "start":    0, "end": 1000, "n_cpgs": 10,
         "meth_diff": 0.3, "pvalue": p_each, "dmr_type": "hyper"},
        {"chrom": "chr1", "start": 1000, "end": 2000, "n_cpgs": 10,
         "meth_diff": 0.3, "pvalue": p_each, "dmr_type": "hyper"},
    ])
    out = _merge_adjacent_tiles(df)
    assert len(out) == 1
    got_p = out["pvalue"][0]
    assert np.isclose(got_p, expected_p, rtol=1e-6)
