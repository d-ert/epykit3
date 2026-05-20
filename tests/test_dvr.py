"""Differentially Variable Regions -- density-based aggregation."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

import epykit as ep
from epykit.dvc import call_dvr_density


def _toy_dvc_table() -> pl.DataFrame:
    """Build a small DVC table with a clear DVC-dense region on chr1.

    chr1: positions 0..49 -- every 5th CpG is a DVC (10/50 = 20% locally).
    chr1: positions 1000..1049 -- 1 DVC out of 50 (background-like).
    chr2: 100 CpGs, 0 DVCs (pure null).
    """
    rows = []
    for i in range(50):
        rows.append({
            "chrom": "chr1", "pos": i,
            "var_log_ratio": 1.0 if i % 5 == 0 else 0.05,
            "is_dvc": (i % 5 == 0),
        })
    for i in range(50):
        rows.append({
            "chrom": "chr1", "pos": 1000 + i,
            "var_log_ratio": -0.5 if i == 0 else 0.05,
            "is_dvc": (i == 0),
        })
    for i in range(100):
        rows.append({
            "chrom": "chr2", "pos": i,
            "var_log_ratio": 0.05,
            "is_dvc": False,
        })
    return pl.DataFrame(rows)


def test_call_dvr_density_finds_enriched_tile():
    dvc = _toy_dvc_table()
    # Use small tiles so the dense region (positions 0..49) lands in
    # one or two tiles and the binomial test has resolution.
    dvr = call_dvr_density(
        dvc, tile_size_bp=50, min_cpgs_per_tile=5, alpha=0.05,
    )
    # Schema sanity.
    for col in (
        "chrom", "start", "end", "n_cpgs", "n_dvc", "frac_dvc",
        "pvalue", "qvalue", "mean_var_log_ratio", "dvr_type", "is_dvr",
    ):
        assert col in dvr.columns
    # Enriched tile: chr1 positions 0..49 with 10 DVCs out of 50.
    enriched = dvr.filter(
        (pl.col("chrom") == "chr1") & (pl.col("start") <= 0)
    )
    assert enriched.height == 1
    row = enriched.row(0, named=True)
    assert row["n_cpgs"] == 50
    assert row["n_dvc"] == 10
    # 10/50 = 20%, genome-wide background = 11/200 ~= 5.5%. Should clear
    # BH-q < 0.05.
    assert row["is_dvr"], (
        f"Enriched tile should be flagged DVR but isn't: q={row['qvalue']}"
    )
    # var_log_ratio is +1.0 on DVCs in this tile -> dvr_type 'var_up'.
    assert row["dvr_type"] == "var_up"


def test_call_dvr_density_empty_input():
    dvc = pl.DataFrame(schema={
        "chrom": pl.Utf8, "pos": pl.Int64,
        "var_log_ratio": pl.Float64, "is_dvc": pl.Boolean,
    })
    dvr = call_dvr_density(dvc)
    assert dvr.is_empty()


def test_call_dvr_density_no_dvcs():
    dvc = pl.DataFrame({
        "chrom": ["chr1"] * 10, "pos": list(range(10)),
        "var_log_ratio": [0.0] * 10, "is_dvc": [False] * 10,
    })
    dvr = call_dvr_density(dvc, tile_size_bp=100, min_cpgs_per_tile=5)
    # No is_dvr should be True since p0 = 0 -> every tile p-value = 1.
    if dvr.height > 0:
        assert not dvr.get_column("is_dvr").any()


def test_call_dvr_density_rejects_missing_columns():
    bad = pl.DataFrame({"chrom": ["chr1"], "pos": [0]})
    with pytest.raises(ValueError, match="missing columns"):
        call_dvr_density(bad)


def test_tl_dvr_orchestrator(synth_md_filtered):
    """End-to-end: tl.dvc -> tl.dvr -> md.uns['dvr'] populated."""
    ep.tl.dvc(synth_md_filtered)
    assert "dvc" in synth_md_filtered.varm
    ep.tl.dvr(synth_md_filtered, tile_size_bp=1000, min_cpgs_per_tile=3)
    assert "dvr" in synth_md_filtered.uns
    dvr_df = synth_md_filtered.uns["dvr"]
    assert isinstance(dvr_df, pl.DataFrame)
    params = synth_md_filtered.uns["dvr_params"]
    assert params["method"] == "density"
    assert params["tile_size_bp"] == 1000
    # The synth fixture isn't designed to produce DVCs (no variance-shift
    # seeds), so we don't require is_dvr.any() -- just that the call ran.


def test_tl_dvr_errors_without_dvc(synth_md_filtered):
    if "dvc" in synth_md_filtered.varm:
        del synth_md_filtered.varm["dvc"]
    with pytest.raises(ValueError, match="dvc"):
        ep.tl.dvr(synth_md_filtered)


# ---- DVC tests (merged from test_dvc.py) ----------------------------------


def test_dvc_writes_expected_schema(synth_md_filtered):
    md = synth_md_filtered
    ep.tl.dvc(md, test="bartlett")
    assert "dvc" in md.varm
    df = md.varm["dvc"]
    for col in (
        "chrom", "pos", "n_treatment", "n_control",
        "var_treatment", "var_control", "var_log_ratio",
        "p_variance", "q_variance", "p_mean", "q_mean", "is_dvc",
    ):
        assert col in df.columns, f"missing column {col}"
    assert df.schema["is_dvc"] == pl.Boolean


@pytest.mark.parametrize("bad_test", ["levene", "brown_forsythe", "f_test"])
def test_dvc_rejects_unsupported_tests(synth_md_filtered, bad_test):
    """Only 'bartlett' is supported; others should raise a clear ValueError."""
    md = synth_md_filtered
    with pytest.raises(ValueError, match="bartlett"):
        ep.tl.dvc(md, test=bad_test)
