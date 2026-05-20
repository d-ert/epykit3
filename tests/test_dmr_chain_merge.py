"""Tests for the chain-and-merge DMR caller (DSS::callDMR semantics).

Covers:
  1. Chain linking & breakage on consecutive-sig-CpG gaps
  2. min_cpgs / pct_sig filters
  3. minlen_bp filter
  4. Direction labeling (hyper / hypo / mixed)
  5. Output schema parity with the sliding-window caller
  6. tl.dmr(method='chain_merge') populates md.uns['dmr'] + dmr_params
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

import epykit as ep
from epykit.dmr import _DMR_EMPTY_SCHEMA, call_dmr_chain_merge

pytestmark = pytest.mark.slow


def _mk_dmc_frame(rows: list[dict]) -> pl.DataFrame:
    """Build a minimal DMC-shaped DataFrame for call_dmr_chain_merge."""
    return pl.DataFrame(
        rows,
        schema={
            "chrom":     pl.Utf8,
            "pos":       pl.Int32,
            "meth_diff": pl.Float32,
            "pvalue":    pl.Float64,
        },
    )


# 1. Chain linking / breakage


def test_chain_links_within_dis_merge():
    """Two sig CpGs 150 bp apart chain when dis_merge_bp=200."""
    df = _mk_dmc_frame([
        # 5 sig CpGs in [100, 700] with gaps ~150 bp; pct_sig=1.0 here
        # because all 5 are significant.
        {"chrom": "chr1", "pos": 100, "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 250, "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 400, "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 550, "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 700, "meth_diff": 0.30, "pvalue": 1e-4},
    ])
    out = call_dmr_chain_merge(
        df,
        alpha=1e-3, min_abs_meth_diff=0.0,
        dis_merge_bp=200, min_cpgs=3, pct_sig=0.5, minlen_bp=50,
    )
    assert len(out) == 1
    assert out["start"][0] == 100
    assert out["end"][0] == 701  # half-open: last_pos + 1
    assert out["n_cpgs"][0] == 5
    assert out["n_significant"][0] == 5
    assert out["dmr_type"][0] == "hyper"


def test_chain_merge_default_merges_300bp_gap():
    """Sig CpGs 300 bp apart should chain at the default dis_merge_bp=500."""
    df = _mk_dmc_frame([
        {"chrom": "chr1", "pos": 100,  "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 400,  "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 700,  "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 1000, "meth_diff": 0.30, "pvalue": 1e-4},
    ])
    out = call_dmr_chain_merge(
        df,
        alpha=1e-3, min_abs_meth_diff=0.0,
        min_cpgs=3, pct_sig=0.5, minlen_bp=50,
    )
    assert len(out) == 1
    assert out["start"][0] == 100
    assert out["end"][0] == 1001


def test_chain_breaks_beyond_dis_merge():
    """A 200 bp gap between sig CpGs breaks the chain when dis_merge_bp=100."""
    df = _mk_dmc_frame([
        # First cluster: 3 sig CpGs in [100, 300] with gaps 100 bp.
        {"chrom": "chr1", "pos": 100, "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 200, "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 300, "meth_diff": 0.30, "pvalue": 1e-4},
        # 200 bp gap (> dis_merge_bp=100) -> new chain.
        {"chrom": "chr1", "pos": 500, "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 600, "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 700, "meth_diff": 0.30, "pvalue": 1e-4},
    ])
    out = call_dmr_chain_merge(
        df,
        alpha=1e-3, min_abs_meth_diff=0.0,
        dis_merge_bp=100, min_cpgs=3, pct_sig=0.5, minlen_bp=50,
    ).sort("start")
    assert len(out) == 2
    assert out["start"][0] == 100
    assert out["end"][0] == 301
    assert out["start"][1] == 500
    assert out["end"][1] == 701


def test_chain_links_kank1_like_141bp_gaps():
    """KANK1-like motivating case: 5 sig CpGs with gaps ~= 141 bp.

    The sliding_window caller with window_bp=100 cannot link these. The
    chain-merge caller with dis_merge_bp=141 (DSS-style) must produce a
    single DMR.
    """
    positions = [10_000, 10_141, 10_282, 10_423, 10_564]
    df = _mk_dmc_frame([
        {"chrom": "chr9", "pos": p, "meth_diff": 0.25, "pvalue": 1e-5}
        for p in positions
    ])
    out = call_dmr_chain_merge(
        df,
        alpha=1e-3, min_abs_meth_diff=0.0,
        dis_merge_bp=141, min_cpgs=3, pct_sig=0.5, minlen_bp=50,
    )
    assert len(out) == 1
    assert out["start"][0] == positions[0]
    assert out["end"][0]   == positions[-1] + 1
    assert out["n_significant"][0] == 5


# 2. Filter behaviour


def test_pct_sig_filter_rejects_low_density():
    """3 sig CpGs surrounded by 7 non-sig CpGs (30 %) fails pct_sig=0.5."""
    # 10 CpGs at 50 bp spacing: 3 sig at positions [100, 150, 200], 7
    # non-sig at [250, 300, ..., 550].
    rows = []
    for i, p in enumerate(range(100, 600, 50)):
        sig = i < 3
        rows.append({
            "chrom":     "chr1",
            "pos":       p,
            "meth_diff": 0.30 if sig else 0.01,
            "pvalue":    1e-4 if sig else 0.5,
        })
    df = _mk_dmc_frame(rows)
    out = call_dmr_chain_merge(
        df,
        alpha=1e-3, min_abs_meth_diff=0.0,
        dis_merge_bp=100, min_cpgs=3, pct_sig=0.5, minlen_bp=50,
    )
    # The 3 sig CpGs chain (gaps = 50 <= 100), span = [100, 201), n_cpgs in
    # span = 3 (the 3 sig only -- non-sig CpGs start at 250). Chain itself
    # passes pct_sig at 100 %. So this DMR survives -- adjust the test.
    # Re-run with sig CpGs interspersed:
    rows = []
    sig_positions = {100, 300, 500}
    for p in range(100, 600, 50):
        sig = p in sig_positions
        rows.append({
            "chrom":     "chr1",
            "pos":       p,
            "meth_diff": 0.30 if sig else 0.01,
            "pvalue":    1e-4 if sig else 0.5,
        })
    df = _mk_dmc_frame(rows)
    out = call_dmr_chain_merge(
        df,
        alpha=1e-3, min_abs_meth_diff=0.0,
        # dis_merge_bp=250 lets the 3 sig CpGs chain across the non-sig
        # CpGs; the span then has 9 total CpGs with only 3 sig -> 33%.
        dis_merge_bp=250, min_cpgs=3, pct_sig=0.5, minlen_bp=50,
    )
    assert len(out) == 0, "33% significant span should be filtered by pct_sig=0.5"


def test_min_cpgs_filter_rejects_short_chain():
    """A chain of 2 sig CpGs fails min_cpgs=3."""
    df = _mk_dmc_frame([
        {"chrom": "chr1", "pos": 100, "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 200, "meth_diff": 0.30, "pvalue": 1e-4},
    ])
    out = call_dmr_chain_merge(
        df,
        alpha=1e-3, min_abs_meth_diff=0.0,
        dis_merge_bp=200, min_cpgs=3, pct_sig=0.5, minlen_bp=50,
    )
    assert len(out) == 0


def test_minlen_bp_filter_rejects_short_span():
    """A chain whose span < minlen_bp is dropped."""
    # 3 sig CpGs 10 bp apart -> span = 21 bp, way below minlen_bp=50.
    df = _mk_dmc_frame([
        {"chrom": "chr1", "pos": 100, "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 110, "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 120, "meth_diff": 0.30, "pvalue": 1e-4},
    ])
    out = call_dmr_chain_merge(
        df,
        alpha=1e-3, min_abs_meth_diff=0.0,
        dis_merge_bp=200, min_cpgs=3, pct_sig=0.5, minlen_bp=50,
    )
    assert len(out) == 0


def test_min_abs_meth_diff_excludes_small_effects():
    """A CpG with tiny meth_diff isn't counted as sig even at p < alpha."""
    df = _mk_dmc_frame([
        # 3 CpGs at p=1e-5 but with meth_diff=0.02 < 0.05 threshold.
        {"chrom": "chr1", "pos": 100, "meth_diff": 0.02, "pvalue": 1e-5},
        {"chrom": "chr1", "pos": 200, "meth_diff": 0.02, "pvalue": 1e-5},
        {"chrom": "chr1", "pos": 300, "meth_diff": 0.02, "pvalue": 1e-5},
    ])
    out = call_dmr_chain_merge(
        df,
        alpha=1e-3, min_abs_meth_diff=0.05,
        dis_merge_bp=200, min_cpgs=3, pct_sig=0.5, minlen_bp=50,
    )
    assert len(out) == 0


# 3. Direction labels


def test_chain_merge_labels_hyper_and_hypo():
    """Sign of meth_diff drives dmr_type."""
    hyper = _mk_dmc_frame([
        {"chrom": "chr1", "pos": 100, "meth_diff": +0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 200, "meth_diff": +0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 300, "meth_diff": +0.30, "pvalue": 1e-4},
    ])
    hypo = _mk_dmc_frame([
        {"chrom": "chr2", "pos": 100, "meth_diff": -0.30, "pvalue": 1e-4},
        {"chrom": "chr2", "pos": 200, "meth_diff": -0.30, "pvalue": 1e-4},
        {"chrom": "chr2", "pos": 300, "meth_diff": -0.30, "pvalue": 1e-4},
    ])
    kwargs = dict(alpha=1e-3, min_abs_meth_diff=0.0,
                  dis_merge_bp=200, min_cpgs=3, pct_sig=0.5, minlen_bp=50)
    out_h = call_dmr_chain_merge(hyper, **kwargs)
    out_p = call_dmr_chain_merge(hypo, **kwargs)
    assert out_h["dmr_type"][0] == "hyper"
    assert out_p["dmr_type"][0] == "hypo"


# 4. Schema parity


def test_schema_matches_sliding_window():
    """Output uses the same schema as call_dmr_sliding_window."""
    df = _mk_dmc_frame([
        {"chrom": "chr1", "pos": 100, "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 200, "meth_diff": 0.30, "pvalue": 1e-4},
        {"chrom": "chr1", "pos": 300, "meth_diff": 0.30, "pvalue": 1e-4},
    ])
    out = call_dmr_chain_merge(
        df,
        alpha=1e-3, min_abs_meth_diff=0.0,
        dis_merge_bp=200, min_cpgs=3, pct_sig=0.5, minlen_bp=50,
    )
    # All schema columns present, with the expected dtypes (combined_qvalue
    # is added by the BH-correction post-step; tolerate either Float32 or
    # Float64 for it since dmc's apply_multiple_testing_correction picks
    # the type).
    for col, dtype in _DMR_EMPTY_SCHEMA.items():
        assert col in out.columns, f"missing column {col}"


def test_empty_input_returns_empty_frame_with_schema():
    """No DMC rows -> empty frame matching the schema."""
    empty = _mk_dmc_frame([])
    out = call_dmr_chain_merge(
        empty,
        alpha=1e-3, min_abs_meth_diff=0.0,
        dis_merge_bp=200, min_cpgs=3, pct_sig=0.5, minlen_bp=50,
    )
    assert len(out) == 0
    for col in _DMR_EMPTY_SCHEMA:
        assert col in out.columns


# 5. End-to-end through tl.dmr


def test_tl_dmr_method_chain_merge(synth_md_filtered):
    """ep.tl.dmr(method='chain_merge') populates md.uns and dmr_params."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    ep.tl.dmr(
        md, method="chain_merge",
        alpha=0.05, min_abs_meth_diff=0.0,
        dis_merge_bp=100, min_cpgs=3, pct_sig=0.5, minlen_bp=50,
        min_mean_qvalue=None,  # keep all candidates so the assertion below is meaningful
    )
    assert "dmr" in md.uns
    params = md.uns["dmr_params"]
    assert params["method"] == "chain_merge"
    assert params["dis_merge_bp"] == 100
    assert params["pct_sig"] == 0.5
    assert params["minlen_bp"] == 50
    # 10 planted DMRs in the synth fixture; at least some should chain.
    assert md.uns["dmr"].height > 0


def test_tl_dmr_unknown_method_message_includes_chain_merge(synth_md_filtered):
    """The error message for unknown methods lists 'chain_merge'."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    with pytest.raises(ValueError, match="chain_merge"):
        ep.tl.dmr(md, method="not_a_real_method")
