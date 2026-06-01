"""HMM-DMR caller tests (backward-compat shim tests).

These tests validate that the deprecated ``epykit.dmr_hmm`` import path
still works (shim preserved for call_dmr_hmm), and that the 1.0 removal
of ``tl.dmr(method='hmm')`` is correctly enforced.

Two-pronged validation:
  1. On a hand-built DMC table with a contiguous hypo run, the shim
     (``call_dmr_hmm`` -> ``call_dmr_rule_segment``) emits a single DMR
     with the right boundaries and a ``dmr_type == "hypo"``.
  2. ``ep.tl.dmr(md, method="hmm")`` raises ``ValueError`` (the
     deprecation shim was removed in 1.0). The supported name is
     ``method="segment"``, which produces a frame schema-compatible
     with the tile DMR engine.
"""

from __future__ import annotations

import warnings

import numpy as np
import polars as pl
import pytest

import epykit as ep


def _import_call_dmr_hmm():
    """Import call_dmr_hmm via the deprecated shim, suppressing the DeprecationWarning."""
    import sys
    # Unload dmr_hmm to force re-evaluation of the shim warning.
    for mod in list(sys.modules.keys()):
        if "dmr_hmm" in mod:
            del sys.modules[mod]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        from epykit.dmr_hmm import call_dmr_hmm  # noqa: PLC0415
        return call_dmr_hmm


def test_hmm_dmr_recovers_planted_hypo_run():
    """A run of negative-meth_diff CpGs becomes one DMR (via deprecated shim)."""
    call_dmr_hmm = _import_call_dmr_hmm()
    n = 200
    positions = (np.arange(n) * 100 + 1000).astype(np.int32)
    # Methylation difference: zero everywhere except a planted hypo run.
    md = np.zeros(n, dtype=np.float64)
    md[50:130] = -0.40  # planted hypo run
    dmc = pl.DataFrame({
        "chrom":     ["chr1"] * n,
        "pos":       positions,
        "meth_diff": md.astype(np.float32),
    })
    dmrs = call_dmr_hmm(dmc, min_cpgs=10, min_abs_meth_diff=0.10, alpha=0.05)
    # Expect >= 1 DMR; at least one of dmr_type "hypo" overlapping our run.
    assert dmrs.height >= 1
    hypo = dmrs.filter(pl.col("dmr_type") == "hypo")
    assert hypo.height >= 1, "no hypo DMR called for planted hypo run"
    row = hypo.row(0, named=True)
    # Loose boundary tolerance (segmentation smoothing margins). Planted: pos 6000-13900.
    assert row["start"] <= positions[60]
    assert row["end"] >= positions[120]
    assert row["meth_diff"] < -0.10


def test_hmm_dmr_no_calls_on_pure_noise():
    """A flat-zero meth_diff signal yields no DMRs (via deprecated shim)."""
    call_dmr_hmm = _import_call_dmr_hmm()
    n = 200
    positions = (np.arange(n) * 100 + 1000).astype(np.int32)
    dmc = pl.DataFrame({
        "chrom":     ["chr1"] * n,
        "pos":       positions,
        "meth_diff": np.zeros(n, dtype=np.float32),
    })
    dmrs = call_dmr_hmm(dmc, min_cpgs=5, min_abs_meth_diff=0.10)
    assert dmrs.height == 0



def test_hmm_dmr_rejects_dmc_missing_columns():
    """Shim function still raises for missing required columns."""
    call_dmr_hmm = _import_call_dmr_hmm()
    bad = pl.DataFrame({"chrom": ["chr1"], "pos": [100]})  # no meth_diff
    with pytest.raises(ValueError, match="missing required columns"):
        call_dmr_hmm(bad)


def test_dmr_method_hmm_raises(synth_md_filtered):
    """method='hmm' was deprecated in 0.7.5 with FutureWarning; removed at 1.0.
    Now raises ValueError via the unknown-method dispatch path."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    with pytest.raises(ValueError, match="Unknown DMR method 'hmm'"):
        ep.tl.dmr(md, method="hmm")


def test_dmr_method_segment_still_works(synth_md_filtered):
    """method='segment' (the proper name) continues to work."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    ep.tl.dmr(md, method="segment")
    assert "dmr" in md.uns
    assert md.uns["dmr_params"]["method"] == "segment"
    # Schema check: the segment engine's frame must be compatible with the tile engine's columns.
    dmr_frame = md.uns["dmr"]
    expected_cols = {"chrom", "start", "end", "n_cpgs", "meth_diff", "dmr_type"}
    actual_cols = set(dmr_frame.columns)
    missing = expected_cols - actual_cols
    assert not missing, f"Segment DMR frame missing expected columns: {missing}"
