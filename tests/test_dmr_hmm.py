"""HMM-DMR caller tests (backward-compat shim tests).

These tests validate that the deprecated ``epykit.dmr_hmm`` import path
and ``tl.dmr(method='hmm')`` continue to work after the rename to
``dmr_segment`` in 0.7.5.

Two-pronged validation:
  1. On a hand-built DMC table with a contiguous hypo run, the shim
     (``call_dmr_hmm`` -> ``call_dmr_rule_segment``) emits a single DMR
     with the right boundaries and a ``dmr_type == "hypo"``.
  2. On the seeded synth fixture, ``ep.tl.dmr(md, method="hmm")``
     produces a frame schema-compatible with the tile DMR engine
     (same columns). The call raises a FutureWarning.
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


def test_hmm_dmr_via_tl_dmr_method_hmm(synth_md_filtered):
    """ep.tl.dmr(md, method='hmm') still works (with FutureWarning) and gives the right schema."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    with pytest.warns(FutureWarning, match="segment"):
        ep.tl.dmr(md, method="hmm", min_cpgs=4, min_abs_meth_diff=0.05)
    assert "dmr" in md.uns
    expected_cols = {"chrom", "start", "end", "n_cpgs", "meth_diff", "dmr_type"}
    actual = set(md.uns["dmr"].columns)
    assert expected_cols.issubset(actual), (
        f"missing DMR columns: {expected_cols - actual}; got {actual}"
    )
    # method='hmm' now maps to 'segment' internally
    assert md.uns["dmr_params"]["method"] == "segment"


def test_hmm_dmr_rejects_dmc_missing_columns():
    """Shim function still raises for missing required columns."""
    call_dmr_hmm = _import_call_dmr_hmm()
    bad = pl.DataFrame({"chrom": ["chr1"], "pos": [100]})  # no meth_diff
    with pytest.raises(ValueError, match="missing required columns"):
        call_dmr_hmm(bad)
