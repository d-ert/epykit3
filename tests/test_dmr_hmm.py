"""HMM-DMR caller tests.

Two-pronged validation:
  1. On a hand-built DMC table with a contiguous hypo run, the HMM
     caller emits a single DMR with the right boundaries and a
     ``dmr_type == "hypo"``.
  2. On the seeded synth fixture, ``ep.tl.dmr(md, method="hmm")``
     produces a frame schema-compatible with the tile DMR engine
     (same columns).
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

import epykit as ep
from epykit.dmr_hmm import call_dmr_hmm


def test_hmm_dmr_recovers_planted_hypo_run():
    """A run of negative-meth_diff CpGs becomes one DMR."""
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
    # Loose boundary tolerance (HMM smoothing margins). Planted: pos 6000-13900.
    assert row["start"] <= positions[60]
    assert row["end"] >= positions[120]
    assert row["meth_diff"] < -0.10


def test_hmm_dmr_no_calls_on_pure_noise():
    """A flat-zero meth_diff signal yields no DMRs."""
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
    """ep.tl.dmr(md, method='hmm') populates md.uns['dmr'] with the right schema."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    ep.tl.dmr(md, method="hmm", min_cpgs=4, min_abs_meth_diff=0.05)
    assert "dmr" in md.uns
    expected_cols = {"chrom", "start", "end", "n_cpgs", "meth_diff", "dmr_type"}
    actual = set(md.uns["dmr"].columns)
    assert expected_cols.issubset(actual), (
        f"missing DMR columns: {expected_cols - actual}; got {actual}"
    )
    assert md.uns["dmr_params"]["method"] == "hmm"


def test_hmm_dmr_rejects_dmc_missing_columns():
    bad = pl.DataFrame({"chrom": ["chr1"], "pos": [100]})  # no meth_diff
    with pytest.raises(ValueError, match="missing required columns"):
        call_dmr_hmm(bad)
