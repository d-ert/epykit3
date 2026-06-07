"""M2: empirical_fdr=True was silently a no-op for chain_merge (the
default), sliding_window, and segment. The calibration-warning note was
also suppressed in those branches, leaving users with no signal that
combined_qvalue is anti-conservative. Batch-1 contract: empirical_fdr=True
must produce columns (tile) OR raise NotImplementedError (others); never
silently no-op."""
from __future__ import annotations

import pytest

import epykit as ep


@pytest.fixture
def md_with_dmc(synth_md_filtered):
    """MethylData with a DMC result populated -- enough for tl.dmr to
    dispatch each method. Uses the canonical session-scoped two-group
    synth bundle from conftest.py."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    return md


@pytest.mark.parametrize("method", ["chain_merge", "sliding_window", "segment"])
def test_non_tile_empirical_fdr_raises_notimplemented(md_with_dmc, method):
    with pytest.raises(NotImplementedError, match="empirical_fdr.*tile"):
        ep.tl.dmr(md_with_dmc, method=method, empirical_fdr=True, n_perm=10)


def test_tile_empirical_fdr_still_works(md_with_dmc):
    ep.tl.dmr(md_with_dmc, method="tile", empirical_fdr=True, n_perm=5, perm_seed=0)
    dmr = md_with_dmc.uns["dmr"]
    assert "empirical_pvalue" in dmr.columns
    assert "empirical_qvalue" in dmr.columns
