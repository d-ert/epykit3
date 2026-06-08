"""M2: empirical_fdr=True was silently a no-op for chain_merge (the
default), sliding_window, and segment. Contract: empirical_fdr=True must
produce columns OR raise NotImplementedError; never silently no-op.

As of the count-ratio follow-up, the API (tl.dmr) supports empirical_fdr for
both ``tile`` and ``chain_merge``; only ``sliding_window`` / ``segment`` still
raise. The CLI (``_cmd_dmr``) still gates everything but ``tile`` (CLI wiring
for chain_merge is a deferred follow-up)."""
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


@pytest.mark.parametrize("method", ["sliding_window", "segment"])
def test_non_tile_empirical_fdr_raises_notimplemented(md_with_dmc, method):
    with pytest.raises(NotImplementedError, match="empirical_fdr.*tile"):
        ep.tl.dmr(md_with_dmc, method=method, empirical_fdr=True, n_perm=10)


def test_tile_empirical_fdr_still_works(md_with_dmc):
    ep.tl.dmr(md_with_dmc, method="tile", empirical_fdr=True, n_perm=5, perm_seed=0)
    dmr = md_with_dmc.uns["dmr"]
    assert "empirical_pvalue" in dmr.columns
    assert "empirical_qvalue" in dmr.columns


@pytest.mark.slow
def test_chain_merge_empirical_fdr_end_to_end(md_with_dmc):
    """API supports empirical_fdr for chain_merge: each shuffle recomputes the
    DMC then chain-merges. Real end-to-end on the synth bundle (n_perm small).
    dmr_params records the run regardless of how many DMRs survive; when DMRs
    exist the count-ratio columns are present."""
    ep.tl.dmr(md_with_dmc, method="chain_merge", empirical_fdr=True,
              n_perm=3, perm_seed=0)
    params = md_with_dmc.uns["dmr_params"]
    assert params["method"] == "chain_merge"
    assert params["empirical_fdr"] is True
    assert params["fdr_method"] == "region"
    dmr = md_with_dmc.uns["dmr"]
    if dmr.height > 0:
        assert "empirical_qvalue" in dmr.columns
        assert "empirical_fdr_set" in dmr.columns


@pytest.mark.parametrize("method", ["chain_merge", "sliding_window", "segment"])
def test_cli_dmr_non_tile_empirical_fdr_raises_notimplemented(method):
    """CLI mirror of the API gate. Pre-fix `_cmd_dmr` accepted --empirical-fdr
    against any --method and silently dropped it on non-tile callers; users
    were left thresholding combined_qvalue as if it were FDR-controlled."""
    import argparse
    from epykit.cli import _cmd_dmr
    args = argparse.Namespace(method=method, empirical_fdr=True)
    with pytest.raises(NotImplementedError, match=r"empirical[-_]fdr.*tile"):
        _cmd_dmr(args)


def test_tile_empirical_fdr_propagates_merge_adjacent_and_backend(
    md_with_dmc, monkeypatch,
):
    """m-perm-2: permutation null must use the same merge_adjacent and
    backend as the observed run; otherwise observed and null distributions
    are computed under different region definitions, producing a distorted
    empirical_pvalue. Mirrors Task 1.3's M3 fix for DMC."""
    import polars as pl
    import epykit.dmr as ep_dmr

    captured: list[dict] = []
    original = ep_dmr.call_dmr_tile_based

    def _fake_call_dmr_tile_based(*args, **kwargs):
        captured.append(dict(kwargs))
        # Use a deterministic non-empty result so the empirical machinery
        # proceeds through to the per-perm path (returning an empty frame
        # would short-circuit it). The observed call gets one real-looking
        # tile; per-perm calls get the same shape.
        return pl.DataFrame({
            "chrom": ["chr1"], "start": [100], "end": [200],
            "pvalue": [0.01], "meth_diff": [0.3],
            "qvalue": [0.01],
        })

    monkeypatch.setattr(ep_dmr, "call_dmr_tile_based", _fake_call_dmr_tile_based)
    # tl.dmr imports call_dmr_tile_based at module load, so we also patch the
    # alias bound in epykit.tl to intercept the observed call.
    import epykit.tl as ep_tl
    monkeypatch.setattr(ep_tl, "call_dmr_tile_based", _fake_call_dmr_tile_based)

    md_with_dmc.uns.pop("dmr", None)
    import epykit as ep
    ep.tl.dmr(
        md_with_dmc, method="tile",
        empirical_fdr=True, n_perm=2, perm_seed=0,
        merge_adjacent=False, backend="sequential",
    )

    # The observed call + per-perm calls all hit the patched function.
    # The OBSERVED call must carry merge_adjacent/backend; every per-perm
    # call must carry the same values.
    assert len(captured) >= 3, (
        f"Expected >=3 patched calls (1 observed + 2 perms); got {len(captured)}"
    )
    for kwargs in captured:
        assert kwargs.get("merge_adjacent") is False, (
            f"merge_adjacent not forwarded: {kwargs}"
        )
        assert kwargs.get("backend") == "sequential", (
            f"backend not forwarded: {kwargs}"
        )
