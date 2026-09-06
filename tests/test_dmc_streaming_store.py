"""Tests for the persistent per-chromosome DMC store cache layer in
``process_chromosomes_dmc``.

The function caches per-chromosome results in
``<methylstore>/.cache/dmc/<test>/`` keyed by an ``input_sig`` hash of
the call parameters. Two cache-hit paths exist:

* **Strict hit** -- manifest's ``input_sig`` matches the current call's
  signature exactly. Serve the cached result verbatim.
* **Weak hit** -- the manifest *predates* the ``input_sig`` field
  (legacy format) but every per-chrom parquet is on disk. Serve the
  cached result and upgrade the manifest in place.

The weak-hit branch must NOT fire when the manifest has an
``input_sig`` that simply differs from the current call -- that is a
real cache-invalidation event and must trigger recompute.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import epykit as ep
from epykit import _cache
from epykit._dmc_store import _MANIFEST_NAME
from epykit.dmc import _resolve_dmc_store_dir, process_chromosomes_dmc

pytestmark = pytest.mark.slow


def _run_dmc(md, treatment, control):
    """Run lr DMC on ``md`` with the given case/control split and return
    the result DataFrame."""
    store_path = Path(md.store)
    return process_chromosomes_dmc(
        methylstore_path=str(store_path),
        samples_treatment=list(treatment),
        samples_control=list(control),
        test="lr",
    )


def test_dmc_recomputes_when_input_sig_mismatches(synth_md_filtered):
    """Swapping samples between treatment and control between two calls
    must trigger a recompute, even if the per-chrom parquets from the
    first call are still on disk. Before the cache-invalidation fix
    the weak-hit branch would serve the stale frame and the two runs
    produced bit-identical results.
    """
    md = synth_md_filtered
    treat = list(md.treatment_ids)
    ctrl = list(md.control_ids)
    assert len(treat) >= 2 and len(ctrl) >= 2, (
        "fixture sanity: need at least 2 samples per group"
    )

    # First call: original labels.
    df_a = _run_dmc(md, treat, ctrl)

    # Sanity: the cache manifest now exists and has a non-trivial input_sig.
    staging = _resolve_dmc_store_dir(Path(md.store), "lr", None)
    manifest = _cache.load_json(staging / _MANIFEST_NAME)
    assert manifest is not None, "DMC cache manifest should be present after the first call"
    assert manifest.get("input_sig"), (
        "cache manifest must carry an input_sig (this regression test relies on it)"
    )

    # Second call: an *asymmetric* relabeling so the result genuinely
    # differs from call A. Move the last original treatment into
    # control; otherwise a full swap is two-sidedly symmetric for LR
    # and would yield bit-identical results post-fix.
    mixed_treat = treat[:-1]
    mixed_ctrl = [*ctrl, treat[-1]]
    df_b = _run_dmc(md, mixed_treat, mixed_ctrl)

    # The signatures differ, so the cache should have been invalidated
    # and recomputation should produce a different result frame.
    assert len(df_a) == len(df_b), "row count is invariant under relabeling"
    pvals_a = df_a["pvalue"].to_numpy()
    pvals_b = df_b["pvalue"].to_numpy()
    assert not np.array_equal(pvals_a, pvals_b), (
        "p-values are bit-identical across two different case/control "
        "splits -- the weak-hit cache branch is serving stale results."
    )

    # The manifest should have been rewritten with the new input_sig.
    manifest_b = _cache.load_json(staging / _MANIFEST_NAME)
    assert manifest_b is not None
    assert manifest_b.get("input_sig") != manifest.get("input_sig"), (
        "input_sig should have changed between the two calls"
    )


def test_dmc_legacy_manifest_without_input_sig_still_weak_hits(synth_md_filtered):
    """The weak-hit branch's documented use case: a cached manifest
    written by an older epykit version that has no ``input_sig``
    field. Serving the cached result and upgrading the manifest in
    place must continue to work.
    """
    md = synth_md_filtered
    treat = list(md.treatment_ids)
    ctrl = list(md.control_ids)

    # Populate the cache normally.
    df_a = _run_dmc(md, treat, ctrl)

    staging = _resolve_dmc_store_dir(Path(md.store), "lr", None)
    manifest_path = staging / _MANIFEST_NAME
    manifest = _cache.load_json(manifest_path)
    assert manifest is not None

    # Strip the input_sig field to simulate a manifest written by a
    # pre-input_sig epykit. Keep the per-chrom parquets intact.
    legacy_manifest = {k: v for k, v in manifest.items() if k != "input_sig"}
    assert "input_sig" not in legacy_manifest
    _cache.write_json(manifest_path, legacy_manifest)

    # Re-run with the *same* labels. The weak-hit branch should fire,
    # serve the cached parquets, and upgrade the manifest in place.
    df_b = _run_dmc(md, treat, ctrl)
    assert len(df_a) == len(df_b)

    # Manifest should have been upgraded with a fresh input_sig.
    upgraded = _cache.load_json(manifest_path)
    assert upgraded is not None
    assert upgraded.get("input_sig"), (
        "weak-hit path must upgrade legacy manifests with a fresh input_sig"
    )


# ---------------------------------------------------------------------------
# materialize= contract (C2): the default assembles the full per-CpG table
# onto md.varm; materialize=False keeps only the streaming DMCStore handle so
# peak memory stays O(largest chromosome) end-to-end.
# ---------------------------------------------------------------------------


def test_materialize_true_default_populates_varm(synth_md_filtered):
    """The default (materialize=True) assembles the full table onto md.varm,
    md.dmc and md.dmc_store both resolve, their sizes agree, and the canonical
    DMC schema is unchanged (1.0 back-compat)."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr", tsv=False)

    assert md.uns["dmc"]["materialized"] is True
    df = md.dmc
    assert df is not None
    store = md.dmc_store
    assert store is not None
    assert len(df) == store.total_sites == md.uns["dmc"]["n_sites"] > 0

    for col in (
        "chrom", "pos", "pvalue", "qvalue", "meth_diff",
        "meth_diff_ci_lo", "meth_diff_ci_hi",
    ):
        assert col in df.columns, f"canonical column {col!r} missing"


def test_materialize_false_keeps_only_store(synth_md_filtered):
    """materialize=False leaves no dmc_* key on md.varm but exposes the
    streaming DMCStore; md.dmc materialises on demand to the same size."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr", materialize=False)

    meta = md.uns["dmc"]
    assert meta["materialized"] is False
    # No eager varm table was assembled.
    assert not any(k.startswith("dmc") for k in md.varm), dict(md.varm)
    # Streaming handle is present and sized from the manifest.
    store = md.dmc_store
    assert store is not None
    assert store.total_sites == meta["n_sites"] > 0
    # On-demand materialisation works and matches the store.
    on_demand = md.dmc
    assert on_demand is not None
    assert len(on_demand) == store.total_sites


def test_materialize_false_dmr_equivalent(synth_md_filtered):
    """sliding-window DMR is driven by the on-disk DMCStore, so it produces
    the same result whether or not the eager table was materialised onto
    md.varm."""
    md = synth_md_filtered

    # Streaming path: no eager varm table.
    ep.tl.dmc(md, test="lr", materialize=False)
    ep.tl.dmr(md, method="sliding_window", window_bp=1000, step_bp=500, min_cpgs=2)
    dmr_stream = md.uns["dmr"]

    # Materialised path on the same store (per-chrom cache hit -> same
    # store_path); md.varm now carries the eager table.
    ep.tl.dmc(md, test="lr", materialize=True, tsv=False)
    assert any(k.startswith("dmc") for k in md.varm)
    ep.tl.dmr(md, method="sliding_window", window_bp=1000, step_bp=500, min_cpgs=2)
    dmr_mat = md.uns["dmr"]

    assert dmr_stream.equals(dmr_mat), (
        "sliding-window DMR differs between materialize=False and "
        "materialize=True paths -- it should be store-driven and identical"
    )


def test_materialize_false_rejects_neighbour_combine(synth_md_filtered):
    """materialize=False cannot run the eager-only post-processors; it raises
    rather than silently producing different output."""
    md = synth_md_filtered
    with pytest.raises(ValueError, match="materialize=False"):
        ep.tl.dmc(md, test="lr", materialize=False, neighbour_combine=True)
