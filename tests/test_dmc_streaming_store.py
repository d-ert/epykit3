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
    mixed_ctrl = ctrl + [treat[-1]]
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
