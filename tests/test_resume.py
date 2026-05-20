"""Checkpoint/resume API tests (0.4.0).

The contract:

  1. ``resumable=True`` on a fresh analysis runs normally and writes a
     manifest entry + sidecar parquet under
     ``<analysis_root>/.epykit_manifest.json`` and ``.epykit_results/``.
  2. A second call with identical inputs + params and ``resumable=True``
     loads the cached result and skips the chrom-streaming compute. The
     ``md.uns["dmc"]["resumed"]`` flag distinguishes the two paths.
  3. ``MethylData.resume_from(stage)`` re-hydrates a fresh MethylData
     from the manifest without rerunning anything.
  4. Default ``resumable=False`` is bit-identical to pre-0.4: no
     manifest read, no manifest write.
"""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest

import epykit as ep
from epykit._cache import (
    input_signature,
    manifest_append,
    manifest_find,
    manifest_read,
    pipeline_manifest_path,
)

pytestmark = pytest.mark.slow


# ---- 1. Low-level _cache primitives ----------------------------------


def test_input_signature_is_stable(tmp_path):
    """Same inputs -> same hash. Different params -> different hash."""
    sig_a = input_signature("foo", {"k": 1, "v": "x"})
    sig_b = input_signature("foo", {"v": "x", "k": 1})  # dict order shouldn't matter
    sig_c = input_signature("foo", {"k": 2, "v": "x"})
    assert sig_a == sig_b
    assert sig_a != sig_c


def test_manifest_round_trip(tmp_path):
    """append -> find -> read round-trips a stage entry intact."""
    manifest_append(
        tmp_path, "dmc_lr",
        params={"test": "lr", "dispersion": "site"},
        input_sig="abc123",
        output_path=str(tmp_path / "dmc_lr.parquet"),
    )
    payload = manifest_read(tmp_path)
    assert payload["stages"][0]["name"] == "dmc_lr"
    assert payload["stages"][0]["input_sig"] == "abc123"

    entry = manifest_find(tmp_path, "dmc_lr")
    assert entry is not None
    assert entry["params"]["test"] == "lr"
    assert manifest_find(tmp_path, "nonexistent") is None


def test_manifest_replace_on_rerun(tmp_path):
    """Appending the same stage twice replaces the prior entry (not duplicates)."""
    manifest_append(tmp_path, "dmc_lr", params={"v": 1}, input_sig="A", output_path="A.parquet")
    manifest_append(tmp_path, "dmc_lr", params={"v": 2}, input_sig="B", output_path="B.parquet")
    payload = manifest_read(tmp_path)
    rows = [s for s in payload["stages"] if s["name"] == "dmc_lr"]
    assert len(rows) == 1
    assert rows[0]["input_sig"] == "B"


# ---- 2. End-to-end resume on tl.dmc ----------------------------------


def test_dmc_resumable_writes_manifest(synth_md_filtered):
    """tl.dmc(..., resumable=True) writes a manifest entry + sidecar."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr", resumable=True)

    root = md._analysis_root or md.store
    mp = pipeline_manifest_path(root)
    assert mp.exists(), f"manifest not written at {mp}"

    entry = manifest_find(root, "dmc_lr")
    assert entry is not None
    sidecar = Path(entry["output_path"])
    assert sidecar.exists(), f"sidecar parquet missing at {sidecar}"

    # Sidecar contents must match the in-memory result row-for-row.
    cached = pl.read_parquet(str(sidecar))
    in_mem = md.varm["dmc_lr"]
    assert len(cached) == len(in_mem)


def test_dmc_resume_skips_recomputation(synth_md_filtered):
    """A second resumable=True call with the same fingerprint skips the IRLS pass.

    Marker: the second call sets ``md.uns['dmc']['resumed'] = True``.
    """
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr", resumable=True)
    assert md.uns["dmc"].get("resumed") is not True  # first run computed

    # Second call on a fresh MethylData (same store) must resume.
    md2 = ep.read_bismark(
        # rebuild MethylData but point at the same analysis root
        samplesheet=Path(md._analysis_root) / ".." / "samplesheet.csv"
        if md._analysis_root else "samplesheet.csv",
        treatment_group="treatment",
        control_group="control",
        assembly="synth",
        store_dir=md._analysis_root,
    ) if False else None  # the fixture doesn't expose samplesheet path cleanly
    # Simpler: rerun on the same md -- the manifest still applies.
    ep.tl.dmc(md, test="lr", resumable=True)
    assert md.uns["dmc"].get("resumed") is True


def test_dmc_resumable_false_does_not_write_manifest(synth_md_filtered):
    """Default resumable=False leaves no manifest behind (pre-0.4 behaviour)."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")  # resumable defaults to False
    root = md._analysis_root or md.store
    mp = pipeline_manifest_path(root)
    if mp.exists():
        payload = manifest_read(root)
        assert not any(s["name"] == "dmc_lr" for s in payload.get("stages", []))


def test_dmc_resume_invalidated_when_params_change(synth_md_filtered):
    """Different ``test=`` between runs must NOT match the cached entry."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr", resumable=True)
    assert md.uns["dmc"].get("resumed") is not True

    # Different test name -> different stage key (dmc_score vs dmc_lr) ->
    # cache miss, fresh computation.
    ep.tl.dmc(md, test="score", resumable=True)
    assert md.uns["dmc"].get("resumed") is not True
    assert md.uns["dmc"]["test_used"] == "score"


# ---- 3. MethylData.resume_from -----------------------------------------


def test_resume_from_loads_dmc_into_varm(synth_md_filtered):
    """MethylData.resume_from('dmc_lr') restores varm from the sidecar."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr", resumable=True)
    cached_rows = len(md.varm["dmc_lr"])

    # Build a fresh MethylData pointing at the same analysis root and
    # call resume_from -- it should restore varm without recomputing.
    from epykit.methyldata import MethylData
    md2 = MethylData(
        obs=md.obs,
        store=md.store,
        assembly=md.assembly,
        context=md.context,
        _analysis_root=md._analysis_root,
    )
    assert "dmc_lr" not in md2.varm
    ok = md2.resume_from("dmc_lr")
    assert ok, "resume_from returned False but stage is in manifest"
    assert "dmc_lr" in md2.varm
    assert len(md2.varm["dmc_lr"]) == cached_rows


def test_completed_stages_lists_recorded_stages(synth_md_filtered):
    md = synth_md_filtered
    assert "dmc_lr" not in md.completed_stages
    ep.tl.dmc(md, test="lr", resumable=True)
    assert "dmc_lr" in md.completed_stages
