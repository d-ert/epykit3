"""Pytest configuration for the benchmark-script test suite.

These tests are independent of the main epykit test suite — run them via
`uv run pytest benchmark/scripts/tests/` from the repo root. They exercise
the simulator, null-calibration runner, and CI helpers without touching
the epykit package internals.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Make the scripts directory importable as a flat package so tests can do
# `from simulate_piao import simulate_dmc` rather than messing with PYTHONPATH.
SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

# Re-export the main tests/ conftest fixtures so benchmark-script
# tests can use synth_md, synth_md_filtered without duplicating fixture code.
# We import the fixture-generating module directly (not via tests.conftest) to
# avoid circular-import issues with pytest's own conftest collection.
_REPO = Path(__file__).resolve().parents[3]
_FIXTURES_DIR = _REPO / "tests" / "fixtures"
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

import importlib.util as _ilu

# Load tests/fixtures/synth.py by file path to avoid package-resolution issues.
_synth_spec = _ilu.spec_from_file_location(
    "tests.fixtures.synth",
    str(_FIXTURES_DIR / "synth.py"),
)
_synth_mod = _ilu.module_from_spec(_synth_spec)
sys.modules.setdefault("tests.fixtures.synth", _synth_mod)
_synth_spec.loader.exec_module(_synth_mod)  # type: ignore[union-attr]

from dataclasses import dataclass as _dataclass
from pathlib import Path as _P
from typing import Optional as _Optional

import numpy as _np
import polars as _pl
import pytest


SimConfig = _synth_mod.SimConfig
_generate = _synth_mod.generate


@_dataclass
class SynthBundle:
    """Bundle of paths + truth table + ids passed around by fixtures."""

    samplesheet: str
    truth: _pl.DataFrame
    store_root: str
    treatment_ids: list
    control_ids: list
    n_dmcs_true: int
    n_dmrs: int
    config: object


@pytest.fixture(scope="session")
def synth_bundle(tmp_path_factory) -> SynthBundle:
    """Generate the Bismark .cov fixture once per session."""
    cfg = SimConfig()
    out_dir = tmp_path_factory.mktemp("synth")
    result = _generate(cfg, out_dir)
    truth = _pl.read_parquet(result["truth"])
    return SynthBundle(
        samplesheet=result["samplesheet"],
        truth=truth,
        store_root=str(out_dir / "methyl_store"),
        treatment_ids=result["treatment_ids"],
        control_ids=result["control_ids"],
        n_dmcs_true=result["n_dmcs_true"],
        n_dmrs=result["n_dmrs"],
        config=cfg,
    )


@pytest.fixture
def synth_md(synth_bundle: SynthBundle, tmp_path):
    """Fresh MethylData pointing at the session methylstore."""
    import epykit as ep
    return ep.read_bismark(
        synth_bundle.samplesheet,
        treatment_group="treatment",
        control_group="control",
        assembly="synth",
        store_dir=str(tmp_path / "store"),
    )


@pytest.fixture
def synth_md_filtered(synth_md):
    """MethylData that has been filter_coverage'd; ready for DMC."""
    import epykit as ep
    ep.pp.filter_coverage(synth_md, lo_count=5, hi_perc=99.9)
    ep.pp.unite(synth_md, type="intersect")
    return synth_md
