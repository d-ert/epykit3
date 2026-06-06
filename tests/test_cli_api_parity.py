"""CLI <-> Python API parity for `dmc` (M-PKG2).

Historically `epykit dmc` called process_chromosomes_dmc without a
`dispersion` argument, inheriting that function's `dispersion="site"`
default, while `ep.tl.dmc` defaulted to `dispersion="eb"`. Identical input
therefore produced different q-values depending on whether the CLI or the
Python API was used -- and the benchmark paper claims around the API
default. These tests pin the two paths together.
"""

from __future__ import annotations

import inspect
import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest

import epykit as ep
from tests.fixtures.synth import SimConfig, generate


def test_tl_dmc_dispersion_default_is_eb():
    """The Python API default the CLI must mirror."""
    default = inspect.signature(ep.tl.dmc).parameters["dispersion"].default
    assert default == "eb"


def test_cli_dmc_qvalues_match_tl_dmc(tmp_path):
    """`epykit dmc` (default flags) and `ep.tl.dmc(test='lr')` produce
    identical q-values on the same store.

    min-samples are passed explicitly-equal on both sides so the test
    isolates the dispersion/reference defaults (the M-PKG2 regression).
    """
    cfg = SimConfig(
        n_per_group=3,
        chromosomes=("chr1",),
        cpgs_per_chrom=200,
        n_dmrs=1,
        n_scattered_dmcs=40,
    )
    res = generate(cfg, tmp_path / "synth")
    sheet = res["samplesheet"]

    # --- Python API path ---
    md_api = ep.read_bismark(
        sheet, treatment_group="treatment", control_group="control",
        store_dir=str(tmp_path / "api_store"),
    )
    ep.tl.dmc(
        md_api, test="lr",
        min_samples_treatment=0, min_samples_control=0,
    )  # dispersion="eb", reference="adaptive" by default
    api = md_api.dmc.select(["chrom", "pos", "qvalue"]).sort(["chrom", "pos"])

    # --- CLI path (fresh store to avoid cache collision) ---
    md_cli = ep.read_bismark(
        sheet, treatment_group="treatment", control_group="control",
        store_dir=str(tmp_path / "cli_store"),
    )
    out = tmp_path / "cli_dmc.parquet"
    proc = subprocess.run(
        [
            sys.executable, "-m", "epykit.cli", "dmc",
            "--methylstore", str(md_cli.store),
            "--samplesheet", sheet,
            "--treatment-group", "treatment",
            "--control-group", "control",
            "--output", str(out),
            "--test", "lr",
            "--min-samples-treatment", "0",
            "--min-samples-control", "0",
            "--no-csv",
        ],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, f"CLI dmc failed:\n{proc.stdout}\n{proc.stderr}"
    assert out.exists(), f"CLI produced no output; stderr:\n{proc.stderr}"

    cli = pl.read_parquet(out).select(["chrom", "pos", "qvalue"]).sort(["chrom", "pos"])

    assert api.height == cli.height, (
        f"row counts differ: api={api.height}, cli={cli.height}"
    )
    a = api["qvalue"].to_numpy()
    c = cli["qvalue"].to_numpy()
    assert np.allclose(a, c, atol=1e-12, equal_nan=True), (
        "CLI and tl.dmc q-values diverge -- CLI/API dispersion default "
        "mismatch (M-PKG2)?"
    )
