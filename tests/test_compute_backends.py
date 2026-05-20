"""Non-regression tests for the _compute.run_chrom_pipeline dispatcher.

The 0.4.0 chrom-loop refactor introduced a backend dispatcher in
``src/epykit/_compute.py``. The contract:

  1. ``backend="sequential"`` produces results bit-identical to the
     pre-0.4 in-line loop (the prior tests in test_accuracy.py /
     test_dmc_ci_and_rename.py / test_dvc.py already cover this -- if
     they pass, the sequential path didn't regress).
  2. ``backend="dask"`` produces results numerically identical to the
     sequential path on the same fixture (just sliced across workers).
  3. Unknown backends raise ValueError. Missing optional deps raise
     ImportError with the install hint.
"""

from __future__ import annotations

import polars as pl
import pytest

import epykit as ep
from epykit._compute import run_chrom_pipeline

pytestmark = pytest.mark.slow


# ---- 1. Direct dispatcher tests --------------------------------------


def test_dispatcher_sequential_identity():
    """run_chrom_pipeline yields exactly what the handler returns, in order."""
    calls: list[str] = []

    def handler(chrom: str) -> pl.DataFrame:
        calls.append(chrom)
        return pl.DataFrame({"chrom": [chrom], "value": [hash(chrom) % 1000]})

    chroms = ["chr1", "chr2", "chr3"]
    results = list(run_chrom_pipeline(chroms, handler, backend="sequential"))

    assert [c for c, _ in results] == chroms
    assert calls == chroms
    assert results[0][1]["value"][0] == hash("chr1") % 1000


def test_dispatcher_skips_empty_and_none():
    """None / empty DataFrame returns are filtered out, not yielded."""

    def handler(chrom: str):
        if chrom == "skip_none":
            return None
        if chrom == "skip_empty":
            return pl.DataFrame({"chrom": [], "value": []})
        return pl.DataFrame({"chrom": [chrom], "value": [1]})

    chroms = ["chr1", "skip_none", "chr2", "skip_empty", "chr3"]
    yielded = [c for c, _ in run_chrom_pipeline(chroms, handler, backend="sequential")]
    assert yielded == ["chr1", "chr2", "chr3"]


def test_dispatcher_unknown_backend_raises():
    with pytest.raises(ValueError, match="Unknown compute backend"):
        list(run_chrom_pipeline(["chr1"], lambda c: None, backend="quantum"))


def test_dispatcher_dask_missing_extra_message():
    """Even on Dask-less environments, the error message names the extra."""
    try:
        import dask.distributed  # noqa: F401
        pytest.skip("dask is installed; cannot test the missing-extra error path")
    except ImportError:
        pass
    with pytest.raises(ImportError, match="epykit\\[distributed\\]"):
        list(run_chrom_pipeline(["chr1"], lambda c: None, backend="dask"))


def test_dispatcher_ray_missing_extra_message():
    try:
        import ray  # noqa: F401
        pytest.skip("ray is installed; cannot test the missing-extra error path")
    except ImportError:
        pass
    with pytest.raises(ImportError, match="epykit\\[ray\\]"):
        list(run_chrom_pipeline(["chr1"], lambda c: None, backend="ray"))


# ---- 2. End-to-end DMC parity (sequential vs dask) -------------------


def _sort_canonical(df: pl.DataFrame) -> pl.DataFrame:
    """Canonicalise row order so identity comparison is robust to
    iteration / submission-order differences."""
    return df.sort(["chrom", "pos"])


def _run_dmc(synth_md_filtered, *, backend: str):
    """Run tl.dmc with the requested backend, return the DMC frame."""
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr", backend=backend, n_workers=2)
    key = md.uns["dmc"]["last_key"]
    return md.varm[key]


def test_dmc_sequential_backend_returns_results(synth_md_filtered):
    """Smoke: the refactored sequential path still produces DMC output."""
    result = _run_dmc(synth_md_filtered, backend="sequential")
    assert len(result) > 0
    assert "pvalue" in result.columns
    assert "meth_diff" in result.columns


def test_dmc_dask_matches_sequential(synth_bundle, tmp_path):
    """Dask backend must produce numerically-identical DMC results.

    Skips when dask is not installed (most dev environments without
    [distributed] extra).
    """
    pytest.importorskip("dask.distributed")

    md_seq = ep.read_bismark(
        synth_bundle.samplesheet, treatment_group="treatment",
        control_group="control", assembly="synth",
        store_dir=str(tmp_path / "seq_store"),
    )
    ep.pp.filter_coverage(md_seq, lo_count=5, hi_perc=99.9)
    ep.pp.unite(md_seq, type="intersect")
    ep.tl.dmc(md_seq, test="lr", backend="sequential")
    seq_df = _sort_canonical(md_seq.varm[md_seq.uns["dmc"]["last_key"]])

    md_dask = ep.read_bismark(
        synth_bundle.samplesheet, treatment_group="treatment",
        control_group="control", assembly="synth",
        store_dir=str(tmp_path / "dask_store"),
    )
    ep.pp.filter_coverage(md_dask, lo_count=5, hi_perc=99.9)
    ep.pp.unite(md_dask, type="intersect")
    ep.tl.dmc(md_dask, test="lr", backend="dask", n_workers=2)
    dask_df = _sort_canonical(md_dask.varm[md_dask.uns["dmc"]["last_key"]])

    assert len(seq_df) == len(dask_df)
    for col in ("pvalue", "meth_diff", "log2_odds_ratio"):
        seq_vals = seq_df[col].to_numpy()
        dask_vals = dask_df[col].to_numpy()
        # Bit-identity for integer / categorical cols; tight tolerance
        # for floats (no algorithmic difference, only iteration order).
        import numpy as np
        np.testing.assert_allclose(
            seq_vals, dask_vals, rtol=1e-12, atol=1e-12, equal_nan=True,
            err_msg=f"DMC column {col!r} diverges between sequential and dask",
        )
