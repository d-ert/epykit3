"""Smoke test for run_epykit_simulator.py.

Runs the simulator runner on a single headline seed (cell minimal) and
asserts the eval_per_seed parquet, IQR aggregate, and MANIFEST are
written with the expected schema. Skipped automatically when the
held-out simulator data has not been generated under
``benchmark/data/study1b_simulator/`` (the case for fresh clones until
``chore(benchmark): generate simulator data for Phase 4`` (commit
d438b81) has been replayed).
"""

from __future__ import annotations

import importlib
from pathlib import Path

import polars as pl
import pytest


HEADLINE_TEST_SEED = 2026000


def _simulator_data_available() -> bool:
    p = (
        Path(__file__).resolve().parents[2]
        / "data" / "study1b_simulator"
        / f"seed={HEADLINE_TEST_SEED}" / "truth.parquet"
    )
    return p.exists()


@pytest.mark.skipif(
    not _simulator_data_available(),
    reason="Held-out simulator data not generated; run simulate_piao.py first.",
)
def test_run_epykit_simulator_single_headline_seed(monkeypatch, tmp_path):
    """Run on a single headline seed and assert output parquets are well-formed.

    Uses the real simulator inputs (they are 100K rows on 1 chromosome --
    each cell takes a few seconds) but redirects all OUTPUTS to tmp_path
    so the real eval_per_seed.parquet is not clobbered. The convert /
    run caches are also redirected; the source AMP txt files under
    ``benchmark/data/study1b_simulator/seed=<S>/`` are read-only.
    """
    runner = importlib.import_module("run_epykit_simulator")
    importlib.reload(runner)

    # Redirect all OUTPUT paths to tmp_path. The SIM_BASE input data
    # under benchmark/data/study1b_simulator/seed=*/ stays at the real
    # location -- we only redirect things the runner writes to.
    tmp_sim_base = tmp_path / "study1b_simulator"
    # SIM_BASE input dir stays the same (it has the truth + AMP txt files);
    # but we redirect the per-output parquet paths into tmp_sim_base.
    monkeypatch.setattr(runner, "CONVERT_CACHE", tmp_path / "_converted_simulator")
    monkeypatch.setattr(runner, "RUN_CACHE", tmp_path / "_runs_simulator")
    monkeypatch.setattr(runner, "EVAL_PER_SEED", tmp_sim_base / "eval_per_seed.parquet")
    monkeypatch.setattr(
        runner, "EVAL_FROZEN_GRID", tmp_sim_base / "eval_frozen_grid.parquet"
    )
    monkeypatch.setattr(runner, "EVAL_SEED_IQR", tmp_sim_base / "eval_seed_iqr.parquet")
    monkeypatch.setattr(runner, "TIMINGS_PATH", tmp_sim_base / "timings_simulator.parquet")
    monkeypatch.setattr(runner, "MANIFEST_PATH", tmp_sim_base / "MANIFEST.txt")
    # SIM_BASE is read for inputs AND for mkdir() in main(); recreating
    # it at the original location is fine since main() only mkdirs (doesn't
    # truncate). Leave it alone.

    rc = runner.main(["--only", f"headline:{HEADLINE_TEST_SEED}", "--skip-ci"])
    assert rc == 0, "runner should exit 0 with no failed cells"

    # eval_per_seed.parquet should exist with one row per engine.
    assert runner.EVAL_PER_SEED.exists()
    df = pl.read_parquet(runner.EVAL_PER_SEED)
    assert df.height == 4, f"expected 4 rows (one per engine), got {df.height}"
    # All four engines present
    tools = set(df["tool"].to_list())
    assert tools == {
        "epykit_lr", "epykit_lrplus", "epykit_welch_t", "epykit_fisher",
    }, f"unexpected tools: {tools}"
    # Schema sanity: required columns
    for col in ("seed", "tool", "tpr", "fpr", "f1", "auroc",
                "tp", "fp", "tn", "fn", "n_dmc_called"):
        assert col in df.columns, f"missing column {col!r}"
    # Every engine got finite TPR / AUROC -- the cell is large enough
    # (100K CpGs, 20K true DMCs) that nothing should degenerate to NaN.
    for r in df.iter_rows(named=True):
        assert r["tpr"] is not None and r["tpr"] >= 0.0
        assert r["auroc"] is not None
        assert r["wall_s"] is not None and r["wall_s"] > 0.0

    # IQR aggregate exists; 4 rows (one per engine).
    assert runner.EVAL_SEED_IQR.exists()
    iqr = pl.read_parquet(runner.EVAL_SEED_IQR)
    assert iqr.height == 4, f"expected 4 IQR rows, got {iqr.height}"
    assert set(iqr["tool"].to_list()) == tools
    for col in ("tpr_median", "tpr_q1", "tpr_q3",
                "fpr_median", "f1_median", "auroc_median"):
        assert col in iqr.columns

    # Manifest written and non-empty.
    assert runner.MANIFEST_PATH.exists()
    manifest = runner.MANIFEST_PATH.read_text(encoding="utf-8")
    assert "study1b_simulator" in manifest
    assert "Engine tag" in manifest
    assert "lr+" in manifest


def test_iqr_handles_empty_input():
    """``compute_seed_iqr`` returns an empty (but well-typed) frame on empty input."""
    runner = importlib.import_module("run_epykit_simulator")

    empty = pl.DataFrame({
        "tool": [], "seed": [], "tpr": [], "fpr": [], "f1": [], "auroc": [],
    }, schema={
        "tool": pl.Utf8, "seed": pl.Int64,
        "tpr": pl.Float64, "fpr": pl.Float64,
        "f1": pl.Float64, "auroc": pl.Float64,
    })
    out = runner.compute_seed_iqr(empty)
    assert out.height == 0
    # Required columns present even when empty
    for col in ("tool", "tpr_median", "tpr_q1", "tpr_q3", "fpr_median",
                "f1_median", "auroc_median"):
        assert col in out.columns
