"""Smoke test for run_epykit_study1.py.

Runs the smallest scenario cell (dmc_replicate=2, which is 1-vs-1 and only
allows the fisher engine) end-to-end and asserts that:
 * the per-cell parquet has the expected DMC schema columns,
 * the timings parquet has a row for it,
 * the reassembled eval_summary_post_phase3.parquet contains
   1) at least one fresh epykit_fisher row for the cell, and
   2) no stale epykit_bb_lr rows (Phase 3 removed the engine).

Skipped automatically if the raw Piao AMP simulator data is not wired in
under benchmark/raw_sim_data/ — that's the case for fresh clones until
``_legacy_benchmark`` is reattached.
"""

from __future__ import annotations

import importlib
from pathlib import Path

import polars as pl
import pytest


def _legacy_amp_available() -> bool:
    p = Path(__file__).resolve().parents[2] / "raw_sim_data" / "simulated_datasets" / "dmc_simulation" / "replicate" / "amp.replicate=2.sample1.txt"
    return p.exists()


@pytest.mark.skipif(
    not _legacy_amp_available(),
    reason="Piao AMP simulator data not wired in; need benchmark/raw_sim_data junction.",
)
def test_run_epykit_study1_dmc_replicate_2(monkeypatch, tmp_path):
    """Run the runner on dmc_replicate=2 only; assert outputs are present."""
    runner = importlib.import_module("run_epykit_study1")
    importlib.reload(runner)  # in case prior test mutated module-level paths

    # Redirect all output paths to tmp_path so we don't trash the real store.
    monkeypatch.setattr(runner, "DATA_STUDY1", tmp_path / "study1")
    monkeypatch.setattr(runner, "TRUTH_DIR", tmp_path / "study1" / "ground_truth")
    monkeypatch.setattr(runner, "OUT_BASE", tmp_path / "study1" / "epykit_post_phase3")
    monkeypatch.setattr(runner, "CONVERT_CACHE", tmp_path / "_converted")
    monkeypatch.setattr(runner, "RUN_CACHE", tmp_path / "_runs")
    monkeypatch.setattr(
        runner, "EVAL_SUMMARY_OLD",
        tmp_path / "study1" / "eval_summary.parquet",
    )
    monkeypatch.setattr(
        runner, "EVAL_SUMMARY_NEW",
        tmp_path / "study1" / "eval_summary_post_phase3.parquet",
    )
    monkeypatch.setattr(
        runner, "TIMINGS_NEW",
        tmp_path / "study1" / "timings_post_phase3.parquet",
    )

    # Stage a fake old eval_summary baseline with a few non-epykit rows + one
    # stale epykit_bb_lr row to verify the filter on reassembly.
    runner.DATA_STUDY1.mkdir(parents=True, exist_ok=True)
    old = pl.DataFrame({
        "tool": ["methylkit", "dss", "epykit_bb_lr"],
        "scenario": ["dmc_replicate", "dmc_replicate", "dmc_replicate"],
        "parameter": ["n_total", "n_total", "n_total"],
        "parameter_value": [2, 2, 2],
        "test": [None, None, "bb_lr"],
        "meth_diff_bin": ["all", "all", "all"],
        "threshold_kind": ["qvalue", "qvalue", "qvalue"],
        "threshold": [0.05, 0.05, 0.05],
        "tpr": [0.5, 0.6, 0.7],
        "fpr": [0.01, 0.01, 0.01],
        "precision": [0.9, 0.9, 0.9],
        "f1": [0.6, 0.7, 0.8],
        "auroc": [None, None, None],
        "tp": [None, None, None],
        "fp": [None, None, None],
        "tn": [None, None, None],
        "fn": [None, None, None],
    })
    old.write_parquet(runner.EVAL_SUMMARY_OLD)

    # Stage a minimal dmc_truth.parquet matching the n=2 replicate data's
    # chrom/pos so the score function has something to join against. The
    # real truth lives at benchmark/data/study1/ground_truth/ but we can
    # copy it over to keep this test self-contained.
    real_truth_dir = (
        Path(__file__).resolve().parents[3]
        / "benchmark" / "data" / "study1" / "ground_truth"
    )
    runner.TRUTH_DIR.mkdir(parents=True, exist_ok=True)
    for name in ("dmc_truth.parquet", "dmc_truth_dmr_sim.parquet", "dmr_truth.parquet"):
        src = real_truth_dir / name
        if src.exists():
            import shutil
            shutil.copyfile(src, runner.TRUTH_DIR / name)

    # Run only the dmc_replicate=2 cell. n_per_group=1 -> tests=('fisher',).
    rc = runner.main(["--only", "dmc_replicate:2"])
    assert rc == 0, "runner should exit 0 with no failed cells"

    # Per-cell parquet should exist with expected columns.
    per_cell = runner.OUT_BASE / "dmc_replicate" / "dmc_rep2_fisher.parquet"
    assert per_cell.exists(), f"expected parquet at {per_cell}"
    df = pl.read_parquet(per_cell)
    for col in ("chrom", "pos", "pvalue", "qvalue", "meth_diff"):
        assert col in df.columns, f"per-cell parquet missing {col!r}"

    # Timings parquet
    assert runner.TIMINGS_NEW.exists()
    timings = pl.read_parquet(runner.TIMINGS_NEW)
    sub = timings.filter(
        (pl.col("scenario") == "dmc_replicate")
        & (pl.col("parameter_value") == 2)
        & (pl.col("tool") == "epykit_fisher")
    )
    assert sub.height == 1, "expected one fisher timing row for dmc_rep2"
    assert sub["ok"][0] is True
    assert sub["wall_s"][0] > 0.0

    # Reassembled eval_summary
    assert runner.EVAL_SUMMARY_NEW.exists()
    summary = pl.read_parquet(runner.EVAL_SUMMARY_NEW)
    # Stale rows dropped
    assert "epykit_bb_lr" not in summary["tool"].unique().to_list(), (
        "epykit_bb_lr rows should be filtered out on reassembly"
    )
    # Non-epykit baseline rows preserved
    assert "methylkit" in summary["tool"].unique().to_list()
    # Fresh epykit_fisher rows added for the cell
    fresh = summary.filter(
        (pl.col("tool") == "epykit_fisher")
        & (pl.col("scenario") == "dmc_replicate")
        & (pl.col("parameter_value") == 2)
    )
    assert fresh.height > 0, (
        "expected fresh epykit_fisher rows for dmc_replicate=2 cell"
    )
