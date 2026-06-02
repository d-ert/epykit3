"""Tests for run_phase4_null_calibration.py.

Heavy methylstore builds + real subprocess invocations of
``run_null_calibration.py`` are skipped in CI -- those are integration runs.
Instead we exercise:

  * The aggregation kernel ``aggregate_summary`` on a synthetic on-disk
    parquet tree (so the per-(dataset, engine) median/IQR/CI math is locked).
  * The spec registry shape (right datasets, right scenarios, right engine
    sets, right k-shuffle counts).
  * The MANIFEST writer (smoke test that it produces a well-formed file).
  * Per-engine parquet path resolution.

The sealed ``run_null_calibration.py`` is not invoked here -- its own
test module covers its kernel.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest


def _write_fake_engine_parquet(
    path: Path, *, engine: str, scenario: str, k: int,
    fdrs: list[float], n_called: list[int], n_total: int,
) -> None:
    """Write a parquet matching run_null_calibration.py's per-shuffle output."""
    rows = []
    for i, (fdr, nc) in enumerate(zip(fdrs, n_called)):
        # Trivial "Wilson" stand-in: lo = fdr/2, hi = min(1, fdr*1.5).
        rows.append({
            "engine": engine,
            "scenario": scenario,
            "k_shuffle": i,
            "n_called": nc,
            "n_total": n_total,
            "observed_fdr": fdr,
            "observed_fdr_ci_lo": fdr / 2,
            "observed_fdr_ci_hi": min(1.0, fdr * 1.5),
        })
    path.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(rows).write_parquet(path)


def test_aggregate_summary_computes_per_cell_median_iqr(tmp_path, monkeypatch):
    """aggregate_summary collapses k rows per cell into one (median, IQR, CI)."""
    import run_phase4_null_calibration as mod

    # Repoint output dirs under tmp.
    monkeypatch.setattr(mod, "NULL_BASE", tmp_path / "null_calibration")
    monkeypatch.setattr(mod, "PIAO_OUT", tmp_path / "null_calibration" /
                        "piao_distributed" / "cov10_3v3")
    monkeypatch.setattr(mod, "SIM_OUT", tmp_path / "null_calibration" /
                        "simulator" / "sim_cov10_3v3")
    monkeypatch.setattr(mod, "GSE_OUT", tmp_path / "null_calibration" /
                        "gse263850")
    monkeypatch.setattr(mod, "SUMMARY_PATH",
                        tmp_path / "null_calibration" / "summary.parquet")

    # Patch spec out_dirs to the redirected paths.
    specs = mod.make_specs()
    specs["piao"].out_dir = mod.PIAO_OUT
    specs["simulator"].out_dir = mod.SIM_OUT
    specs["gse"].out_dir = mod.GSE_OUT
    # Narrow engine sets to keep the synthetic fixture small.
    specs["piao"].engines = ("lr", "fisher")
    specs["simulator"].engines = ("lr",)
    specs["gse"].engines = ("lr_plus",)

    # Piao/lr: 5 shuffles with observed_fdr [0.01..0.05]; expect median=0.03.
    _write_fake_engine_parquet(
        mod.PIAO_OUT / "lr.parquet",
        engine="lr", scenario="cov10_3v3", k=5,
        fdrs=[0.01, 0.02, 0.03, 0.04, 0.05],
        n_called=[10, 20, 30, 40, 50], n_total=1000,
    )
    _write_fake_engine_parquet(
        mod.PIAO_OUT / "fisher.parquet",
        engine="fisher", scenario="cov10_3v3", k=3,
        fdrs=[0.10, 0.20, 0.30],
        n_called=[100, 200, 300], n_total=1000,
    )
    _write_fake_engine_parquet(
        mod.SIM_OUT / "lr.parquet",
        engine="lr", scenario="sim_cov10_3v3", k=4,
        fdrs=[0.04, 0.05, 0.06, 0.07],
        n_called=[40, 50, 60, 70], n_total=1000,
    )
    _write_fake_engine_parquet(
        mod.GSE_OUT / "lr_plus.parquet",
        engine="lr_plus", scenario="gse263850", k=2,
        fdrs=[0.001, 0.002], n_called=[1, 2], n_total=1000,
    )

    summary = mod.aggregate_summary(specs)
    assert summary.height == 4
    # Expected schema (order doesn't matter; presence does).
    expected_cols = {
        "engine", "dataset", "scenario", "k_shuffles",
        "observed_fdr_median", "observed_fdr_q1", "observed_fdr_q3",
        "observed_fdr_ci_lo", "observed_fdr_ci_hi",
        "n_sites_called_median", "n_sites_total",
    }
    assert expected_cols.issubset(set(summary.columns))

    # Piao/lr median = 0.03 (middle of 0.01..0.05).
    piao_lr = summary.filter(
        (pl.col("dataset") == "piao_distributed") & (pl.col("engine") == "lr")
    ).row(0, named=True)
    assert abs(piao_lr["observed_fdr_median"] - 0.03) < 1e-9
    # q1/q3 in linear interpolation.
    assert abs(piao_lr["observed_fdr_q1"] - 0.02) < 1e-9
    assert abs(piao_lr["observed_fdr_q3"] - 0.04) < 1e-9
    # CI envelope spans the per-row Wilson lows/highs.
    assert piao_lr["observed_fdr_ci_lo"] == pytest.approx(0.01 / 2)
    assert piao_lr["observed_fdr_ci_hi"] == pytest.approx(min(1.0, 0.05 * 1.5))
    assert piao_lr["n_sites_total"] == 1000
    assert piao_lr["n_sites_called_median"] == 30.0


def test_aggregate_summary_skips_missing_parquets(tmp_path, monkeypatch):
    """If a per-engine parquet is missing, that row is dropped (not crashed)."""
    import run_phase4_null_calibration as mod

    monkeypatch.setattr(mod, "PIAO_OUT", tmp_path / "piao")
    monkeypatch.setattr(mod, "SIM_OUT", tmp_path / "sim")
    monkeypatch.setattr(mod, "GSE_OUT", tmp_path / "gse")

    specs = mod.make_specs()
    specs["piao"].out_dir = mod.PIAO_OUT
    specs["simulator"].out_dir = mod.SIM_OUT
    specs["gse"].out_dir = mod.GSE_OUT
    specs["piao"].engines = ("lr",)
    specs["simulator"].engines = ("lr",)
    specs["gse"].engines = ("lr",)

    # Only piao/lr has a parquet.
    _write_fake_engine_parquet(
        mod.PIAO_OUT / "lr.parquet",
        engine="lr", scenario="cov10_3v3", k=2,
        fdrs=[0.01, 0.02], n_called=[10, 20], n_total=100,
    )
    summary = mod.aggregate_summary(specs)
    assert summary.height == 1
    assert summary["dataset"][0] == "piao_distributed"


def test_make_specs_matches_task7_plan():
    """Sanity: spec registry matches the plan's dataset/engine/k matrix."""
    import run_phase4_null_calibration as mod

    specs = mod.make_specs()
    assert set(specs) == {"piao", "simulator", "gse"}
    assert specs["piao"].scenario == "cov10_3v3"
    assert specs["simulator"].scenario == "sim_cov10_3v3"
    assert specs["gse"].scenario == "gse263850"
    assert specs["piao"].k_shuffles == 20
    assert specs["simulator"].k_shuffles == 20
    assert specs["gse"].k_shuffles == 10
    assert specs["piao"].engines == ("lr", "lr_plus", "welch_t", "fisher")
    assert specs["simulator"].engines == ("lr", "lr_plus", "welch_t", "fisher")
    assert specs["gse"].engines == ("lr", "lr_plus", "welch_t", "fisher", "glm")


def test_engine_parquet_path_layout(tmp_path):
    """Per-engine parquets land at <out_dir>/<engine>.parquet."""
    import run_phase4_null_calibration as mod
    out = tmp_path / "x"
    assert mod._engine_parquet(out, "lr") == out / "lr.parquet"
    assert mod._engine_parquet(out, "lr_plus") == out / "lr_plus.parquet"


def test_write_manifest_smoke(tmp_path, monkeypatch):
    """Manifest writer produces a non-empty file with expected sections."""
    import run_phase4_null_calibration as mod

    monkeypatch.setattr(mod, "NULL_BASE", tmp_path / "nc")
    monkeypatch.setattr(mod, "MANIFEST_PATH", tmp_path / "nc" / "MANIFEST.txt")

    specs = mod.make_specs()
    timings = [
        {"dataset": "piao_distributed", "engine": "lr",
         "scenario": "cov10_3v3", "k_shuffles": 20,
         "wall_s": 12.3, "ok": True, "error": None},
        {"dataset": "gse263850", "engine": "glm",
         "scenario": "gse263850", "k_shuffles": 10,
         "wall_s": 99.9, "ok": False, "error": "boom"},
    ]
    summary = pl.DataFrame([
        {"engine": "lr", "dataset": "piao_distributed",
         "scenario": "cov10_3v3", "k_shuffles": 20,
         "observed_fdr_median": 0.04, "observed_fdr_q1": 0.03,
         "observed_fdr_q3": 0.05,
         "observed_fdr_ci_lo": 0.01, "observed_fdr_ci_hi": 0.07,
         "n_sites_called_median": 40.0, "n_sites_total": 1000},
    ])
    mod.write_manifest(specs, timings, summary, q_thresh=0.05)
    text = mod.MANIFEST_PATH.read_text(encoding="utf-8")
    assert "Phase 4 Task 7" in text
    assert "epykit version" in text
    assert "Datasets:" in text
    assert "Summary table" in text
    assert "median_FDR=0.0400" in text
    assert "FAIL" in text
    assert "boom" in text
