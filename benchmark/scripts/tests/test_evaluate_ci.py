"""evaluate.py --ci-only appends Wilson + bootstrap CI columns."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest

SCRIPT = Path(__file__).parents[1] / "evaluate.py"


def _make_synth_eval(path: Path):
    """Write a minimal eval_summary.parquet for testing."""
    pl.DataFrame({
        "tool":      ["epykit_lr"] * 5,
        "scenario":  [f"cov10_n3_seed{i}" for i in range(5)],
        "test":      ["lr"] * 5,
        "tp":        [90, 85, 70, 60, 95],
        "fp":        [10,  5, 20, 30,  2],
        "tn":        [890, 905, 880, 870, 898],
        "fn":        [10,  5, 30, 40,  5],
        "tpr":       [0.9,  0.944, 0.7,  0.6,  0.95],
        "fpr":       [0.011, 0.0055, 0.022, 0.033, 0.0022],
        "f1":        [0.9,  0.94,  0.74, 0.63, 0.96],
        "auroc":     [0.97, 0.95,  0.88, 0.80, 0.99],
    }).write_parquet(str(path))


def test_evaluate_ci_adds_required_columns(tmp_path):
    """--ci-only must add tpr_ci_lo, tpr_ci_hi, fpr_ci_lo, fpr_ci_hi."""
    eval_in = tmp_path / "eval_summary.parquet"
    _make_synth_eval(eval_in)
    result = subprocess.run(
        [sys.executable, str(SCRIPT),
         "--ci-only", "--eval-summary", str(eval_in)],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, (
        f"evaluate.py --ci-only failed:\n{result.stderr}"
    )
    out = pl.read_parquet(str(eval_in))
    for col in ("tpr_ci_lo", "tpr_ci_hi", "fpr_ci_lo", "fpr_ci_hi"):
        assert col in out.columns, f"Missing column {col}; got {out.columns}"


def test_evaluate_ci_brackets_point_estimate(tmp_path):
    """CI bounds must bracket the point estimate for all rows."""
    eval_in = tmp_path / "eval_summary.parquet"
    _make_synth_eval(eval_in)
    subprocess.run(
        [sys.executable, str(SCRIPT),
         "--ci-only", "--eval-summary", str(eval_in)],
        check=True, capture_output=True, text=True,
    )
    out = pl.read_parquet(str(eval_in))
    lo = out["tpr_ci_lo"].to_numpy()
    hi = out["tpr_ci_hi"].to_numpy()
    p  = out["tpr"].to_numpy()
    valid = np.isfinite(lo) & np.isfinite(hi) & np.isfinite(p)
    assert valid.all(), "All CI values should be finite"
    assert ((lo[valid] <= p[valid]) & (p[valid] <= hi[valid])).all(), (
        "CI must bracket the point estimate"
    )
