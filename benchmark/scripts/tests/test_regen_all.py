"""regen_all.py --verify reads claims.yaml and asserts each claim
matches its source parquet to the printed precision."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

SCRIPT = Path(__file__).parents[1] / "regen_all.py"


def _make_parquet(path: Path, auroc: float, tpr: float):
    pl.DataFrame({
        "tool":     ["epykit_lr"],
        "scenario": ["cov10_3v3"],
        "auroc":    [auroc],
        "tpr":      [tpr],
    }).write_parquet(str(path))


def test_verify_passes_on_matching_claim(tmp_path):
    parquet = tmp_path / "vals.parquet"
    _make_parquet(parquet, auroc=0.987, tpr=0.95)
    claims = tmp_path / "claims.yaml"
    claims.write_text(
        f"- claim_id: study1_auroc\n"
        f"  parquet: {parquet}\n"
        f"  column: auroc\n"
        f"  filter:\n"
        f"    tool: epykit_lr\n"
        f"    scenario: cov10_3v3\n"
        f"  expected: 0.987\n"
        f"  precision: 0.001\n"
    )
    paper = tmp_path / "paper.md"
    paper.write_text(
        "AUROC was 0.987 <!-- claim: study1_auroc --> in Study 1.\n"
    )
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--verify",
         "--claims", str(claims), "--paper", str(paper)],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, (
        f"verify should pass on matching claim; got:\n{result.stdout}\n{result.stderr}"
    )


def test_verify_fails_on_off_by_precision(tmp_path):
    parquet = tmp_path / "vals.parquet"
    _make_parquet(parquet, auroc=0.987, tpr=0.95)
    claims = tmp_path / "claims.yaml"
    # Expected is 0.887 but parquet has 0.987 — off by 0.1 > precision 0.001.
    claims.write_text(
        f"- claim_id: study1_auroc\n"
        f"  parquet: {parquet}\n"
        f"  column: auroc\n"
        f"  filter:\n"
        f"    tool: epykit_lr\n"
        f"    scenario: cov10_3v3\n"
        f"  expected: 0.887\n"
        f"  precision: 0.001\n"
    )
    paper = tmp_path / "paper.md"
    paper.write_text("AUROC 0.887 <!-- claim: study1_auroc -->\n")
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--verify",
         "--claims", str(claims), "--paper", str(paper)],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode != 0, "verify should fail on precision mismatch"
    assert "study1_auroc" in (result.stdout + result.stderr)
