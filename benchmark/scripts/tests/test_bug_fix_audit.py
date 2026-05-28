"""bug_fix_audit.py: pre/post per-cell delta with commit-message
Affects:-trailer attribution."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

SCRIPT = Path(__file__).parents[1] / "bug_fix_audit.py"


def _write_parquets(pre_path: Path, post_path: Path):
    pl.DataFrame({
        "tool":     ["epykit_lr", "epykit_lr"],
        "scenario": ["cov10_3v3", "cov15_3v3"],
        "tpr":      [0.90, 0.85],
        "fpr":      [0.05, 0.04],
    }).write_parquet(str(pre_path))
    pl.DataFrame({
        "tool":     ["epykit_lr", "epykit_lr"],
        "scenario": ["cov10_3v3", "cov15_3v3"],
        "tpr":      [0.92, 0.85],   # cov10 changed; cov15 unchanged
        "fpr":      [0.07, 0.04],
    }).write_parquet(str(post_path))


def _write_commits(commits_path: Path, commits: list[dict]):
    commits_path.write_text(json.dumps(commits))


def test_attribution_success_exits_zero(tmp_path):
    """Changed cells that are attributed to a fix -> exit 0."""
    pre = tmp_path / "pre.parquet"
    post = tmp_path / "post.parquet"
    commits = tmp_path / "commits.json"
    out = tmp_path / "audit.parquet"
    _write_parquets(pre, post)
    _write_commits(commits, [
        {"subject": "fix(dmc) P1-1: Fisher mid-p", "body": "Affects: lr@cov10_3v3"},
    ])
    result = subprocess.run(
        [sys.executable, str(SCRIPT),
         "--pre", str(pre), "--post", str(post),
         "--commits-json", str(commits), "--out", str(out)],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, (
        f"audit should exit 0 when all changed cells attributed:\n"
        f"stdout={result.stdout}\nstderr={result.stderr}"
    )
    audit = pl.read_parquet(str(out))
    changed = audit.filter((pl.col("tool") == "epykit_lr") & (pl.col("scenario") == "cov10_3v3"))
    assert changed.height >= 1
    assert (changed["fix_id"] == "P1-1").all()


def test_unattributed_cell_exits_nonzero(tmp_path):
    """Changed cells with no matching Affects: trailer -> exit non-zero."""
    pre = tmp_path / "pre.parquet"
    post = tmp_path / "post.parquet"
    commits = tmp_path / "commits.json"
    out = tmp_path / "audit.parquet"
    _write_parquets(pre, post)
    _write_commits(commits, [])  # No commits -> no attribution
    result = subprocess.run(
        [sys.executable, str(SCRIPT),
         "--pre", str(pre), "--post", str(post),
         "--commits-json", str(commits), "--out", str(out)],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode != 0, "audit should exit non-zero when cells are UNATTRIBUTED"
    assert "UNATTRIBUTED" in (result.stdout + result.stderr)
