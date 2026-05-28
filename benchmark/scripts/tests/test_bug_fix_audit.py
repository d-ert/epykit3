"""bug_fix_audit.py: pre/post per-cell delta with commit-message
Affects:-trailer attribution."""
from __future__ import annotations

import json
import logging
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

SCRIPT = Path(__file__).parents[1] / "bug_fix_audit.py"

# Import the audit module directly for in-process unit tests that exercise
# join-key resolution / NaN-skip / coalesce branches without the subprocess
# round-trip used by the smoke tests above.
sys.path.insert(0, str(SCRIPT.parent))
import bug_fix_audit as bfa  # noqa: E402


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


def test_audit_handles_extra_join_keys():
    """When `threshold` is in both frames, it joins on it so the outer join
    stays 1:1 per (tool, scenario, threshold) cell instead of cross-producting.

    Without threshold in the join key, two pre rows x two post rows per
    (tool, scenario) would explode to 4 rows. With threshold in the join,
    the audit emits exactly 2 changed rows (one per real cell).
    """
    pre = pl.DataFrame({
        "tool":      ["epykit_lr", "epykit_lr"],
        "scenario":  ["cov10_3v3", "cov10_3v3"],
        "threshold": [0.1, 0.2],
        "tpr":       [0.50, 0.60],
        "fpr":       [0.05, 0.05],
    })
    post = pl.DataFrame({
        "tool":      ["epykit_lr", "epykit_lr"],
        "scenario":  ["cov10_3v3", "cov10_3v3"],
        "threshold": [0.1, 0.2],
        "tpr":       [0.55, 0.65],
        "fpr":       [0.05, 0.05],   # unchanged - should be dropped
    })
    attribution = {("lr", "cov10_3v3"): "P1-1"}
    audit_df, n_unattr = bfa.audit(pre, post, attribution)
    # 2 real cells, both changed in tpr only; fpr unchanged -> dropped.
    assert audit_df.height == 2, (
        f"expected 2 rows (one per threshold), got {audit_df.height}; "
        f"a cross-product would give 4. rows={audit_df.to_dicts()}"
    )
    assert n_unattr == 0
    assert set(audit_df["metric"].to_list()) == {"tpr"}
    assert (audit_df["fix_id"] == "P1-1").all()


def test_audit_skips_nan_pre():
    """A NaN pre value can't be diffed against a finite post value -- the
    audit silently drops that (tool, scenario, metric) row. This documents
    the current intentional behaviour: comparing NaN -> X has no defined
    delta sign, so the cell is excluded rather than emitted as UNATTRIBUTED.
    """
    pre = pl.DataFrame({
        "tool":     ["epykit_lr", "epykit_lr"],
        "scenario": ["cov10_3v3", "cov15_3v3"],
        "tpr":      [float("nan"), 0.80],
        "fpr":      [0.05, 0.04],
    })
    post = pl.DataFrame({
        "tool":     ["epykit_lr", "epykit_lr"],
        "scenario": ["cov10_3v3", "cov15_3v3"],
        "tpr":      [0.90, 0.80],   # cov10: NaN -> 0.90 should be skipped; cov15 unchanged
        "fpr":      [0.05, 0.04],   # both unchanged
    })
    attribution = {("lr", "cov10_3v3"): "P1-1", ("lr", "cov15_3v3"): "P1-1"}
    audit_df, n_unattr = bfa.audit(pre, post, attribution)
    # The cov10 tpr row is NaN-skipped; cov15 + all fpr rows are unchanged
    # and so they're dropped by the |delta| < 1e-9 filter -> empty audit.
    rows = audit_df.to_dicts()
    matching = [
        r for r in rows
        if r["tool"] == "epykit_lr" and r["scenario"] == "cov10_3v3" and r["metric"] == "tpr"
    ]
    assert matching == [], (
        f"NaN pre value should be silently skipped, but got: {matching}"
    )
    assert audit_df.height == 0, f"expected no rows, got {rows}"
    assert n_unattr == 0


def test_audit_coalesces_post_only_rows():
    """A (tool, scenario) that exists only in `post` (e.g. methylkit_tuned
    introduced in Phase 4) comes back through the outer join with the join
    keys under `tool_post` / `scenario_post`. The audit coalesces and
    attributes the row to the right fix_id.
    """
    pre = pl.DataFrame({
        "tool":     ["epykit_lr"],
        "scenario": ["cov10_3v3"],
        "tpr":      [0.80],
        "fpr":      [0.05],
    })
    post = pl.DataFrame({
        "tool":     ["epykit_lr", "methylkit_tuned"],
        "scenario": ["cov10_3v3", "cov10_3v3"],
        "tpr":      [0.85, 0.70],
        "fpr":      [0.05, 0.06],
    })
    attribution = {
        ("lr", "cov10_3v3"):              "P1-1",
        ("methylkit_tuned", "cov10_3v3"): "P0-9",
    }
    audit_df, n_unattr = bfa.audit(pre, post, attribution)
    # epykit_lr tpr: 0.80 -> 0.85 (changed); fpr unchanged.
    # methylkit_tuned: pre rows missing -> pre_v is None -> skipped (None branch).
    # So the only emitted row should be epykit_lr / tpr.
    rows = audit_df.to_dicts()
    assert len(rows) == 1, f"expected only epykit_lr tpr row, got {rows}"
    assert rows[0]["tool"] == "epykit_lr"
    assert rows[0]["scenario"] == "cov10_3v3"
    assert rows[0]["metric"] == "tpr"
    assert rows[0]["fix_id"] == "P1-1"
    assert n_unattr == 0


def test_audit_warns_on_asymmetric_extra_key(caplog):
    """If `threshold` is in pre but not post, _resolve_join_keys silently
    drops it from the join. This is a footgun (cross-product on that axis),
    so the resolver emits a logger.warning that surfaces in CI logs.
    """
    pre = pl.DataFrame({
        "tool":      ["epykit_lr"],
        "scenario":  ["cov10_3v3"],
        "threshold": [0.1],
        "tpr":       [0.50],
    })
    post = pl.DataFrame({
        "tool":     ["epykit_lr"],
        "scenario": ["cov10_3v3"],
        "tpr":      [0.60],
    })
    with caplog.at_level(logging.WARNING, logger=bfa.logger.name):
        keys = bfa._resolve_join_keys(pre, post)
    assert "threshold" not in keys, (
        "asymmetric column must NOT enter the join key set"
    )
    warnings = [
        rec for rec in caplog.records
        if rec.levelno >= logging.WARNING and "threshold" in rec.getMessage()
    ]
    assert warnings, (
        "expected a WARNING about asymmetric `threshold` column, "
        f"got records: {[(r.levelname, r.getMessage()) for r in caplog.records]}"
    )
