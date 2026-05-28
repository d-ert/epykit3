"""bug_fix_audit.py -- pre/post-fix per-cell delta with commit attribution.

Diffs pre-fix and post-fix eval_summary.parquet per (tool, scenario,
metric), attributes each delta to a P0/P1 fix via ``Affects: engine@scenario``
trailers parsed from commit messages, and exits non-zero if any changed
cell is unattributed.

Usage:
    git log --format='{"subject":"%s","body":"%b"}%n---' v0.7.2..v0.7.5 | \\
        python commits_to_json.py > commits.json
    python bug_fix_audit.py \\
        --pre  benchmark/data/study1/eval_summary.parquet \\
        --post benchmark/data/study1/eval_summary_post.parquet \\
        --commits-json commits.json \\
        --out  benchmark/data/audit/bug_fix_deltas.parquet

For Phase 3, test with fixture parquets and a hand-written commits.json.
Phase 4 runs this with the real pre-Phase-1 baseline vs post-Phase-3 outputs.

``Affects: engine@scenario`` trailer format:
    One or more ``engine@scenario`` tokens, comma-separated.
    ``engine`` matches the ``tool`` column after stripping ``epykit_`` prefix.
    ``scenario`` matches the ``scenario`` column literally.
    Example: ``Affects: lr@cov10_3v3, glm@cov10_3v3``
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import polars as pl

AFFECTS_RE  = re.compile(r"^Affects:\s*(.+)$", re.MULTILINE)
CELL_RE     = re.compile(r"([A-Za-z0-9_+\-]+)@([A-Za-z0-9_]+)")
SUBJECT_PID = re.compile(r"\b(P[01]-\d+[a-z]?)\b")


METRICS = ("tpr", "fpr", "f1", "auroc", "precision")


def _parse_commits(commits: list[dict]) -> dict[tuple[str, str], str]:
    """Map (engine, scenario) -> most-recent fix_id with an Affects: trailer."""
    attribution: dict[tuple[str, str], str] = {}
    for c in commits:
        body    = c.get("body", "")
        subject = c.get("subject", "")
        pid_m   = SUBJECT_PID.search(subject)
        fix_id  = pid_m.group(1) if pid_m else subject[:40]
        for m in AFFECTS_RE.finditer(body):
            for cm in CELL_RE.finditer(m.group(1)):
                key = (cm.group(1), cm.group(2))
                attribution[key] = fix_id  # later commit wins
    return attribution


def _engine_from_tool(tool: str) -> str:
    """Strip 'epykit_' prefix. 'epykit_lr' -> 'lr'."""
    return tool[len("epykit_"):] if tool.startswith("epykit_") else tool


def audit(
    pre_df: pl.DataFrame,
    post_df: pl.DataFrame,
    attribution: dict[tuple[str, str], str],
    metrics: tuple[str, ...] = METRICS,
) -> tuple[pl.DataFrame, int]:
    """Return (audit_df, n_unattributed)."""
    join_cols = ["tool", "scenario"]
    joined = pre_df.join(
        post_df.select(join_cols + [m for m in metrics if m in post_df.columns]),
        on=join_cols, how="outer", suffix="_post",
    )
    rows: list[dict] = []
    n_unattributed = 0
    for r in joined.iter_rows(named=True):
        for m in metrics:
            pre_v  = r.get(m)
            post_v = r.get(f"{m}_post")
            if pre_v is None or post_v is None:
                continue
            try:
                delta = float(post_v) - float(pre_v)
            except (TypeError, ValueError):
                continue
            if abs(delta) < 1e-9:
                continue  # unchanged
            engine = _engine_from_tool(r.get("tool", ""))
            scen   = r.get("scenario", "")
            fix_id = attribution.get((engine, scen), "UNATTRIBUTED")
            if fix_id == "UNATTRIBUTED":
                n_unattributed += 1
            rows.append({
                "tool": r.get("tool", ""), "scenario": scen,
                "metric": m,
                "pre_value": float(pre_v), "post_value": float(post_v),
                "delta": delta, "fix_id": fix_id,
            })
    if not rows:
        schema = {
            "tool": pl.Utf8, "scenario": pl.Utf8, "metric": pl.Utf8,
            "pre_value": pl.Float64, "post_value": pl.Float64,
            "delta": pl.Float64, "fix_id": pl.Utf8,
        }
        return pl.DataFrame(schema=schema), n_unattributed
    return pl.DataFrame(rows), n_unattributed


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pre",          required=True, type=Path)
    parser.add_argument("--post",         required=True, type=Path)
    parser.add_argument("--commits-json", required=True, type=Path)
    parser.add_argument("--out",          required=True, type=Path)
    args = parser.parse_args(argv)

    pre_df   = pl.read_parquet(str(args.pre))
    post_df  = pl.read_parquet(str(args.post))
    commits  = json.loads(args.commits_json.read_text())
    attr     = _parse_commits(commits)

    audit_df, n_unattr = audit(pre_df, post_df, attr)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    audit_df.write_parquet(str(args.out))
    print(f"wrote {args.out} ({audit_df.height} rows; {n_unattr} UNATTRIBUTED)")
    if n_unattr > 0:
        print(
            f"FAIL: {n_unattr} changed cells have no Affects: attribution.\n"
            "Add an 'Affects: engine@scenario' trailer to the responsible "
            "commit or document the delta in Limitations S10.5."
        )
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
