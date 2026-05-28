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

# Columns that, when present in both pre and post frames, further sub-divide a
# (tool, scenario) cell into individual rows. Phase 4 eval_summary parquets
# carry these threshold/parameter axes; the Phase 3 fixture parquets do not.
# The audit auto-detects which of these are present in BOTH frames and adds
# them to the join key so the outer join doesn't collapse to a cross-product.
EXTRA_KEY_CANDIDATES = (
    "test",
    "parameter",
    "parameter_value",
    "threshold",
    "threshold_kind",
    "meth_diff_bin",
)


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


def _resolve_join_keys(pre_df: pl.DataFrame, post_df: pl.DataFrame) -> list[str]:
    """Pick join columns: always (tool, scenario), plus any EXTRA_KEY_CANDIDATES
    that are present in BOTH frames."""
    keys = ["tool", "scenario"]
    for c in EXTRA_KEY_CANDIDATES:
        if c in pre_df.columns and c in post_df.columns:
            keys.append(c)
    return keys


def audit(
    pre_df: pl.DataFrame,
    post_df: pl.DataFrame,
    attribution: dict[tuple[str, str], str],
    metrics: tuple[str, ...] = METRICS,
) -> tuple[pl.DataFrame, int]:
    """Return (audit_df, n_unattributed).

    Joins pre and post on (tool, scenario) plus any EXTRA_KEY_CANDIDATES
    columns present in both frames, so threshold/parameter axes in Phase 4
    eval_summary parquets don't collapse the outer join to a cross-product.
    """
    join_cols = _resolve_join_keys(pre_df, post_df)
    metric_cols = [m for m in metrics if m in post_df.columns]
    joined = pre_df.join(
        post_df.select(join_cols + metric_cols),
        on=join_cols, how="outer", suffix="_post",
    )
    rows: list[dict] = []
    n_unattributed = 0
    for r in joined.iter_rows(named=True):
        # Outer join with `coalesce=False` (Polars default for "outer" pre-1.0)
        # leaves right-side keys under e.g. "tool_post" when only the right
        # row exists. Coalesce manually so we always have a tool/scenario.
        tool = r.get("tool") or r.get("tool_post") or ""
        scen = r.get("scenario") or r.get("scenario_post") or ""
        for m in metrics:
            pre_v  = r.get(m)
            post_v = r.get(f"{m}_post")
            if pre_v is None or post_v is None:
                continue
            try:
                pre_f, post_f = float(pre_v), float(post_v)
            except (TypeError, ValueError):
                continue
            # NaN propagation: if either side is NaN, skip (can't compare).
            if pre_f != pre_f or post_f != post_f:
                continue
            delta = post_f - pre_f
            if abs(delta) < 1e-9:
                continue  # unchanged
            engine = _engine_from_tool(tool)
            fix_id = attribution.get((engine, scen), "UNATTRIBUTED")
            if fix_id == "UNATTRIBUTED":
                n_unattributed += 1
            rows.append({
                "tool": tool, "scenario": scen,
                "metric": m,
                "pre_value": pre_f, "post_value": post_f,
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
