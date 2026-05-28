"""extract_affects_commits.py -- safe git-log -> commits.json extractor.

The Phase 4 audit (`bug_fix_audit.py`) consumes a JSON array of
``{"subject": ..., "body": ...}`` objects. Building this JSON via inline
``git log --format='{"subject": %s ...}'`` breaks on multi-line bodies, embedded
quotes, and backslashes. This helper uses ASCII US (``\\x1f``) as the field
separator and ASCII RS (``\\x1e``) as the record separator, neither of which
appears in commit messages, then constructs dicts in Python so ``json.dumps``
handles escaping correctly.

Usage:
    python extract_affects_commits.py <git-revision-range> [--out path]

Example (Phase 4):
    python benchmark/scripts/extract_affects_commits.py \\
        main..v0.7.5-phase3-engines-frozen \\
        --out benchmark/data/audit/commits.json

Why this range for Phase 4: ``v0.7.2`` referenced in the original plan does not
exist as a tag (the tag set is ``v0.7.3-p0-complete``, ``v0.7.4-phase2-scripts``,
``v0.7.5-phase3-engines-frozen``). The merge-base of ``p0-fixes`` with ``main``
is the natural pre-Phase-1 boundary, and ``v0.7.5-phase3-engines-frozen`` is the
post-Phase-3 freeze point that produced ``eval_summary_post_phase3.parquet``.
"""
from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

logger = logging.getLogger(__name__)

FIELD_SEP = "\x1f"  # ASCII Unit Separator
REC_SEP = "\x1e"    # ASCII Record Separator


def extract(revision_range: str, repo_dir: Path | None = None) -> list[dict]:
    """Run ``git log`` over ``revision_range`` and return a list of
    ``{"hash": str, "subject": str, "body": str}`` dicts.

    Uses ASCII US/RS separators that cannot appear in commit messages, so
    multi-line bodies with quotes/backslashes round-trip safely.
    """
    fmt = f"%H{FIELD_SEP}%s{FIELD_SEP}%b{REC_SEP}"
    cmd = ["git", "log", f"--format={fmt}", revision_range]
    result = subprocess.run(
        cmd, cwd=str(repo_dir) if repo_dir else None,
        capture_output=True, text=True, encoding="utf-8", check=True,
    )
    out: list[dict] = []
    for record in result.stdout.split(REC_SEP):
        record = record.strip("\n\r")
        if not record:
            continue
        parts = record.split(FIELD_SEP)
        if len(parts) < 3:
            logger.warning("skipping malformed record: %r", record[:80])
            continue
        commit_hash, subject, body = parts[0], parts[1], FIELD_SEP.join(parts[2:])
        out.append({
            "hash": commit_hash.strip(),
            "subject": subject.strip(),
            "body": body.strip("\n"),
        })
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "revision_range",
        help="git revision range, e.g. 'main..v0.7.5-phase3-engines-frozen'",
    )
    parser.add_argument(
        "--out", type=Path, default=None,
        help="output JSON path; default writes to stdout",
    )
    parser.add_argument(
        "--repo", type=Path, default=None,
        help="git repository directory (default: current working directory)",
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    commits = extract(args.revision_range, repo_dir=args.repo)
    payload = json.dumps(commits, ensure_ascii=False, indent=2)
    if args.out is None:
        print(payload)
    else:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(payload, encoding="utf-8")
        print(f"wrote {args.out} ({len(commits)} commits)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
