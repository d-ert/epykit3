"""regen_all.py — acceptance gate for paper claims.

Modes:
  --verify       (default in CI): asserts every claim in claims.yaml
                 matches its cited parquet to the printed precision.
                 Exits non-zero on any mismatch.
  --run-cheap    Phase 4 expansion target. Exits with message.
  --run-all      Phase 4 expansion target. Exits with message.

claims.yaml schema (YAML list, one entry per claim):
  - claim_id:   stable identifier
    parquet:    absolute or repo-relative path to the source parquet
    column:     column name in the parquet to read
    filter:     mapping of column -> value (selects exactly one row)
    expected:   numeric expected value
    precision:  absolute tolerance (|actual - expected| <= precision to PASS)

Paper markup: place ``<!-- claim: <claim_id> -->`` adjacent to any
numeric value in paper.md. --verify scans for these comments and
checks that the cited claim_id exists in claims.yaml.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import polars as pl
import yaml

CLAIM_COMMENT_RE = re.compile(r"<!--\s*claim:\s*([A-Za-z0-9_\-]+)\s*-->")


def _load_claims(claims_path: Path) -> dict:
    text = claims_path.read_text(encoding="utf-8")
    data = yaml.safe_load(text)
    if not data:
        return {}
    return {c["claim_id"]: c for c in data}


def _parse_paper_claims(paper_path: Path) -> set[str]:
    if not paper_path.exists():
        return set()
    text = paper_path.read_text(encoding="utf-8")
    return set(CLAIM_COMMENT_RE.findall(text))


def _read_claim_value(claim: dict) -> float:
    parquet_path = Path(str(claim["parquet"]))
    if not parquet_path.is_absolute():
        # Resolve relative to repo root (two levels above benchmark/scripts/).
        repo = Path(__file__).resolve().parents[2]
        parquet_path = repo / parquet_path
    df = pl.read_parquet(str(parquet_path))
    for col, val in (claim.get("filter") or {}).items():
        df = df.filter(pl.col(col) == val)
    if df.height != 1:
        raise ValueError(
            f"claim '{claim['claim_id']}': filter selected {df.height} rows "
            f"(expected 1) from {parquet_path}"
        )
    return float(df[claim["column"]][0])


def verify(claims_path: Path, paper_path: Path) -> int:
    claims    = _load_claims(claims_path)
    referenced = _parse_paper_claims(paper_path)

    # Claims in the paper but not in claims.yaml → fail.
    missing = referenced - set(claims)
    if missing:
        for cid in sorted(missing):
            print(f"FAIL: paper references unknown claim '{cid}'")
        return 1

    if not referenced:
        print("OK: no claims referenced in paper (empty seed manifest)")
        return 0

    failures = 0
    for cid in sorted(referenced):
        claim = claims[cid]
        try:
            actual = _read_claim_value(claim)
        except Exception as exc:
            print(f"FAIL: {cid}: could not read value — {exc}")
            failures += 1
            continue
        expected  = float(claim["expected"])
        precision = float(claim.get("precision", 0.0))
        if abs(actual - expected) > precision:
            print(
                f"FAIL: {cid}: parquet={actual:.6f} expected={expected:.6f} "
                f"diff={actual - expected:+.6f} tolerance={precision}"
            )
            failures += 1
        else:
            print(f"OK:   {cid}: {actual:.6f} (expected {expected:.6f} ± {precision})")
    return 0 if failures == 0 else 1


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true",
                        help="Verify all claims against parquets (default mode).")
    parser.add_argument("--run-cheap", action="store_true")
    parser.add_argument("--run-all",   action="store_true")
    parser.add_argument("--claims", default="benchmark/scripts/claims.yaml", type=Path)
    parser.add_argument("--paper",  default="benchmark/paper/paper.md",    type=Path)
    args = parser.parse_args(argv)

    if args.run_cheap or args.run_all:
        print("--run-cheap / --run-all are Phase 4 expansion targets; "
              "not implemented in 0.7.5.")
        sys.exit(2)

    # --verify is the default action.
    sys.exit(verify(args.claims, args.paper))


if __name__ == "__main__":
    main()
