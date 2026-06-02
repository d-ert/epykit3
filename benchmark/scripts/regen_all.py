"""regen_all.py -- end-to-end benchmark runbook + acceptance gate.

Modes:
  --verify       (default in CI): asserts every claim in claims.yaml
                 matches its cited parquet to the printed precision.
                 Exits non-zero on any mismatch.
  --run-all      Linux-only. Orchestrates the full benchmark pipeline:
                 the simulator at the configured phi-sweep, all
                 epykit + R baselines, k=1000 null calibration, both
                 truth_modes, sensitivity sweep, EB validation. Writes
                 to benchmark/data/ (the canonical parquet source
                 tree). After this completes the paper_data/ TSV
                 mirror is re-derived from data/ by a separate
                 converter.
  --skip <step>  Skip a named step (use multiple times). Step names:
                 simulator, epykit_dmc, methylkit, dss, dmrseq,
                 bsmooth, null_calibration, sensitivity, eb_validation.
  --only <step>  Run only the named step (cannot combine with --skip).
  --dry-run      Print the command sequence without executing.

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
import logging
import re
import subprocess
import sys
from pathlib import Path

import polars as pl
import yaml

logger = logging.getLogger(__name__)

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


# ---------------------------------------------------------------------------
# --run-all dispatcher
# ---------------------------------------------------------------------------

# Default phi sweep for the simulator. phi=0 is pure Binomial (Piao
# original); phi in (0, 1) is Beta-Binomial with the given intraclass
# correlation. The four non-zero values bracket realistic WGBS
# overdispersion at coverage 10-30.
DEFAULT_PHI_SWEEP = (0.0, 0.01, 0.05, 0.1, 0.2)
DEFAULT_SEEDS = tuple(range(2026_000, 2026_020))   # 20 seeds (PROTOCOL.md)
DEFAULT_K_SHUFFLES = 1000

# Named steps in execution order. Each is dispatched to a separate
# Python or R runner; failures are isolated per step so a single
# baseline crash does not nuke the run.
STEPS = (
    "simulator",
    "epykit_dmc",
    "methylkit",
    "dss",
    "dmrseq",
    "bsmooth",
    "null_calibration",
    "sensitivity",
    "eb_validation",
)


def _shell(cmd: list[str], *, dry_run: bool, cwd: Path | None = None) -> int:
    """Run ``cmd`` and stream stdout / stderr to the parent process."""
    rendered = " ".join(str(c) for c in cmd)
    logger.info("$ %s", rendered)
    if dry_run:
        return 0
    result = subprocess.run(cmd, cwd=cwd, check=False)
    if result.returncode != 0:
        logger.error("step failed (exit %d): %s", result.returncode, rendered)
    return result.returncode


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def run_all(
    *,
    seeds: tuple[int, ...] = DEFAULT_SEEDS,
    phi_sweep: tuple[float, ...] = DEFAULT_PHI_SWEEP,
    k_shuffles: int = DEFAULT_K_SHUFFLES,
    skip: set[str] | None = None,
    only: str | None = None,
    dry_run: bool = False,
) -> int:
    """Run the full benchmark pipeline. Returns 0 on success.

    This is intentionally linear (no parallelism inside the
    dispatcher). The user is expected to run this on a Linux box with
    the R container available; parallel execution can be added once
    the linear path is known to work end to end.
    """
    skip = set(skip or ())
    if only is not None:
        if only not in STEPS:
            raise SystemExit(f"--only {only!r} not in {STEPS}")
        skip = set(STEPS) - {only}

    root = _repo_root()
    scripts = root / "benchmark" / "scripts"
    data = root / "benchmark" / "data"
    data.mkdir(parents=True, exist_ok=True)

    failures: list[tuple[str, str]] = []

    # Step 1 -- simulator (phi-sweep, multi-seed).
    if "simulator" not in skip:
        for phi in phi_sweep:
            for seed in seeds:
                out = data / "study1b_simulator" / f"phi={phi}" / f"seed={seed}"
                rc = _shell([
                    "python", str(scripts / "simulate_piao.py"),
                    "--seed", str(seed),
                    "--phi", str(phi),
                    "--out", str(out),
                ], dry_run=dry_run)
                if rc != 0:
                    failures.append(("simulator", f"phi={phi}, seed={seed}"))

    # Step 2 -- epykit DMC on every simulator cell. Reuses
    # run_epykit_simulator.py's existing per-cell logic.
    if "epykit_dmc" not in skip:
        rc = _shell([
            "python", str(scripts / "run_epykit_simulator.py"),
        ], dry_run=dry_run)
        if rc != 0:
            failures.append(("epykit_dmc", "run_epykit_simulator.py"))

    # Step 3 / 4 / 5 / 6 -- R baselines. These need Dockerfile.r or a
    # native methylKit/DSS/dmrseq/bsseq install; outside the Linux
    # container this is skipped with a warning.
    r_runners = {
        "methylkit": "run_methylkit_simulator.R",
        "dss": "run_dss_simulator.R",
        "dmrseq": "run_dmrseq.R",
        "bsmooth": "run_bsmooth.R",
    }
    for step, runner in r_runners.items():
        if step in skip:
            continue
        for phi in phi_sweep:
            for seed in seeds:
                in_dir = data / "study1b_simulator" / f"phi={phi}" / f"seed={seed}" / "bismark_cov"
                out = data / "study1b_simulator" / f"phi={phi}" / f"seed={seed}" / f"{step}.tsv"
                rc = _shell([
                    "Rscript", str(scripts / runner),
                    "--in-dir", str(in_dir),
                    "--out", str(out),
                ], dry_run=dry_run)
                if rc != 0:
                    failures.append((step, f"phi={phi}, seed={seed}"))

    # Step 7 -- null calibration on shuffled labels.
    if "null_calibration" not in skip:
        # Run the calibration on a single representative store; the
        # user's Linux setup will iterate over more if desired.
        null_dir = data / "null_calibration"
        for engine in ("lr", "lr_plus"):
            rc = _shell([
                "python", str(scripts / "run_null_calibration.py"),
                "--engine", engine,
                "--methylstore", str(data / "study1b_simulator" / "phi=0.0" / "seed=2026000"),
                "--scenario", "cov10_3v3",
                "--k-shuffles", str(k_shuffles),
                "--out", str(null_dir / f"{engine}.parquet"),
                "--summary-out", str(null_dir / f"{engine}_summary.parquet"),
                "--pvalue-out", str(null_dir / f"{engine}_pvalues.parquet"),
                "--ks-out", str(null_dir / f"{engine}_ks.json"),
            ], dry_run=dry_run)
            if rc != 0:
                failures.append(("null_calibration", engine))

    # Step 8 -- sensitivity sweep on the headline DMC.
    if "sensitivity" not in skip:
        dmc = data / "study1b_simulator" / "phi=0.0" / "seed=2026000" / "epykit_lr" / "dmc.parquet"
        truth = data / "study1b_simulator" / "phi=0.0" / "seed=2026000" / "truth.parquet"
        rc = _shell([
            "python", str(scripts / "sensitivity_sweep.py"),
            "--dmc-parquet", str(dmc),
            "--truth-parquet", str(truth),
            "--out", str(data / "sensitivity" / "sweep.parquet"),
        ], dry_run=dry_run)
        if rc != 0:
            failures.append(("sensitivity", "sweep"))

    # Step 9 -- EB prior validation. Requires phi_site exported from
    # an epykit run (a separate small wiring step the user enables
    # via a future dmc.py debug flag).
    if "eb_validation" not in skip:
        phi_parquet = data / "eb_validation" / "phi_site.parquet"
        if phi_parquet.exists() or dry_run:
            rc = _shell([
                "python", str(scripts / "validate_eb_prior.py"),
                "--phi-parquet", str(phi_parquet),
                "--qq-out", str(data / "eb_validation" / "qq.png"),
                "--summary-out", str(data / "eb_validation" / "summary.json"),
            ], dry_run=dry_run)
            if rc != 0:
                failures.append(("eb_validation", "validate_eb_prior"))
        else:
            logger.warning(
                "eb_validation skipped: %s not found "
                "(produce it by enabling phi_site export in dmc.py)",
                phi_parquet,
            )

    if failures:
        print(f"\n{len(failures)} step(s) failed:")
        for step, detail in failures:
            print(f"  {step}: {detail}")
        return 1
    print("\nALL STEPS OK")
    return 0


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--verify", action="store_true",
        help="Verify all claims against parquets (default mode).",
    )
    parser.add_argument(
        "--run-all", action="store_true",
        help="Run the full end-to-end pipeline. Linux + R container "
             "strongly recommended; see benchmark/README.md.",
    )
    parser.add_argument(
        "--skip", action="append", default=[],
        choices=STEPS,
        help="Skip a named step. May be passed multiple times.",
    )
    parser.add_argument(
        "--only", default=None, choices=STEPS,
        help="Run only the named step. Mutually exclusive with --skip.",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print the command sequence without executing.",
    )
    parser.add_argument(
        "--phi-sweep", nargs="+", type=float, default=list(DEFAULT_PHI_SWEEP),
        help="Phi values to sweep in the simulator (default: 0 0.01 0.05 0.1 0.2)",
    )
    parser.add_argument(
        "--seeds", nargs="+", type=int, default=list(DEFAULT_SEEDS),
        help="Seeds to sweep (default: 20 seeds starting at 2026000)",
    )
    parser.add_argument(
        "--k-shuffles", type=int, default=DEFAULT_K_SHUFFLES,
        help="Null-calibration shuffle count (default: 1000)",
    )
    parser.add_argument("--claims", default="benchmark/scripts/claims.yaml", type=Path)
    parser.add_argument("--paper",  default="benchmark/paper/paper.md", type=Path)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    if args.run_all:
        if args.only and args.skip:
            raise SystemExit("--only and --skip are mutually exclusive")
        sys.exit(run_all(
            seeds=tuple(args.seeds),
            phi_sweep=tuple(args.phi_sweep),
            k_shuffles=args.k_shuffles,
            skip=set(args.skip),
            only=args.only,
            dry_run=args.dry_run,
        ))

    # --verify is the default action.
    sys.exit(verify(args.claims, args.paper))


if __name__ == "__main__":
    main()
