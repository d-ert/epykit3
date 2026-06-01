"""run_external_simulator_sweep.py -- Phase 4 Task 5 extension.

Runs methylKit + DSS (smoothing=TRUE and FALSE) on every available
simulator seed, with skip-if-exists logic so this is resumable across
interruptions.

For each seed in ``benchmark/data/study1b_simulator/seed=*/``:
  1. Build ``bismark_cov/*.cov.gz`` from the AMP files (if missing).
  2. Run methylKit -> ``methylkit.tsv`` (if missing).
  3. Run DSS smoothing=TRUE -> ``dss.tsv`` (if missing).
  4. Run DSS smoothing=FALSE -> ``dss_nosmooth.tsv`` (if missing).

Per-seed walltime on chr1 / 100k CpGs is ~10 min total (methylKit
~6.5 min + DSS smoothed ~1.2 min + DSS unsmoothed ~1.1 min + R startup
overhead). 20 seeds × ~10 min = ~3 h end-to-end if run sequentially.

After this script finishes, run::

    uv run python benchmark/scripts/eval_simulator_intrinsic.py --all-seeds

to score everything against truth and produce the multi-seed
intrinsic-truth parallel column.

Usage::

    uv run python benchmark/scripts/run_external_simulator_sweep.py
        [--seeds 2026000 2026001 ...]
        [--coverage 10]
        [--skip-methylkit]
        [--skip-dss]
        [--skip-dss-nosmooth]
"""

from __future__ import annotations

import argparse
import gzip
import logging
import shutil
import subprocess
import sys
import time
from pathlib import Path

import polars as pl

logger = logging.getLogger("run_external_simulator_sweep")

ROOT = Path(__file__).resolve().parents[2]
SIM_ROOT = ROOT / "benchmark" / "data" / "study1b_simulator"
SCRIPTS_DIR = Path(__file__).resolve().parent

METHYLKIT_RUNNER = SCRIPTS_DIR / "run_methylkit_simulator.R"
DSS_RUNNER = SCRIPTS_DIR / "run_dss_simulator.R"


def convert_amp_to_cov(seed_dir: Path, coverage: int) -> int:
    """Build bismark .cov.gz files from AMP for one seed. Returns count."""
    dst = seed_dir / "bismark_cov"
    dst.mkdir(exist_ok=True)
    amp_files = sorted(seed_dir.glob(f"amp.coverage={coverage}.sample*.txt"))
    if not amp_files:
        raise FileNotFoundError(
            f"no AMP files at coverage={coverage} in {seed_dir}"
        )
    n_converted = 0
    for amp in amp_files:
        gz_path = dst / amp.name.replace(".txt", ".cov.gz")
        if gz_path.exists():
            continue
        df = pl.read_csv(amp, separator="\t")
        out = (
            df.with_columns(
                (pl.col("freqC") / 100.0 * pl.col("coverage"))
                .round().cast(pl.Int64).alias("count_M")
            )
            .with_columns((pl.col("coverage") - pl.col("count_M")).alias("count_U"))
            .select([
                pl.col("chr"),
                pl.col("base").alias("start"),
                pl.col("base").alias("end"),
                pl.col("freqC").alias("beta"),
                pl.col("count_M"),
                pl.col("count_U"),
            ])
        )
        raw = dst / amp.name.replace(".txt", ".cov")
        out.write_csv(raw, separator="\t", include_header=False)
        with open(raw, "rb") as fi, gzip.open(gz_path, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        raw.unlink()
        n_converted += 1
    return n_converted


def run_r_runner(
    runner: Path, in_dir: Path, out: Path,
    extra_args: list[str] | None = None,
) -> dict:
    """Call an R runner via Rscript; return a small status dict."""
    cmd = [
        "Rscript", str(runner),
        "--in-dir", str(in_dir),
        "--out", str(out),
    ]
    if extra_args:
        cmd.extend(extra_args)
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    if proc.returncode != 0:
        logger.error("%s failed (rc=%d): %s",
                     runner.name, proc.returncode, proc.stderr[-400:])
        return {"ok": False, "wall_s": wall, "error": proc.stderr[-400:]}
    return {"ok": True, "wall_s": wall, "stdout_tail": proc.stdout[-200:]}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--seeds", nargs="+", type=int, default=None,
        help="Seeds to process (default: every seed=NNN dir in study1b_simulator/)",
    )
    parser.add_argument(
        "--coverage", type=int, default=10,
        help="AMP coverage level to convert (default: 10)",
    )
    parser.add_argument("--skip-methylkit", action="store_true")
    parser.add_argument("--skip-dss", action="store_true")
    parser.add_argument("--skip-dss-nosmooth", action="store_true")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        force=True,
    )

    if args.seeds is None:
        seeds = sorted(
            int(d.name.split("=")[1])
            for d in SIM_ROOT.glob("seed=*")
            if d.is_dir() and d.name.split("=")[1].isdigit()
        )
    else:
        seeds = sorted(args.seeds)
    logger.info("processing %d seeds: %s", len(seeds), seeds)

    overall_t0 = time.perf_counter()
    per_seed_status: list[dict] = []

    for seed in seeds:
        seed_dir = SIM_ROOT / f"seed={seed}"
        if not seed_dir.exists():
            logger.warning("seed=%d: dir not found, skipping", seed)
            continue
        logger.info("=== seed=%d ===", seed)
        status = {"seed": seed}

        # 1. AMP -> .cov.gz
        t0 = time.perf_counter()
        n_new = convert_amp_to_cov(seed_dir, args.coverage)
        status["convert_s"] = time.perf_counter() - t0
        status["convert_n_new"] = n_new
        logger.info(
            "[seed=%d] convert: %d new .cov.gz files in %.1fs",
            seed, n_new, status["convert_s"],
        )

        in_dir = seed_dir / "bismark_cov"

        # 2. methylKit
        if not args.skip_methylkit:
            mk_out = seed_dir / "methylkit.tsv"
            if mk_out.exists():
                logger.info("[seed=%d] methylkit: %s exists, skipping",
                            seed, mk_out.name)
                status["methylkit"] = {"ok": True, "wall_s": 0.0, "skipped": True}
            else:
                logger.info("[seed=%d] methylkit: launching...", seed)
                status["methylkit"] = run_r_runner(METHYLKIT_RUNNER, in_dir, mk_out)
                logger.info("[seed=%d] methylkit: %s in %.1fs",
                            "ok" if status["methylkit"]["ok"] else "FAIL",
                            seed, status["methylkit"]["wall_s"])

        # 3. DSS smoothing=TRUE
        if not args.skip_dss:
            dss_out = seed_dir / "dss.tsv"
            if dss_out.exists():
                logger.info("[seed=%d] dss: %s exists, skipping",
                            seed, dss_out.name)
                status["dss"] = {"ok": True, "wall_s": 0.0, "skipped": True}
            else:
                logger.info("[seed=%d] dss (smoothing=TRUE): launching...", seed)
                status["dss"] = run_r_runner(
                    DSS_RUNNER, in_dir, dss_out,
                    extra_args=["--smoothing", "TRUE"],
                )
                logger.info("[seed=%d] dss: %s in %.1fs",
                            "ok" if status["dss"]["ok"] else "FAIL",
                            seed, status["dss"]["wall_s"])

        # 4. DSS smoothing=FALSE
        if not args.skip_dss_nosmooth:
            dss_ns_out = seed_dir / "dss_nosmooth.tsv"
            if dss_ns_out.exists():
                logger.info("[seed=%d] dss_nosmooth: %s exists, skipping",
                            seed, dss_ns_out.name)
                status["dss_nosmooth"] = {"ok": True, "wall_s": 0.0, "skipped": True}
            else:
                logger.info("[seed=%d] dss (smoothing=FALSE): launching...", seed)
                status["dss_nosmooth"] = run_r_runner(
                    DSS_RUNNER, in_dir, dss_ns_out,
                    extra_args=["--smoothing", "FALSE"],
                )
                logger.info("[seed=%d] dss_nosmooth: %s in %.1fs",
                            "ok" if status["dss_nosmooth"]["ok"] else "FAIL",
                            seed, status["dss_nosmooth"]["wall_s"])

        per_seed_status.append(status)
        logger.info("[seed=%d] done", seed)

    overall_wall = time.perf_counter() - overall_t0
    logger.info("=== SWEEP COMPLETE: %d seeds in %.1fs (%.1f min) ===",
                len(per_seed_status), overall_wall, overall_wall / 60.0)

    # Aggregate the status into a small log for the user.
    status_path = SIM_ROOT / "external_sweep_status.txt"
    with open(status_path, "w", encoding="utf-8") as f:
        f.write(f"run_external_simulator_sweep.py status\n")
        f.write(f"date_finished: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"seeds_processed: {len(per_seed_status)}\n")
        f.write(f"total_wallclock_s: {overall_wall:.1f}\n\n")
        for s in per_seed_status:
            f.write(f"seed={s['seed']}: convert={s.get('convert_n_new', 0)} new files\n")
            for tool in ("methylkit", "dss", "dss_nosmooth"):
                if tool in s:
                    r = s[tool]
                    flag = "skipped" if r.get("skipped") else (
                        "ok" if r.get("ok") else "FAIL"
                    )
                    f.write(f"  {tool}: {flag} ({r.get('wall_s', 0):.1f}s)\n")
    logger.info("wrote status: %s", status_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
