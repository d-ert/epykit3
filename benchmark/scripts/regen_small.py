"""regen_small.py -- CI-gateable engine-regression slice.

Phase 1.4 of the GB resubmission plan asks for a <10-minute slice of
``regen_all.py`` that produces a handful of reference hashes the CI
can diff against. The goal is to catch engine regressions (changes to
``dmc.py``, ``_glm.py``, ``dmr.py`` that move the per-CpG numerics)
in pull requests, BEFORE they land in the headline benchmark numbers.

Approach: run a small deterministic DMC sweep against a fixed-seed
Piao simulator instance and hash the engine's output columns. The
hash file (``benchmark/scripts/regen_small_hashes.json``) is committed
to git; CI re-runs the same slice and asserts the hashes match.

The slice is much smaller than the headline benchmark (10k CpGs, one
seed, two engines) so it finishes in well under 10 minutes on the CI
runners (Linux, py3.12). It exercises:

  - The Piao simulator (with phi=0, the binomial default).
  - epykit's bare ``lr`` engine (the recommended default).
  - epykit's ``power_stack="lr+"`` (the opt-in tunable; we hash it
    here so engine refactors of any of its four components also
    show up).

What it does NOT do (deliberately):
  - Run methylKit, DSS, dmrseq, BSmooth -- those need the R container
    and are too slow for a per-PR CI lane.
  - Score against truth -- the goal here is engine-output stability,
    not benchmark metric stability.

Usage:
  python regen_small.py             # run the slice + diff hashes
  python regen_small.py --update    # re-write the hash file (use after
                                    # an intentional engine change)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import polars as pl

logger = logging.getLogger(__name__)

_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

REFERENCE_PATH = _SCRIPTS_DIR / "regen_small_hashes.json"

# Hashed columns from the DMC parquet. We hash the rounded-to-8-decimal
# numerical columns so platform-specific float jitter at the 1e-9 level
# does not flake the CI. Engine refactors that move output by >1e-8
# (the meaningful regression signal) still trigger a hash mismatch.
HASHED_COLUMNS = ("chrom", "pos", "pvalue", "qvalue", "meth_diff")
NUMERIC_ROUND_DIGITS = 8


def _slice_config() -> dict:
    """Frozen config for the regression slice -- DO NOT change without
    re-snapshotting the hashes."""
    return {
        "n_cpgs": 10_000,
        "n_per_group": 3,
        "coverage": 10,
        "seed": 2026_000,
        "dmc_fraction": 0.20,
        "phi": 0.0,
    }


def _hash_parquet(parquet_path: Path) -> str:
    """SHA-256 of selected DMC columns at fixed numeric precision.

    Returns a hex digest. Robust to row-order shuffles because we
    sort by (chrom, pos) before hashing.
    """
    df = pl.read_parquet(str(parquet_path))
    keep = [c for c in HASHED_COLUMNS if c in df.columns]
    if "chrom" not in keep or "pos" not in keep:
        raise RuntimeError(
            f"parquet {parquet_path} missing chrom/pos columns "
            f"required for stable hashing; got {df.columns}"
        )
    df = df.sort(["chrom", "pos"]).select(keep)
    # Round numeric columns to a stable precision; cast all to string
    # to get a canonical byte representation regardless of column dtype.
    numeric_cols = [c for c in keep if c not in ("chrom", "pos")]
    df = df.with_columns([
        pl.col(c).round(NUMERIC_ROUND_DIGITS).cast(pl.Utf8)
        for c in numeric_cols
    ])
    df = df.with_columns([
        pl.col("chrom").cast(pl.Utf8),
        pl.col("pos").cast(pl.Utf8),
    ])
    buf = df.write_csv(separator="\t", include_header=True).encode("utf-8")
    return hashlib.sha256(buf).hexdigest()


def _build_methyldata_from_amp(amp_files, sample_ids, groups, store_dir):
    """Build a MethylData from AMP-format simulator output.

    The simulator writes AMP-format files (Piao layout). We convert
    each to 6-col bismark .cov.gz, drop them in `store_dir`, and run
    ``ep.read_bismark`` to build the methylstore. Mirrors
    ``run_epykit_simulator.py`` but stripped to the minimum needed
    for the small slice.
    """
    import epykit as ep

    cov_dir = Path(store_dir) / "bismark_cov"
    cov_dir.mkdir(parents=True, exist_ok=True)

    # Re-export each AMP file as a 6-col bismark .cov.gz that
    # ep.read_bismark can ingest.
    cov_paths = []
    for amp_path, sid in zip(amp_files, sample_ids):
        amp = pl.read_csv(str(amp_path), separator="\t")
        # AMP columns: chrBase, chr, base, strand, coverage, freqC, freqT
        cov = amp.select(
            chrom=pl.col("chr"),
            start=pl.col("base"),
            end=pl.col("base"),
            beta=pl.col("freqC"),
            count_M=(pl.col("freqC") / 100.0 * pl.col("coverage")).round().cast(pl.Int64),
            count_U=(pl.col("freqT") / 100.0 * pl.col("coverage")).round().cast(pl.Int64),
        )
        # Plain .cov (not .cov.gz) -- polars rejects the .gz extension
        # on uncompressed writes, and ep.read_bismark auto-detects the
        # bismark layout from the column count, not the extension.
        out = cov_dir / f"{sid}.cov"
        cov.write_csv(str(out), separator="\t", include_header=False)
        cov_paths.append(out)

    # Build samplesheet.
    sheet = pl.DataFrame({
        "sample_id": sample_ids,
        "group": groups,
        "path": [str(p) for p in cov_paths],
    })
    sheet_path = Path(store_dir) / "samples.csv"
    sheet.write_csv(str(sheet_path))

    md = ep.read_bismark(
        str(sheet_path),
        store_dir=str(store_dir),
        treatment_group="treat",
        control_group="ctrl",
    )
    return md


def regen_small(workdir: Path) -> dict:
    """Run the slice and return ``{engine_label: hash}``."""
    from simulate_piao import simulate_dmc
    import epykit as ep

    cfg = _slice_config()
    sim_dir = workdir / "sim"
    sim_dir.mkdir(parents=True, exist_ok=True)
    sim = simulate_dmc(
        n_cpgs=cfg["n_cpgs"],
        n_per_group=cfg["n_per_group"],
        coverage=cfg["coverage"],
        seed=cfg["seed"],
        dmc_fraction=cfg["dmc_fraction"],
        phi=cfg["phi"],
        out_dir=sim_dir,
    )
    logger.info("simulated %d cpgs, %d files", cfg["n_cpgs"], len(sim["amp_files"]))

    sample_ids = [f"treat_{i}" for i in range(1, cfg["n_per_group"] + 1)] + \
                 [f"ctrl_{i}"  for i in range(1, cfg["n_per_group"] + 1)]
    groups = ["treat"] * cfg["n_per_group"] + ["ctrl"] * cfg["n_per_group"]

    store_dir = workdir / "store"
    md = _build_methyldata_from_amp(
        sim["amp_files"], sample_ids, groups, store_dir,
    )

    hashes = {}

    # Engine 1: bare lr (recommended default).
    ep.tl.dmc(md, test="lr", power_stack="off")
    out_lr = workdir / "dmc_lr.parquet"
    md.dmc.write_parquet(str(out_lr))
    hashes["lr"] = _hash_parquet(out_lr)

    # Engine 2: lr+ (opt-in research tunable).
    ep.tl.dmc(md, test="lr", power_stack="lr+")
    out_lrplus = workdir / "dmc_lrplus.parquet"
    md.dmc.write_parquet(str(out_lrplus))
    hashes["lr_plus"] = _hash_parquet(out_lrplus)

    return hashes


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--update", action="store_true",
        help="Re-write the reference hash file. Use after an "
             "intentional engine change.",
    )
    parser.add_argument(
        "--keep-workdir", action="store_true",
        help="Leave the temp work dir on disk for debugging.",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    workdir = Path(tempfile.mkdtemp(prefix="regen_small_"))
    try:
        hashes = regen_small(workdir)
    finally:
        if not args.keep_workdir:
            shutil.rmtree(workdir, ignore_errors=True)

    print("computed hashes:")
    for k, v in sorted(hashes.items()):
        print(f"  {k:10s} {v}")

    if args.update:
        REFERENCE_PATH.write_text(json.dumps(hashes, indent=2) + "\n")
        print(f"updated reference: {REFERENCE_PATH}")
        return 0

    if not REFERENCE_PATH.exists():
        print(
            f"\nNo reference hash file at {REFERENCE_PATH}.\n"
            f"Run with --update once on Linux to snapshot the hashes;\n"
            f"commit regen_small_hashes.json so CI has a target."
        )
        return 1

    reference = json.loads(REFERENCE_PATH.read_text())
    mismatches = []
    for k, v in hashes.items():
        if reference.get(k) != v:
            mismatches.append((k, reference.get(k, "<missing>"), v))
    if mismatches:
        print("\nHASH MISMATCH (engine output regression detected):")
        for k, ref, got in mismatches:
            print(f"  {k}: reference={ref} got={got}")
        print(
            "\nIf this change is intentional, re-run with --update and "
            "commit the new hash file alongside the engine change. "
            "If not, investigate before merging."
        )
        return 1

    print(f"\nOK: all {len(hashes)} hashes match {REFERENCE_PATH.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
