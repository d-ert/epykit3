"""run_epykit_study1.py -- Phase 4 Task 2.

Re-runs epykit on the Piao-as-distributed Study 1+2 simulation grid (the
same raw inputs that produced ``eval_summary.parquet`` pre-Phase-1) using
the post-Phase-3 surviving engine surface. Writes:

* per-cell DMC/DMR parquets under
  ``benchmark/data/study1/epykit_post_phase3/<scenario>/<cell>_<test>.parquet``
* ``benchmark/data/study1/timings_post_phase3.parquet`` (per-cell wallclock)
* ``benchmark/data/study1/eval_summary_post_phase3.parquet`` (rebuilt by
  keeping the non-epykit baseline rows from the pre-existing
  ``eval_summary.parquet`` and concatenating fresh epykit rows)

The runner is self-contained: it converts the AMP-format simulated files
to gzipped 6-col Bismark ``.cov.gz`` (epykit's native ingestion path), runs
``ep.tl.dmc`` / ``ep.tl.dmr`` for each surviving engine, scores against the
canonical ground-truth parquets, and reassembles the summary table.

Engines:
* DMC: lr (default), lr+ (power_stack=True explicit kwargs), welch_t,
  fisher. glm is skipped: this dataset has no covariates and reduces to lr.
* DMR: tile, chain_merge, sliding_window, segment.

The ``epykit_bb_lr`` rows from the old ``eval_summary.parquet`` are
filtered out on reassembly: the bb_lr engine no longer exists post-Phase-3
(``ValueError`` with a migration hint to ``test='lr'``).

Usage:
    uv run python benchmark/scripts/run_epykit_study1.py
    uv run python benchmark/scripts/run_epykit_study1.py --only dmc_replicate:2
    uv run python benchmark/scripts/run_epykit_study1.py --skip-eval

Exit non-zero on any per-cell failure.
"""

from __future__ import annotations

import argparse
import gc
import gzip
import json
import logging
import os
import shutil
import sys
import time
from pathlib import Path

import polars as pl

# Make sibling scripts importable when this file is invoked as a script
# from the repo root (not as a package).
_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

# Re-exported scoring helpers + constants live in _epykit_scoring.py
# (single source of truth across Study 1 + Study 1b). The bare names
# stay accessible on this module for back-compat with the existing test
# suite -- F401-noqa marks each one as a deliberate re-export, not dead
# code. Ruff is configured to F-only in this repo, so a per-import
# noqa keeps the lint baseline green without disabling the broader
# F-rule.
from _epykit_scoring import (  # noqa: E402, F401
    DMR_OVERLAP_THRESHOLD,
    ENGINE_EXCEPTIONS,
    METH_DIFF_BINS,
    P_THRESHOLDS,
    Q_THRESHOLD,
    STALE_EPYKIT_TOOLS,
    _auroc,
    _confusion,
    _dmc_kwargs,
    _join_with_truth,
    _merge_intervals,
    score_dmc_parquet,
    score_dmr_parquet,
    split_for_ci,
)
from _epykit_scoring import reassemble_eval_summary as _reassemble_eval_summary  # noqa: E402


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent  # benchmark/
RAW_SIM = ROOT / "raw_sim_data" / "simulated_datasets"
DATA_STUDY1 = ROOT / "data" / "study1"
TRUTH_DIR = DATA_STUDY1 / "ground_truth"
OUT_BASE = DATA_STUDY1 / "epykit_post_phase3"
CONVERT_CACHE = ROOT / "_converted_post_phase3"
RUN_CACHE = ROOT / "_runs_post_phase3"

EVAL_SUMMARY_OLD = DATA_STUDY1 / "eval_summary.parquet"
EVAL_SUMMARY_NEW = DATA_STUDY1 / "eval_summary_post_phase3.parquet"
TIMINGS_NEW = DATA_STUDY1 / "timings_post_phase3.parquet"

logger = logging.getLogger("run_epykit_study1")


# ---------------------------------------------------------------------------
# Cell grid
# ---------------------------------------------------------------------------

COVERAGES_DMC = (5, 10, 15, 20, 25)
REPLICATES = (2, 4, 6, 8, 10)
COVERAGES_DMR = (5, 10, 15, 20, 25)

DMC_TESTS = ("lr", "lr+", "welch_t", "fisher")
DMR_METHODS = ("tile", "chain_merge", "sliding_window", "segment")

# DMC threshold grid + DMR overlap threshold now live in _epykit_scoring
# (single source of truth across Study 1 + Study 1b). They are re-exported
# at the top of this module so callers can keep using the bare names.


# ---------------------------------------------------------------------------
# AMP -> Bismark .cov.gz conversion
# ---------------------------------------------------------------------------

_AMP_SCHEMA = {
    "chrBase": pl.Utf8,
    "chr": pl.Utf8,
    "base": pl.Int64,
    "strand": pl.Utf8,
    "coverage": pl.Int64,
    "freqC": pl.Float64,
    "freqT": pl.Float64,
}


def amp_to_bismark_cov(src: Path, dst: Path) -> Path:
    """Convert one AMP-format ``.txt`` to gzipped 6-col Bismark ``.cov``.

    Re-derives integer M/U counts from ``round(coverage * freqC / 100)`` so
    downstream binomial / quasi-binomial tests see clean integers. Idempotent.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        return dst

    df = pl.read_csv(src, separator="\t", schema=_AMP_SCHEMA, has_header=True)
    df = (
        df.with_columns(
            count_M=(pl.col("coverage") * pl.col("freqC") / 100.0)
            .round()
            .cast(pl.Int64),
        )
        .with_columns(
            count_U=(pl.col("coverage") - pl.col("count_M")).clip(lower_bound=0),
        )
        .select(
            chrom=pl.col("chr"),
            start=pl.col("base"),
            end=pl.col("base"),
            meth_pct=pl.col("freqC"),
            count_M=pl.col("count_M"),
            count_U=pl.col("count_U"),
        )
    )

    tmp = dst.with_suffix(dst.suffix + ".tmp")
    df.write_csv(tmp, separator="\t", include_header=False)
    with open(tmp, "rb") as fin, gzip.open(dst, "wb", compresslevel=4) as fout:
        shutil.copyfileobj(fin, fout)
    tmp.unlink()
    return dst


def _convert_scenario(src_dir: Path, names: list[str], out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    for src_name, dst_name in names:
        amp_to_bismark_cov(src_dir / src_name, out_dir / dst_name)
    return out_dir


def convert_dmc_coverage(cov: int) -> Path:
    src = RAW_SIM / "dmc_simulation" / "coverage"
    out = CONVERT_CACHE / f"dmc_coverage_{cov}"
    names = [
        (f"amp.coverage={cov}.sample{i}.txt", f"sample{i}.cov.gz") for i in range(1, 7)
    ]
    return _convert_scenario(src, names, out)


def convert_dmc_replicate(n_total: int) -> Path:
    src = RAW_SIM / "dmc_simulation" / "replicate"
    out = CONVERT_CACHE / f"dmc_replicate_{n_total}"
    names = [
        (f"amp.replicate={n_total}.sample{i}.txt", f"sample{i}.cov.gz")
        for i in range(1, n_total + 1)
    ]
    return _convert_scenario(src, names, out)


def convert_dmr_coverage(cov: int) -> Path:
    src = RAW_SIM / "dmr_simulation" / "coverage"
    out = CONVERT_CACHE / f"dmr_coverage_{cov}"
    names = [
        (f"amp.coverage={cov}.sample{i}.txt", f"sample{i}.cov.gz") for i in range(1, 7)
    ]
    return _convert_scenario(src, names, out)


# ---------------------------------------------------------------------------
# Samplesheet builder
# ---------------------------------------------------------------------------


def write_samplesheet(sample_dir: Path, n_per_group: int, sheet_path: Path) -> Path:
    """Build samplesheet: first ``n_per_group`` are 'treat', next are 'ctrl'."""
    rows: list[tuple[str, str, str]] = []
    for i in range(1, n_per_group + 1):
        rows.append((f"treat_{i}", "treat", str(sample_dir / f"sample{i}.cov.gz")))
    for j in range(1, n_per_group + 1):
        i = n_per_group + j
        rows.append((f"ctrl_{j}", "ctrl", str(sample_dir / f"sample{i}.cov.gz")))
    df = pl.DataFrame(rows, schema=["sample_id", "group", "path"], orient="row")
    sheet_path.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(sheet_path)
    return sheet_path


# ---------------------------------------------------------------------------
# Per-scenario epykit drivers
# ---------------------------------------------------------------------------


def run_dmc_cell(
    sample_dir: Path,
    samplesheet: Path,
    label: str,
    scenario: str,
    parameter: str,
    parameter_value: int,
    tests: tuple[str, ...],
    n_per_group: int,
    out_dir: Path,
) -> list[dict]:
    """Run epykit DMC for one (scenario, parameter_value) cell across `tests`.

    Returns list of timing rows. Writes per-test parquets to ``out_dir``.
    """
    import epykit as ep

    store_dir = RUN_CACHE / scenario / f"cell_{parameter_value}"
    if store_dir.exists():
        shutil.rmtree(store_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    md = ep.read_bismark(
        str(samplesheet),
        treatment_group="treat",
        control_group="ctrl",
        assembly="hg18",
        store_dir=str(store_dir),
    )
    ep.pp.unite(md, type="intersect")

    allow_n1 = n_per_group < 2
    timings: list[dict] = []

    for test in tests:
        backend_test, kwargs = _dmc_kwargs(test, allow_n1=allow_n1)
        fname_test = test.replace("+", "plus")
        out_parquet = out_dir / f"{label}_{fname_test}.parquet"

        t0 = time.perf_counter()
        try:
            ep.tl.dmc(md, **kwargs)
        except ENGINE_EXCEPTIONS as exc:
            elapsed = time.perf_counter() - t0
            logger.error("[%s] DMC test=%s FAILED in %.1fs: %r",
                         label, test, elapsed, exc)
            timings.append({
                "scenario": scenario, "parameter": parameter,
                "parameter_value": parameter_value,
                "tool": f"epykit_{fname_test}", "test": test,
                "wall_s": elapsed, "n_dmc_called": None,
                "ok": False, "error": str(exc),
            })
            continue
        elapsed = time.perf_counter() - t0

        df = md.get_dmc(test=backend_test)
        df.write_parquet(out_parquet)
        n_sig = (
            df.filter(pl.col("qvalue") < Q_THRESHOLD).height
            if "qvalue" in df.columns
            else None
        )
        timings.append({
            "scenario": scenario, "parameter": parameter,
            "parameter_value": parameter_value,
            "tool": f"epykit_{fname_test}", "test": test,
            "wall_s": elapsed, "n_dmc_called": n_sig,
            "ok": True, "error": None,
        })
        logger.info("[%s] DMC test=%s: %s sites, %s sig (q<%.2f), %.1fs",
                    label, test, f"{df.height:,}",
                    f"{n_sig:,}" if n_sig is not None else "?",
                    Q_THRESHOLD, elapsed)

    del md
    gc.collect()
    return timings


def run_dmr_cell(
    sample_dir: Path,
    samplesheet: Path,
    label: str,
    scenario: str,
    parameter: str,
    parameter_value: int,
    methods: tuple[str, ...],
    n_per_group: int,
    out_dir: Path,
) -> list[dict]:
    """Run epykit DMR for one cell. DMC (lr) is run once first so the DMC
    store is present for sliding_window / chain_merge / segment.
    """
    import epykit as ep

    store_dir = RUN_CACHE / scenario / f"cell_{parameter_value}"
    if store_dir.exists():
        shutil.rmtree(store_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    md = ep.read_bismark(
        str(samplesheet),
        treatment_group="treat",
        control_group="ctrl",
        assembly="hg18",
        store_dir=str(store_dir),
    )
    ep.pp.unite(md, type="intersect")
    ep.tl.dmc(md, test="lr", allow_n1=n_per_group < 2)

    timings: list[dict] = []
    for method in methods:
        out_parquet = out_dir / f"{label}_dmr_{method}.parquet"

        t0 = time.perf_counter()
        try:
            ep.tl.dmr(md, method=method)
        except ENGINE_EXCEPTIONS as exc:
            elapsed = time.perf_counter() - t0
            logger.error("[%s] DMR method=%s FAILED in %.1fs: %r",
                         label, method, elapsed, exc)
            timings.append({
                "scenario": scenario, "parameter": parameter,
                "parameter_value": parameter_value,
                "tool": f"epykit_dmr_{method}", "test": method,
                "wall_s": elapsed, "n_dmc_called": None,
                "ok": False, "error": str(exc),
            })
            continue
        elapsed = time.perf_counter() - t0

        dmr_df = md.uns.get("dmr")
        if isinstance(dmr_df, dict):
            dmr_df = dmr_df.get("frame") or next(iter(dmr_df.values()))
        dmr_df.write_parquet(out_parquet)
        n_sig = (
            dmr_df.filter(pl.col("qvalue") < Q_THRESHOLD).height
            if "qvalue" in dmr_df.columns
            else dmr_df.height
        )
        timings.append({
            "scenario": scenario, "parameter": parameter,
            "parameter_value": parameter_value,
            "tool": f"epykit_dmr_{method}", "test": method,
            "wall_s": elapsed, "n_dmc_called": n_sig,
            "ok": True, "error": None,
        })
        logger.info("[%s] DMR method=%s: %s regions, %s sig, %.1fs",
                    label, method, f"{dmr_df.height:,}",
                    f"{n_sig:,}", elapsed)

    del md
    gc.collect()
    return timings


# ---------------------------------------------------------------------------
# Evaluation + reassembly
#
# Scoring primitives (_confusion / _auroc / _join_with_truth /
# score_dmc_parquet / score_dmr_parquet / _merge_intervals) and the
# stale-engine filter (STALE_EPYKIT_TOOLS / reassemble_eval_summary) now
# live in _epykit_scoring.py and are imported at the top of this module.
# This module passes EVAL_SUMMARY_OLD as the baseline path explicitly so
# the helper stays pure (no module-level globals).
# ---------------------------------------------------------------------------


def reassemble_eval_summary(new_rows: list[dict]) -> pl.DataFrame:
    """Thin wrapper around ``_epykit_scoring.reassemble_eval_summary`` that
    binds ``EVAL_SUMMARY_OLD`` for back-compat with the test suite.
    """
    return _reassemble_eval_summary(new_rows, EVAL_SUMMARY_OLD)


# Back-compat alias for the test suite (the old name was a module global).
_STALE_EPYKIT_TOOLS = STALE_EPYKIT_TOOLS


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def _cells_for(only: str | None) -> list[tuple[str, int]]:
    """Resolve --only filter; returns (scenario, parameter_value) list."""
    full = [("dmc_coverage", c) for c in COVERAGES_DMC]
    full += [("dmc_replicate", r) for r in REPLICATES]
    full += [("dmr_coverage", c) for c in COVERAGES_DMR]
    if only is None:
        return full
    if ":" not in only:
        raise SystemExit(f"--only expects 'scenario:value', got {only!r}")
    scen, val = only.split(":", 1)
    return [(scen, int(val))]


def _git_head() -> str:
    try:
        import subprocess
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT.parent, text=True
        )
        return out.strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def _write_manifest(timings: pl.DataFrame, new_row_count: int) -> None:
    import epykit as ep
    OUT_BASE.mkdir(parents=True, exist_ok=True)
    manifest_path = OUT_BASE / "MANIFEST.txt"

    total_wall = float(timings["wall_s"].sum()) if timings.height else 0.0

    per_engine = (
        timings.filter(pl.col("ok"))
        .group_by("tool").agg(
            n_cells=pl.len(),
            total_s=pl.col("wall_s").sum(),
            median_s=pl.col("wall_s").median(),
        )
        .sort("tool")
        if timings.height else None
    )

    lines = [
        "epykit_post_phase3 -- Phase 4 Task 2 (run_epykit_study1.py)",
        "",
        f"Date           : {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"epykit version : {ep.__version__}",
        "Engine tag     : v0.7.5-phase3-engines-frozen",
        f"Git HEAD       : {_git_head()}",
        "",
        "DMC engines run: lr, lr+ (power stack), welch_t, fisher",
        "DMR methods run: tile, chain_merge, sliding_window, segment",
        "Skipped        : glm (no covariate columns in obs; reduces to lr).",
        "",
        "lr+ kwargs (explicit on top of test='lr'):",
        "    neighbour_combine=True, neighbour_bp=500, sep_fallback=True,",
        "    sep_threshold=0.9, fdr_method='fdr_tsbh', dispersion='eb'",
        "",
        f"Total wallclock: {total_wall:.1f}s",
        f"eval_summary_post_phase3.parquet rows added: {new_row_count}",
        "",
        "Per-engine wallclock summary:",
    ]
    if per_engine is not None and per_engine.height:
        for r in per_engine.iter_rows(named=True):
            lines.append(
                f"  {r['tool']:<28}  n={r['n_cells']:<2}  "
                f"total={r['total_s']:.1f}s  median={r['median_s']:.1f}s"
            )

    lines += [
        "",
        "Reassembly:",
        "    - Read benchmark/data/study1/eval_summary.parquet (pre-fix baseline).",
        "    - Dropped all rows whose tool starts with 'epykit_' (stale pre-Phase-3).",
        "      Explicitly stale: epykit_bb_lr (engine removed in 0.7.5).",
        "    - Concatenated fresh post-Phase-3 epykit rows.",
        "    - Wrote benchmark/data/study1/eval_summary_post_phase3.parquet.",
        "    - Original eval_summary.parquet kept untouched as pre-fix baseline.",
        "",
        "Per-cell intermediate parquets live under this directory:",
        "    epykit_post_phase3/dmc_coverage/<label>_<test>.parquet",
        "    epykit_post_phase3/dmc_replicate/<label>_<test>.parquet",
        "    epykit_post_phase3/dmr_coverage/<label>_dmr_<method>.parquet",
        "These are gitignored bulk artefacts; regenerable via run_epykit_study1.py.",
        "",
    ]
    manifest_path.write_text("\n".join(lines), encoding="utf-8")
    logger.info("wrote %s", manifest_path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--only", default=None,
        help="Run a single cell, e.g. --only dmc_replicate:2",
    )
    parser.add_argument(
        "--skip-run", action="store_true",
        help="Skip the runner stage; re-score / reassemble from cached parquets.",
    )
    parser.add_argument(
        "--skip-eval", action="store_true",
        help="Skip scoring + reassembly; only run epykit and write parquets.",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="DEBUG-level logging.",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        force=True,
    )

    cells = _cells_for(args.only)
    logger.info("running %d cell(s)", len(cells))

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    all_timings: list[dict] = []

    if not args.skip_run:
        for scenario, val in cells:
            cell_dir = OUT_BASE / scenario
            cell_dir.mkdir(parents=True, exist_ok=True)

            if scenario == "dmc_coverage":
                cov = val
                sample_dir = convert_dmc_coverage(cov)
                sheet = CONVERT_CACHE / f"dmc_coverage_{cov}" / "samplesheet.csv"
                write_samplesheet(sample_dir, n_per_group=3, sheet_path=sheet)
                label = f"dmc_cov{cov}"
                logger.info("=== %s (cov=%d) ===", scenario, cov)
                tims = run_dmc_cell(
                    sample_dir, sheet, label, scenario, "coverage", cov,
                    tests=DMC_TESTS, n_per_group=3, out_dir=cell_dir,
                )
                all_timings.extend(tims)

            elif scenario == "dmc_replicate":
                n_total = val
                n_per_group = n_total // 2
                sample_dir = convert_dmc_replicate(n_total)
                sheet = CONVERT_CACHE / f"dmc_replicate_{n_total}" / "samplesheet.csv"
                write_samplesheet(
                    sample_dir, n_per_group=n_per_group, sheet_path=sheet,
                )
                label = f"dmc_rep{n_total}"
                # n_per_group=1 (n_total=2) -- only fisher is valid.
                tests = DMC_TESTS if n_per_group >= 2 else ("fisher",)
                logger.info("=== %s (n=%d, %dv%d) ===",
                            scenario, n_total, n_per_group, n_per_group)
                tims = run_dmc_cell(
                    sample_dir, sheet, label, scenario, "n_total", n_total,
                    tests=tests, n_per_group=n_per_group, out_dir=cell_dir,
                )
                all_timings.extend(tims)

            elif scenario == "dmr_coverage":
                cov = val
                sample_dir = convert_dmr_coverage(cov)
                sheet = CONVERT_CACHE / f"dmr_coverage_{cov}" / "samplesheet.csv"
                write_samplesheet(sample_dir, n_per_group=3, sheet_path=sheet)
                label = f"dmr_cov{cov}"
                logger.info("=== %s (cov=%d) ===", scenario, cov)
                tims = run_dmr_cell(
                    sample_dir, sheet, label, scenario, "coverage", cov,
                    methods=DMR_METHODS, n_per_group=3, out_dir=cell_dir,
                )
                all_timings.extend(tims)

            else:
                raise SystemExit(f"unknown scenario {scenario!r}")

        if all_timings:
            timings_df = pl.DataFrame(all_timings)
            # Merge with any existing timings parquet if --only was used.
            if TIMINGS_NEW.exists() and args.only is not None:
                old_t = pl.read_parquet(TIMINGS_NEW)
                key_cols = ["scenario", "parameter_value", "tool", "test"]
                # Drop rows in old_t that we just regenerated, then concat.
                pairs = (
                    timings_df.select(key_cols).unique().to_dicts()
                )
                mask = pl.lit(False)
                for p in pairs:
                    cond = pl.lit(True)
                    for k in key_cols:
                        cond = cond & (pl.col(k) == p[k])
                    mask = mask | cond
                old_t = old_t.filter(~mask)
                timings_df = pl.concat([old_t, timings_df], how="diagonal_relaxed")
            TIMINGS_NEW.parent.mkdir(parents=True, exist_ok=True)
            timings_df.write_parquet(TIMINGS_NEW)
            logger.info("wrote %s (%d rows)", TIMINGS_NEW, timings_df.height)

    n_fail = sum(1 for t in all_timings if not t.get("ok", True))
    if n_fail:
        logger.error("%d cell-test combinations FAILED", n_fail)

    if args.skip_eval:
        logger.info("--skip-eval set; not rebuilding eval_summary")
        return 1 if n_fail else 0

    # --- Score + reassemble ------------------------------------------------
    dmc_truth = pl.read_parquet(TRUTH_DIR / "dmc_truth.parquet")
    dmc_truth_dmr_sim = pl.read_parquet(TRUTH_DIR / "dmc_truth_dmr_sim.parquet")
    dmr_truth = pl.read_parquet(TRUTH_DIR / "dmr_truth.parquet")

    all_rows: list[dict] = []
    # DMC coverage
    for parquet in sorted((OUT_BASE / "dmc_coverage").glob("dmc_cov*_*.parquet")):
        stem = parquet.stem  # dmc_cov10_lr  /  dmc_cov10_lrplus
        parts = stem.split("_")
        cov = int(parts[1].removeprefix("cov"))
        test = "_".join(parts[2:])
        tool_test = test if test != "lrplus" else "lrplus"
        all_rows.extend(score_dmc_parquet(
            parquet, dmc_truth,
            tool=f"epykit_{tool_test}",
            scenario="dmc_coverage", parameter="coverage",
            parameter_value=cov, test=test,
        ))

    # DMC replicate
    for parquet in sorted((OUT_BASE / "dmc_replicate").glob("dmc_rep*_*.parquet")):
        stem = parquet.stem  # dmc_rep4_lr
        parts = stem.split("_")
        n_total = int(parts[1].removeprefix("rep"))
        test = "_".join(parts[2:])
        all_rows.extend(score_dmc_parquet(
            parquet, dmc_truth,
            tool=f"epykit_{test}",
            scenario="dmc_replicate", parameter="n_total",
            parameter_value=n_total, test=test,
        ))

    # DMR coverage
    for parquet in sorted((OUT_BASE / "dmr_coverage").glob("dmr_cov*_dmr_*.parquet")):
        stem = parquet.stem  # dmr_cov10_dmr_tile  /  dmr_cov10_dmr_chain_merge
        # split off the cov prefix; everything after '_dmr_' is the method.
        before, _, after = stem.partition("_dmr_")
        cov = int(before.split("_", 1)[1].removeprefix("cov"))
        method = after
        all_rows.extend(score_dmr_parquet(
            parquet, dmr_truth,
            tool=f"epykit_dmr_{method}",
            scenario="dmr_coverage", parameter_value=cov, method=method,
        ))

    combined = reassemble_eval_summary(all_rows)
    EVAL_SUMMARY_NEW.parent.mkdir(parents=True, exist_ok=True)
    combined.write_parquet(EVAL_SUMMARY_NEW)
    logger.info("wrote %s (%d rows; %d new epykit rows)",
                EVAL_SUMMARY_NEW, combined.height, len(all_rows))

    # --- CI columns ------------------------------------------------------
    # evaluate.py --ci-only expects every row to have valid tp/fp/tn/fn
    # (so Wilson can binomtest on integer counts). Our DMR rows + the
    # transcribed Piao baseline rows have NULL tp/fn/fp/tn -- those are
    # split off here (via _epykit_scoring.split_for_ci), the DMC subset
    # gets Wilson + bootstrap CIs added, and the no-counts subset is
    # concatenated back unchanged.
    has_counts, no_counts = split_for_ci(combined)
    logger.info("CI split: %d rows with counts, %d without",
                has_counts.height, no_counts.height)

    if has_counts.height:
        tmp_path = EVAL_SUMMARY_NEW.parent / ".eval_summary_ci_tmp.parquet"
        has_counts.write_parquet(tmp_path)
        import subprocess
        # Use the Phase-3-sealed evaluate.py without modification.
        evaluate_script = HERE / "evaluate.py"
        proc = subprocess.run(
            [sys.executable, str(evaluate_script),
             "--ci-only", "--eval-summary", str(tmp_path)],
            check=False, capture_output=True, text=True,
        )
        if proc.returncode != 0:
            logger.error("evaluate.py --ci-only FAILED:\n%s", proc.stderr)
            return 1
        logger.info("evaluate.py: %s", proc.stdout.strip())
        with_ci = pl.read_parquet(tmp_path)
        tmp_path.unlink(missing_ok=True)

        # Align schemas of no_counts to with_ci (add NaN CI cols).
        for c in ("tpr_ci_lo", "tpr_ci_hi", "fpr_ci_lo", "fpr_ci_hi",
                  "auroc_ci_lo", "auroc_ci_hi", "f1_ci_lo", "f1_ci_hi"):
            if c not in no_counts.columns:
                no_counts = no_counts.with_columns(pl.lit(float("nan")).alias(c))
        # Reorder cols so the concat is clean.
        all_cols = with_ci.columns
        no_counts = no_counts.select(all_cols)
        final = pl.concat([with_ci, no_counts], how="vertical_relaxed")
        final.write_parquet(EVAL_SUMMARY_NEW)
        logger.info("wrote %s (%d rows; CI cols added)",
                    EVAL_SUMMARY_NEW, final.height)

    # Manifest
    if TIMINGS_NEW.exists():
        timings_df = pl.read_parquet(TIMINGS_NEW)
    else:
        timings_df = pl.DataFrame()
    _write_manifest(timings_df, new_row_count=len(all_rows))

    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
