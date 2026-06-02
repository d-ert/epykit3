"""run_phase4_null_calibration.py -- Phase 4 Task 7.

Orchestrates null-calibration runs across (engine x dataset) cells for the
post-Phase-3 surviving epykit engine surface. For each dataset:

  1. Build a per-dataset methylstore (Bismark .cov.gz or combined-strand BED).
  2. For each engine, shell out to ``run_null_calibration.py`` to run
     ``--k-shuffles`` label-permutations and write a per-engine parquet.
  3. Clean up the methylstore cache.

After all datasets run, aggregate every per-engine parquet into
``benchmark/data/null_calibration/summary.parquet`` with one row per
(dataset, engine) pair carrying median observed_fdr + IQR + bootstrap CIs.

Datasets
--------

* **piao_distributed**   : Piao-as-distributed cov10_3v3 -- 6 AMP files at
                           ``benchmark/raw_sim_data/simulated_datasets/
                           dmc_simulation/coverage/amp.coverage=10.sample{1..6}.txt``.
                           Engines: lr, lr_plus, welch_t, fisher. 20 shuffles.
* **simulator**          : Held-out simulator seed=2026000 cov10_3v3 -- 6
                           AMP files at ``benchmark/data/study1b_simulator/
                           seed=2026000/amp.coverage=10.sample{1..6}.txt``.
                           Engines: lr, lr_plus, welch_t, fisher. 20 shuffles.
* **gse263850**          : Real cohort (3v3) at
                           ``../epykit2/GSE263850_RAW/*.bed.gz``. Only
                           ``C(6,3)/2 = 10`` unique 3v3 assignments exist
                           so k_shuffles=10. Engines: lr, lr_plus, welch_t,
                           fisher, glm (no covariates -- glm reduces to
                           single-treatment-coef lr, but we report it for
                           the Table-S-Calib row count).

Sealed Phase-3 scripts (``run_null_calibration.py``, ``_null_engines.py``,
``evaluate.py``, ``wilson_bootstrap_ci.py``) are invoked as-is; this
wrapper never imports their internals.

Usage
-----
    uv run python benchmark/scripts/run_phase4_null_calibration.py
    uv run python benchmark/scripts/run_phase4_null_calibration.py --dataset piao
    uv run python benchmark/scripts/run_phase4_null_calibration.py --skip-gse
    uv run python benchmark/scripts/run_phase4_null_calibration.py \\
        --dataset gse --engines lr lr_plus

Exit non-zero on any per-engine failure.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

# Reuse the AMP -> Bismark conversion + samplesheet writer from Task 2.
# These are pure helpers (no module-level state) so import is cheap.
from run_epykit_study1 import (  # noqa: E402
    amp_to_bismark_cov,
    write_samplesheet,
)
from run_epykit_gse import _resolve_samplesheet  # noqa: E402


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                          # benchmark/
REPO = ROOT.parent                          # epykit3/
RUN_NULL_SCRIPT = HERE / "run_null_calibration.py"

PIAO_AMP_DIR = ROOT / "raw_sim_data" / "simulated_datasets" / "dmc_simulation" / "coverage"
SIMULATOR_AMP_DIR = ROOT / "data" / "study1b_simulator" / "seed=2026000"
GSE_SAMPLESHEET = ROOT / "data" / "study3" / "samplesheet_gse263850.csv"

NULL_BASE = ROOT / "data" / "null_calibration"
PIAO_OUT = NULL_BASE / "piao_distributed" / "cov10_3v3"
SIM_OUT = NULL_BASE / "simulator" / "sim_cov10_3v3"
GSE_OUT = NULL_BASE / "gse263850"
MANIFEST_PATH = NULL_BASE / "MANIFEST.txt"
SUMMARY_PATH = NULL_BASE / "summary.parquet"

CACHE_BASE = ROOT / "_runs_null_calibration"

logger = logging.getLogger("run_phase4_null_calibration")


# ---------------------------------------------------------------------------
# Dataset specs
# ---------------------------------------------------------------------------


@dataclass
class DatasetSpec:
    """One dataset's null-calibration run plan."""

    key: str
    scenario: str
    out_dir: Path
    k_shuffles: int
    engines: tuple[str, ...]
    n_per_group: int = 3
    # Methylstore build callback: returns the absolute store_dir path.
    builder: callable = field(default=lambda: None)


def _build_piao_store() -> Path:
    """Build methylstore for Piao cov=10 cell (samples 1-3 treat, 4-6 ctrl).

    Returns the saved-snapshot path: ``MethylData.load()`` requires a
    ``methyldata.json`` snapshot, so after ingest + unite we ``md.save()``
    to a sibling directory and point the sealed runner at that.
    """
    cov = 10
    cache_root = CACHE_BASE / "piao_cov10"
    sample_dir = cache_root / "samples"
    raw_store = cache_root / "ingest_store"
    snapshot = cache_root / "snapshot"
    if raw_store.exists():
        shutil.rmtree(raw_store)
    if snapshot.exists():
        shutil.rmtree(snapshot)
    sample_dir.mkdir(parents=True, exist_ok=True)
    for i in range(1, 7):
        src = PIAO_AMP_DIR / f"amp.coverage={cov}.sample{i}.txt"
        dst = sample_dir / f"sample{i}.cov.gz"
        amp_to_bismark_cov(src, dst)
    sheet = cache_root / "samplesheet.csv"
    write_samplesheet(sample_dir, n_per_group=3, sheet_path=sheet)
    import epykit as ep  # noqa: PLC0415
    md = ep.read_bismark(
        str(sheet),
        treatment_group="treat",
        control_group="ctrl",
        assembly="hg18",
        store_dir=str(raw_store),
    )
    ep.pp.unite(md, type="intersect")
    md.save(str(snapshot))
    logger.info("piao methylstore built: %d samples, snapshot=%s",
                md.obs.height, snapshot)
    return snapshot


def _build_simulator_store() -> Path:
    """Build methylstore for simulator seed=2026000 cov=10 cell.

    Returns the saved-snapshot path (see ``_build_piao_store`` for the
    save/load rationale).
    """
    cov = 10
    cache_root = CACHE_BASE / "simulator_seed2026000_cov10"
    sample_dir = cache_root / "samples"
    raw_store = cache_root / "ingest_store"
    snapshot = cache_root / "snapshot"
    if raw_store.exists():
        shutil.rmtree(raw_store)
    if snapshot.exists():
        shutil.rmtree(snapshot)
    sample_dir.mkdir(parents=True, exist_ok=True)
    for i in range(1, 7):
        src = SIMULATOR_AMP_DIR / f"amp.coverage={cov}.sample{i}.txt"
        dst = sample_dir / f"sample{i}.cov.gz"
        amp_to_bismark_cov(src, dst)
    sheet = cache_root / "samplesheet.csv"
    write_samplesheet(sample_dir, n_per_group=3, sheet_path=sheet)
    import epykit as ep  # noqa: PLC0415
    md = ep.read_bismark(
        str(sheet),
        treatment_group="treat",
        control_group="ctrl",
        assembly="hg38",
        store_dir=str(raw_store),
    )
    ep.pp.unite(md, type="intersect")
    md.save(str(snapshot))
    logger.info("simulator methylstore built: %d samples, snapshot=%s",
                md.obs.height, snapshot)
    return snapshot


def _build_gse_store() -> Path:
    """Build methylstore for GSE263850 (3 sbp009 vs 3 clone).

    Returns the saved-snapshot path (see ``_build_piao_store``).
    """
    cache_root = CACHE_BASE / "gse263850"
    raw_store = cache_root / "ingest_store"
    snapshot = cache_root / "snapshot"
    if raw_store.exists():
        shutil.rmtree(raw_store)
    if snapshot.exists():
        shutil.rmtree(snapshot)
    cache_root.mkdir(parents=True, exist_ok=True)
    resolved_sheet = cache_root / "samplesheet_resolved.csv"
    sheet = _resolve_samplesheet(GSE_SAMPLESHEET, resolved_sheet)
    import epykit as ep  # noqa: PLC0415
    md = ep.read_combined_strand_bed(
        str(sheet),
        treatment_group="clone",
        control_group="sbp009",
        assembly="hg38",
        store_dir=str(raw_store),
    )
    ep.pp.unite(md, type="intersect")
    md.save(str(snapshot))
    logger.info("gse methylstore built: %d samples, snapshot=%s",
                md.obs.height, snapshot)
    return snapshot


def make_specs() -> dict[str, DatasetSpec]:
    """Construct the full (dataset -> spec) registry."""
    return {
        "piao": DatasetSpec(
            key="piao_distributed",
            scenario="cov10_3v3",
            out_dir=PIAO_OUT,
            k_shuffles=20,
            engines=("lr", "lr_plus", "welch_t", "fisher"),
            builder=_build_piao_store,
        ),
        "simulator": DatasetSpec(
            key="simulator",
            scenario="sim_cov10_3v3",
            out_dir=SIM_OUT,
            k_shuffles=20,
            engines=("lr", "lr_plus", "welch_t", "fisher"),
            builder=_build_simulator_store,
        ),
        "gse": DatasetSpec(
            key="gse263850",
            scenario="gse263850",
            out_dir=GSE_OUT,
            k_shuffles=10,    # C(6,3)/2 = 10 unique 3v3 assignments
            engines=("lr", "lr_plus", "welch_t", "fisher", "glm"),
            builder=_build_gse_store,
        ),
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def _engine_parquet(out_dir: Path, engine: str) -> Path:
    """Where one engine's per-shuffle parquet lives."""
    return out_dir / f"{engine}.parquet"


def run_engine(
    engine: str, store_dir: Path, spec: DatasetSpec,
    seed: int, q_thresh: float,
) -> dict:
    """Run one (engine, dataset) cell via the sealed runner; return timing row."""
    out_pq = _engine_parquet(spec.out_dir, engine)
    spec.out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable, str(RUN_NULL_SCRIPT),
        "--engine", engine,
        "--methylstore", str(store_dir),
        "--scenario", spec.scenario,
        "--k-shuffles", str(spec.k_shuffles),
        "--seed", str(seed),
        "--q-thresh", str(q_thresh),
        "--out", str(out_pq),
    ]
    logger.info("[%s/%s] launching: %s", spec.key, engine, " ".join(cmd))
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    elapsed = time.perf_counter() - t0
    if proc.returncode != 0:
        logger.error("[%s/%s] FAILED in %.1fs:\nSTDOUT:\n%s\nSTDERR:\n%s",
                     spec.key, engine, elapsed, proc.stdout, proc.stderr)
        return {
            "dataset": spec.key, "engine": engine, "scenario": spec.scenario,
            "k_shuffles": spec.k_shuffles, "wall_s": elapsed,
            "ok": False, "error": proc.stderr.strip().splitlines()[-1]
            if proc.stderr.strip() else "non-zero exit",
        }
    logger.info("[%s/%s] ok in %.1fs -- %s",
                spec.key, engine, elapsed, proc.stdout.strip())
    return {
        "dataset": spec.key, "engine": engine, "scenario": spec.scenario,
        "k_shuffles": spec.k_shuffles, "wall_s": elapsed,
        "ok": True, "error": None,
    }


def run_dataset(
    spec: DatasetSpec, engines: tuple[str, ...] | None,
    seed: int, q_thresh: float, keep_store: bool,
) -> list[dict]:
    """Build the methylstore once, then run every engine against it."""
    engines = engines or spec.engines
    logger.info("=== dataset=%s  scenario=%s  engines=%s  k=%d ===",
                spec.key, spec.scenario, engines, spec.k_shuffles)
    t0 = time.perf_counter()
    store_dir = spec.builder()
    logger.info("[%s] methylstore ready in %.1fs at %s",
                spec.key, time.perf_counter() - t0, store_dir)

    timings: list[dict] = []
    for engine in engines:
        timings.append(run_engine(engine, store_dir, spec, seed, q_thresh))

    if not keep_store:
        # store_dir is the snapshot path under CACHE_BASE/<dataset>/snapshot;
        # we delete the dataset's full cell root (snapshot + ingest_store +
        # samples + samplesheet).
        cell_root = store_dir.parent if store_dir.name == "snapshot" \
            else store_dir
        if cell_root.exists() and CACHE_BASE in cell_root.parents:
            try:
                shutil.rmtree(cell_root)
                logger.info("[%s] cleaned methylstore cache %s",
                            spec.key, cell_root)
            except OSError as exc:
                logger.warning("[%s] could not clean %s: %r",
                               spec.key, cell_root, exc)
    return timings


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------


def aggregate_summary(specs: dict[str, DatasetSpec]) -> pl.DataFrame:
    """Aggregate every per-engine parquet into the summary table.

    Schema (one row per (dataset, scenario, engine)):
        engine, dataset, scenario, k_shuffles,
        observed_fdr_median, observed_fdr_q1, observed_fdr_q3,
        observed_fdr_ci_lo, observed_fdr_ci_hi,
        n_sites_called_median, n_sites_total.

    The per-shuffle parquet already carries Wilson CI bounds per row; the
    summary's ``observed_fdr_ci_lo``/``observed_fdr_ci_hi`` are the
    across-shuffle aggregates (min lo / max hi) -- a conservative envelope
    that bounds the observed FDR over the shuffle set.
    """
    rows: list[dict] = []
    for spec_key, spec in specs.items():
        if not spec.out_dir.exists():
            continue
        for engine in spec.engines:
            pq = _engine_parquet(spec.out_dir, engine)
            if not pq.exists():
                logger.warning("missing parquet for %s/%s: %s",
                               spec.key, engine, pq)
                continue
            df = pl.read_parquet(pq)
            if df.height == 0:
                logger.warning("empty parquet for %s/%s", spec.key, engine)
                continue
            fdr = df["observed_fdr"].to_numpy()
            ci_lo = df["observed_fdr_ci_lo"].to_numpy()
            ci_hi = df["observed_fdr_ci_hi"].to_numpy()
            n_called = df["n_called"].to_numpy()
            n_total = df["n_total"].to_numpy()
            # Use numpy via polars Series for percentiles to handle NaN safely.
            rows.append({
                "engine": engine,
                "dataset": spec.key,
                "scenario": spec.scenario,
                "k_shuffles": int(spec.k_shuffles),
                "observed_fdr_median": float(pl.Series(fdr).median()),
                "observed_fdr_q1": float(pl.Series(fdr).quantile(0.25, "linear")),
                "observed_fdr_q3": float(pl.Series(fdr).quantile(0.75, "linear")),
                "observed_fdr_ci_lo": float(pl.Series(ci_lo).min()),
                "observed_fdr_ci_hi": float(pl.Series(ci_hi).max()),
                "n_sites_called_median": float(pl.Series(n_called).median()),
                "n_sites_total": int(pl.Series(n_total).max()),
            })
    if not rows:
        return pl.DataFrame(schema={
            "engine": pl.Utf8, "dataset": pl.Utf8, "scenario": pl.Utf8,
            "k_shuffles": pl.Int64,
            "observed_fdr_median": pl.Float64,
            "observed_fdr_q1": pl.Float64,
            "observed_fdr_q3": pl.Float64,
            "observed_fdr_ci_lo": pl.Float64,
            "observed_fdr_ci_hi": pl.Float64,
            "n_sites_called_median": pl.Float64,
            "n_sites_total": pl.Int64,
        })
    return pl.DataFrame(rows).sort(["dataset", "engine"])


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------


def _git_head() -> str:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True,
        )
        return out.strip()
    except (subprocess.SubprocessError, OSError):
        return "unknown"


def write_manifest(
    specs: dict[str, DatasetSpec], timings: list[dict],
    summary: pl.DataFrame, q_thresh: float,
) -> None:
    import epykit as ep  # noqa: PLC0415
    NULL_BASE.mkdir(parents=True, exist_ok=True)
    per_ds_wall: dict[str, float] = {}
    for t in timings:
        per_ds_wall[t["dataset"]] = per_ds_wall.get(t["dataset"], 0.0) \
            + float(t.get("wall_s", 0.0))
    total_wall = sum(per_ds_wall.values())
    lines = [
        "null_calibration -- Phase 4 Task 7 (run_phase4_null_calibration.py)",
        "",
        f"Date           : {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"epykit version : {ep.__version__}",
        "Engine tag     : v0.7.5-phase3-engines-frozen",
        f"Git HEAD       : {_git_head()}",
        f"q-threshold    : {q_thresh}",
        "",
        "Datasets:",
    ]
    for key, spec in specs.items():
        wall = per_ds_wall.get(spec.key, 0.0)
        lines.append(
            f"  {spec.key:<22}  scenario={spec.scenario:<16}  "
            f"engines={','.join(spec.engines):<35}  k={spec.k_shuffles}  "
            f"wall={wall:.1f}s"
        )
    lines += [
        "",
        f"Total wallclock: {total_wall:.1f}s",
        "",
        "Per-(dataset, engine) outcome:",
    ]
    for t in timings:
        status = "ok" if t.get("ok") else "FAIL"
        lines.append(
            f"  {t['dataset']:<22} {t['engine']:<10} {status:<4} "
            f"wall={t.get('wall_s', 0.0):.1f}s"
        )
        if not t.get("ok"):
            lines.append(f"      error: {t.get('error')}")
    lines += [
        "",
        "Summary table (benchmark/data/null_calibration/summary.parquet):",
        f"  rows={summary.height}",
    ]
    for r in summary.iter_rows(named=True):
        lines.append(
            f"  {r['dataset']:<22} {r['engine']:<10} "
            f"median_FDR={r['observed_fdr_median']:.4f}  "
            f"IQR=[{r['observed_fdr_q1']:.4f}, {r['observed_fdr_q3']:.4f}]  "
            f"n_called_med={r['n_sites_called_median']:.0f}/"
            f"{r['n_sites_total']}"
        )
    lines += [
        "",
        "Per-engine bulk parquets:",
        "  benchmark/data/null_calibration/piao_distributed/cov10_3v3/<engine>.parquet",
        "  benchmark/data/null_calibration/simulator/sim_cov10_3v3/<engine>.parquet",
        "  benchmark/data/null_calibration/gse263850/<engine>.parquet",
        "These are gitignored as bulk artefacts; regenerable via",
        "run_phase4_null_calibration.py. summary.parquet + this MANIFEST are tracked.",
        "",
        "Sealed Phase-3 helpers invoked as subprocess (never modified):",
        "  benchmark/scripts/run_null_calibration.py",
        "  benchmark/scripts/_null_engines.py",
        "  benchmark/scripts/wilson_bootstrap_ci.py",
        "",
        "Methodology:",
        "  For each (dataset, engine) cell we relabel the 6-sample cohort",
        "  via k uniformly-random 3v3 assignments (rng=np.random.default_rng(seed)),",
        "  rerun the engine, and record the proportion of CpGs called",
        "  significant at q<{q}. Under a calibrated test with no true",
        "  DMCs in the shuffled design, this proportion ~= q. <<{q} means",
        "  the test is conservative; >>{q} means anti-conservative.".format(q=q_thresh),
        "",
    ]
    MANIFEST_PATH.write_text("\n".join(lines), encoding="utf-8")
    logger.info("wrote %s", MANIFEST_PATH)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset", choices=("piao", "simulator", "gse", "all"),
        default="all",
        help="Restrict to one dataset (default: all). Composes with --skip-*.",
    )
    parser.add_argument("--skip-piao", action="store_true")
    parser.add_argument("--skip-simulator", action="store_true")
    parser.add_argument("--skip-gse", action="store_true")
    parser.add_argument(
        "--engines", nargs="+", default=None,
        help="Override engine list (default: per-dataset default).",
    )
    parser.add_argument(
        "--seed", type=int, default=0,
        help="Seed for the label-shuffle RNG (default 0).",
    )
    parser.add_argument(
        "--q-thresh", type=float, default=0.05,
        help="Nominal q threshold for 'called' (default 0.05).",
    )
    parser.add_argument(
        "--keep-store", action="store_true",
        help="Don't delete methylstore caches after the run "
             "(useful for incremental engine reruns).",
    )
    parser.add_argument(
        "--skip-aggregate", action="store_true",
        help="Don't (re)build summary.parquet -- only run engines.",
    )
    parser.add_argument(
        "--only-aggregate", action="store_true",
        help="Skip all engine runs; just (re)aggregate existing parquets.",
    )
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="DEBUG-level logging.")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        force=True,
    )

    specs = make_specs()
    selected: list[str] = []
    if args.dataset == "all":
        selected = list(specs.keys())
    else:
        selected = [args.dataset]
    if args.skip_piao and "piao" in selected:
        selected.remove("piao")
    if args.skip_simulator and "simulator" in selected:
        selected.remove("simulator")
    if args.skip_gse and "gse" in selected:
        selected.remove("gse")

    timings: list[dict] = []
    if not args.only_aggregate:
        for key in selected:
            spec = specs[key]
            engines = tuple(args.engines) if args.engines else spec.engines
            timings.extend(run_dataset(
                spec, engines=engines,
                seed=args.seed, q_thresh=args.q_thresh,
                keep_store=args.keep_store,
            ))

    n_fail = sum(1 for t in timings if not t.get("ok", True))
    if n_fail:
        logger.error("%d (dataset, engine) cells FAILED", n_fail)

    if not args.skip_aggregate:
        summary = aggregate_summary(specs)
        NULL_BASE.mkdir(parents=True, exist_ok=True)
        SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)
        summary.write_parquet(SUMMARY_PATH)
        logger.info("wrote %s (%d rows)", SUMMARY_PATH, summary.height)
        write_manifest(specs, timings, summary, q_thresh=args.q_thresh)

    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
