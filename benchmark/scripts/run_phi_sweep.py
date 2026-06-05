"""run_phi_sweep.py -- dispersion (phi) sweep on the intrinsic-truth simulator.

Closes review issue M1 (the phi-sweep that was scripted but never run) and
M2 (dmrseq / BSmooth locally on the simulator), and captures **per-tool peak
RSS** for every tool (review/user request -- the simulator R runners only
recorded wall/CPU, never memory).

For each (rho, seed) cell it:
  1. simulates a Beta-Binomial cell at intraclass correlation ``rho`` (the
     simulator's ``--phi``; rho=0 is pure Binomial / Piao original),
  2. converts AMP -> bismark .cov.gz (shared with the R tools),
  3. runs each tool as a *monitored subprocess* (peak RSS + wall + CPU):
       epykit_lr, epykit_lrplus  (per-CpG DMC)
       methylkit (mc.cores=8), dss (smooth), dss_nosmooth  (per-CpG DMC)
       dmrseq, bsmooth  (region/DMR callers)
  4. scores every tool against the per-CpG intrinsic truth:
       - DMC tools via the canonical score_dmc_parquet (q<0.05, all-bins),
       - DMR tools by predicting every CpG inside a significant called region
         as positive (region -> member-CpG), then a per-CpG confusion matrix.

Each cell emits one row per tool. We report both ``rho`` (the simulator ICC)
and the implied Pearson overdispersion ``phi_pearson = 1 + (coverage-1)*rho``
at the test coverage, resolving the rho-vs-Pearson-phi unit ambiguity the
review flagged.

Outputs (under benchmark/data/study1b_simulator/):
  - eval_phi_sweep_per_cell.parquet  (one row per rho x seed x tool)
  - eval_phi_sweep_iqr.parquet       (median + IQR across seeds per rho x tool)
  - phi_sweep/phi=<rho>/seed=<seed>/ ... (gitignored bulk cell artefacts)

Usage:
  python run_phi_sweep.py                      # default 6 rho x 10 seeds
  python run_phi_sweep.py --rhos 0 0.1 --seeds 2026000 2026001
  python run_phi_sweep.py --skip-existing      # resume (skip tools w/ output)
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from _resource_monitor import run_subprocess_monitored  # noqa: E402
from _epykit_scoring import Q_THRESHOLD, score_dmc_parquet  # noqa: E402
from eval_simulator_intrinsic import _load_methylkit, _load_dss  # noqa: E402
from run_external_simulator_sweep import convert_amp_to_cov  # noqa: E402

logger = logging.getLogger("run_phi_sweep")

ROOT = Path(__file__).resolve().parents[2]
SIM_ROOT = ROOT / "benchmark" / "data" / "study1b_simulator"
SWEEP_ROOT = SIM_ROOT / "phi_sweep"
CACHE_ROOT = ROOT / "benchmark" / "_phi_sweep_cache"

DEFAULT_RHOS = (0.0, 0.05, 0.1, 0.2, 0.3, 0.44)
DEFAULT_SEEDS = tuple(range(2026_000, 2026_010))  # 10 seeds
DMC_TOOLS = ("epykit_lr", "epykit_lrplus", "methylkit", "dss", "dss_nosmooth")
DMR_TOOLS = ("dmrseq", "bsmooth")
COV_GLOB = "amp.coverage=*.sample*.cov.gz"


# --------------------------------------------------------------------------- #
# Per-tool runners (each returns the resource dict from the monitor)
# --------------------------------------------------------------------------- #

def _py() -> str:
    return sys.executable


def _run_tool(tool: str, cell: Path, coverage: int, interval: float) -> dict:
    """Run one tool on one cell as a monitored subprocess. Returns the
    resource dict augmented with the output path under key 'out'."""
    in_dir = cell / "bismark_cov"
    if tool == "epykit_lr":
        out = cell / "epykit_lr.parquet"
        cmd = [_py(), str(_HERE / "run_epykit_cell.py"), "--in-dir", in_dir,
               "--engine", "lr", "--out", out,
               "--store-dir", CACHE_ROOT / cell.parent.name / cell.name / "store_lr",
               "--glob", COV_GLOB]
    elif tool == "epykit_lrplus":
        out = cell / "epykit_lrplus.parquet"
        cmd = [_py(), str(_HERE / "run_epykit_cell.py"), "--in-dir", in_dir,
               "--engine", "lr+", "--out", out,
               "--store-dir", CACHE_ROOT / cell.parent.name / cell.name / "store_lrplus",
               "--glob", COV_GLOB]
    elif tool == "methylkit":
        out = cell / "methylkit.tsv"
        cmd = ["Rscript", str(_HERE / "run_methylkit_simulator.R"),
               "--in-dir", in_dir, "--out", out, "--cores", "8"]
    elif tool == "dss":
        out = cell / "dss.tsv"
        cmd = ["Rscript", str(_HERE / "run_dss_simulator.R"),
               "--in-dir", in_dir, "--out", out, "--smoothing", "TRUE"]
    elif tool == "dss_nosmooth":
        out = cell / "dss_nosmooth.tsv"
        cmd = ["Rscript", str(_HERE / "run_dss_simulator.R"),
               "--in-dir", in_dir, "--out", out, "--smoothing", "FALSE"]
    elif tool == "dmrseq":
        out = cell / "dmrseq.tsv"
        cmd = ["Rscript", str(_HERE / "run_dmrseq.R"), "--in-dir", in_dir, "--out", out]
    elif tool == "bsmooth":
        out = cell / "bsmooth.tsv"
        cmd = ["Rscript", str(_HERE / "run_bsmooth.R"), "--in-dir", in_dir, "--out", out]
    else:
        raise ValueError(f"unknown tool {tool!r}")
    res = run_subprocess_monitored([str(c) for c in cmd], interval=interval)
    res["out"] = out
    return res


# --------------------------------------------------------------------------- #
# Scoring
# --------------------------------------------------------------------------- #

def _score_dmc(out_pq_or_df, truth: pl.DataFrame, tool: str, coverage: int) -> dict:
    """Score a per-CpG DMC tool at the headline cell (q<0.05, all bins)."""
    import tempfile
    if isinstance(out_pq_or_df, pl.DataFrame):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td) / "x.parquet"
            out_pq_or_df.write_parquet(tmp)
            rows = score_dmc_parquet(tmp, truth, tool=tool, scenario="phi_sweep",
                                     parameter="rho", parameter_value=coverage,
                                     test=tool.replace("epykit_", ""))
    else:
        rows = score_dmc_parquet(out_pq_or_df, truth, tool=tool, scenario="phi_sweep",
                                 parameter="rho", parameter_value=coverage,
                                 test=tool.replace("epykit_", ""))
    for r in rows:
        if (r["threshold_kind"] == "qvalue" and abs(r["threshold"] - Q_THRESHOLD) < 1e-9
                and r["meth_diff_bin"] == "all"):
            return r
    return {}


def _score_dmr_as_percpg(tsv: Path, truth: pl.DataFrame, tool: str,
                         q_threshold: float = Q_THRESHOLD) -> dict:
    """Score a region caller against per-CpG truth: every CpG inside a called
    region is predicted positive, then a per-CpG confusion matrix vs is_dmc.

    NOTE: the intrinsic simulator scatters DMCs per-CpG; it has no contiguous
    reference DMRs. Region callers are therefore evaluated here on a per-CpG
    basis with that caveat.

    Significance convention differs by tool, deliberately:
      * dmrseq emits a real BH ``qvalue`` -> filter qvalue < q_threshold.
      * BSmooth (dmrFinder) returns regions already thresholded on |t-stat| and
        has *no native FDR*; its ``qvalue`` is only a rank-of-|t-stat| surrogate
        (run_bsmooth.R). Filtering it by q<0.05 would keep just the single
        top-ranked region. We therefore treat every BSmooth output region as a
        call -- and flag in the manuscript that BSmooth's calls are t-stat
        thresholded, not FDR-controlled (review M2 caveat).
    """
    if not tsv.exists():
        return {}
    dmr = pl.read_csv(tsv, separator="\t")
    if tool == "dmrseq" and "qvalue" in dmr.columns:
        sig = dmr.filter(pl.col("qvalue") < q_threshold)
    else:  # bsmooth: dmrFinder output is already the call set (no native FDR)
        sig = dmr
    t = truth.select(["chrom", "pos", "is_dmc"])
    covered = np.zeros(t.height, dtype=bool)
    pos = t["pos"].to_numpy()
    chrom = t["chrom"].to_numpy()
    n_sig = sig.height
    for row in sig.iter_rows(named=True):
        c = str(row["chr"]) if "chr" in row else str(row.get("chrom"))
        s, e = int(row["start"]), int(row["end"])
        covered |= (chrom == c) & (pos >= s) & (pos <= e)
    is_dmc = t["is_dmc"].to_numpy().astype(bool)
    tp = int(np.sum(covered & is_dmc))
    fp = int(np.sum(covered & ~is_dmc))
    fn = int(np.sum(~covered & is_dmc))
    tn = int(np.sum(~covered & ~is_dmc))
    tpr = tp / (tp + fn) if (tp + fn) else float("nan")
    fpr = fp / (fp + tn) if (fp + tn) else float("nan")
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    f1 = (2 * precision * tpr / (precision + tpr)
          if precision and tpr and not np.isnan(precision) and not np.isnan(tpr) else float("nan"))
    return dict(tool=tool, tp=tp, fp=fp, tn=tn, fn=fn, tpr=tpr, fpr=fpr,
                precision=precision, f1=f1, auroc=float("nan"),
                n_called=tp + fp, n_sig_regions=n_sig, threshold=q_threshold,
                threshold_kind="dmr_qvalue_to_percpg", meth_diff_bin="all")


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #

def run_cell(rho: float, seed: int, coverage: int, tools: list[str],
             interval: float, skip_existing: bool) -> list[dict]:
    cell = SWEEP_ROOT / f"phi={rho}" / f"seed={seed}"
    cell.mkdir(parents=True, exist_ok=True)
    truth_pq = cell / "truth.parquet"

    # 1. simulate (Beta-Binomial at rho)
    if not truth_pq.exists() or not list(cell.glob(f"amp.coverage={coverage}.sample*.txt")):
        rc = run_subprocess_monitored(
            [_py(), str(_HERE / "simulate_piao.py"), "--seed", str(seed),
             "--phi", str(rho), "--coverage", str(coverage), "--out", str(cell)],
            interval=interval, stdout=None, stderr=None)["returncode"]
        if rc != 0:
            logger.error("simulate failed rho=%s seed=%s", rho, seed)
            return []
    truth = pl.read_parquet(truth_pq)

    # 2. convert AMP -> bismark_cov
    convert_amp_to_cov(cell, coverage)

    phi_pearson = 1.0 + (coverage - 1) * rho
    rows: list[dict] = []
    for tool in tools:
        # resume: skip if output already present
        out_name = {"epykit_lr": "epykit_lr.parquet",
                    "epykit_lrplus": "epykit_lrplus.parquet",
                    "methylkit": "methylkit.tsv", "dss": "dss.tsv",
                    "dss_nosmooth": "dss_nosmooth.tsv", "dmrseq": "dmrseq.tsv",
                    "bsmooth": "bsmooth.tsv"}[tool]
        out_path = cell / out_name
        res_path = cell / f"{tool}.resources.json"
        if skip_existing and out_path.exists() and res_path.exists():
            res = json.loads(res_path.read_text())
            logger.info("[rho=%s seed=%s] %s: cached", rho, seed, tool)
        else:
            t0 = time.time()
            res = _run_tool(tool, cell, coverage, interval)
            logger.info("[rho=%s seed=%s] %s: rc=%s wall=%.1fs rss_peak=%sMB (%.0fs)",
                        rho, seed, tool, res["returncode"], res["wall_s"],
                        res.get("rss_peak_mb"), time.time() - t0)
            res.pop("samples", None)
            res["out"] = str(out_path)
            res_path.write_text(json.dumps(res, indent=2))
            if res["returncode"] != 0:
                logger.error("[rho=%s seed=%s] %s FAILED (rc=%s)",
                             rho, seed, tool, res["returncode"])

        # 3. score
        score: dict = {}
        try:
            if tool in ("epykit_lr", "epykit_lrplus"):
                if out_path.exists():
                    score = _score_dmc(out_path, truth, tool, coverage)
            elif tool == "methylkit":
                if out_path.exists():
                    score = _score_dmc(_load_methylkit(out_path), truth, "methylkit", coverage)
            elif tool in ("dss", "dss_nosmooth"):
                if out_path.exists():
                    score = _score_dmc(_load_dss(out_path), truth, tool, coverage)
            elif tool in DMR_TOOLS:
                score = _score_dmr_as_percpg(out_path, truth, tool)
        except Exception as exc:  # scoring must never abort the sweep
            logger.exception("[rho=%s seed=%s] scoring %s failed: %r", rho, seed, tool, exc)

        tp_, fp_ = score.get("tp"), score.get("fp")
        n_called = score.get("n_called")
        if n_called is None and tp_ is not None and fp_ is not None:
            n_called = tp_ + fp_  # DMC tools: score_dmc_parquet emits tp/fp, not n_called
        rows.append({
            "rho": rho, "phi_pearson": round(phi_pearson, 3), "seed": seed,
            "coverage": coverage, "tool": tool,
            "kind": "dmr" if tool in DMR_TOOLS else "dmc",
            "tpr": score.get("tpr"), "fpr": score.get("fpr"),
            "precision": score.get("precision"), "f1": score.get("f1"),
            "auroc": score.get("auroc"),
            "tp": tp_, "fp": fp_,
            "tn": score.get("tn"), "fn": score.get("fn"),
            "n_called": n_called,
            "wall_s": res.get("wall_s"), "cpu_s": res.get("cpu_s"),
            "rss_peak_mb": res.get("rss_peak_mb"), "rss_mean_mb": res.get("rss_mean_mb"),
            "uss_peak_mb": res.get("uss_peak_mb"),
            "cpu_percent_peak": res.get("cpu_percent_peak"),
            "num_processes_peak": res.get("num_processes_peak"),
            "returncode": res.get("returncode"),
        })
    # free the per-cell epykit methylstore scratch (bulk, regenerable)
    shutil.rmtree(CACHE_ROOT / cell.parent.name / cell.name, ignore_errors=True)
    return rows


def aggregate(per_cell: pl.DataFrame) -> pl.DataFrame:
    num = ["tpr", "fpr", "precision", "f1", "auroc", "rss_peak_mb", "wall_s"]
    aggs = [pl.len().alias("n_seeds")]
    for c in num:
        aggs += [pl.col(c).median().alias(f"{c}_median"),
                 pl.col(c).quantile(0.25).alias(f"{c}_q1"),
                 pl.col(c).quantile(0.75).alias(f"{c}_q3")]
    return (per_cell.group_by(["tool", "rho", "phi_pearson", "kind"])
            .agg(aggs).sort(["tool", "rho"]))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rhos", nargs="+", type=float, default=list(DEFAULT_RHOS))
    ap.add_argument("--seeds", nargs="+", type=int, default=list(DEFAULT_SEEDS))
    ap.add_argument("--coverage", type=int, default=10)
    ap.add_argument("--tools", nargs="+", default=list(DMC_TOOLS + DMR_TOOLS))
    ap.add_argument("--interval", type=float, default=0.2, help="RSS sample interval (s)")
    ap.add_argument("--skip-existing", action="store_true",
                    help="Skip a tool if its output + resources.json already exist")
    ap.add_argument("--verbose", "-v", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s", force=True)

    SWEEP_ROOT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    all_rows: list[dict] = []
    n_cells = len(args.rhos) * len(args.seeds)
    k = 0
    for rho in args.rhos:
        for seed in args.seeds:
            k += 1
            logger.info("=== cell %d/%d: rho=%s seed=%s ===", k, n_cells, rho, seed)
            all_rows.extend(run_cell(rho, seed, args.coverage, args.tools,
                                     args.interval, args.skip_existing))
            # checkpoint after every cell so a crash keeps progress
            pl.DataFrame(all_rows).write_parquet(SIM_ROOT / "eval_phi_sweep_per_cell.parquet")

    per_cell = pl.DataFrame(all_rows)
    per_cell.write_parquet(SIM_ROOT / "eval_phi_sweep_per_cell.parquet")
    iqr = aggregate(per_cell.filter(pl.col("returncode") == 0))
    iqr.write_parquet(SIM_ROOT / "eval_phi_sweep_iqr.parquet")
    logger.info("wrote eval_phi_sweep_{per_cell,iqr}.parquet (%d cells, %d rows, %.0f min)",
                n_cells, per_cell.height, (time.time() - t0) / 60)
    n_fail = int((per_cell["returncode"] != 0).sum())
    if n_fail:
        logger.warning("%d tool runs returned non-zero", n_fail)
    return 0


if __name__ == "__main__":
    sys.exit(main())
