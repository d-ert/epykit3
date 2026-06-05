"""F7 — Resource comparison across DMR callers on GSE263850 (Linux host).

Bars for wall time, total CPU time, peak RSS. Four callers:
  methylKit-tile, epykit-tile, epykit-chain_merge, DSS-from-scratch.

All numbers are the Linux host (pivoine, 24 logical cores) measurements
that match paper Table 5b:
  - methylKit-tile and epykit-tile come from their step_benchmarks.csv
    files (psutil + /proc VmHWM, so Linux); methylKit ran single-threaded
    (mc.cores = 1, explicit on Linux for a fair single-core comparison).
  - epykit-chain_merge is the cached-store DMC+DMR re-call (run_log on
    pivoine; RSS shares epykit-tile's backing store).
  - DSS-from-scratch is the committed Linux rerun: dss/resources.json
    (peak RSS 14.3 GB) + step_timings_resume.tsv (resume wall/cpu) plus
    the ~2,044 s initial DMLfit+DMLtest (dss/summary.md).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _style import (DATA_DIR, THREE_WAY, CALLER_COLOR, setup, save_dual)

MK_BENCH  = Path(r"D:/Coding/Projeler/methyl_lib/methylkıt_realResults/"
                 r"scripts_and_results/methylkit_results/benchmark/"
                 r"step_benchmarks.csv")
EK_BENCH  = Path(r"D:/Coding/Projeler/methyl_lib/benchmarkin_merges/"
                 r"epykit_vs_methylkit(GSE263850)/"
                 r"epykit_results/benchmark/step_benchmarks.csv")
DSS_JSON  = DATA_DIR / "dss" / "resources.json"
DSS_STEPS = DATA_DIR / "dss" / "step_timings_resume.tsv"

# Linux-rerun chain_merge (host pivoine): cached-store DMC+DMR re-call.
# Wall from run_log brackets (~92 s); DMC 72.5 s per dss summary. RSS not
# psutil-profiled — shares epykit-tile's backing store, used as estimate.
CM_WALL_S_LINUX = 92.0
CM_CPU_S_LINUX  = 261.9
# Initial DSS DMLfit+DMLtest wall on pivoine (dss/summary.md); single-
# threaded, so CPU ~= wall for this leg.
DSS_INITIAL_WALL_S = 2044.0


def main() -> None:
    setup()

    mk = pd.read_csv(MK_BENCH)
    ek = pd.read_csv(EK_BENCH)
    dss = json.loads(DSS_JSON.read_text(encoding="utf-8"))
    dss_steps = pd.read_csv(DSS_STEPS, sep="\t")

    # DSS Linux total = initial DMLfit/DMLtest + the resume steps.
    dss_resume_wall = float(dss_steps["wall_seconds"].sum())
    dss_resume_cpu  = float(dss_steps["total_cpu_seconds"].sum())
    dss_wall = DSS_INITIAL_WALL_S + dss_resume_wall
    dss_cpu  = DSS_INITIAL_WALL_S + dss_resume_cpu
    dss_rss  = float(dss["resources"]["rss_peak_mb"])

    # ---- per-caller summary (Linux host pivoine, matching Table 5b) ----
    callers = [
        ("methylKit-tile", {
            "wall_s": float(mk["wall_seconds"].sum()),
            "cpu_s":  float((mk["user_cpu_seconds"]+mk["sys_cpu_seconds"]).sum()),
            "rss_mb": float(mk["vmhwm_mb"].max()),
        }),
        ("epykit-tile", {
            "wall_s": float(ek["wall_seconds"].sum()),
            "cpu_s":  float((ek["user_cpu_seconds"].fillna(0)
                              + ek["sys_cpu_seconds"].fillna(0)).sum()),
            "rss_mb": float(ek["vmhwm_mb"].max()),
        }),
        ("epykit-chain_merge-100", {
            "wall_s": CM_WALL_S_LINUX,
            "cpu_s":  CM_CPU_S_LINUX,
            # No direct RSS profiling — same backing store + DMC backend
            # as epykit-tile, so we use that as a faithful estimate.
            "rss_mb": float(ek["vmhwm_mb"].max()),
        }),
        ("DSS-from-scratch", {
            "wall_s": dss_wall,
            "cpu_s":  dss_cpu,
            "rss_mb": dss_rss,
        }),
    ]
    df_summary = pd.DataFrame([dict(caller=c, **vals) for c, vals in callers])
    df_summary.to_csv(THREE_WAY / "F7_resources_data.csv", index=False)

    labels = [c[0] for c in callers]
    walls  = [c[1]["wall_s"]  for c in callers]
    cpus   = [c[1]["cpu_s"]   for c in callers]
    rsses  = [c[1]["rss_mb"]  for c in callers]
    colors = [CALLER_COLOR[c] for c in labels]

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2), constrained_layout=True)

    # Wall time (log scale)
    ax = axes[0]
    bars = ax.bar(labels, walls, color=colors, edgecolor="black", linewidth=0.6)
    ax.set_yscale("log")
    ax.set_ylabel("Wall time (s, log scale)")
    ax.set_title("A · Pipeline wall time")
    for b, v in zip(bars, walls):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.08,
                 f"{v:,.0f} s", ha="center", va="bottom", fontsize=9)
    ax.set_xticklabels(labels, rotation=18, ha="right", fontsize=9)
    ax.grid(axis="y", alpha=0.2, which="both")

    # CPU time
    ax = axes[1]
    bars = ax.bar(labels, cpus, color=colors, edgecolor="black", linewidth=0.6)
    ax.set_yscale("log")
    ax.set_ylabel("Total CPU time (s, log scale)")
    ax.set_title("B · Total CPU time")
    for b, v in zip(bars, cpus):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.08,
                 f"{v:,.0f} s", ha="center", va="bottom", fontsize=9)
    ax.set_xticklabels(labels, rotation=18, ha="right", fontsize=9)
    ax.grid(axis="y", alpha=0.2, which="both")

    # Peak RSS (linear, GB)
    ax = axes[2]
    bars = ax.bar(labels, [r / 1024 for r in rsses],
                   color=colors, edgecolor="black", linewidth=0.6)
    ax.set_ylabel("Peak RSS (GB)")
    ax.set_title("C · Peak RAM (RSS)")
    for b, v in zip(bars, rsses):
        ax.text(b.get_x() + b.get_width() / 2, (v / 1024) * 1.02,
                 f"{v / 1024:.1f} GB", ha="center", va="bottom", fontsize=9)
    ax.set_xticklabels(labels, rotation=18, ha="right", fontsize=9)
    ax.grid(axis="y", alpha=0.2)

    fig.suptitle("F7 · Resource cost across DMR callers on GSE263850\n"
                  "(15.6 M CpGs, n=6, hg38, Linux host pivoine; methylKit mc.cores=1, "
                  "DSS single-thread)",
                  fontsize=11, y=1.04)
    save_dual(fig, THREE_WAY / "F7_resources")
    plt.close(fig)


if __name__ == "__main__":
    main()
