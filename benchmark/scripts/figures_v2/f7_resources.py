"""F7 — Resource comparison across DMR callers on GSE263850.

Bars for wall time, total CPU time, peak RSS. Four callers:
  methylKit-tile, epykit-tile, epykit-chain_merge, DSS-from-scratch.

This figure reproduces paper Table 5b, which is the **Windows-host**
comparison (the only platform on which all four callers were profiled:
methylKit's `mc.cores` is a no-op on Windows, so it ran single-threaded).
methylKit-tile and epykit-tile come from the Windows step_benchmarks.csv
files; epykit-chain_merge and DSS-from-scratch use the Windows-run
values that Table 5b reports. (The committed Linux-rerun
`dss/resources.json` is a *separate* platform — host pivoine — and is
NOT what Table 5b describes; see paper §4.3 for the Windows-vs-Linux
timing split.)
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

# Windows-host values reported in paper Table 5b for the two callers
# without a committed Windows resources file (the committed DSS
# resources.json is the Linux-pivoine rerun, a different platform).
# Keeping these fixed makes the figure reproduce Table 5b exactly.
DSS_WALL_S_WIN = 2820.0
DSS_CPU_S_WIN  = 2756.0
DSS_RSS_MB_WIN = 9.3 * 1024            # 9.3 GB peak (Windows run)
CM_WALL_S_WIN  = 443.0                 # chain_merge wall, Table 5b
CM_CPU_S_WIN   = 261.9                 # chain_merge CPU, Table 5b


def main() -> None:
    setup()

    mk = pd.read_csv(MK_BENCH)
    ek = pd.read_csv(EK_BENCH)

    # ---- per-caller summary (Windows host, matching Table 5b) ----------
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
            "wall_s": CM_WALL_S_WIN,
            "cpu_s":  CM_CPU_S_WIN,
            # No direct RSS profiling — same backing store + DMC backend
            # as epykit-tile, so we use that as a faithful estimate.
            "rss_mb": float(ek["vmhwm_mb"].max()),
        }),
        ("DSS-from-scratch", {
            "wall_s": DSS_WALL_S_WIN,
            "cpu_s":  DSS_CPU_S_WIN,
            "rss_mb": DSS_RSS_MB_WIN,
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
                  "(15.6 M CpGs, n=6, hg38, Windows host, single-process)",
                  fontsize=11, y=1.04)
    save_dual(fig, THREE_WAY / "F7_resources")
    plt.close(fig)


if __name__ == "__main__":
    main()
