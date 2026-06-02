"""F7 — Resource comparison across DMR callers on GSE263850.

Bars for wall time, total CPU time, peak RSS. Four callers:
  methylKit-tile, epykit-tile, epykit-chain_merge, DSS-from-scratch.

Numbers from each pipeline's benchmark CSV / resources.json. The
epykit-chain_merge run wasn't profiled with psutil 1-Hz sampling — RSS
estimated from epykit-tile (same store + DMC backend).
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
CM_LOG    = DATA_DIR / "chain_merge" / "run_log.txt"


def parse_cm_wall_cpu() -> tuple[float, float]:
    """Total wall + summed CPU from the chain_merge log."""
    text = CM_LOG.read_text(encoding="utf-8").splitlines()
    # Use ISO timestamps to bracket the run
    first = last = None
    for line in text:
        if line[:4].isdigit() and " " in line:
            ts = line.split(",", 1)[0]
            if first is None:
                first = ts
            last = ts
    from datetime import datetime
    fmt = "%Y-%m-%d %H:%M:%S"
    dt = (datetime.strptime(last, fmt) -
          datetime.strptime(first, fmt)).total_seconds()
    # CPU we don't have directly; use the DMC + DMR + annotation snippets
    cpu_s = 219.8 + 20.6 + 21.5   # DMC + DMR + (annotate ~21s from log)
    return dt, cpu_s


def main() -> None:
    setup()

    mk = pd.read_csv(MK_BENCH)
    ek = pd.read_csv(EK_BENCH)
    dss = json.loads(DSS_JSON.read_text(encoding="utf-8"))
    cm_wall, cm_cpu = parse_cm_wall_cpu()

    # ---- per-caller summary --------------------------------------------
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
            "wall_s": cm_wall,
            "cpu_s":  cm_cpu,
            # No direct RSS profiling — same backing store + DMC backend
            # as epykit-tile, so we use that as a faithful estimate.
            "rss_mb": float(ek["vmhwm_mb"].max()),
        }),
        ("DSS-from-scratch", {
            # Full run = initial DMLfit+DMLtest (2044+145+75+51+5 ≈ 2320s wall)
            # + resume (499s). Use the sum.
            "wall_s": 75.4 + 50.6 + 2044.5 + 145.4 + 4.7
                       + dss["resources"]["wall_seconds_resume_only"],
            "cpu_s":  67.5 + 49.4 + 2029.0 + 144.3 + 4.5
                       + 14.0 + 4.8 + 86.8 + 361.1,  # from step_timings parts
            "rss_mb": float(dss["resources"]["rss_peak_mb"]),
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
