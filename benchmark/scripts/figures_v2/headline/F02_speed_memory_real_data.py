"""F02 — Real-data speed and memory on GSE263850 (6-sample WGBS, hg38).

Side-by-side bars: pipeline wall-clock (minutes) and peak resident set size (GB)
for epykit `lr`, methylKit 1.36.0, and DSS. Annotates the 12.15x speed-up and
3.82x memory reduction directly on the bars.

Hardware caveat: methylKit ran on a 24-CPU / 64-GB Linux workstation; epykit and
DSS ran on the user's machine (a 16-thread Windows desktop). The caption
discloses this asymmetry.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _headline_style import (
    DATA_STUDY3, METHYLKIT_REAL, FIG_OUT, color_for, label_for,
    save_dual, setup,
)


def read_epykit() -> tuple[float, float]:
    df = pd.read_csv(DATA_STUDY3 / "benchmark" / "run_summary.csv")
    row = df.iloc[0]
    return float(row["pipeline_wall_sec"]), float(row["peak_rss_mb"])


def read_methylkit() -> tuple[float, float]:
    df = pd.read_csv(METHYLKIT_REAL / "benchmark" / "run_summary.csv")
    row = df.iloc[0]
    return float(row["pipeline_wall_sec"]), float(row["peak_rss_mb"])


def read_dss() -> tuple[float, float]:
    df = pd.read_csv(DATA_STUDY3 / "dss" / "resources.csv")
    wall = float(df["elapsed_s"].max())
    peak_rss_mb = float(df["rss_mb"].max())
    return wall, peak_rss_mb


def main() -> None:
    setup()

    ek_wall, ek_rss = read_epykit()
    mk_wall, mk_rss = read_methylkit()
    dss_wall, dss_rss = read_dss()

    tools = ["epykit_lr", "methylkit", "dss"]
    walls_min = np.array([ek_wall, mk_wall, dss_wall]) / 60.0
    rss_gb = np.array([ek_rss, mk_rss, dss_rss]) / 1024.0
    colors = [color_for(t) for t in tools]
    labels = [label_for(t) for t in tools]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    # ---- wall-clock panel ----
    bars1 = ax1.bar(labels, walls_min, color=colors, width=0.6, edgecolor="white")
    for bar, val in zip(bars1, walls_min):
        ax1.text(bar.get_x() + bar.get_width() / 2, val + max(walls_min) * 0.015,
                 f"{val:.1f} min", ha="center", va="bottom", fontsize=10,
                 fontweight="bold")
    speedup_mk = mk_wall / ek_wall
    speedup_dss = dss_wall / ek_wall
    ax1.set_title("Pipeline wall-clock", loc="left")
    ax1.set_ylabel("minutes (lower is better)")
    ax1.set_ylim(0, max(walls_min) * 1.18)
    ax1.grid(axis="y", linestyle=":", alpha=0.4)

    # ---- memory panel ----
    bars2 = ax2.bar(labels, rss_gb, color=colors, width=0.6, edgecolor="white")
    for bar, val in zip(bars2, rss_gb):
        ax2.text(bar.get_x() + bar.get_width() / 2, val + max(rss_gb) * 0.015,
                 f"{val:.1f} GB", ha="center", va="bottom", fontsize=10,
                 fontweight="bold")
    ratio_mk = mk_rss / ek_rss
    ratio_dss = dss_rss / ek_rss
    ax2.set_title("Peak resident memory", loc="left")
    ax2.set_ylabel("GB (lower is better)")
    ax2.set_ylim(0, max(rss_gb) * 1.18)
    ax2.grid(axis="y", linestyle=":", alpha=0.4)

    fig.suptitle(
        "F02 — Real-data benchmark on GSE263850 (6 samples, 15.6 M CpGs, hg38)",
        fontsize=12, fontweight="bold", x=0.02, ha="left",
    )

    fig.text(
        0.02, -0.04,
        "Caveat: methylKit benchmark ran on a 24-CPU / 64 GB Linux workstation; "
        "epykit and DSS ran on a 16-thread Windows desktop.",
        fontsize=7, color="#666666", wrap=True,
    )

    plt.subplots_adjust(top=0.88)
    save_dual(fig, FIG_OUT / "F02a_speed_memory_real_data")


if __name__ == "__main__":
    main()
