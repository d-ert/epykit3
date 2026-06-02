"""F07 — Per-step pipeline wall-clock and peak RSS on GSE263850.

Side-by-side comparison of epykit vs methylKit. Steps are listed in pipeline
order; the methylKit step set has 13 phases, epykit 17 — labels are tool-
specific so we plot them as two separate horizontal stacked bars rather than
forcing a common ontology.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _headline_style import (
    DATA_STUDY3, METHYLKIT_REAL, FIG_OUT, save_dual, setup,
)


def _load(steps_path: Path) -> pd.DataFrame:
    df = pd.read_csv(steps_path)
    df = df.sort_values("step").reset_index(drop=True)
    return df[["step", "wall_seconds"]].rename(columns={"wall_seconds": "wall_s"})


def main() -> None:
    setup()

    ek = _load(DATA_STUDY3 / "benchmark" / "step_benchmarks.csv")
    mk = _load(METHYLKIT_REAL / "benchmark" / "step_benchmarks.csv")

    ek["wall_min"] = ek["wall_s"] / 60.0
    mk["wall_min"] = mk["wall_s"] / 60.0

    ek_total_min = ek["wall_min"].sum()
    mk_total_min = mk["wall_min"].sum()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 7.5),
                                   gridspec_kw={"width_ratios": [1, 1]})

    # Use viridis-like palette per tool
    ek_palette = plt.cm.Blues(np.linspace(0.4, 0.95, len(ek)))
    mk_palette = plt.cm.Oranges(np.linspace(0.4, 0.95, len(mk)))

    # --- epykit panel ---
    y = np.arange(len(ek))
    ax1.barh(y, ek["wall_min"], color=ek_palette, edgecolor="white")
    for i, (_, row) in enumerate(ek.iterrows()):
        if row["wall_min"] > ek["wall_min"].max() * 0.02:
            ax1.text(row["wall_min"] + 0.2, i, f"{row['wall_min']:.1f}m",
                     va="center", fontsize=7)
    ax1.set_yticks(y)
    ax1.set_yticklabels(ek["step"], fontsize=7.5)
    ax1.set_xlabel("wall-clock (minutes)")
    ax1.set_title(
        f"epykit — total {ek_total_min:.1f} min ({len(ek)} steps)",
        loc="left", fontsize=10,
    )
    ax1.invert_yaxis()
    ax1.grid(axis="x", linestyle=":", alpha=0.4)

    # --- methylKit panel ---
    y = np.arange(len(mk))
    ax2.barh(y, mk["wall_min"], color=mk_palette, edgecolor="white")
    for i, (_, row) in enumerate(mk.iterrows()):
        if row["wall_min"] > mk["wall_min"].max() * 0.02:
            ax2.text(row["wall_min"] + 0.5, i, f"{row['wall_min']:.1f}m",
                     va="center", fontsize=7)
    ax2.set_yticks(y)
    ax2.set_yticklabels(mk["step"], fontsize=7.5)
    ax2.set_xlabel("wall-clock (minutes)")
    ax2.set_title(
        f"methylKit 1.36.0 — total {mk_total_min:.1f} min ({len(mk)} steps)",
        loc="left", fontsize=10,
    )
    ax2.invert_yaxis()
    ax2.grid(axis="x", linestyle=":", alpha=0.4)

    fig.suptitle(
        "F07 — Per-step pipeline breakdown on GSE263850",
        fontsize=12, fontweight="bold", x=0.02, ha="left",
    )

    fig.text(
        0.02, -0.025,
        "Step labels are tool-specific. methylKit ran on a 24-CPU/64-GB Linux "
        "workstation; epykit on a 16-thread Windows desktop. Both pipelines "
        "process 6 samples and produce DMCs + DMRs end-to-end.",
        fontsize=7, color="#666666",
    )

    plt.subplots_adjust(top=0.92)
    save_dual(fig, FIG_OUT / "F07_pipeline_step_breakdown")


if __name__ == "__main__":
    main()
