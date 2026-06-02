"""F2 — dis.merge sweep curves.

Plots how dis.merge ∈ {100, 150, 200, 250, 500} affects:
- recall vs paper (any-bp, J>=0.25, J>=0.5, J>=0.75)
- precision vs paper
- median DMR length
- panel-E gene capture

Source: chain_merge_dis_merge_sweep/sweep_summary.csv
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _style import (DATA_DIR, THREE_WAY, setup, save_dual)


def main() -> None:
    setup()
    df = pd.read_csv(DATA_DIR / "chain_merge_dis_merge_sweep" /
                     "sweep_summary.csv")

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), constrained_layout=True)

    # ---- Panel A: recall stratified by Jaccard ---------------------------
    ax = axes[0]
    ax.plot(df["dis_merge_bp"], df["recall_anybp"] * 100,
            marker="o", lw=2, label="any-bp overlap", color="#3498db")
    ax.plot(df["dis_merge_bp"], df["recall_J_0_25"] * 100,
            marker="s", lw=2, label="J ≥ 0.25", color="#2980b9")
    ax.plot(df["dis_merge_bp"], df["recall_J_0_5"] * 100,
            marker="^", lw=2, label="J ≥ 0.5", color="#1f618d")
    ax.plot(df["dis_merge_bp"], df["recall_J_0_75"] * 100,
            marker="D", lw=2, label="J ≥ 0.75", color="#154360")
    ax.axvline(100, color="#7f8c8d", ls="--", alpha=0.6, lw=1)
    ax.annotate("paper\ndis.merge=100", xy=(100, 5), xytext=(105, 5),
                fontsize=8, color="#7f8c8d")
    ax.axvline(250, color="#27ae60", ls="--", alpha=0.6, lw=1)
    ax.annotate("morphology-\nmatched (250)", xy=(250, 5), xytext=(255, 5),
                fontsize=8, color="#27ae60")
    ax.set_xlabel("dis.merge (bp)")
    ax.set_ylabel("Paper-DMR recall (%)")
    ax.set_title("A · Recall vs paper Supp Table 5\n(813 DMRs, stratified by Jaccard)")
    ax.set_ylim(0, 80)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.grid(alpha=0.2)

    # ---- Panel B: precision + DMR count + median length ------------------
    ax = axes[1]
    ax2 = ax.twinx()
    ax.plot(df["dis_merge_bp"], df["precision_anybp"] * 100,
            marker="o", color="#e74c3c", lw=2, label="precision (any bp)")
    ax.set_xlabel("dis.merge (bp)")
    ax.set_ylabel("Precision (%)", color="#e74c3c")
    ax.tick_params(axis="y", labelcolor="#e74c3c")
    ax.set_ylim(40, 80)
    ax2.plot(df["dis_merge_bp"], df["n_dmr"],
             marker="s", color="#34495e", lw=2, label="DMR count", ls="--")
    ax2.set_ylabel("DMR count", color="#34495e")
    ax2.tick_params(axis="y", labelcolor="#34495e")
    ax2.axhline(813, color="#7f8c8d", ls=":", alpha=0.5)
    ax2.annotate("paper 813", xy=(500, 813), xytext=(380, 830),
                 fontsize=8, color="#7f8c8d")
    ax.set_title("B · Precision and DMR count\n(precision–recall trade-off)")
    ax.grid(alpha=0.2)

    # ---- Panel C: morphology + gene captures -----------------------------
    ax = axes[2]
    ax2 = ax.twinx()
    ax.plot(df["dis_merge_bp"], df["median_length_bp"],
            marker="o", color="#9b59b6", lw=2, label="median length")
    ax.set_xlabel("dis.merge (bp)")
    ax.set_ylabel("Median DMR length (bp)", color="#9b59b6")
    ax.tick_params(axis="y", labelcolor="#9b59b6")
    ax.axhline(240, color="#7f8c8d", ls=":", alpha=0.5)
    ax.annotate("paper 240 bp", xy=(500, 240), xytext=(360, 250),
                fontsize=8, color="#7f8c8d")
    ax2.plot(df["dis_merge_bp"], df["panel_e_recall_nearest_tss"] * 100,
             marker="s", color="#16a085", lw=2,
             label="Panel-E gene capture", ls="--")
    ax2.set_ylabel("Panel-E gene capture (%)", color="#16a085")
    ax2.tick_params(axis="y", labelcolor="#16a085")
    ax2.set_ylim(55, 75)
    ax.set_title("C · DMR morphology and Panel-E capture\n(Supp Table 8, 46 genes)")
    ax.grid(alpha=0.2)

    fig.suptitle("F2 · dis.merge sweep — paper $D_{is.merge}$=100 vs morphology-matched 250",
                 fontsize=12, y=1.04)
    save_dual(fig, THREE_WAY / "F2_dis_merge_sweep")
    plt.close(fig)


if __name__ == "__main__":
    main()
