"""F5 — HOMER annotation feature pie (paper Fig 3C reproduction).

Side-by-side pie charts (or stacked-bar) of feature_type distribution
across paper / methylKit-tile / ek-tile / ek-chain_merge-100 /
ek-chain_merge-250 / DSS.

Source: comparisons/annotation_distribution.csv (long form).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _style import (DATA_DIR, THREE_WAY, FEAT_COLOR, setup, save_dual)

FEAT_ORDER = ["promoter-TSS", "5UTR", "exon", "intron", "3UTR",
              "TTS", "non-coding", "intergenic"]


def main() -> None:
    setup()
    df = pd.read_csv(DATA_DIR / "comparisons" / "annotation_distribution.csv")
    callers = [
        "paper-DSS (Supp Table 5)",
        "methylKit-tile",
        "epykit-tile",
        "epykit-chain_merge-100",
        "epykit-chain_merge-250",
        "DSS-from-scratch",
    ]
    short = {
        "paper-DSS (Supp Table 5)": "paper (DSS) — Fig 3C",
        "methylKit-tile":            "methylKit-tile",
        "epykit-tile":               "epykit-tile",
        "epykit-chain_merge-100":    "epykit-chain_merge-100",
        "epykit-chain_merge-250":    "epykit-chain_merge-250",
        "DSS-from-scratch":          "DSS-from-scratch",
    }

    # ---- pie panel ------------------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 8.0), constrained_layout=True)
    axes = axes.ravel()
    for ax, caller in zip(axes, callers):
        sub = df[df["caller"] == caller].set_index("feature_type")
        n = int(sub["n_total"].iloc[0])
        fracs = [sub.loc[f, "fraction"] if f in sub.index else 0.0
                 for f in FEAT_ORDER]
        colors = [FEAT_COLOR[f] for f in FEAT_ORDER]
        # Only show wedges > 0
        used = [(f, c, fr) for f, c, fr in zip(FEAT_ORDER, colors, fracs) if fr > 0]
        wedges, _ = ax.pie(
            [fr for _, _, fr in used],
            colors=[c for _, c, _ in used],
            startangle=90, counterclock=False,
            wedgeprops=dict(edgecolor="white", linewidth=1.0),
        )
        ax.set_title(f"{short[caller]}\n(n = {n:,})", fontsize=10)

    # Single legend at bottom
    legend_handles = [plt.Rectangle((0, 0), 1, 1, color=FEAT_COLOR[f])
                       for f in FEAT_ORDER]
    fig.legend(legend_handles, FEAT_ORDER, loc="lower center",
                ncol=8, frameon=False, fontsize=9,
                bbox_to_anchor=(0.5, -0.04))
    fig.suptitle("F5 · HOMER genomic-feature distribution across callers (paper Fig 3C repro)",
                 fontsize=12, y=1.02)
    save_dual(fig, THREE_WAY / "F5_annotation_pie")
    plt.close(fig)

    # ---- stacked-bar companion -----------------------------------------
    fig2, ax = plt.subplots(figsize=(10.5, 4.5), constrained_layout=True)
    width = 0.7
    bottoms = np.zeros(len(callers))
    pivot = df.pivot(index="caller", columns="feature_type", values="fraction")
    pivot = pivot.reindex(callers).reindex(columns=FEAT_ORDER).fillna(0)
    for f in FEAT_ORDER:
        ax.bar(callers, pivot[f] * 100, bottom=bottoms * 100, label=f,
                color=FEAT_COLOR[f], edgecolor="white", linewidth=0.5,
                width=width)
        bottoms += pivot[f].values
    ax.set_xticks(range(len(callers)))
    ax.set_xticklabels([short[c] for c in callers], rotation=20, ha="right",
                         fontsize=9)
    ax.set_ylabel("% of DMRs in feature type")
    ax.set_ylim(0, 100)
    ax.set_title("F5b · stacked bar (same data)")
    ax.legend(fontsize=8, frameon=False, loc="center left",
               bbox_to_anchor=(1.0, 0.5), ncol=1)
    save_dual(fig2, THREE_WAY / "F5b_annotation_stacked_bar")
    plt.close(fig2)


if __name__ == "__main__":
    main()
