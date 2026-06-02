"""F6 — Per-sample methylation heatmap of paper-named top-20 DMR genes.

Reproduces (in spirit) the paper's Fig 3B: rows = paper top-10 hyper +
top-10 hypo genes, columns = 6 samples (3 KO + 3 WT). Each cell is the
mean methylation β across the DMR for that sample. The DMR
coordinates come from paper Supp Table 5 columns
SBP009_*_mean / Het-AKAP11-KO-Clone*_mean.

This is "paper-side" data — every cell is from the paper's own table.
That's the cleanest reproduction; we can compare *our* (epykit / DSS)
per-sample betas in a separate figure if needed by reading our per-DMR
beta columns at the same coordinates.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _style import (DATA_DIR, THREE_WAY, setup, save_dual)

PAPER_T5 = Path(r"D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW/"
                r"Paper resources/DMR_total_list.xlsx")


def main() -> None:
    setup()
    paper = pd.read_excel(PAPER_T5).rename(columns={"chr": "chrom"})
    paper["abs_dB"] = paper["diff.meth_mean"].abs()
    paper["gene"] = paper["Gene.Name"].astype(str)

    # Top 10 hyper + top 10 hypo by |diff.meth_mean|, drop dup genes
    hyper = (paper[paper["diff.meth_mean"] > 0]
              .sort_values("abs_dB", ascending=False)
              .drop_duplicates("gene").head(10))
    hypo  = (paper[paper["diff.meth_mean"] < 0]
              .sort_values("abs_dB", ascending=False)
              .drop_duplicates("gene").head(10))
    selection = pd.concat([hyper, hypo], ignore_index=True)

    sample_cols = ["SBP009_1_mean", "SBP009_2_mean", "SBP009_3_mean",
                    "Het-AKAP11-KO-Clone16_mean",
                    "Het-AKAP11-KO-Clone20_mean",
                    "Het-AKAP11-KO-Clone21_mean"]
    short_cols = ["WT_Rep1", "WT_Rep2", "WT_Rep3",
                   "Clone16", "Clone20", "Clone21"]
    selection["paper_coord"] = (
        selection["chrom"] + ":" + selection["start"].astype(str)
        + "-" + selection["end"].astype(str))

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 7.0),
                              constrained_layout=True,
                              gridspec_kw={"width_ratios": [6, 1]})

    # Row labels with direction marker
    row_labels = []
    for _, r in selection.iterrows():
        marker = "▲" if r["diff.meth_mean"] > 0 else "▼"
        row_labels.append(f"{marker} {r['gene']}  ({r['paper_coord']})")

    # Methylation heatmap
    ax = axes[0]
    mat = selection[sample_cols].values.astype(float)
    im = ax.imshow(mat, aspect="auto", cmap="RdYlBu_r",
                     vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(len(short_cols)))
    ax.set_xticklabels(short_cols, rotation=0, fontsize=9)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8.5)
    ax.set_xlabel("Sample")
    ax.set_title("Per-sample methylation β within each DMR\n"
                  "(paper Supp Table 5)", fontsize=10)
    # Group separator
    ax.axvline(2.5, color="black", lw=1.2)
    # Diff color bar
    cax = fig.add_axes([0.12, 0.04, 0.4, 0.018])
    cb = plt.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("Methylation β (0 = unmethylated, 1 = methylated)",
                  fontsize=9)
    cb.ax.tick_params(labelsize=8)

    # diff.meth_mean bar
    ax2 = axes[1]
    diffs = selection["diff.meth_mean"].values
    colors = ["#c0392b" if d > 0 else "#2980b9" for d in diffs]
    bars = ax2.barh(range(len(diffs)), diffs, color=colors,
                       edgecolor="black", linewidth=0.4)
    ax2.axvline(0, color="black", lw=0.8)
    ax2.set_yticks([]); ax2.set_xlabel("Δβ (KO − WT)")
    ax2.set_title("Δβ", fontsize=10)
    ax2.invert_yaxis()
    ax2.set_xlim(-1, 1)
    for spine in ("top", "right"):
        ax2.spines[spine].set_visible(False)
    ax.invert_yaxis()

    fig.suptitle("F6 · Paper Fig 3B reproduction — methylation heatmap of "
                  "top-10 hyper ▲ + top-10 hypo ▼ DMR-associated genes",
                  fontsize=11, y=1.02)
    save_dual(fig, THREE_WAY / "F6_methylation_heatmap")
    plt.close(fig)


if __name__ == "__main__":
    main()
