"""F3 — DMR length distributions across all callers.

Violin + box, one row per caller (paper, methylKit-tile, ek-tile,
ek-chain_merge-100, ek-chain_merge-250, DSS). Median lines + count
annotations.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _style import (DATA_DIR, SWEEP_DIR, THREE_WAY, CALLER_COLOR, PAPER_T5_XLSX,
                    MK_TILE_DIR, setup, save_dual)

RAW = PAPER_T5_XLSX.parent.parent
MK  = MK_TILE_DIR / "dmr_significant_lenient.csv"


def length_df(df: pd.DataFrame, label: str) -> pd.DataFrame:
    df = df.copy()
    df["length"] = df["end"].astype(int) - df["start"].astype(int)
    df["caller"] = label
    return df[["length", "caller"]]


def main() -> None:
    setup()
    rows = []
    # paper
    paper = pd.read_excel(RAW / "Paper resources" / "DMR_total_list.xlsx")
    paper = paper.rename(columns={"chr": "chrom"})
    rows.append(length_df(paper, "paper-DSS (Supp Table 5)"))
    # methylKit-tile
    rows.append(length_df(pd.read_csv(MK), "methylKit-tile"))
    # epykit-tile
    rows.append(length_df(
        pd.read_csv(DATA_DIR / "dmr_significant_lenient.csv"),
        "epykit-tile"))
    # ek-chain_merge-100
    ek100 = pl.read_parquet(DATA_DIR / "chain_merge" /
                            "dmr_chain_merge.parquet").to_pandas()
    rows.append(length_df(ek100, "epykit-chain_merge-100"))
    # ek-chain_merge-250
    ek250 = pl.read_parquet(SWEEP_DIR / "dis_merge_250" / "dmr.parquet").to_pandas()
    rows.append(length_df(ek250, "epykit-chain_merge-250"))
    # DSS
    rows.append(length_df(
        pd.read_csv(DATA_DIR / "dss" / "dmr_dss.csv"),
        "DSS-from-scratch"))

    df = pd.concat(rows, ignore_index=True)
    caller_order = [
        "paper-DSS (Supp Table 5)",
        "methylKit-tile",
        "epykit-tile",
        "epykit-chain_merge-100",
        "epykit-chain_merge-250",
        "DSS-from-scratch",
    ]
    df["caller"] = pd.Categorical(df["caller"], categories=caller_order,
                                   ordered=True)

    fig, ax = plt.subplots(figsize=(9.0, 5.5), constrained_layout=True)

    # Cap at 2000 bp for visibility; tail extends to ~2700
    cap = 2000
    df_plot = df.copy()
    df_plot["length_capped"] = df_plot["length"].clip(upper=cap)

    parts = ax.violinplot(
        [df_plot[df_plot["caller"] == c]["length_capped"].values
         for c in caller_order],
        positions=range(len(caller_order)),
        widths=0.8, showmeans=False, showmedians=False,
        showextrema=False,
    )
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(CALLER_COLOR[caller_order[i]])
        pc.set_alpha(0.55)
        pc.set_edgecolor("black")
        pc.set_linewidth(0.7)

    # Box overlay
    bp = ax.boxplot(
        [df_plot[df_plot["caller"] == c]["length_capped"].values
         for c in caller_order],
        positions=range(len(caller_order)),
        widths=0.18, patch_artist=True, showfliers=False,
    )
    for p, c in zip(bp["boxes"], caller_order):
        p.set_facecolor("white"); p.set_edgecolor("black"); p.set_linewidth(0.9)
    for med in bp["medians"]:
        med.set_color("black"); med.set_linewidth(1.4)

    # Count + median annotations
    for i, c in enumerate(caller_order):
        sub = df[df["caller"] == c]["length"]
        ax.text(i, -120, f"n={len(sub):,}", ha="center", va="top",
                fontsize=8.5, color="#34495e")
        ax.text(i, -240, f"med {int(sub.median())} bp", ha="center",
                va="top", fontsize=8, color="#7f8c8d")

    ax.set_xticks(range(len(caller_order)))
    ax.set_xticklabels(["paper\n(DSS)", "methylKit\ntile", "epykit\ntile",
                         "epykit\nchain_merge\n100", "epykit\nchain_merge\n250",
                         "DSS\n(local)"],
                        fontsize=9)
    ax.set_ylabel("DMR length (bp; capped at 2000 for plot)")
    ax.set_ylim(-300, cap + 100)
    ax.set_title("F3 · DMR length distributions across callers\n"
                  "(paper Fig 3A morphology — paper median 239 bp)")
    ax.axhline(239, color="#2c3e50", ls=":", alpha=0.6, lw=1)
    ax.text(5.5, 250, "paper median", ha="right", va="bottom",
            fontsize=8, color="#2c3e50", alpha=0.7)
    ax.grid(axis="y", alpha=0.2)
    save_dual(fig, THREE_WAY / "F3_length_distributions")
    plt.close(fig)


if __name__ == "__main__":
    main()
