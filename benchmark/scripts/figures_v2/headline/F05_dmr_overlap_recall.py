"""F05 — DMR recovery on GSE263850.

Three-tool comparison: how well do epykit chain_merge / methylKit / DSS each
recover the DMRs called by (a) the original paper's published call set and
(b) each other. Direction-agreement on matched DMRs is 100% across all
comparisons — the disagreement is in WHICH DMRs are called, not in the SIGN
of the methylation change.

Source: data/study3/comparisons_post_phase3/dmr_iou.parquet
        data/study3/comparisons/chain_merge_vs_paper/headline.json
        data/study3/comparisons/epykit_vs_dss/headline.json
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from _headline_style import DATA_STUDY3, FIG_OUT, save_dual, setup


def main() -> None:
    setup()

    cmp_dir = DATA_STUDY3 / "comparisons"

    chain_vs_paper = json.loads(
        (cmp_dir / "chain_merge_vs_paper" / "headline.json").read_text()
    )
    ek_vs_dss = json.loads(
        (cmp_dir / "epykit_vs_dss" / "headline.json").read_text()
    )
    iou = pl.read_parquet(
        DATA_STUDY3 / "comparisons_post_phase3" / "dmr_iou.parquet"
    ).to_pandas()

    # ----- Panel A: recall of paper's published 813 DMRs by each tool -----
    paper_n = chain_vs_paper["input_counts"]["paper_dmrs"]
    paper_recovered = {
        "epykit chain_merge (100 bp)": chain_vs_paper["coord_overlap"]["paper_DMRs_with_any_overlap"],
    }
    # Direction-agreement fraction on matched DMRs (always 1.0 for chain_merge vs paper)
    direction_agree = chain_vs_paper["coord_overlap"]["direction_agreement_on_matched"]

    # ----- Panel B: pairwise recall + precision + Jaccard from dmr_iou -----
    pair_rows = []
    for _, row in iou.iterrows():
        ta, tb = row["tool_a"], row["tool_b"]
        n_a, n_b = int(row["n_a"]), int(row["n_b"])
        n_pairs = int(row["n_overlapping_pairs"])
        # symmetric pair Jaccard from dmr_iou is bp-level, not DMR-count level
        jacc = float(row["jaccard"])
        dir_frac = float(row["direction_agree_frac"])
        # For DMR-count recall: n_pairs / n_a means "% of A matched to ≥ 1 B"
        recall_a = n_pairs / n_a if n_a else 0.0
        recall_b = n_pairs / n_b if n_b else 0.0
        pair_rows.append({
            "pair":   f"{_pretty(ta)} ↔ {_pretty(tb)}",
            "n_a":    n_a,
            "n_b":    n_b,
            "matched_pairs": n_pairs,
            "matched_a_frac": recall_a,
            "matched_b_frac": recall_b,
            "jaccard_bp":     jacc,
            "direction_agree": dir_frac,
        })

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5),
                                   gridspec_kw={"width_ratios": [0.9, 1.4]})

    # ---- Panel A: paper recall bar
    pct = paper_recovered["epykit chain_merge (100 bp)"] / paper_n * 100
    bars = ax1.bar(
        ["paper-DSS (Supp Table 5)\n813 published DMRs",
         "recovered by epykit\nchain_merge (100 bp)"],
        [paper_n, paper_recovered["epykit chain_merge (100 bp)"]],
        color=["#2c3e50", "#0F4C81"],
        width=0.55,
    )
    for bar, val in zip(bars, [paper_n, paper_recovered["epykit chain_merge (100 bp)"]]):
        ax1.text(bar.get_x() + bar.get_width() / 2,
                 val + paper_n * 0.015,
                 f"n = {val}", ha="center", va="bottom", fontsize=10,
                 fontweight="bold")
    ax1.set_ylim(0, paper_n * 1.18)
    ax1.set_ylabel("DMRs")
    ax1.set_title(
        f"Recall of Farhangdoost et al. 2024 paper DMRs\n"
        f"any-bp recall = {pct:.1f}%, direction agreement = {direction_agree:.0%}",
        loc="left", fontsize=10,
    )

    # ---- Panel B: pairwise overlap stats
    ax2.set_title(
        "Pairwise DMR-set overlap (this work)\n"
        "matched = % of one tool's DMRs with ≥ 1 bp overlap in the other",
        loc="left", fontsize=10,
    )
    pairs = [r["pair"] for r in pair_rows]
    y = np.arange(len(pairs))
    width = 0.36
    matched_a = [r["matched_a_frac"] for r in pair_rows]
    matched_b = [r["matched_b_frac"] for r in pair_rows]

    bars_a = ax2.barh(y - width / 2, matched_a, height=width,
                      color="#0F4C81", label="matched (A side)")
    bars_b = ax2.barh(y + width / 2, matched_b, height=width,
                      color="#E07B39", label="matched (B side)")
    for bar, val in zip(bars_a, matched_a):
        ax2.text(val + 0.01, bar.get_y() + bar.get_height() / 2,
                 f"{val:.1%}", va="center", fontsize=8, color="#0F4C81")
    for bar, val in zip(bars_b, matched_b):
        ax2.text(val + 0.01, bar.get_y() + bar.get_height() / 2,
                 f"{val:.1%}", va="center", fontsize=8, color="#E07B39")

    ax2.set_yticks(y)
    ax2.set_yticklabels(pairs, fontsize=9)
    ax2.set_xlim(0, 1.1)
    ax2.set_xlabel("fraction of DMRs matched to the other tool")
    ax2.legend(loc="lower right", frameon=False, fontsize=8)
    ax2.grid(axis="x", linestyle=":", alpha=0.4)

    fig.suptitle(
        "F05 — DMR set agreement across tools on GSE263850",
        fontsize=12, fontweight="bold", x=0.02, ha="left",
    )

    fig.text(
        0.02, -0.04,
        "Tools agree on DIRECTION 100% of the time on matched DMRs — disagreement "
        "is over which intervals to call, not the sign of the methylation change. "
        "methylKit emits ~50× more DMRs than DSS or epykit chain_merge here, "
        "which inflates its 'matched A' but dilutes its 'matched B' fraction.",
        fontsize=7, color="#666666", wrap=True,
    )

    plt.subplots_adjust(top=0.86, wspace=0.45)
    save_dual(fig, FIG_OUT / "F05_dmr_overlap_recall")


def _pretty(t: str) -> str:
    return {
        "epykit_chain_merge": "epykit chain_merge",
        "methylkit":          "methylKit",
        "dss":                "DSS",
    }.get(t, t)


if __name__ == "__main__":
    main()
