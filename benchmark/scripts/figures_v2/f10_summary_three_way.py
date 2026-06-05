"""F10 — Summary composite figure (replaces S4_dmr_engine_choice).

6 panels in a 2×3 grid:
  A. Recall vs paper Supp Table 5 (any-bp + J>=0.5), 4 callers as bars
  B. dis.merge sweep curves (recall stratified by Jaccard)
  C. DMR length violin distribution across callers
  D. Paper-named gene hits (top 20) — bar chart, count of any-bp hits
  E. HOMER feature distribution stacked-bar
  F. Resource bars (wall time, log-y)

All on one PNG/SVG for deck use.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _style import (DATA_DIR, SWEEP_DIR, THREE_WAY, CALLER_COLOR, FEAT_COLOR,
                     PAPER_T5_XLSX, MK_TILE_DIR, setup, save_dual)

PAPER_T5 = PAPER_T5_XLSX
MK_TILE  = MK_TILE_DIR / "dmr_significant_lenient.csv"

FEAT_ORDER = ["promoter-TSS", "5UTR", "exon", "intron", "3UTR",
              "TTS", "non-coding", "intergenic"]


def main() -> None:
    setup()

    fig, axes = plt.subplots(2, 3, figsize=(16.0, 9.5), constrained_layout=True)

    # ---- A: Recall + precision bars (vs DSS-from-scratch, 922) ----------
    ax = axes[0, 0]
    callers = ["methylKit-tile", "epykit-chain_merge-100",
                "epykit-chain_merge-250", "DSS-from-scratch"]
    short_l = ["mk-tile", "ek-cm-100", "ek-cm-250", "DSS"]
    # Post-rerun coordinate concordance vs the locally-rerun DSS-922 call
    # set (headline.json + dis_merge_vs_dss_sensitivity.csv). DSS is the
    # reference set, so its self-recall/precision is 1.0.
    recall_anybp    = [0.087, 0.638, 0.773, 1.0]
    precision_anybp = [0.030, 0.744, 0.630, 1.0]
    j5_recall       = [0.0,   0.345, 0.642, 1.0]
    width = 0.27
    xs = np.arange(len(callers))
    bars1 = ax.bar(xs - width, [r * 100 for r in recall_anybp],
                    width=width, label="recall any-bp",
                    color="#3498db", edgecolor="black", linewidth=0.5)
    bars2 = ax.bar(xs,         [r * 100 for r in precision_anybp],
                    width=width, label="precision any-bp",
                    color="#e74c3c", edgecolor="black", linewidth=0.5)
    bars3 = ax.bar(xs + width, [r * 100 for r in j5_recall],
                    width=width, label="recall J ≥ 0.5",
                    color="#2ecc71", edgecolor="black", linewidth=0.5)
    ax.set_xticks(xs); ax.set_xticklabels(short_l, fontsize=9)
    ax.set_ylabel("%"); ax.set_ylim(0, 105)
    ax.set_title("A · Coordinate concordance vs DSS-922 (DSS = reference)")
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.grid(axis="y", alpha=0.2)

    # ---- B: dis.merge sweep ---------------------------------------------
    ax = axes[0, 1]
    df = pd.read_csv(SWEEP_DIR / "sweep_summary.csv")
    ax.plot(df["dis_merge_bp"], df["recall_anybp"] * 100,
             marker="o", lw=2, label="any-bp", color="#3498db")
    ax.plot(df["dis_merge_bp"], df["recall_J_0_5"] * 100,
             marker="^", lw=2, label="J ≥ 0.5", color="#1f618d")
    ax.plot(df["dis_merge_bp"], df["precision_anybp"] * 100,
             marker="s", lw=2, label="precision", color="#e74c3c")
    ax.axvline(100, color="#7f8c8d", ls="--", alpha=0.5)
    ax.axvline(250, color="#27ae60", ls="--", alpha=0.5)
    ax.set_xlabel("dis.merge (bp)"); ax.set_ylabel("%")
    ax.set_title("B · dis.merge sweep")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.2)

    # ---- C: length violins ---------------------------------------------
    ax = axes[0, 2]
    paper = pd.read_excel(PAPER_T5).rename(columns={"chr": "chrom"})
    paper["len"] = paper["end"] - paper["start"]
    mk = pd.read_csv(MK_TILE)
    mk["len"] = mk["end"] - mk["start"]
    ek100 = pl.read_parquet(DATA_DIR / "chain_merge" /
                              "dmr_chain_merge.parquet").to_pandas()
    ek100["len"] = ek100["end"] - ek100["start"]
    ek250 = pl.read_parquet(SWEEP_DIR / "dis_merge_250" / "dmr.parquet").to_pandas()
    ek250["len"] = ek250["end"] - ek250["start"]
    dss = pd.read_csv(DATA_DIR / "dss" / "dmr_dss.csv")
    dss["len"] = dss["end"] - dss["start"]

    series = [
        ("paper", paper["len"], CALLER_COLOR["paper"]),
        ("mk-tile", mk["len"], CALLER_COLOR["methylKit-tile"]),
        ("ek-cm-100", ek100["len"], CALLER_COLOR["epykit-chain_merge-100"]),
        ("ek-cm-250", ek250["len"], CALLER_COLOR["epykit-chain_merge-250"]),
        ("DSS", dss["len"], CALLER_COLOR["DSS-from-scratch"]),
    ]
    cap = 1500
    parts = ax.violinplot(
        [s[1].clip(upper=cap).values for s in series],
        positions=range(len(series)),
        widths=0.85, showmeans=False, showmedians=False, showextrema=False)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(series[i][2]); pc.set_alpha(0.6)
        pc.set_edgecolor("black"); pc.set_linewidth(0.6)
    bp = ax.boxplot(
        [s[1].clip(upper=cap).values for s in series],
        positions=range(len(series)),
        widths=0.15, patch_artist=True, showfliers=False)
    for p in bp["boxes"]:
        p.set_facecolor("white")
    for med in bp["medians"]:
        med.set_color("black"); med.set_linewidth(1.2)
    ax.set_xticks(range(len(series)))
    ax.set_xticklabels([s[0] for s in series], fontsize=9, rotation=0)
    ax.set_ylabel("DMR length (bp; capped 1500)")
    ax.set_title("C · DMR length distributions")
    ax.axhline(239, color="black", ls=":", alpha=0.5)
    ax.text(4.6, 250, "paper med", ha="right", fontsize=7,
              color="#2c3e50", alpha=0.8)

    # ---- D: paper-named gene hits ---------------------------------------
    ax = axes[1, 0]
    f4 = pd.read_csv(THREE_WAY / "F4_top_named_gene_hits_data.csv")
    callers_d = [("mk-tile", "jaccard_methylKit_tile"),
                  ("ek-cm-100", "jaccard_ek_chain_merge_100"),
                  ("ek-cm-250", "jaccard_ek_chain_merge_250"),
                  ("DSS", "jaccard_DSS_from_scratch")]
    n_anybp = [(f4[c[1]] > 0).sum() for c in callers_d]
    n_strong = [(f4[c[1]] >= 0.5).sum() for c in callers_d]
    width = 0.35; xs = np.arange(len(callers_d))
    ax.bar(xs - width / 2, n_anybp, width=width, label="any-bp",
            color="#3498db", edgecolor="black", linewidth=0.5)
    ax.bar(xs + width / 2, n_strong, width=width, label="J ≥ 0.5",
            color="#1f618d", edgecolor="black", linewidth=0.5)
    ax.set_xticks(xs); ax.set_xticklabels([c[0] for c in callers_d],
                                              fontsize=9)
    ax.set_ylabel("# of 20 paper-named genes")
    ax.set_ylim(0, 21)
    ax.set_title("D · Paper Fig 3B gene hits (top 10 hyper + 10 hypo)")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(axis="y", alpha=0.2)

    # ---- E: feature distribution stacked-bar ----------------------------
    ax = axes[1, 1]
    af = pd.read_csv(DATA_DIR / "comparisons" / "annotation_distribution.csv")
    callers_e = ["paper-DSS (Supp Table 5)", "methylKit-tile",
                  "epykit-chain_merge-100", "epykit-chain_merge-250",
                  "DSS-from-scratch"]
    short_e = ["paper", "mk-tile", "ek-cm-100", "ek-cm-250", "DSS"]
    pivot = af.pivot(index="caller", columns="feature_type",
                       values="fraction").reindex(callers_e) \
                .reindex(columns=FEAT_ORDER).fillna(0)
    bottoms = np.zeros(len(callers_e))
    for f in FEAT_ORDER:
        ax.bar(range(len(callers_e)), pivot[f] * 100, bottom=bottoms * 100,
                color=FEAT_COLOR[f], edgecolor="white", linewidth=0.4,
                width=0.7, label=f)
        bottoms += pivot[f].values
    ax.set_xticks(range(len(callers_e)))
    ax.set_xticklabels(short_e, fontsize=9)
    ax.set_ylabel("% of DMRs"); ax.set_ylim(0, 100)
    ax.set_title("E · HOMER feature distribution")
    ax.legend(fontsize=7, frameon=False, ncol=2, loc="lower left",
                bbox_to_anchor=(1.02, 0.0))

    # ---- F: resource bars -----------------------------------------------
    ax = axes[1, 2]
    rdata = pd.read_csv(THREE_WAY / "F7_resources_data.csv")
    label_map = {"methylKit-tile": "mk-tile",
                  "epykit-tile": "ek-tile",
                  "epykit-chain_merge-100": "ek-cm-100",
                  "DSS-from-scratch": "DSS"}
    rdata["short"] = rdata["caller"].map(label_map)
    rdata = rdata[rdata["short"].notna()]
    colors_f = [CALLER_COLOR[c] for c in rdata["caller"]]
    bars = ax.bar(rdata["short"], rdata["wall_s"],
                    color=colors_f, edgecolor="black", linewidth=0.6)
    ax.set_yscale("log"); ax.set_ylabel("Wall time (s, log)")
    ax.set_title("F · Pipeline wall time")
    for b, v in zip(bars, rdata["wall_s"]):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.06,
                  f"{v:,.0f}s", ha="center", va="bottom", fontsize=9)
    ax.grid(axis="y", alpha=0.2, which="both")

    fig.suptitle("F10 · Study 3 — three-way DMR-caller comparison "
                  "(GSE263850, AKAP11 KO vs WT)",
                  fontsize=14, fontweight="bold", y=1.02)
    save_dual(fig, THREE_WAY / "F10_summary_three_way")
    plt.close(fig)


if __name__ == "__main__":
    main()
