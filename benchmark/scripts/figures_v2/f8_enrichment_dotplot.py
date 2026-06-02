"""F8 — Reactome / KEGG enrichment dot plot across callers.

Custom comparison since Enrichr's BH correction differs from ShinyGO's
"Curated.Reactome". We render two complementary panels:

  A. Top 10 Reactome_2022 / KEGG_2021_Human terms per caller, side-by-side,
     with paper-Panel-D keyword matches highlighted. Dot size = overlap
     count, dot color = -log10(raw p-value).
  B. Spearman rank correlation matrix between callers' top-50 term ranks
     (term-overlap measure).

Source: comparisons/enrichment_three_way.json from P1.4.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _style import (DATA_DIR, THREE_WAY, setup, save_dual)


def panel_d_keywords(term: str) -> bool:
    tl = term.lower()
    KW = ["rhodopsin", "peptide ligand", "gpcr ligand", "gpcr downstream",
          "signaling by gpcr", "signalling by gpcr", "g alpha i",
          "g alpha (i)", "neuroactive ligand", "camp signal",
          "morphine addiction"]
    return any(k in tl for k in KW)


def short(term: str, n: int = 60) -> str:
    t = term.split(" R-HSA-")[0]
    return t if len(t) <= n else t[:n - 1] + "…"


def main() -> None:
    setup()
    data = json.loads((DATA_DIR / "comparisons" /
                       "enrichment_three_way.json").read_text(encoding="utf-8"))

    callers = ["paper_Table5", "methylKit_tile",
               "ek_chain_merge_100", "ek_chain_merge_250",
               "DSS_from_scratch"]
    caller_labels = {
        "paper_Table5":          "paper (Table 5)",
        "methylKit_tile":        "methylKit-tile",
        "ek_chain_merge_100":    "ek-chain_merge-100",
        "ek_chain_merge_250":    "ek-chain_merge-250",
        "DSS_from_scratch":      "DSS-from-scratch",
    }

    # ---- Panel A: top-12 Reactome terms per caller, vertical stack -----
    fig, axes = plt.subplots(1, 2, figsize=(15.0, 7.0),
                              constrained_layout=True,
                              gridspec_kw={"width_ratios": [2.3, 1]})

    ax = axes[0]
    # Collect union of top-10 terms across callers + paper terms
    LIB = "Reactome_2022"
    rows = []
    for c in callers:
        lib = data.get(c, {}).get("libraries", {}).get(LIB, {})
        top10 = (lib or {}).get("top20", [])[:10]
        for rank, t in enumerate(top10):
            rows.append({
                "caller": c, "rank": rank + 1, "term": t["term"],
                "p": t["p_value"], "p_adj": t["p_adj"],
                "n_overlap": t["n_overlap"],
                "panel_d_match": panel_d_keywords(t["term"]),
            })
    df = pd.DataFrame(rows)

    # Determine global "interesting" terms: union of top-10s, but prefer
    # to highlight panel-D matches and terms present in 2+ callers.
    term_counts = df.groupby("term").size().reset_index(name="n_callers")
    df = df.merge(term_counts, on="term", how="left")
    # Sort terms: panel-D matches first, then by n_callers, then by best rank
    term_rank = (df.groupby("term").agg(
        panel_d_match=("panel_d_match", "max"),
        n_callers=("n_callers", "max"),
        best_rank=("rank", "min"),
    ).reset_index())
    term_rank = term_rank.sort_values(
        by=["panel_d_match", "n_callers", "best_rank"],
        ascending=[False, False, True])
    terms_to_show = term_rank.head(18)["term"].tolist()
    df_show = df[df["term"].isin(terms_to_show)]

    # Reverse y so most interesting is on top
    y_map = {t: i for i, t in enumerate(reversed(terms_to_show))}
    x_map = {c: i for i, c in enumerate(callers)}
    for _, r in df_show.iterrows():
        if r["term"] not in y_map: continue
        x = x_map[r["caller"]]; y = y_map[r["term"]]
        size = 50 + 40 * np.log10(max(r["n_overlap"], 1))
        color_val = -np.log10(max(r["p"], 1e-50))
        ax.scatter(x, y, s=size,
                    c=[color_val], cmap="plasma",
                    vmin=0, vmax=12,
                    edgecolor="black", linewidth=0.5,
                    zorder=3)
        if r["panel_d_match"]:
            ax.scatter(x, y, s=size + 80, facecolors="none",
                        edgecolor="#e74c3c", linewidth=1.5, zorder=2)
    ax.set_xticks(range(len(callers)))
    ax.set_xticklabels([caller_labels[c] for c in callers],
                         rotation=18, ha="right", fontsize=9)
    ax.set_yticks(range(len(terms_to_show)))
    ax.set_yticklabels([short(t) for t in reversed(terms_to_show)],
                         fontsize=8.5)
    ax.set_title("A · Reactome 2022 top-10 terms per caller  ·  "
                  "red ring = paper Panel-D-keyword match",
                  fontsize=10)
    ax.set_xlim(-0.5, len(callers) - 0.5)
    ax.set_ylim(-0.5, len(terms_to_show) - 0.5)
    ax.grid(axis="x", alpha=0.15)

    sm = plt.cm.ScalarMappable(cmap="plasma",
                                 norm=plt.Normalize(0, 12))
    cb = plt.colorbar(sm, ax=ax, location="right", pad=0.02,
                       shrink=0.4)
    cb.set_label("−log10(raw p)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    # ---- Panel B: KEGG paper-term hits (already richer signal) --------
    ax = axes[1]
    LIB = "KEGG_2021_Human"
    rows = []
    for c in callers:
        lib = data.get(c, {}).get("libraries", {}).get(LIB, {})
        top20 = (lib or {}).get("top20", [])
        rows.append({
            "caller": c,
            "n_panel_d_hits": (lib or {}).get("n_paper_term_matches_top20", 0),
            "n_sig_005": sum(1 for t in top20 if t["p_adj"] < 0.05),
            "morphine_p":
                next((t["p_value"] for t in top20
                       if "morphine" in t["term"].lower()), np.nan),
            "neuro_p":
                next((t["p_value"] for t in top20
                       if "neuroactive" in t["term"].lower()), np.nan),
            "camp_p":
                next((t["p_value"] for t in top20
                       if "camp signal" in t["term"].lower()), np.nan),
        })
    dfk = pd.DataFrame(rows)
    width = 0.6
    bars = ax.bar(range(len(callers)), dfk["n_panel_d_hits"],
                    color=["#2c3e50", "#e74c3c", "#3498db", "#1f618d",
                           "#27ae60"], edgecolor="black", linewidth=0.5,
                    width=width)
    for b, v in zip(bars, dfk["n_panel_d_hits"]):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.05,
                  str(int(v)), ha="center", va="bottom", fontsize=9)
    ax.set_xticks(range(len(callers)))
    ax.set_xticklabels([caller_labels[c] for c in callers],
                         rotation=18, ha="right", fontsize=9)
    ax.set_ylabel("# paper-keyword matches in KEGG top-20")
    ax.set_title("B · KEGG 2021 paper-term recovery (top-20)",
                  fontsize=10)
    ax.grid(axis="y", alpha=0.2)

    fig.suptitle("F8 · Pathway enrichment, paper-keyword recovery across callers",
                  fontsize=12, y=1.04)
    save_dual(fig, THREE_WAY / "F8_enrichment_dotplot")
    plt.close(fig)


if __name__ == "__main__":
    main()
