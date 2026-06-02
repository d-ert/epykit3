"""F4 — Top-10 hyper + top-10 hypo paper-named gene hits across callers.

The paper's Fig 3B labels the top 10 hyper- and top 10 hypo-methylated
DMR-associated genes. We:
  1. Pull those names from Supp Table 5 (sort by |diff.meth_mean|).
  2. For each named gene, locate the widest paper-DMR row for that gene.
  3. Check overlap (any bp) against each caller's DMR set.

Render: heatmap. Rows = genes; columns = {paper (always ✓), methylKit,
ek-100, ek-250, DSS}. Cells: dot if overlap, blank if not. Annotated
with Jaccard.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _style import (DATA_DIR, THREE_WAY, CALLER_COLOR, setup, save_dual)

PAPER_T5 = Path(r"D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW/"
                r"Paper resources/DMR_total_list.xlsx")
MK_TILE  = Path(r"D:/Coding/Projeler/methyl_lib/methylkıt_realResults/"
                r"scripts_and_results/methylkit_results/dmr_significant_lenient.csv")


def jaccard(a_s, a_e, b_s, b_e) -> float:
    inter = max(0, min(a_e, b_e) - max(a_s, b_s))
    union = max(a_e, b_e) - min(a_s, b_s)
    return inter / max(1, union)


def index_dmrs(df: pd.DataFrame) -> dict[str, list]:
    by = defaultdict(list)
    for _, r in df.iterrows():
        by[str(r["chrom"])].append((int(r["start"]), int(r["end"])))
    for ch in by:
        by[ch].sort(key=lambda t: t[0])
    return by


def best_overlap(chrom: str, s: int, e: int, idx: dict) -> float:
    cands = idx.get(chrom, [])
    best = 0.0
    for ts, te in cands:
        if te < s: continue
        if ts > e: break
        j = jaccard(s, e, ts, te)
        if j > best: best = j
    return best


def main() -> None:
    setup()
    paper = pd.read_excel(PAPER_T5).rename(columns={"chr": "chrom"})
    paper["chrom"] = paper["chrom"].astype(str)
    paper["abs_dB"] = paper["diff.meth_mean"].abs()
    paper["length"] = paper["end"].astype(int) - paper["start"].astype(int)

    # Per Fig 3B caption: "top 10 hyper- and 10 hypo-methylated
    # DMR-associated genes". Sort by signed diff.meth_mean and take
    # top 10 in each direction by |diff.meth_mean|.
    hyper = paper[paper["diff.meth_mean"] > 0]
    hypo  = paper[paper["diff.meth_mean"] < 0]
    top_hyper = (hyper.sort_values("abs_dB", ascending=False)
                       .drop_duplicates("Gene.Name").head(10))
    top_hypo  = (hypo.sort_values("abs_dB", ascending=False)
                       .drop_duplicates("Gene.Name").head(10))
    top = pd.concat([top_hyper.assign(group="hyper"),
                      top_hypo.assign(group="hypo")], ignore_index=True)
    top["gene"] = top["Gene.Name"].astype(str)
    top["paper_dmr_str"] = (
        top["chrom"] + ":" + top["start"].astype(str) + "-"
        + top["end"].astype(str))

    # Caller DMR indices
    mk_idx = index_dmrs(pd.read_csv(MK_TILE))
    ek100_idx = index_dmrs(
        pl.read_parquet(DATA_DIR / "chain_merge" /
                         "dmr_chain_merge.parquet").to_pandas())
    ek250_idx = index_dmrs(
        pl.read_parquet(DATA_DIR / "chain_merge_dis_merge_sweep" /
                         "dis_merge_250" / "dmr.parquet").to_pandas())
    dss_idx = index_dmrs(pd.read_csv(DATA_DIR / "dss" / "dmr_dss.csv"))

    results = []
    for _, r in top.iterrows():
        results.append({
            "gene": r["gene"],
            "group": r["group"],
            "diff_mean": float(r["diff.meth_mean"]),
            "paper_dmr": r["paper_dmr_str"],
            "paper_length": int(r["length"]),
            "jaccard_methylKit_tile":
                best_overlap(r["chrom"], int(r["start"]),
                             int(r["end"]), mk_idx),
            "jaccard_ek_chain_merge_100":
                best_overlap(r["chrom"], int(r["start"]),
                             int(r["end"]), ek100_idx),
            "jaccard_ek_chain_merge_250":
                best_overlap(r["chrom"], int(r["start"]),
                             int(r["end"]), ek250_idx),
            "jaccard_DSS_from_scratch":
                best_overlap(r["chrom"], int(r["start"]),
                             int(r["end"]), dss_idx),
        })
    res = pd.DataFrame(results)
    res.to_csv(THREE_WAY / "F4_top_named_gene_hits_data.csv", index=False)

    # ---- Heatmap render ------------------------------------------------
    callers = [
        ("methylKit-tile",          "jaccard_methylKit_tile"),
        ("epykit-chain_merge-100",  "jaccard_ek_chain_merge_100"),
        ("epykit-chain_merge-250",  "jaccard_ek_chain_merge_250"),
        ("DSS-from-scratch",        "jaccard_DSS_from_scratch"),
    ]
    res_sorted = pd.concat([
        res[res["group"] == "hyper"].sort_values("diff_mean", ascending=False),
        res[res["group"] == "hypo"].sort_values("diff_mean", ascending=True),
    ], ignore_index=True)

    n = len(res_sorted)
    fig, ax = plt.subplots(figsize=(8.0, max(5.0, n * 0.32)),
                             constrained_layout=True)

    mat = res_sorted[[c[1] for c in callers]].values.astype(float)
    n_caller = len(callers)
    # Background grid
    for i in range(n):
        for j in range(n_caller):
            j_val = mat[i, j]
            facecolor = "#f5f5f5" if j_val == 0 else None
            ax.add_patch(plt.Rectangle((j - 0.45, i - 0.45), 0.9, 0.9,
                                         facecolor=facecolor or "none",
                                         edgecolor="#ecf0f1", linewidth=0.5))
            if j_val > 0:
                # Color by Jaccard (deeper = higher)
                col = plt.cm.viridis(j_val * 0.85 + 0.15)
                ax.add_patch(plt.Rectangle((j - 0.45, i - 0.45), 0.9, 0.9,
                                             facecolor=col, edgecolor="white",
                                             linewidth=0.6))
                # Print Jaccard
                txt_color = "white" if j_val > 0.4 else "#1c1c1c"
                ax.text(j, i, f"{j_val:.2f}", ha="center", va="center",
                          color=txt_color, fontsize=8, fontweight="bold")

    ax.set_xticks(range(n_caller))
    ax.set_xticklabels([c[0] for c in callers], rotation=18, ha="right",
                         fontsize=9)
    # Y labels: gene + direction marker
    yticks = []
    for _, r in res_sorted.iterrows():
        marker = "▲" if r["group"] == "hyper" else "▼"
        ymask = "#c0392b" if r["group"] == "hyper" else "#2980b9"
        yticks.append(f"{marker} {r['gene']}")
    ax.set_yticks(range(n))
    ax.set_yticklabels(yticks, fontsize=9)
    ax.set_xlim(-0.6, n_caller - 0.4)
    ax.set_ylim(-0.6, n - 0.4)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    for spine in ("top", "right", "bottom", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(left=False, bottom=False)
    ax.set_title("F4 · Paper Fig 3B labeled genes — coordinate-overlap "
                  "hits across callers\n(top 10 hyper ▲ + top 10 hypo ▼ "
                  "DMR-associated genes; values = Jaccard of best overlap)",
                  fontsize=11, pad=12)

    # Colorbar (Jaccard)
    cax = fig.add_axes([0.13, 0.04, 0.45, 0.018])
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis,
                                 norm=plt.Normalize(0.15, 1.0))
    cb = plt.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("Jaccard of best overlap (blank = no overlap)",
                  fontsize=8)
    cb.ax.tick_params(labelsize=8)

    save_dual(fig, THREE_WAY / "F4_top_named_gene_hits")
    plt.close(fig)

    # Per-caller hit count summary
    print("\nHit count per caller (out of 20 named genes):")
    for col_label, col in callers:
        n_hit = int((res_sorted[col] > 0).sum())
        n_strong = int((res_sorted[col] >= 0.5).sum())
        print(f"  {col_label}: {n_hit}/20 any-bp; {n_strong}/20 J>=0.5")


if __name__ == "__main__":
    main()
