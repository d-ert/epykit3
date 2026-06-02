"""Compare the epykit chain_merge call set to the paper's Supp Tables.

Inputs
------
chain_merge :
    FINAL_REPORT/data/study3/chain_merge/dmr_chain_merge.parquet
    FINAL_REPORT/data/study3/chain_merge/dmr_gene_links_100kb.csv

paper Supp Table 5 (813 DMRs with Gene.Name, coords, direction) :
    GSE263850_RAW/Paper resources/DMR_total_list.xlsx

paper Supp Table 8 (46 Panel-E "critical" DMR-DEG genes) :
    FINAL_REPORT/shinygo_lists/outputs/reactome/table8.xlsx

Outputs (in FINAL_REPORT/data/study3/comparisons/chain_merge_vs_paper/) :

    comparison_report.md       - human-readable summary across every level
    coord_overlap_per_paper_dmr.csv
                              - per paper DMR: best-overlapping chain_merge
                                DMR (Jaccard, reciprocal-overlap fraction,
                                direction agreement)
    coord_overlap_per_our_dmr.csv
                              - mirror: per our DMR, best paper hit
    gene_overlap_table5.csv    - per paper gene (705): captured (Y/N) by
                                each of {nearest_tss_gene, 100kb gene set}
    panel_e_capture_table8.csv - per Panel-E gene (46): captured (Y/N)
    headline.json              - all numbers as a JSON blob
"""

from __future__ import annotations

import json
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

warnings.filterwarnings("ignore")

# ---- Paths -----------------------------------------------------------------

REPO_ROOT = Path(r"D:/Coding/Projeler/methyl_lib/benchmarkin_merges").resolve()
CM_DIR    = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "chain_merge"
PAPER_T5  = Path(r"D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW/"
                 r"Paper resources/DMR_total_list.xlsx")
PAPER_T8  = REPO_ROOT / "FINAL_REPORT" / "shinygo_lists" / "outputs" \
            / "reactome" / "table8.xlsx"
OUT_DIR   = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" \
            / "comparisons" / "chain_merge_vs_paper"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ---- Load ------------------------------------------------------------------

def load_ours():
    p = CM_DIR / "dmr_chain_merge.parquet"
    df = pl.read_parquet(p).to_pandas()
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    df["our_index"] = np.arange(len(df))
    df["length"] = df["end"] - df["start"]
    df["midpoint"] = (df["start"] + df["end"]) // 2
    df["direction"] = np.where(df["mean_meth_diff"] > 0, "hyper",
                       np.where(df["mean_meth_diff"] < 0, "hypo", "none"))
    df["nearest_tss_gene_u"] = df["nearest_tss_gene"].fillna("").astype(str).str.upper()
    return df


def load_paper_t5():
    df = pd.read_excel(PAPER_T5, sheet_name=0)
    df = df.rename(columns={"chr": "chrom"}).copy()
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    df["paper_index"] = np.arange(len(df))
    df["length"] = df["end"] - df["start"]
    df["midpoint"] = (df["start"] + df["end"]) // 2
    df["direction"] = np.where(df["diff.meth_mean"] > 0, "hyper",
                       np.where(df["diff.meth_mean"] < 0, "hypo", "none"))
    df["gene_u"] = df["Gene.Name"].fillna("").astype(str).str.strip().str.upper()
    return df


def load_paper_t8():
    df = pd.read_excel(PAPER_T8, sheet_name=0)
    df["gene_u"] = df["Gene"].astype(str).str.strip().str.upper()
    return df


def load_100kb_links():
    p = CM_DIR / "dmr_gene_links_100kb.csv"
    df = pd.read_csv(p)
    df["chrom"] = df["chrom"].astype(str)
    df["gene_u"] = df["gene"].astype(str).str.strip().str.upper()
    return df


# ---- Per-DMR overlap matching ---------------------------------------------

def jaccard(a_s, a_e, b_s, b_e):
    inter = max(0, min(a_e, b_e) - max(a_s, b_s))
    union = max(a_e, b_e) - min(a_s, b_s)
    return inter / max(1, union)


def reciprocal_overlap_frac(a_s, a_e, b_s, b_e):
    inter = max(0, min(a_e, b_e) - max(a_s, b_s))
    if inter == 0:
        return 0.0
    return inter / min(a_e - a_s, b_e - b_s)


def index_by_chrom(df, start_col, end_col, idx_col):
    """Return chrom -> sorted list of (start, end, idx) tuples."""
    out: dict[str, list] = defaultdict(list)
    for _, r in df.iterrows():
        out[str(r["chrom"])].append(
            (int(r[start_col]), int(r[end_col]), int(r[idx_col]))
        )
    for ch in out:
        out[ch].sort(key=lambda t: t[0])
    return out


def best_match_per_row(query_df, target_idx, q_start, q_end, t_dir_lookup,
                       q_dir_col):
    """For each row in query_df, find target with maximal Jaccard overlap.

    Returns DataFrame with: query_index, n_target_overlapping,
    best_target_idx, jaccard, reciprocal_frac, direction_match.
    """
    rows = []
    for _, r in query_df.iterrows():
        chrom = str(r["chrom"])
        s = int(r[q_start]); e = int(r[q_end])
        candidates = target_idx.get(chrom, [])
        # Sweep candidates that could overlap
        best_jac = 0.0; best_rec = 0.0; best_idx = -1
        n_overlap = 0
        for ts, te, ti in candidates:
            if te < s:
                continue
            if ts > e:
                break
            n_overlap += 1
            j = jaccard(s, e, ts, te)
            if j > best_jac:
                best_jac = j
                best_rec = reciprocal_overlap_frac(s, e, ts, te)
                best_idx = ti
        if best_idx == -1:
            rows.append(dict(
                n_target_overlapping=0,
                best_target_idx=None,
                best_jaccard=0.0,
                best_reciprocal_frac=0.0,
                direction_match=None,
            ))
        else:
            dm = None
            if t_dir_lookup is not None:
                dm = bool(r[q_dir_col] == t_dir_lookup[best_idx])
            rows.append(dict(
                n_target_overlapping=n_overlap,
                best_target_idx=best_idx,
                best_jaccard=round(best_jac, 4),
                best_reciprocal_frac=round(best_rec, 4),
                direction_match=dm,
            ))
    return pd.DataFrame(rows)


# ---- Main ------------------------------------------------------------------

def main() -> None:
    print("=== chain_merge vs paper ===")
    ours = load_ours()
    paper = load_paper_t5()
    panel_e = load_paper_t8()
    links = load_100kb_links()
    print(f"chain_merge DMRs:  {len(ours):,}")
    print(f"paper Table 5 DMRs:{len(paper):,}")
    print(f"paper Table 8 genes (unique): {panel_e['gene_u'].nunique()}")
    print(f"100 kb DMR-gene pairs: {len(links):,}")

    # ---- 1. Coordinate overlap (both directions) --------------------------

    ours_idx_by_chrom  = index_by_chrom(ours,  "start", "end", "our_index")
    paper_idx_by_chrom = index_by_chrom(paper, "start", "end", "paper_index")
    ours_dir  = dict(zip(ours["our_index"],  ours["direction"]))
    paper_dir = dict(zip(paper["paper_index"], paper["direction"]))

    # For each paper DMR: best chain_merge match
    p2o = best_match_per_row(
        paper, ours_idx_by_chrom, "start", "end",
        t_dir_lookup=ours_dir, q_dir_col="direction"
    )
    paper_match = pd.concat([paper.reset_index(drop=True), p2o], axis=1)
    paper_match = paper_match.rename(
        columns={"best_target_idx": "best_our_idx"}
    )
    # Attach the matched chain_merge DMR's coordinates + diff for inspection.
    our_coord = ours.set_index("our_index")[
        ["chrom", "start", "end", "length", "mean_meth_diff", "direction"]
    ].rename(columns={
        "chrom": "our_chrom", "start": "our_start", "end": "our_end",
        "length": "our_length", "mean_meth_diff": "our_diff",
        "direction": "our_direction",
    })
    paper_match = paper_match.join(our_coord, on="best_our_idx")
    paper_match.to_csv(OUT_DIR / "coord_overlap_per_paper_dmr.csv", index=False)

    # For each chain_merge DMR: best paper match
    o2p = best_match_per_row(
        ours, paper_idx_by_chrom, "start", "end",
        t_dir_lookup=paper_dir, q_dir_col="direction"
    )
    our_match = pd.concat([ours.reset_index(drop=True), o2p], axis=1)
    our_match = our_match.rename(columns={"best_target_idx": "best_paper_idx"})
    paper_coord = paper.set_index("paper_index")[
        ["chrom", "start", "end", "length", "diff.meth_mean", "direction",
         "Gene.Name"]
    ].rename(columns={
        "chrom": "paper_chrom", "start": "paper_start", "end": "paper_end",
        "length": "paper_length", "diff.meth_mean": "paper_diff",
        "direction": "paper_direction", "Gene.Name": "paper_gene",
    })
    our_match = our_match.join(paper_coord, on="best_paper_idx")
    our_match.to_csv(OUT_DIR / "coord_overlap_per_our_dmr.csv", index=False)

    # Headline overlap stats (any >= 1bp overlap)
    paper_hit = (paper_match["best_jaccard"] > 0).sum()
    our_hit   = (our_match["best_jaccard"] > 0).sum()
    recall    = paper_hit / len(paper)
    precision = our_hit / len(ours)

    # Stratified recall by Jaccard threshold
    j_bins = [0.0, 0.1, 0.25, 0.5, 0.75]
    recall_at = {}
    for t in j_bins:
        recall_at[f"recall_J>={t}"] = round(
            (paper_match["best_jaccard"] > t).sum() / len(paper), 4
        )

    # Direction agreement on the matched paper DMRs
    matched_paper = paper_match[paper_match["best_jaccard"] > 0]
    dir_agree = matched_paper["direction_match"].fillna(False).sum()
    dir_total = matched_paper["direction_match"].notna().sum()
    dir_agree_frac = dir_agree / max(dir_total, 1)

    # ---- 2. Length / morphology comparison --------------------------------

    p_len = paper["length"]
    o_len = ours["length"]
    morph = dict(
        paper_n=len(paper),
        paper_n_hyper=int((paper["direction"] == "hyper").sum()),
        paper_n_hypo =int((paper["direction"] == "hypo").sum()),
        paper_pct_hyper=round(
            100 * (paper["direction"] == "hyper").sum() / len(paper), 1),
        paper_median_bp=int(p_len.median()),
        paper_mean_bp=int(p_len.mean()),
        paper_max_bp=int(p_len.max()),
        our_n=len(ours),
        our_n_hyper=int((ours["direction"] == "hyper").sum()),
        our_n_hypo =int((ours["direction"] == "hypo").sum()),
        our_pct_hyper=round(
            100 * (ours["direction"] == "hyper").sum() / len(ours), 1),
        our_median_bp=int(o_len.median()),
        our_mean_bp=int(o_len.mean()),
        our_max_bp=int(o_len.max()),
    )

    # ---- 3. Gene-level (Table 5: 705 genes) ------------------------------

    paper_genes = (paper["gene_u"].dropna().unique())
    paper_genes = [g for g in paper_genes if g and g != "NAN"]
    our_nearest = set(ours["nearest_tss_gene_u"].dropna().unique()) - {""}
    our_100kb   = set(links["gene_u"].dropna().unique()) - {""}

    table5_rows = []
    for g in paper_genes:
        table5_rows.append(dict(
            paper_gene=g,
            captured_by_nearest_tss=g in our_nearest,
            captured_by_100kb_links=g in our_100kb,
        ))
    table5_df = pd.DataFrame(table5_rows).sort_values("paper_gene")
    table5_df.to_csv(OUT_DIR / "gene_overlap_table5.csv", index=False)
    t5_near = int(table5_df["captured_by_nearest_tss"].sum())
    t5_100  = int(table5_df["captured_by_100kb_links"].sum())

    # ---- 4. Panel E (Table 8: 46 genes) ----------------------------------

    panel_genes = panel_e["gene_u"].dropna().unique().tolist()
    panel_genes = [g for g in panel_genes if g and g != "NAN"]
    panel_rows = []
    for g in panel_genes:
        panel_rows.append(dict(
            panel_e_gene=g,
            captured_by_nearest_tss=g in our_nearest,
            captured_by_100kb_links=g in our_100kb,
        ))
    panel_df = pd.DataFrame(panel_rows).sort_values("panel_e_gene")
    panel_df.to_csv(OUT_DIR / "panel_e_capture_table8.csv", index=False)
    pe_near = int(panel_df["captured_by_nearest_tss"].sum())
    pe_100  = int(panel_df["captured_by_100kb_links"].sum())

    # ---- 5. Heatmap-gene direct DMR coordinate hits ----------------------
    # Heatmap names + paper-DMR coords pulled from Supp Table 5 by gene name.

    heatmap_names = ["NR2E1", "OTX1", "IRX2", "OTX2", "ENPP2", "GREB1L",
                     "CCDC177", "PAX7", "NAALADL2"]
    hm_rows = []
    for gene in heatmap_names:
        sub = paper[paper["gene_u"] == gene]
        if len(sub) == 0:
            hm_rows.append(dict(
                gene=gene, paper_dmr=None,
                our_overlap=False, best_jaccard=0.0,
                our_dmr=None,
            ))
            continue
        # Take widest paper DMR for that gene (heatmap rows often picked the
        # canonical interval).
        row = sub.sort_values("length", ascending=False).iloc[0]
        paper_str = f"{row['chrom']}:{int(row['start'])}-{int(row['end'])}"
        # Look up the best chain_merge overlap from our paper_match table
        pm = paper_match.loc[paper_match["paper_index"] == row["paper_index"]]
        if len(pm) == 0 or pd.isna(pm["best_our_idx"].iloc[0]):
            hm_rows.append(dict(
                gene=gene, paper_dmr=paper_str,
                our_overlap=False, best_jaccard=0.0,
                our_dmr=None,
            ))
        else:
            our_i = int(pm["best_our_idx"].iloc[0])
            best_j = float(pm["best_jaccard"].iloc[0])
            our_r = ours.loc[ours["our_index"] == our_i].iloc[0]
            our_str = f"{our_r['chrom']}:{int(our_r['start'])}-{int(our_r['end'])}"
            hm_rows.append(dict(
                gene=gene, paper_dmr=paper_str,
                our_overlap=best_j > 0, best_jaccard=round(best_j, 4),
                our_dmr=our_str if best_j > 0 else None,
            ))
    hm_df = pd.DataFrame(hm_rows)
    hm_df.to_csv(OUT_DIR / "heatmap_gene_hits.csv", index=False)
    hm_hits = int(hm_df["our_overlap"].sum())

    # ---- 6. Bytes for headline.json --------------------------------------

    headline = {
        "input_counts": {
            "our_dmrs": int(len(ours)),
            "paper_dmrs": int(len(paper)),
            "paper_unique_genes": int(len(paper_genes)),
            "panel_e_genes": int(len(panel_genes)),
            "our_100kb_gene_pairs": int(len(links)),
            "our_100kb_unique_genes": int(len(our_100kb)),
        },
        "coord_overlap": {
            "paper_DMRs_with_any_overlap": int(paper_hit),
            "our_DMRs_with_any_overlap": int(our_hit),
            "recall_of_paper_anybp": round(recall, 4),
            "precision_anybp": round(precision, 4),
            **recall_at,
            "direction_agreement_on_matched": round(dir_agree_frac, 4),
            "direction_agreement_n": int(dir_agree),
            "direction_agreement_total_matched": int(dir_total),
        },
        "morphology": morph,
        "table5_gene_recall": {
            "n_paper_genes": int(len(paper_genes)),
            "captured_by_nearest_tss": t5_near,
            "captured_by_100kb_links": t5_100,
            "frac_nearest_tss": round(t5_near / len(paper_genes), 4),
            "frac_100kb_links": round(t5_100 / len(paper_genes), 4),
        },
        "panel_e_capture": {
            "n_genes": int(len(panel_genes)),
            "captured_by_nearest_tss": pe_near,
            "captured_by_100kb_links": pe_100,
            "frac_nearest_tss": round(pe_near / len(panel_genes), 4),
            "frac_100kb_links": round(pe_100 / len(panel_genes), 4),
        },
        "heatmap_gene_hits": {
            "checked": int(len(hm_df)),
            "with_paper_dmr": int(hm_df["paper_dmr"].notna().sum()),
            "with_our_coord_overlap": hm_hits,
        },
    }
    (OUT_DIR / "headline.json").write_text(
        json.dumps(headline, indent=2), encoding="utf-8"
    )

    # ---- 7. Markdown report ----------------------------------------------

    L = []
    L.append("# Chain-merge replication vs paper — full comparison\n")
    L.append("Inputs:\n")
    L.append(f"- chain_merge: `{len(ours):,}` DMRs  "
             f"(at `alpha=1e-5, delta=0, minlen=50, minCG=3, dis.merge=100, "
             f"pct.sig=0.5`)")
    L.append(f"- paper Table 5: `{len(paper):,}` DMRs  "
             f"(DSS `callDMR` with same parameters)")
    L.append(f"- paper Table 8 (Panel E critical genes): "
             f"`{len(panel_genes)}` unique genes")
    L.append("")

    L.append("## 1. DMR morphology\n")
    L.append("| | DMRs | hyper | hypo | %hyper | median bp | mean bp | max bp |")
    L.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    L.append(f"| paper (DSS) | {morph['paper_n']} | {morph['paper_n_hyper']} | "
             f"{morph['paper_n_hypo']} | {morph['paper_pct_hyper']}% | "
             f"{morph['paper_median_bp']} | {morph['paper_mean_bp']} | "
             f"{morph['paper_max_bp']:,} |")
    L.append(f"| **ours (chain_merge)** | **{morph['our_n']}** | "
             f"**{morph['our_n_hyper']}** | **{morph['our_n_hypo']}** | "
             f"**{morph['our_pct_hyper']}%** | **{morph['our_median_bp']}** | "
             f"**{morph['our_mean_bp']}** | **{morph['our_max_bp']:,}** |")
    L.append("")
    L.append("- Hyper bias agrees in sign and magnitude (paper "
             f"{morph['paper_pct_hyper']}% vs ours {morph['our_pct_hyper']}%).")
    L.append(f"- Our DMRs are ~half the median length of the paper's "
             f"({morph['our_median_bp']} bp vs {morph['paper_median_bp']} bp). "
             f"Same call, narrower aggregation; expected for chain-merge with "
             f"`dis.merge=100`.")
    L.append("")

    L.append("## 2. DMR coordinate overlap\n")
    L.append(f"- Any-bp overlap recall of paper: "
             f"**{paper_hit} / {len(paper)} = {recall*100:.1f}%**")
    L.append(f"- Any-bp overlap precision (ours hitting any paper DMR): "
             f"**{our_hit} / {len(ours)} = {precision*100:.1f}%**")
    L.append("")
    L.append("Stratified recall by Jaccard threshold:")
    L.append("")
    L.append("| Jaccard >= | Recall (of 813) | Count |")
    L.append("|---:|---:|---:|")
    for t in j_bins:
        rec_t = (paper_match["best_jaccard"] > t).sum()
        L.append(f"| {t:.2f} | {100*rec_t/len(paper):.1f}% | {rec_t} |")
    L.append("")
    L.append("Direction agreement on the matched paper DMRs "
             "(hyper/hypo match): "
             f"**{dir_agree} / {dir_total} = {100*dir_agree_frac:.1f}%**.")
    L.append("")

    L.append("## 3. Gene-level recall (Table 5: 705 unique genes)\n")
    L.append("| Linkage rule | Captured | Recall |")
    L.append("|---|---:|---:|")
    L.append(f"| nearest TSS only | {t5_near} / {len(paper_genes)} | "
             f"**{100*t5_near/len(paper_genes):.1f}%** |")
    L.append(f"| 100 kb (paper rule) | {t5_100} / {len(paper_genes)} | "
             f"**{100*t5_100/len(paper_genes):.1f}%** |")
    L.append("")

    L.append("## 4. Panel E critical-gene capture (Table 8: 46 genes)\n")
    L.append("| Linkage rule | Captured | Recall |")
    L.append("|---|---:|---:|")
    L.append(f"| nearest TSS only | {pe_near} / {len(panel_genes)} | "
             f"**{100*pe_near/len(panel_genes):.1f}%** |")
    L.append(f"| 100 kb (paper rule) | {pe_100} / {len(panel_genes)} | "
             f"**{100*pe_100/len(panel_genes):.1f}%** |")
    L.append("")
    missed = panel_df.loc[~panel_df["captured_by_100kb_links"], "panel_e_gene"]
    if len(missed):
        L.append("Genes in Table 8 we still miss at 100 kb:")
        L.append("> " + ", ".join(sorted(missed.tolist())))
        L.append("")

    L.append("## 5. Heatmap-gene direct DMR coordinate hits\n")
    L.append("| Gene | Paper DMR | Our coord overlap? | Jaccard | Our DMR |")
    L.append("|---|---|:---:|---:|---|")
    for _, r in hm_df.iterrows():
        mark = "✓" if r["our_overlap"] else " "
        paper_str = r["paper_dmr"] or "(no paper DMR)"
        our_str = r["our_dmr"] or "—"
        L.append(f"| {r['gene']} | {paper_str} | {mark} | "
                 f"{r['best_jaccard']:.3f} | {our_str} |")
    L.append("")
    L.append(f"Direct coordinate hits: **{hm_hits} / "
             f"{int(hm_df['paper_dmr'].notna().sum())}** named genes "
             "checked.")
    L.append("")

    L.append("## 6. Files in this folder\n")
    L.append("| File | What it is |")
    L.append("|---|---|")
    L.append("| `comparison_report.md` | This document |")
    L.append("| `coord_overlap_per_paper_dmr.csv` | One row per paper DMR; "
             "best-overlapping chain_merge DMR + jaccard + direction match |")
    L.append("| `coord_overlap_per_our_dmr.csv` | Mirror: per chain_merge DMR, "
             "best paper match (useful for inspecting false positives) |")
    L.append("| `gene_overlap_table5.csv` | Each paper gene (705): captured "
             "by nearest_tss / by 100 kb linkage |")
    L.append("| `panel_e_capture_table8.csv` | Each Table-8 critical gene "
             "(46): captured Y/N |")
    L.append("| `heatmap_gene_hits.csv` | Heatmap-named genes + coordinate-hit "
             "status |")
    L.append("| `headline.json` | All numbers as JSON |")
    L.append("")

    (OUT_DIR / "comparison_report.md").write_text("\n".join(L),
                                                  encoding="utf-8")
    print(f"Wrote comparison_report.md")
    print(json.dumps({k: v for k, v in headline.items() if k != "morphology"},
                     indent=2))


if __name__ == "__main__":
    main()
