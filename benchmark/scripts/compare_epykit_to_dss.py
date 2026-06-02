"""Compare epykit chain_merge call set against DSS-from-scratch.

Mirrors compare_chain_merge_to_paper.py but uses DSS-922 as the target
instead of the paper's 813. This lets us separate "epykit missed paper"
from "epykit missed DSS" — i.e., distinguish test-statistic divergence
from DSS-vs-paper version drift.

Inputs
------
- FINAL_REPORT/data/study3/chain_merge/dmr_chain_merge.parquet (702 DMRs)
- FINAL_REPORT/data/study3/chain_merge_dis_merge_sweep/dis_merge_250/dmr.parquet
  (940 DMRs at dis.merge=250 — also compared for the "morphology-matched"
  point)
- FINAL_REPORT/data/study3/dss/dmr_dss.csv (922 DMRs)

Outputs (FINAL_REPORT/data/study3/comparisons/epykit_vs_dss/)
- coord_overlap_per_dss_dmr.csv   (per DSS, best epykit match)
- coord_overlap_per_our_dmr.csv   (per epykit, best DSS match)
- coord_overlap_per_dss_dmr_dm250.csv  (same, but using ek-250)
- gene_overlap.csv                (DSS gene set captured by epykit)
- panel_e_capture_dss.csv         (DSS's recall of Table 8 panel-E genes)
- heatmap_gene_hits_dss.csv       (DSS's direct hits on paper Fig 3B labels)
- headline.json                   (all numbers)
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

REPO_ROOT = Path(r"D:/Coding/Projeler/methyl_lib/benchmarkin_merges").resolve()
CM_DIR    = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "chain_merge"
DM250_PQ  = (REPO_ROOT / "FINAL_REPORT" / "data" / "study3"
             / "chain_merge_dis_merge_sweep" / "dis_merge_250" / "dmr.parquet")
DSS_CSV   = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "dss" / "dmr_dss.csv"
PAPER_T5  = Path(r"D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW/"
                 r"Paper resources/DMR_total_list.xlsx")
PAPER_T8  = REPO_ROOT / "FINAL_REPORT" / "shinygo_lists" / "outputs" \
            / "reactome" / "table8.xlsx"
OUT_DIR   = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" \
            / "comparisons" / "epykit_vs_dss"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ---- Helpers (same algorithm as compare_chain_merge_to_paper.py) -----------

def _interval_index(df: pd.DataFrame, start_col="start", end_col="end",
                    idx_col="idx") -> dict[str, list]:
    out: dict[str, list] = defaultdict(list)
    for _, r in df.iterrows():
        out[str(r["chrom"])].append(
            (int(r[start_col]), int(r[end_col]), int(r[idx_col]))
        )
    for ch in out:
        out[ch].sort(key=lambda t: t[0])
    return out


def jaccard(a_s, a_e, b_s, b_e) -> float:
    inter = max(0, min(a_e, b_e) - max(a_s, b_s))
    union = max(a_e, b_e) - min(a_s, b_s)
    return inter / max(1, union)


def reciprocal_overlap(a_s, a_e, b_s, b_e) -> float:
    inter = max(0, min(a_e, b_e) - max(a_s, b_s))
    if inter == 0:
        return 0.0
    return inter / min(a_e - a_s, b_e - b_s)


def best_match_per_row(query_df: pd.DataFrame, target_idx: dict,
                       q_start: str, q_end: str,
                       q_dir_col: str | None,
                       t_dir_lookup: dict | None) -> pd.DataFrame:
    rows = []
    for _, r in query_df.iterrows():
        chrom = str(r["chrom"])
        s = int(r[q_start]); e = int(r[q_end])
        cands = target_idx.get(chrom, [])
        best_j = 0.0; best_rec = 0.0; best_idx = -1
        n_overlap = 0
        for ts, te, ti in cands:
            if te < s:
                continue
            if ts > e:
                break
            n_overlap += 1
            j = jaccard(s, e, ts, te)
            if j > best_j:
                best_j = j
                best_rec = reciprocal_overlap(s, e, ts, te)
                best_idx = ti
        dir_m = None
        if best_idx != -1 and q_dir_col is not None and t_dir_lookup is not None:
            dir_m = bool(r[q_dir_col] == t_dir_lookup[best_idx])
        rows.append(dict(
            n_target_overlapping=n_overlap,
            best_target_idx=best_idx if best_idx != -1 else None,
            best_jaccard=round(best_j, 4),
            best_reciprocal_frac=round(best_rec, 4),
            direction_match=dir_m,
        ))
    return pd.DataFrame(rows)


# ---- Loaders ---------------------------------------------------------------

def load_chain_merge(parquet_path: Path) -> pd.DataFrame:
    df = pl.read_parquet(parquet_path).to_pandas()
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    df["idx"]   = np.arange(len(df))
    df["length"]    = df["end"] - df["start"]
    df["direction"] = np.where(df["mean_meth_diff"] > 0, "hyper",
                       np.where(df["mean_meth_diff"] < 0, "hypo", "none"))
    df["nearest_tss_gene_u"] = (
        df["nearest_tss_gene"].fillna("").astype(str).str.upper()
    )
    return df


def load_dss(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    df["idx"]   = np.arange(len(df))
    df["length"] = df["end"] - df["start"]
    # DSS uses 'dmr_type' column (hyper/hypo) and 'diff_Methy_fromCounts'
    df["direction"] = df["dmr_type"].astype(str)
    df["nearest_tss_gene_u"] = (
        df["nearest_tss_gene"].fillna("").astype(str).str.upper()
    )
    return df


# ---- Two-way comparison ----------------------------------------------------

def two_way_compare(query_df: pd.DataFrame, target_df: pd.DataFrame,
                    query_name: str, target_name: str) -> dict:
    """Compute recall (of target by query) and precision (of query against
    target). 'query' is the one we walk; 'target' is the index we search."""
    q_idx = _interval_index(query_df, "start", "end", "idx")
    t_idx = _interval_index(target_df, "start", "end", "idx")
    q_dir = dict(zip(query_df["idx"], query_df["direction"]))
    t_dir = dict(zip(target_df["idx"], target_df["direction"]))

    # For each target row: best query match (recall, direction match)
    t2q = best_match_per_row(target_df, q_idx, "start", "end",
                              "direction", q_dir)
    target_match = pd.concat([target_df.reset_index(drop=True), t2q], axis=1)
    target_match = target_match.rename(columns={"best_target_idx":
                                                 f"best_{query_name}_idx"})

    # For each query row: best target match (precision)
    q2t = best_match_per_row(query_df, t_idx, "start", "end",
                              "direction", t_dir)
    query_match = pd.concat([query_df.reset_index(drop=True), q2t], axis=1)
    query_match = query_match.rename(columns={"best_target_idx":
                                               f"best_{target_name}_idx"})

    return dict(target_match=target_match, query_match=query_match)


# ---- Headline metrics ------------------------------------------------------

def headline_block(query_df, target_df, target_match: pd.DataFrame,
                    query_match: pd.DataFrame, label: str) -> dict:
    n_q = len(query_df); n_t = len(target_df)
    target_jacs = target_match["best_jaccard"].astype(float).values
    target_hit = (target_jacs > 0).sum()
    query_hit = (query_match["best_jaccard"].astype(float).values > 0).sum()
    matched = target_match[target_match["best_jaccard"] > 0]
    dir_total = matched["direction_match"].notna().sum()
    dir_agree = matched["direction_match"].fillna(False).sum()
    return {
        f"{label}_n_query": int(n_q),
        f"{label}_n_target": int(n_t),
        f"{label}_target_hit_anybp": int(target_hit),
        f"{label}_query_hit_anybp": int(query_hit),
        f"{label}_recall_anybp": round(float(target_hit) / max(n_t, 1), 4),
        f"{label}_precision_anybp": round(float(query_hit) / max(n_q, 1), 4),
        f"{label}_recall_J_0_25": round(float((target_jacs >= 0.25).mean()), 4),
        f"{label}_recall_J_0_5": round(float((target_jacs >= 0.5).mean()), 4),
        f"{label}_recall_J_0_75": round(float((target_jacs >= 0.75).mean()), 4),
        f"{label}_direction_agree_n": int(dir_agree),
        f"{label}_direction_agree_total": int(dir_total),
        f"{label}_direction_agree_frac":
            round(float(dir_agree) / max(dir_total, 1), 4),
    }


# ---- Gene-level + panel-E + heatmap helpers --------------------------------

def gene_overlap_table(query_df, target_df, target_name: str) -> pd.DataFrame:
    """How well does query's nearest-TSS gene set recover target's?"""
    q_genes = set(query_df["nearest_tss_gene_u"].dropna().unique()) - {""}
    t_genes = set(target_df["nearest_tss_gene_u"].dropna().unique()) - {""}
    rows = []
    for g in sorted(t_genes):
        rows.append(dict(
            target_gene=g, captured_by_query=g in q_genes,
        ))
    df = pd.DataFrame(rows)
    return df


# ---- Main ------------------------------------------------------------------

def main() -> None:
    print("Loading chain_merge dis.merge=100 …")
    ek100 = load_chain_merge(CM_DIR / "dmr_chain_merge.parquet")
    print(f"  {len(ek100)} DMRs")

    print("Loading chain_merge dis.merge=250 …")
    ek250 = load_chain_merge(DM250_PQ)
    print(f"  {len(ek250)} DMRs")

    print("Loading DSS …")
    dss = load_dss(DSS_CSV)
    print(f"  {len(dss)} DMRs")

    # ek-100 vs DSS
    print("\n=== ek-chain_merge-100 vs DSS-922 ===")
    cmp100 = two_way_compare(ek100, dss, "ek", "dss")
    h100 = headline_block(ek100, dss,
                          cmp100["target_match"], cmp100["query_match"],
                          "ek100_vs_dss")
    print(json.dumps(h100, indent=2))
    cmp100["target_match"].to_csv(OUT_DIR / "coord_overlap_per_dss_dmr.csv",
                                   index=False)
    cmp100["query_match"].to_csv(OUT_DIR / "coord_overlap_per_our_dmr.csv",
                                  index=False)

    # ek-250 vs DSS (the morphology-matched operating point)
    print("\n=== ek-chain_merge-250 vs DSS-922 ===")
    cmp250 = two_way_compare(ek250, dss, "ek250", "dss")
    h250 = headline_block(ek250, dss,
                          cmp250["target_match"], cmp250["query_match"],
                          "ek250_vs_dss")
    print(json.dumps(h250, indent=2))
    cmp250["target_match"].to_csv(OUT_DIR / "coord_overlap_per_dss_dmr_dm250.csv",
                                   index=False)
    cmp250["query_match"].to_csv(OUT_DIR / "coord_overlap_per_ek250_dmr_dm250.csv",
                                  index=False)

    # Gene overlap (nearest-TSS gene set, DSS as target)
    g_df = gene_overlap_table(ek100, dss, "dss")
    g_df.to_csv(OUT_DIR / "gene_overlap_ek100_vs_dss.csv", index=False)
    g250_df = gene_overlap_table(ek250, dss, "dss")
    g250_df.to_csv(OUT_DIR / "gene_overlap_ek250_vs_dss.csv", index=False)

    # DSS's recall of paper Table 8 (panel-E critical genes)
    panel_e = pd.read_excel(PAPER_T8, sheet_name=0)
    panel_genes = (panel_e["Gene"].astype(str).str.strip().str.upper()
                                    .dropna().unique().tolist())
    panel_genes = [g for g in panel_genes if g and g != "NAN"]
    dss_genes = set(dss["nearest_tss_gene_u"].dropna().unique()) - {""}
    pe_rows = []
    for g in panel_genes:
        pe_rows.append(dict(
            panel_e_gene=g,
            captured_by_dss_nearest_tss=g in dss_genes,
        ))
    pe_df = pd.DataFrame(pe_rows)
    pe_df.to_csv(OUT_DIR / "panel_e_capture_dss.csv", index=False)
    pe_cap = int(pe_df["captured_by_dss_nearest_tss"].sum())

    # DSS's heatmap-gene direct coordinate hits (paper Fig 3B top 10 hyper/hypo)
    paper = pd.read_excel(PAPER_T5, sheet_name=0).rename(columns={"chr": "chrom"})
    paper["chrom"] = paper["chrom"].astype(str)
    paper["gene_u"] = paper["Gene.Name"].fillna("").astype(str).str.upper()
    # Heatmap labels per paper Fig 3B caption:
    heatmap_names = ["NR2E1", "OTX1", "IRX2", "OTX2", "ENPP2", "GREB1L",
                     "CCDC177", "PAX7", "NAALADL2", "GNG11"]
    dss_idx_by_chrom = _interval_index(dss, "start", "end", "idx")
    hm_rows = []
    for gene in heatmap_names:
        sub = paper[paper["gene_u"] == gene]
        if len(sub) == 0:
            hm_rows.append(dict(gene=gene, paper_dmr=None,
                                 dss_overlap=False, best_jaccard=0.0,
                                 dss_dmr=None))
            continue
        row = sub.sort_values("length" if "length" in sub.columns
                              else "diff.meth_mean", ascending=False).iloc[0]
        paper_str = f"{row['chrom']}:{int(row['start'])}-{int(row['end'])}"
        chrom = str(row["chrom"])
        cands = dss_idx_by_chrom.get(chrom, [])
        best_j = 0.0; best_i = -1
        for ts, te, ti in cands:
            if te < int(row["start"]): continue
            if ts > int(row["end"]): break
            j = jaccard(int(row["start"]), int(row["end"]), ts, te)
            if j > best_j:
                best_j = j; best_i = ti
        if best_i == -1:
            hm_rows.append(dict(gene=gene, paper_dmr=paper_str,
                                 dss_overlap=False, best_jaccard=0.0,
                                 dss_dmr=None))
        else:
            d = dss[dss["idx"] == best_i].iloc[0]
            dss_str = f"{d['chrom']}:{int(d['start'])}-{int(d['end'])}"
            hm_rows.append(dict(gene=gene, paper_dmr=paper_str,
                                 dss_overlap=best_j > 0,
                                 best_jaccard=round(best_j, 4),
                                 dss_dmr=dss_str))
    hm_df = pd.DataFrame(hm_rows)
    hm_df.to_csv(OUT_DIR / "heatmap_gene_hits_dss.csv", index=False)
    hm_hit = int(hm_df["dss_overlap"].sum())

    # ---- Bundle ----------------------------------------------------------

    headline = {
        "ek100_vs_dss": h100,
        "ek250_vs_dss": h250,
        "ek100_gene_recall_of_dss": round(
            int(g_df["captured_by_query"].sum()) / max(len(g_df), 1), 4),
        "ek250_gene_recall_of_dss": round(
            int(g250_df["captured_by_query"].sum()) / max(len(g250_df), 1), 4),
        "dss_panel_e_capture": {
            "n_panel_e": int(len(panel_genes)),
            "captured": pe_cap,
            "recall": round(pe_cap / max(len(panel_genes), 1), 4),
        },
        "dss_heatmap_gene_hits": {
            "checked": int(len(hm_df)),
            "with_paper_dmr": int(hm_df["paper_dmr"].notna().sum()),
            "with_dss_coord_overlap": hm_hit,
        },
    }
    (OUT_DIR / "headline.json").write_text(
        json.dumps(headline, indent=2), encoding="utf-8"
    )
    print("\n--- HEADLINE ---")
    print(json.dumps(headline, indent=2))


if __name__ == "__main__":
    main()
