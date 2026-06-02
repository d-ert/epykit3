"""Three-way Reactome+KEGG enrichment for the paper-faithful DMR-gene sets.

Compares pathway enrichment across five gene lists vs the paper's
reported terms (Fig 3D Curated.Reactome top 20):
  1. paper-Table-5  (the 705 unique paper DMR-associated genes)
  2. methylKit-tile  (HOMER-annotated nearest-TSS genes from 500 bp tile DMRs)
  3. epykit-chain_merge-100  (genes within 100 kb of any chain_merge DMR)
  4. epykit-chain_merge-250  (same, dis.merge=250)
  5. DSS-from-scratch  (genes within 100 kb of any DSS DMR)

Backend: Enrichr REST API (no auth) against:
  - Reactome_2022
  - KEGG_2021_Human
  - GO_Molecular_Function_2023  (for paper Panel E reference)

Outputs (FINAL_REPORT/data/study3/comparisons/):
  enrichment_three_way.json  — per gene list × library, top-20 terms
                                with p-value, adj p-value, fold enrichment,
                                gene overlap; paper-term match flags
  enrichment_three_way_summary.md — human-readable comparison

Note: we don't try to reproduce the paper's exact "Curated.Reactome"
(ShinyGO-specific) database because the universe / scoring differ.
Enrichr's Reactome_2022 is the closest portable replacement and
preserves term names.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import pandas as pd
import polars as pl
import requests

REPO_ROOT = Path(r"D:/Coding/Projeler/methyl_lib/benchmarkin_merges").resolve()
RAW_DIR   = Path(r"D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW")

CM_LINKS_100   = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "chain_merge" \
                 / "dmr_gene_links_100kb.csv"
CM_LINKS_250   = (REPO_ROOT / "FINAL_REPORT" / "data" / "study3"
                  / "chain_merge_dis_merge_sweep" / "dis_merge_250"
                  / "dmr_gene_links_100kb.csv")
DSS_LINKS      = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "dss" \
                 / "dmr_gene_links_100kb.csv"
MK_TILE        = Path(r"D:/Coding/Projeler/methyl_lib/methylkıt_realResults/"
                      r"scripts_and_results/methylkit_results/"
                      r"dmr_significant_lenient.csv")
PAPER_T5       = RAW_DIR / "Paper resources" / "DMR_total_list.xlsx"

OUT_DIR        = REPO_ROOT / "FINAL_REPORT" / "data" / "study3" / "comparisons"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ENRICHR = "https://maayanlab.cloud/Enrichr"
LIBRARIES = ["Reactome_2022", "KEGG_2021_Human", "GO_Molecular_Function_2023"]

# Paper-reported terms to flag (from Fig 3D + 3E captions, methods text)
PAPER_TERM_KEYWORDS = [
    # Panel D Reactome
    "rhodopsin", "peptide ligand", "gpcr ligand", "gpcr downstream",
    "signaling by gpcr", "signalling by gpcr", "g alpha i", "g alpha (i)",
    "g protein", "extracellular matrix organization", "non-integrin",
    # Panel E GO MF
    "sequence-specific dna binding", "dna-binding transcription factor",
    "transcription factor activity", "rna polymerase ii", "regulatory region",
    "transcription activator",
    # KEGG analogues
    "neuroactive ligand", "camp signaling", "morphine addiction",
    "pluripotency",
]


# ---- Gene list builders ----------------------------------------------------

def genes_paper_table5() -> list[str]:
    df = pd.read_excel(PAPER_T5, sheet_name=0)
    s = (df["Gene.Name"].dropna().astype(str).str.strip().str.upper()
                          .unique().tolist())
    return sorted(set(g for g in s if g and g != "NAN"))


def genes_from_links(csv_path: Path) -> list[str]:
    df = pd.read_csv(csv_path)
    s = (df["gene"].dropna().astype(str).str.strip().str.upper()
                    .unique().tolist())
    return sorted(set(g for g in s if g and g != "NAN"))


def genes_methylkit_tile() -> list[str]:
    """methylKit DMRs annotated to genes — use HOMER-equivalent nearest TSS."""
    # We already produced HOMER-annotated methylKit in P1.3:
    mk_ann = OUT_DIR / "methylkit_dmrs_annotated.csv"
    if not mk_ann.exists():
        # Fall back: annotate now
        from annotation_distribution import (   # type: ignore
            load_refgene, annotate_via_refgene, REFGENE,
        )
        by_chr = load_refgene(REFGENE)
        df = pd.read_csv(MK_TILE)[["chrom", "start", "end"]]
        df = annotate_via_refgene(df, by_chr)
    else:
        df = pd.read_csv(mk_ann)
    # mk_dmrs_annotated.csv has feature_type but not gene name.
    # Re-derive nearest gene now using the same refgene index.
    from annotation_distribution import load_refgene as _lr, REFGENE  # type: ignore
    by_chr = _lr(REFGENE)
    # Build TSS-sorted catalog
    import bisect
    rg_by_chrom: dict[str, list] = {}
    for ch, lst in by_chr.items():
        sorted_g = sorted(lst, key=lambda g: g["tss"])
        rg_by_chrom[ch] = sorted_g
    genes: set[str] = set()
    for _, r in df.iterrows():
        ch = str(r["chrom"])
        center = (int(r["start"]) + int(r["end"])) // 2
        cands = rg_by_chrom.get(ch, [])
        if not cands:
            continue
        positions = [g["tss"] for g in cands]
        i = bisect.bisect_left(positions, center)
        pick = []
        if i > 0:               pick.append(cands[i - 1])
        if i < len(cands):      pick.append(cands[i])
        if not pick:
            continue
        best = min(pick, key=lambda g: abs(g["tss"] - center))
        genes.add(best["gene"].strip().upper())
    return sorted(g for g in genes if g)


# ---- Enrichr submission ----------------------------------------------------

def enrichr_submit(gene_list: list[str], description: str) -> int:
    payload = {
        "list": (None, "\n".join(gene_list)),
        "description": (None, description),
    }
    r = requests.post(f"{ENRICHR}/addList", files=payload, timeout=60)
    r.raise_for_status()
    return r.json()["userListId"]


def enrichr_enrich(user_id: int, library: str) -> list[dict]:
    r = requests.get(f"{ENRICHR}/enrich",
                     params={"userListId": user_id, "backgroundType": library},
                     timeout=120)
    r.raise_for_status()
    payload = r.json()
    rows = payload.get(library, [])
    out = []
    for row in rows:
        # Enrichr schema: [rank, term, p, z, combined, overlap_genes, p_adj, _, _]
        term = row[1]; pval = row[2]; combined = row[4]
        overlap_genes = row[5]; p_adj = row[6]
        out.append({
            "term": term,
            "p_value": float(pval),
            "p_adj": float(p_adj),
            "combined_score": float(combined),
            "n_overlap": len(overlap_genes),
            "overlap_genes": overlap_genes,
        })
    return out


def flag_paper_term(term: str) -> str | None:
    tl = term.lower()
    for kw in PAPER_TERM_KEYWORDS:
        if kw in tl:
            return kw
    return None


# ---- Main ------------------------------------------------------------------

def run_one_list(label: str, genes: list[str], result: dict) -> None:
    if not genes:
        result[label] = {"error": "empty gene list"}
        return
    print(f"\n=== {label} ({len(genes)} genes) ===", flush=True)
    try:
        uid = enrichr_submit(genes, f"three_way::{label}")
    except Exception as e:
        result[label] = {"error": str(e)}
        return
    result[label] = {"n_genes": len(genes), "userListId": uid, "libraries": {}}
    for lib in LIBRARIES:
        time.sleep(1)
        try:
            rows = enrichr_enrich(uid, lib)
        except Exception as e:
            result[label]["libraries"][lib] = {"error": str(e)}
            continue
        top20 = sorted(rows, key=lambda r: r["p_value"])[:20]
        for r in top20:
            r["paper_term_match"] = flag_paper_term(r["term"])
        n_sig = sum(1 for r in top20 if r["p_adj"] < 0.05)
        n_paper_hits = sum(1 for r in top20 if r["paper_term_match"])
        result[label]["libraries"][lib] = {
            "n_total": len(rows),
            "top20": top20,
            "n_sig_in_top20": n_sig,
            "n_paper_term_matches_top20": n_paper_hits,
        }
        print(f"  {lib}: top-20 sig={n_sig}, paper-term hits={n_paper_hits}",
              flush=True)
        if n_paper_hits:
            for r in top20[:8]:
                if r["paper_term_match"]:
                    print(f"    ✓ {r['term'][:70]:70}  "
                          f"p_adj={r['p_adj']:.2e}  "
                          f"n={r['n_overlap']}",
                          flush=True)


def main() -> None:
    result: dict = {}

    # ---- Build gene lists ----
    print("Building gene lists …", flush=True)
    g_paper = genes_paper_table5()
    print(f"  paper-Table5: {len(g_paper)} genes", flush=True)
    g_ek100 = genes_from_links(CM_LINKS_100)
    print(f"  ek-chain_merge-100 (100kb): {len(g_ek100)} genes", flush=True)
    g_ek250 = genes_from_links(CM_LINKS_250)
    print(f"  ek-chain_merge-250 (100kb): {len(g_ek250)} genes", flush=True)
    g_dss = genes_from_links(DSS_LINKS)
    print(f"  DSS (100kb): {len(g_dss)} genes", flush=True)
    g_mk = genes_methylkit_tile()
    print(f"  methylKit-tile (nearest-TSS): {len(g_mk)} genes", flush=True)

    # Save them to shinygo_lists/ for posterity
    out_list_dir = REPO_ROOT / "FINAL_REPORT" / "shinygo_lists"
    (out_list_dir / "dss_100kb_genes.txt").write_text(
        "\n".join(g_dss), encoding="utf-8")
    (out_list_dir / "methylkit_tile_nearestTSS_genes.txt").write_text(
        "\n".join(g_mk), encoding="utf-8")

    # ---- Run enrichments ----
    run_one_list("paper_Table5", g_paper, result)
    run_one_list("methylKit_tile", g_mk, result)
    run_one_list("ek_chain_merge_100", g_ek100, result)
    run_one_list("ek_chain_merge_250", g_ek250, result)
    run_one_list("DSS_from_scratch", g_dss, result)

    out_json = OUT_DIR / "enrichment_three_way.json"
    out_json.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"\nWrote {out_json.name}", flush=True)

    # ---- Markdown summary ----
    L = ["# Three-way pathway enrichment summary\n",
         "Gene lists from paper Supp Table 5 (705 paper-named genes), our DMR "
         "callers (100 kb linkage for chain_merge/DSS, nearest-TSS for "
         "methylKit-tile). Enrichr REST API; Reactome_2022 + KEGG_2021_Human "
         "+ GO_MF_2023; top 20 by p-value per library.\n",
         "\n## n_paper_term_matches in top-20 (per library, per caller)\n",
         "| Caller | Reactome | KEGG | GO MF | n_genes |",
         "|---|---:|---:|---:|---:|"]
    for caller in ("paper_Table5", "methylKit_tile",
                    "ek_chain_merge_100", "ek_chain_merge_250",
                    "DSS_from_scratch"):
        entry = result.get(caller, {})
        libs = entry.get("libraries", {})
        rc = libs.get("Reactome_2022", {}).get("n_paper_term_matches_top20", 0)
        kg = libs.get("KEGG_2021_Human", {}).get("n_paper_term_matches_top20", 0)
        gm = libs.get("GO_Molecular_Function_2023", {}).get(
            "n_paper_term_matches_top20", 0)
        L.append(f"| {caller} | {rc} | {kg} | {gm} | {entry.get('n_genes', 0)} |")
    L.append("")
    L.append("(See enrichment_three_way.json for the full top-20 + overlap "
             "genes per cell.)")
    (OUT_DIR / "enrichment_three_way_summary.md").write_text("\n".join(L),
                                                              encoding="utf-8")


if __name__ == "__main__":
    main()
