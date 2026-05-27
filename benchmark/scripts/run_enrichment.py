"""Pathway enrichment for our top DMRs/DMCs vs the original AKAP11 paper.

Uses Enrichr REST API (no auth) to enrich gene lists against:
  - Reactome_2022 (compare to paper panel D: GPCR/Class A receptor pathways)
  - GO_Molecular_Function_2023 (compare to panel E: TF binding / DNA binding)

We submit two gene lists separately:
  - "DMR-gene" set: genes within ±50 kb of any significant epykit DMR (top-N tiles)
  - "DMC-gene" set: genes annotated to the top-N significant epykit DMCs by q

The paper's headline pathways are:
  Panel D: Class A/1 Rhodopsin-like receptors · Peptide ligand-binding receptors ·
           GPCR ligand binding · GPCR downstream signalling · Signalling by GPCR ·
           G alpha i signalling events
  Panel E: Sequence-specific DNA binding · DNA-binding TF activity · RNA Pol II
           cis-regulatory region sequence-specific DNA binding · TF binding ·
           DNA-binding transcription activator activity
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT.parent / "epykit_vs_methylkit(GSE263850)"

ENRICHR = "https://maayanlab.cloud/Enrichr"

# Paper-reported pathway terms to match against
PAPER_PANEL_D_KEYWORDS = [
    "rhodopsin", "peptide ligand", "gpcr ligand", "gpcr downstream",
    "signalling by gpcr", "g alpha i", "g-protein", "g protein",
]
PAPER_PANEL_E_KEYWORDS = [
    "sequence-specific dna binding", "transcription factor",
    "rna polymerase ii cis", "transcription activator",
    "dna-binding transcription",
]


def load_dmc_gene_table() -> pd.DataFrame:
    """Return epykit significant DMCs + annotated gene + sortable q-value."""
    sig = pd.read_csv(SOURCE / "epykit_results" / "dmc_significant_qval05.csv")
    hyper_subset = sig[sig["meth_diff"] > 0].reset_index(drop=True)
    hypo_subset = sig[sig["meth_diff"] < 0].reset_index(drop=True)

    rows = []
    for fn, subset in [
        ("tss_distance_hyper.csv", hyper_subset),
        ("tss_distance_hypo.csv", hypo_subset),
    ]:
        t = pd.read_csv(SOURCE / "epykit_results" / fn)
        for _, r in t.iterrows():
            ridx = int(r["target.row"]) - 1
            if 0 <= ridx < len(subset):
                site = subset.iloc[ridx]
                gene = r["feature.name"]
                if pd.isna(gene) or str(gene).lower() == "nan":
                    continue
                rows.append({
                    "chrom": site["chrom"],
                    "pos_1based": int(site["pos"]) + 1,
                    "meth_diff": float(site["meth_diff"]),
                    "qvalue": float(site["qvalue"]),
                    "gene": str(gene),
                    "dist_to_tss": int(r["dist.to.feature"]),
                })
    return pd.DataFrame(rows)


def genes_from_top_dmcs(dmc_ann: pd.DataFrame, k: int, max_dist_bp: int = 50_000) -> list[str]:
    """Take top-k DMCs by q-value (within ±max_dist of a gene)."""
    df = dmc_ann[dmc_ann["dist_to_tss"].abs() <= max_dist_bp]
    df = df.sort_values("qvalue").head(k)
    return sorted(set(df["gene"].astype(str)))


def genes_from_top_dmrs(dmc_ann: pd.DataFrame, dmr_df: pd.DataFrame, k: int) -> list[str]:
    """For each of the top-k DMRs (by q-value), collect genes within tile."""
    dmr_df = dmr_df.sort_values("qvalue").head(k).copy()
    dmr_df["start_1based"] = dmr_df["start"].astype(int) + 1
    gene_set = set()
    for _, dmr in dmr_df.iterrows():
        chrom = dmr["chrom"]
        s = int(dmr["start_1based"])
        e = int(dmr["end"])
        mask = (
            (dmc_ann["chrom"] == chrom)
            & (dmc_ann["pos_1based"] >= s)
            & (dmc_ann["pos_1based"] <= e)
        )
        for g in dmc_ann.loc[mask, "gene"]:
            if pd.notna(g) and str(g).lower() != "nan":
                gene_set.add(str(g))
    return sorted(gene_set)


def enrichr_upload(gene_list: list[str], label: str) -> int:
    """Upload a gene list to Enrichr, return userListId."""
    payload = {
        "list": (None, "\n".join(gene_list)),
        "description": (None, label),
    }
    r = requests.post(f"{ENRICHR}/addList", files=payload, timeout=30)
    r.raise_for_status()
    return r.json()["userListId"]


def enrichr_enrich(user_list_id: int, library: str) -> list[list]:
    """Fetch enrichment results for a library."""
    r = requests.get(
        f"{ENRICHR}/enrich",
        params={"userListId": user_list_id, "backgroundType": library},
        timeout=60,
    )
    r.raise_for_status()
    return r.json().get(library, [])


def fmt_enrichment(rows: list[list], topn: int = 10) -> pd.DataFrame:
    """Enrichr row schema: [rank, term, pvalue, zscore, comb_score, overlap_genes, adj_p, old_p, old_adj_p]."""
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(
        rows,
        columns=[
            "rank", "term", "pvalue", "zscore", "comb_score",
            "overlap_genes", "adj_p", "old_p", "old_adj_p",
        ],
    )
    df["n_overlap"] = df["overlap_genes"].apply(len)
    df["overlap_str"] = df["overlap_genes"].apply(
        lambda gs: ", ".join(gs[:8]) + (f" (+{len(gs) - 8})" if len(gs) > 8 else "")
    )
    return df.head(topn)


def keyword_match(term: str, keywords: list[str]) -> bool:
    t = term.lower()
    return any(kw in t for kw in keywords)


def main():
    print("Loading DMC annotation...")
    dmc_ann = load_dmc_gene_table()
    print(f"  annotated DMCs: {len(dmc_ann):,}, unique genes: {dmc_ann['gene'].nunique():,}")

    print("\nLoading DMRs...")
    ep_dmr = pd.read_csv(ROOT / "data" / "study3" / "dmr_significant_lenient.csv")
    print(f"  epykit DMRs (lenient): {len(ep_dmr):,}")

    # Build gene sets
    print("\nBuilding gene sets...")
    # Use top-500 DMRs (closer to paper's 813) and the genes within those tiles
    dmr_genes_500 = genes_from_top_dmrs(dmc_ann, ep_dmr, k=500)
    dmr_genes_813 = genes_from_top_dmrs(dmc_ann, ep_dmr, k=813)
    dmc_genes_top1000 = genes_from_top_dmcs(dmc_ann, k=1000, max_dist_bp=50_000)
    print(f"  Genes in top-500 DMRs:  {len(dmr_genes_500):,}")
    print(f"  Genes in top-813 DMRs:  {len(dmr_genes_813):,}")
    print(f"  Genes near top-1000 DMCs (<50 kb): {len(dmc_genes_top1000):,}")

    libraries = [
        ("Reactome_2022", "Panel D (Reactome)", PAPER_PANEL_D_KEYWORDS),
        ("GO_Molecular_Function_2023", "Panel E (GO MF)", PAPER_PANEL_E_KEYWORDS),
    ]

    inputs = [
        ("epykit_top500_DMRs", dmr_genes_500),
        ("epykit_top813_DMRs", dmr_genes_813),
        ("epykit_top1000_DMCs_genes", dmc_genes_top1000),
    ]

    lines = []
    lines.append("# Pathway enrichment — epykit DMRs vs the AKAP11 paper\n")
    lines.append(
        "Submitted gene sets to **Enrichr** (Maayan lab) and pulled the top 10 "
        "enriched terms from **Reactome 2022** (vs paper panel D) and "
        "**GO Molecular Function 2023** (vs paper panel E). Paper-reported "
        "headline terms are flagged with a ★ when a related term appears in our "
        "top 10.\n"
    )
    lines.append("## Paper's reported pathways (target)\n")
    lines.append("**Panel D (Reactome):** Class A/1 Rhodopsin-like receptors · Peptide ligand-binding receptors · GPCR ligand binding · GPCR downstream signalling · Signalling by GPCR · G alpha i signalling events\n")
    lines.append("**Panel E (GO MF):** Sequence-specific DNA binding · DNA-binding TF activity · RNA Pol II cis-regulatory region DNA binding · TF binding · DNA-binding TF activator activity\n")

    all_hits = {}

    for label, genes in inputs:
        print(f"\n--- Running enrichment for {label} ({len(genes):,} genes) ---")
        if len(genes) < 10:
            print("  (too few genes, skipping)")
            continue
        try:
            uid = enrichr_upload(genes, label)
            time.sleep(0.6)  # be polite
        except Exception as e:
            print(f"  upload failed: {e}")
            continue

        lines.append(f"\n## Input gene set: `{label}` ({len(genes):,} genes)\n")

        for lib, lib_label, keywords in libraries:
            try:
                rows = enrichr_enrich(uid, lib)
                time.sleep(0.6)
            except Exception as e:
                print(f"  {lib} failed: {e}")
                continue
            df = fmt_enrichment(rows, topn=15)
            if df.empty:
                lines.append(f"### {lib_label} — no results\n")
                continue

            n_match = sum(keyword_match(t, keywords) for t in df["term"])
            all_hits[(label, lib_label)] = n_match
            lines.append(f"### {lib_label} — top 15 (★ = matches paper)\n")
            lines.append("| Rank | Term | adj p | n / overlap | Sample overlap genes |")
            lines.append("|---:|---|---:|---:|---|")
            for _, r in df.iterrows():
                star = "★ " if keyword_match(r["term"], keywords) else ""
                lines.append(
                    f"| {int(r['rank'])} | {star}{r['term']} | "
                    f"{r['adj_p']:.2e} | {r['n_overlap']} | {r['overlap_str']} |"
                )
            lines.append("")
            print(f"  {lib_label}: {n_match} paper-matching terms in top 15")

    # Summary
    lines.append("## Summary table\n")
    lines.append("| Input gene set | Paper-matching terms in Reactome top 15 | Paper-matching in GO MF top 15 |")
    lines.append("|---|---:|---:|")
    for label, _ in inputs:
        n_d = all_hits.get((label, "Panel D (Reactome)"), "—")
        n_e = all_hits.get((label, "Panel E (GO MF)"), "—")
        lines.append(f"| `{label}` | {n_d} | {n_e} |")

    out = ROOT / "enrichment_vs_paper.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
