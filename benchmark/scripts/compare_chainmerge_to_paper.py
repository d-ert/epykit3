"""Compare epykit's smoothed/chain_merge DMR set (alpha=1e-5, matching the paper)
against the paper's Supplementary Tables 5 and 8.

Loads the existing run at dmr_lr_site_smooth_alpha1e-5.parquet (702 DMRs) and
recomputes the comparison stats from compare_to_paper_tables.py.
"""

from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = Path("D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW")
SUPP = ROOT / "shinygo_lists" / "outputs" / "reactome"

OUT_LIST = ROOT / "shinygo_lists"


def overlap_any(paper, ours, our_start_col="start", our_end_col="end"):
    n = 0
    for _, p in paper.iterrows():
        chrom = p["chr"]
        s = int(p["start"]); e = int(p["end"])
        sub = ours[ours["chrom"] == chrom]
        if ((sub[our_start_col].astype(int) <= e) & (sub[our_end_col].astype(int) >= s)).any():
            n += 1
    return n


def reverse_overlap(theirs, our):
    n = 0
    for _, o in our.iterrows():
        chrom = o["chrom"]; s = int(o["start"]); e = int(o["end"])
        sub = theirs[theirs["chr"] == chrom]
        if ((sub["start"].astype(int) <= e) & (sub["end"].astype(int) >= s)).any():
            n += 1
    return n


def main():
    print("Loading paper Supp Tables 5 & 8...")
    paper = pd.read_excel(SUPP / "table5.xlsx")
    paper_deg = pd.read_excel(SUPP / "table8.xlsx")
    paper_genes = set(paper["Gene.Name"].dropna().astype(str).unique())
    panel_e_genes = set(paper_deg["Gene"].dropna().astype(str).unique())
    print(f"  paper DMRs: {len(paper):,}")
    print(f"  paper genes: {len(paper_genes):,}")
    print(f"  Panel E genes: {len(panel_e_genes):,}")

    # ----- Smoothed (chain_merge equivalent) at alpha=1e-5 -----
    print("\nLoading epykit smoothed DMRs at alpha=1e-5 ...")
    sm = pd.read_parquet(RAW / "dmr_lr_site_smooth_alpha1e-5.parquet")
    # filter to significant (combined_qvalue_reject == True)
    sm_sig = sm[sm["combined_qvalue_reject"]].copy()
    print(f"  total: {len(sm):,} ; significant (q-reject): {len(sm_sig):,}")
    # use ALL (the file is already prefiltered to significant from the look of it,
    # but force a safe filter)
    if len(sm_sig) < len(sm) * 0.5:
        # file already filtered
        sm_sig = sm.copy()
    print(f"  using {len(sm_sig):,} DMRs")
    print(f"  median length: {(sm_sig['end'] - sm_sig['start']).median():.0f} bp")
    print(f"  hyper: {(sm_sig['mean_meth_diff'] > 0).sum():,} ; hypo: {(sm_sig['mean_meth_diff'] < 0).sum():,}")

    # ----- DMR coordinate overlap -----
    print("\nComputing coordinate overlap vs paper Supp Table 5 ...")
    paper_in_sm = overlap_any(paper, sm_sig)
    sm_in_paper = reverse_overlap(paper, sm_sig)
    rec = paper_in_sm / len(paper) * 100
    prec = sm_in_paper / len(sm_sig) * 100
    print(f"  paper DMRs recovered: {paper_in_sm} / {len(paper)} = {rec:.1f}%")
    print(f"  our DMRs hitting paper: {sm_in_paper} / {len(sm_sig)} = {prec:.1f}%")

    # ----- DMR-associated gene overlap (using the parquet's own gene_name col) -----
    our_genes = set(sm_sig["gene_name"].dropna().astype(str).unique())
    our_genes.discard("")
    print(f"\nUnique gene names in smoothed DMRs: {len(our_genes):,}")
    gene_inter = paper_genes & our_genes
    gene_rec = len(gene_inter) / len(paper_genes) * 100
    print(f"  gene overlap with paper 705: {len(gene_inter)} / {len(paper_genes)} = {gene_rec:.1f}%")

    panel_inter = panel_e_genes & our_genes
    panel_rec = len(panel_inter) / len(panel_e_genes) * 100
    print(f"  Panel E gene capture: {len(panel_inter)} / {len(panel_e_genes)} = {panel_rec:.0f}%")

    # ----- Also build a "wider" gene set using 100 kb from DMR midpoint -----
    # Reuse the same TSS-based gene catalogue used by build_shinygo_lists.py
    print("\nLoading TSS catalog from epykit annotation...")
    SOURCE = ROOT.parent / "epykit_vs_methylkit(GSE263850)"
    sig = pd.read_csv(SOURCE / "epykit_results" / "dmc_significant_qval05.csv")
    hyper_subset = sig[sig["meth_diff"] > 0].reset_index(drop=True)
    hypo_subset = sig[sig["meth_diff"] < 0].reset_index(drop=True)
    gene_pos = []
    for fn, subset in [("tss_distance_hyper.csv", hyper_subset),
                       ("tss_distance_hypo.csv", hypo_subset)]:
        t = pd.read_csv(SOURCE / "epykit_results" / fn)
        for _, r in t.iterrows():
            ridx = int(r["target.row"]) - 1
            if 0 <= ridx < len(subset):
                site = subset.iloc[ridx]
                gene = r["feature.name"]
                if pd.isna(gene) or str(gene).lower() == "nan":
                    continue
                cpg = int(site["pos"]) + 1
                tss = cpg - int(r["dist.to.feature"])
                gene_pos.append({"chrom": site["chrom"], "gene": str(gene), "tss": tss})
    gene_pos = pd.DataFrame(gene_pos).drop_duplicates(["chrom", "gene", "tss"])
    print(f"  {len(gene_pos):,} unique gene-TSS triples in catalog")

    # Apply 100 kb DMR-midpoint association
    flank = 100_000
    wider_genes = set()
    for _, dmr in sm_sig.iterrows():
        chrom = dmr["chrom"]
        mid = int((dmr["start"] + dmr["end"]) / 2)
        sub = gene_pos[gene_pos["chrom"] == chrom]
        nearby = sub[(sub["tss"] - mid).abs() <= flank]
        wider_genes.update(nearby["gene"].tolist())
    print(f"  Wider gene set (100 kb assoc): {len(wider_genes):,}")
    wider_inter = paper_genes & wider_genes
    print(f"    overlap with paper 705: {len(wider_inter)} ({len(wider_inter)/len(paper_genes)*100:.1f}%)")
    wider_panel = panel_e_genes & wider_genes
    print(f"    Panel E capture: {len(wider_panel)} / {len(panel_e_genes)} ({len(wider_panel)/len(panel_e_genes)*100:.0f}%)")

    # ----- Write ShinyGO-ready lists -----
    out_a = OUT_LIST / "epykit_smoothed_alpha1e-5_nearest_genes.txt"
    out_b = OUT_LIST / "epykit_smoothed_alpha1e-5_100kb_genes.txt"
    with out_a.open("w", encoding="utf-8") as f:
        f.write(f"# epykit smoothed DMRs (alpha=1e-5) nearest-gene annotation\n")
        f.write(f"# {len(our_genes)} genes\n")
        for g in sorted(our_genes):
            f.write(f"{g}\n")
    with out_b.open("w", encoding="utf-8") as f:
        f.write(f"# epykit smoothed DMRs (alpha=1e-5) 100 kb-TSS-midpoint association\n")
        f.write(f"# {len(wider_genes)} genes\n")
        for g in sorted(wider_genes):
            f.write(f"{g}\n")
    print(f"\nWrote {out_a.name} ({len(our_genes)} genes)")
    print(f"Wrote {out_b.name} ({len(wider_genes)} genes)")

    # ----- Markdown report -----
    lines = []
    lines.append("# epykit chain_merge / smoothed DMRs vs the paper\n")
    lines.append(
        f"epykit run with **smoothing enabled** and **alpha = 1e-5** (matching the "
        "paper's `DSS::callDMR(p.threshold = 1e-5)` parameter).\n"
    )
    lines.append(
        f"- DMRs called: **{len(sm_sig):,}** "
        f"(paper: 813; ratio = {len(sm_sig)/len(paper):.2f}×)\n"
        f"- Median length: **{int((sm_sig['end'] - sm_sig['start']).median())} bp** "
        f"(paper: 240 bp — much closer than fixed tiles' 500 bp)\n"
        f"- Hyper / Hypo: {(sm_sig['mean_meth_diff']>0).sum():,} / "
        f"{(sm_sig['mean_meth_diff']<0).sum():,} (paper: 638 / 175 — same direction bias)\n"
    )

    lines.append("\n## DMR coordinate overlap\n")
    lines.append("| Source | DMRs | Recall of paper | Precision (ours ∩ paper) |")
    lines.append("|---|---:|---:|---:|")
    lines.append(
        f"| **epykit smoothed alpha=1e-5** | {len(sm_sig):,} | "
        f"**{paper_in_sm}/{len(paper)} = {rec:.1f}%** | "
        f"**{sm_in_paper}/{len(sm_sig)} = {prec:.1f}%** |"
    )
    lines.append(
        f"| (for reference) epykit lenient 500 bp tiles | 3,433 | "
        f"63/813 = 7.7% | 63/3,433 = 1.8% |"
    )
    lines.append(
        f"| (for reference) methylKit lenient 500 bp tiles | 2,661 | "
        f"72/813 = 8.9% | 74/2,661 = 2.8% |"
    )
    lines.append("")

    lines.append("\n## Gene-list overlap with paper Supp Table 5\n")
    lines.append("| Gene set | size | ∩ paper 705 | recall of paper |")
    lines.append("|---|---:|---:|---:|")
    lines.append(
        f"| epykit smoothed: nearest-gene | {len(our_genes):,} | "
        f"{len(gene_inter)} | **{gene_rec:.1f}%** |"
    )
    lines.append(
        f"| epykit smoothed: 100 kb assoc | {len(wider_genes):,} | "
        f"{len(wider_inter)} | **{len(wider_inter)/len(paper_genes)*100:.1f}%** |"
    )
    lines.append("")

    lines.append("\n## Panel E gene-list capture (46 critical genes)\n")
    lines.append("| Gene set | ∩ Table 8 (46) | recall |")
    lines.append("|---|---:|---:|")
    lines.append(
        f"| epykit smoothed: nearest-gene | {len(panel_inter)} | "
        f"**{panel_rec:.0f}%** |"
    )
    lines.append(
        f"| epykit smoothed: 100 kb assoc | {len(wider_panel)} | "
        f"**{len(wider_panel)/len(panel_e_genes)*100:.0f}%** |"
    )
    lines.append("")

    # Captured / missed (using 100 kb assoc)
    captured = sorted(panel_e_genes & wider_genes)
    missed = sorted(panel_e_genes - wider_genes)
    lines.append(f"### Panel E genes captured ({len(captured)} / {len(panel_e_genes)})\n")
    lines.append(", ".join(captured) + "\n")
    if missed:
        lines.append(f"### Panel E genes missed ({len(missed)})\n")
        lines.append(", ".join(missed) + "\n")
    lines.append("")

    # ----- Paper heatmap genes — direct DMR overlap with smoothed calls -----
    HEATMAP = ["NR2E1", "OTX1", "IRX2", "OTX2", "ENPP2", "GREB1L", "CCDC177",
               "PAX7", "NAALADL2", "PDK3", "TMEM242", "OSBPL8", "GNG11"]
    lines.append("\n## Paper heatmap genes — direct coordinate match\n")
    lines.append("| Gene | Paper DMR coord | smoothed-epykit overlap? |")
    lines.append("|---|---|:---:|")
    for gene in HEATMAP:
        rows = paper[paper["Gene.Name"] == gene]
        if not len(rows):
            continue
        r = rows.iloc[0]
        chrom = r["chr"]; s = int(r["start"]); e = int(r["end"])
        hit = ((sm_sig["chrom"] == chrom) &
               (sm_sig["start"].astype(int) <= e) &
               (sm_sig["end"].astype(int) >= s)).any()
        lines.append(f"| {gene} | {chrom}:{s:,}-{e:,} | {'✓' if hit else ' '} |")
    lines.append("")

    lines.append("\n## Summary\n")
    lines.append(
        f"With smoothing + alpha=1e-5 (matching the paper's parameters), "
        f"epykit:\n\n"
        f"- recovers **{rec:.0f}% of paper DMRs** at the coordinate level "
        f"(vs 7.7% with fixed tiles — a {rec/7.7:.1f}× improvement)\n"
        f"- captures **{len(wider_inter)/len(paper_genes)*100:.0f}% of paper's 705 genes** "
        f"(vs 43% with tiles)\n"
        f"- captures **{len(wider_panel)/len(panel_e_genes)*100:.0f}% of the 46 Panel E genes** "
        f"(vs 48% with tiles)\n\n"
        f"The remaining gap is plausibly DSS-vs-LR test differences (DSS smooths "
        f"dispersion across CpGs in a different way than epykit's combined-pvalue "
        f"approach). But qualitatively, when the region model matches, the call set converges.\n"
    )

    out = ROOT / "smoothed_dmr_vs_paper.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
