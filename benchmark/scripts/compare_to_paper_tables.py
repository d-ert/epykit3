"""Compare our DMR calls and enrichments against the paper's supplementary tables.

Tests we run:
  1. DMR coordinate overlap: paper's 813 DMRs vs our epykit / methylKit DMRs.
  2. DMR-associated gene overlap: paper's 705 unique gene annotations vs ours.
  3. Panel E gene-list capture: does our DMR-near-gene set contain the 46
     genes the paper used for the Panel E (TF binding) enrichment?
  4. For each paper-named heatmap gene: do we have a DMR within 1 kb of the
     paper's exact DMR coordinate?
"""

from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT.parent / "epykit_vs_methylkit(GSE263850)"
MK_SOURCE = ROOT.parent.parent / "methylkıt_realResults" / "scripts_and_results" / "methylkit_results"
SUPP = ROOT / "shinygo_lists" / "outputs" / "reactome"


def load_paper_dmrs():
    t5 = pd.read_excel(SUPP / "table5.xlsx")
    return t5


def load_paper_dmr_deg():
    t8 = pd.read_excel(SUPP / "table8.xlsx")
    return t8


def load_our_dmrs():
    ep_str = pd.read_csv(ROOT / "data" / "study3" / "dmr_significant_strict.csv")
    ep_str["start_1based"] = ep_str["start"].astype(int) + 1
    ep_str = ep_str.sort_values("qvalue").reset_index(drop=True)
    ep_len = pd.read_csv(ROOT / "data" / "study3" / "dmr_significant_lenient.csv")
    ep_len["start_1based"] = ep_len["start"].astype(int) + 1
    ep_len = ep_len.sort_values("qvalue").reset_index(drop=True)

    mk_str = pd.read_csv(MK_SOURCE / "dmr_significant_strict.csv")
    mk_str["start_1based"] = mk_str["start"].astype(int)
    mk_str = mk_str.sort_values("qvalue").reset_index(drop=True)
    mk_len = pd.read_csv(MK_SOURCE / "dmr_significant_lenient.csv")
    mk_len["start_1based"] = mk_len["start"].astype(int)
    mk_len = mk_len.sort_values("qvalue").reset_index(drop=True)
    return ep_str, ep_len, mk_str, mk_len


def overlap_any(paper, ours):
    """For each paper DMR, does any of ours overlap it (positions in 1-based, inclusive)?"""
    n_paper_covered = 0
    for _, p in paper.iterrows():
        chrom = p["chr"]
        s = int(p["start"])
        e = int(p["end"])
        sub = ours[ours["chrom"] == chrom]
        # interval intersect: our.start <= p.end and our.end >= p.start
        hit = (
            (sub["start_1based"].astype(int) <= e)
            & (sub["end"].astype(int) >= s)
        ).any()
        if hit:
            n_paper_covered += 1
    return n_paper_covered


def main():
    paper = load_paper_dmrs()
    paper_deg = load_paper_dmr_deg()
    ep_str, ep_len, mk_str, mk_len = load_our_dmrs()

    lines = []
    lines.append("# Comparison against the paper's supplementary tables\n")
    lines.append(
        f"Loaded **Supp Table 5** (n = {len(paper):,} DMRs, "
        f"{paper['Gene.Name'].nunique():,} unique genes) and "
        f"**Supp Table 8** (n = {len(paper_deg):,} DMR-DEG rows, "
        f"{paper_deg['Gene'].nunique():,} unique genes).\n"
    )

    # ---- DMR size distribution -----------------------------------------
    lines.append("## DMR size distribution\n")
    lines.append("| Source | n DMRs | median len | mean len | min | max |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for name, df in [
        ("Paper (DSS callDMR)", paper),
        ("epykit lenient (500 bp tiles)", ep_len),
        ("epykit strict", ep_str),
        ("methylKit lenient", mk_len),
        ("methylKit strict", mk_str),
    ]:
        if "length" in df.columns:
            lens = df["length"]
        else:
            lens = df["end"].astype(int) - df["start_1based"].astype(int) + 1
        lines.append(
            f"| {name} | {len(df):,} | {int(lens.median()):,} | "
            f"{int(lens.mean()):,} | {int(lens.min()):,} | {int(lens.max()):,} |"
        )
    lines.append("")

    # ---- DMR coordinate overlap -----------------------------------------
    lines.append("## DMR coordinate overlap\n")
    lines.append(
        "For each paper DMR (Supp Table 5), is there any overlap with our "
        "called DMRs (interval intersection)?\n"
    )

    paper_n = len(paper)
    paper_in_ep_strict = overlap_any(paper, ep_str)
    paper_in_ep_lenient = overlap_any(paper, ep_len)
    paper_in_mk_strict = overlap_any(paper, mk_str)
    paper_in_mk_lenient = overlap_any(paper, mk_len)

    def reverse_overlap(theirs_df, our_df):
        n_covered = 0
        for _, o in our_df.iterrows():
            chrom = o["chrom"]
            s = int(o["start_1based"])
            e = int(o["end"])
            sub = theirs_df[theirs_df["chr"] == chrom]
            if ((sub["start"].astype(int) <= e) & (sub["end"].astype(int) >= s)).any():
                n_covered += 1
        return n_covered

    ep_str_overlap_paper = reverse_overlap(paper, ep_str)
    ep_len_overlap_paper = reverse_overlap(paper, ep_len)
    mk_str_overlap_paper = reverse_overlap(paper, mk_str)
    mk_len_overlap_paper = reverse_overlap(paper, mk_len)

    lines.append("| Our call set | n DMRs | recall of paper (paper ⊆ ours) | precision (ours overlap paper) |")
    lines.append("|---|---:|---:|---:|")
    lines.append(
        f"| epykit lenient | {len(ep_len):,} | "
        f"{paper_in_ep_lenient:,} / {paper_n:,} = **{paper_in_ep_lenient/paper_n*100:.1f}%** | "
        f"{ep_len_overlap_paper:,} / {len(ep_len):,} = {ep_len_overlap_paper/len(ep_len)*100:.1f}% |"
    )
    lines.append(
        f"| epykit strict | {len(ep_str):,} | "
        f"{paper_in_ep_strict:,} / {paper_n:,} = **{paper_in_ep_strict/paper_n*100:.1f}%** | "
        f"{ep_str_overlap_paper:,} / {len(ep_str):,} = {ep_str_overlap_paper/len(ep_str)*100:.1f}% |"
    )
    lines.append(
        f"| methylKit lenient | {len(mk_len):,} | "
        f"{paper_in_mk_lenient:,} / {paper_n:,} = **{paper_in_mk_lenient/paper_n*100:.1f}%** | "
        f"{mk_len_overlap_paper:,} / {len(mk_len):,} = {mk_len_overlap_paper/len(mk_len)*100:.1f}% |"
    )
    lines.append(
        f"| methylKit strict | {len(mk_str):,} | "
        f"{paper_in_mk_strict:,} / {paper_n:,} = **{paper_in_mk_strict/paper_n*100:.1f}%** | "
        f"{mk_str_overlap_paper:,} / {len(mk_str):,} = {mk_str_overlap_paper/len(mk_str)*100:.1f}% |"
    )
    lines.append("")

    # ---- DMR-associated gene overlap -----------------------------------
    paper_genes = set(paper["Gene.Name"].dropna().astype(str).unique())
    # Our gene lists (built earlier in shinygo_lists)
    list_dir = ROOT / "shinygo_lists"
    def load_list(fn):
        out = set()
        for line in (list_dir / fn).read_text(encoding="utf-8").splitlines():
            if line.startswith("#") or not line.strip():
                continue
            out.add(line.strip())
        return out

    our_lists = {
        "epykit strict (1,079 genes)": load_list("epykit_dmr_strict_genes.txt"),
        "epykit top-813 (2,493 genes)": load_list("epykit_dmr_top813_genes.txt"),
        "epykit top-500 (1,840 genes)": load_list("epykit_dmr_top500_genes.txt"),
        "methylKit strict (639 genes)": load_list("methylkit_dmr_strict_genes.txt"),
        "methylKit top-813 (2,319 genes)": load_list("methylkit_dmr_top813_genes.txt"),
    }

    lines.append("## DMR-associated gene-list overlap\n")
    lines.append(
        f"Paper Supp Table 5 has **{len(paper_genes):,}** unique DMR-associated genes "
        "(HOMER nearest gene). How many of these appear in our DMR-associated gene sets "
        "(built using 100 kb-TSS-to-DMR-midpoint as in the paper)?\n"
    )
    lines.append("| Our gene set | size | & paper 705 | recall of paper |")
    lines.append("|---|---:|---:|---:|")
    for name, our in our_lists.items():
        inter = paper_genes & our
        rec = len(inter) / len(paper_genes)
        lines.append(f"| {name} | {len(our):,} | {len(inter):,} | **{rec*100:.1f}%** |")
    lines.append("")

    # ---- Panel E gene list (Table 8) check -----------------------------
    panel_e_genes = set(paper_deg["Gene"].dropna().astype(str).unique())
    lines.append("## Panel E gene-list capture (the 46 critical genes)\n")
    lines.append(
        f"Panel E (TF binding enrichment) was computed on **{len(panel_e_genes):,} "
        "DMR-near-DEG genes** (Supp Table 8). Do we have these genes in our calls?\n"
    )
    lines.append("| Our gene set | & Table 8 (46 genes) | recall |")
    lines.append("|---|---:|---:|")
    for name, our in our_lists.items():
        inter = panel_e_genes & our
        rec = len(inter) / len(panel_e_genes)
        lines.append(f"| {name} | {len(inter)} / {len(panel_e_genes)} | **{rec*100:.0f}%** |")
    lines.append("")

    # Show which genes we miss
    best_our = our_lists["epykit top-813 (2,493 genes)"]
    missed = sorted(panel_e_genes - best_our)
    captured = sorted(panel_e_genes & best_our)
    lines.append(f"### Captured by epykit top-813 ({len(captured)} / {len(panel_e_genes)})\n")
    lines.append(", ".join(captured) + "\n")
    if missed:
        lines.append(f"### Missed by epykit top-813 ({len(missed)})\n")
        lines.append(", ".join(missed) + "\n")
    lines.append("")

    # ---- Paper-named heatmap genes — direct DMR coord match ------------
    HEATMAP_GENES = [
        "NR2E1", "OTX1", "IRX2", "OTX2", "ENPP2", "GREB1L", "CCDC177",
        "PAX7", "NAALADL2", "PDK3", "TMEM242", "OSBPL8", "GNG11",
        "KC6", "RPLP0P2", "LOC100506858", "LOC100131655",
    ]
    lines.append("## Paper heatmap genes — direct DMR coordinate check\n")
    lines.append(
        "For each gene named in the paper's Figure 6B heatmap, look up its DMR "
        "coordinate in Supp Table 5 and check whether we called a DMR at that "
        "exact region.\n"
    )
    lines.append("| Gene | Paper DMR coord | length | epykit DMR (any q) overlap? | methylKit DMR (any q) overlap? |")
    lines.append("|---|---|---:|:---:|:---:|")
    for gene in HEATMAP_GENES:
        rows = paper[paper["Gene.Name"] == gene]
        if not len(rows):
            lines.append(f"| {gene} | not in Table 5 | — | — | — |")
            continue
        # take strongest by area stat
        r = rows.iloc[0]
        chrom = r["chr"]; s = int(r["start"]); e = int(r["end"]); L = int(r["length"])
        in_ep = (
            (ep_len["chrom"] == chrom)
            & (ep_len["start_1based"].astype(int) <= e)
            & (ep_len["end"].astype(int) >= s)
        ).any()
        in_mk = (
            (mk_len["chrom"] == chrom)
            & (mk_len["start_1based"].astype(int) <= e)
            & (mk_len["end"].astype(int) >= s)
        ).any()
        lines.append(
            f"| **{gene}** | {chrom}:{s:,}–{e:,} | {L:,} bp | "
            f"{'✓' if in_ep else ' '} | {'✓' if in_mk else ' '} |"
        )
    lines.append("")

    # ---- Summary --------------------------------------------------------
    lines.append("## Headline numbers\n")
    rec_ep = paper_in_ep_lenient / paper_n * 100
    rec_mk = paper_in_mk_lenient / paper_n * 100
    rec_genes_ep = len(paper_genes & our_lists["epykit top-813 (2,493 genes)"]) / len(paper_genes) * 100
    rec_panel_e_ep = len(panel_e_genes & our_lists["epykit top-813 (2,493 genes)"]) / len(panel_e_genes) * 100
    lines.append(
        f"- **epykit lenient recalls {rec_ep:.0f}% of paper DMRs** (interval overlap)\n"
        f"- **methylKit lenient recalls {rec_mk:.0f}% of paper DMRs**\n"
        f"- **epykit top-813 captures {rec_genes_ep:.0f}% of the paper's 705 DMR-associated genes**\n"
        f"- **epykit top-813 captures {rec_panel_e_ep:.0f}% of the 46 Panel E genes**\n"
    )

    out = ROOT / "paper_table_comparison.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out}")

    # console summary
    print(f"\nDMR overlap (interval):")
    print(f"  epykit lenient covers {paper_in_ep_lenient}/{paper_n} = {rec_ep:.1f}% of paper DMRs")
    print(f"  methylKit lenient covers {paper_in_mk_lenient}/{paper_n} = {rec_mk:.1f}% of paper DMRs")
    print(f"\nGene overlap:")
    print(f"  epykit top-813 & paper 705 genes: {len(paper_genes & our_lists['epykit top-813 (2,493 genes)'])}/{len(paper_genes)} = {rec_genes_ep:.1f}%")
    print(f"  epykit top-813 & Panel E 46 genes: {len(panel_e_genes & our_lists['epykit top-813 (2,493 genes)'])}/{len(panel_e_genes)} = {rec_panel_e_ep:.1f}%")


if __name__ == "__main__":
    main()
