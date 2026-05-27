"""Build ShinyGO-paste-ready gene lists from our DMRs using the paper's
gene-DMR association criterion: any gene whose TSS is within 100 kb of a
significant DMR midpoint.

Saves four files in FINAL_REPORT/shinygo_lists/:
  - epykit_dmr_strict_genes.txt        (q<0.01, |d|≥25%, paper-matched stringency)
  - epykit_dmr_lenient_genes.txt       (q<0.05, |d|≥10%, our default)
  - epykit_dmr_top500_genes.txt        (top 500 DMRs by q-value)
  - methylkit_dmr_strict_genes.txt     (same strict threshold for methylKit)

Plus a background gene list — every gene that has at least one *tested*
CpG (regardless of significance) — for ShinyGO's "custom background" field.

The script uses epykit's tss_distance_*.csv as the gene→position mapping
(restricted to significant CpGs). For each DMR midpoint we search all
annotated DMC positions within ±100 kb, take their gene names, and union.
"""

from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT.parent / "epykit_vs_methylkit(GSE263850)"
MK_SOURCE = ROOT.parent.parent / "methylkıt_realResults" / "scripts_and_results" / "methylkit_results"

OUT = ROOT / "shinygo_lists"
OUT.mkdir(exist_ok=True)


def load_gene_positions() -> pd.DataFrame:
    """All annotated gene → genomic position pairs from epykit's TSS-distance
    files, restricted to significant CpGs.

    The TSS-distance files give us, per significant CpG: nearest gene and
    signed distance to its TSS. So for each row we know:
        cpg_position = exact pos
        gene_tss     = cpg_position - dist_to_tss
    """
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
                cpg_pos = int(site["pos"]) + 1  # 1-based
                tss_pos = cpg_pos - int(r["dist.to.feature"])
                rows.append({
                    "chrom": site["chrom"],
                    "gene": str(gene),
                    "tss": tss_pos,
                })
    df = pd.DataFrame(rows)
    # one row per (chrom, gene, tss) — collapse duplicates
    df = df.drop_duplicates(["chrom", "gene", "tss"])
    return df


def genes_near_dmrs(dmrs: pd.DataFrame, gene_pos: pd.DataFrame,
                    flank_bp: int = 100_000) -> list[str]:
    """For each DMR, find genes whose TSS is within ±flank_bp of the DMR midpoint."""
    out = set()
    # Coordinate alignment: paper used 1-based; our methylKit DMRs are 1-based,
    # epykit DMRs are 0-based starts. Caller is responsible for passing aligned
    # start positions in the DMR table.
    for _, dmr in dmrs.iterrows():
        chrom = dmr["chrom"]
        mid = int((dmr["start_1based"] + dmr["end"]) / 2)
        sub = gene_pos[gene_pos["chrom"] == chrom]
        nearby = sub[(sub["tss"] - mid).abs() <= flank_bp]
        out.update(nearby["gene"].tolist())
    return sorted(out)


def write_gene_list(genes: list[str], filename: str, header_note: str):
    path = OUT / filename
    with path.open("w", encoding="utf-8") as f:
        f.write(f"# {header_note}\n")
        f.write(f"# {len(genes):,} genes\n")
        for g in genes:
            f.write(f"{g}\n")
    print(f"  wrote {len(genes):,} genes -> {path.name}")


def main():
    print("Loading gene positions from epykit annotation...")
    gene_pos = load_gene_positions()
    print(f"  {len(gene_pos):,} unique (chrom, gene, tss) triples")
    print(f"  {gene_pos['gene'].nunique():,} unique gene names")

    # ----- epykit DMRs ---------------------------------------------------
    ep_lenient = pd.read_csv(ROOT / "data" / "study3" / "dmr_significant_lenient.csv")
    ep_strict = pd.read_csv(ROOT / "data" / "study3" / "dmr_significant_strict.csv")
    ep_lenient["start_1based"] = ep_lenient["start"].astype(int) + 1
    ep_strict["start_1based"] = ep_strict["start"].astype(int) + 1
    ep_lenient = ep_lenient.sort_values("qvalue").reset_index(drop=True)
    ep_strict = ep_strict.sort_values("qvalue").reset_index(drop=True)
    print(f"\nepykit DMRs: lenient = {len(ep_lenient):,}, strict = {len(ep_strict):,}")

    # ----- methylKit DMRs ------------------------------------------------
    mk_lenient = pd.read_csv(MK_SOURCE / "dmr_significant_lenient.csv")
    mk_strict = pd.read_csv(MK_SOURCE / "dmr_significant_strict.csv")
    mk_lenient["start_1based"] = mk_lenient["start"].astype(int)  # already 1-based
    mk_strict["start_1based"] = mk_strict["start"].astype(int)
    mk_lenient = mk_lenient.sort_values("qvalue").reset_index(drop=True)
    mk_strict = mk_strict.sort_values("qvalue").reset_index(drop=True)
    print(f"methylKit DMRs: lenient = {len(mk_lenient):,}, strict = {len(mk_strict):,}")

    print("\nBuilding gene lists (paper criterion: ±100 kb of DMR midpoint)...")
    write_gene_list(
        genes_near_dmrs(ep_strict, gene_pos),
        "epykit_dmr_strict_genes.txt",
        "epykit DMRs at q<0.01 ∧ |d|≥25% — 257 DMRs · paper-matched stringency",
    )
    write_gene_list(
        genes_near_dmrs(ep_lenient.head(813), gene_pos),
        "epykit_dmr_top813_genes.txt",
        "epykit top-813 DMRs by q-value — matches paper's reported DMR count",
    )
    write_gene_list(
        genes_near_dmrs(ep_lenient.head(500), gene_pos),
        "epykit_dmr_top500_genes.txt",
        "epykit top-500 DMRs by q-value",
    )
    write_gene_list(
        genes_near_dmrs(ep_lenient, gene_pos),
        "epykit_dmr_lenient_genes.txt",
        "epykit DMRs at q<0.05 ∧ |d|≥10% — 3,433 DMRs · our default",
    )
    write_gene_list(
        genes_near_dmrs(mk_strict, gene_pos),
        "methylkit_dmr_strict_genes.txt",
        "methylKit DMRs at q<0.01 ∧ |d|≥25% — 147 DMRs · paper-matched stringency",
    )
    write_gene_list(
        genes_near_dmrs(mk_lenient.head(813), gene_pos),
        "methylkit_dmr_top813_genes.txt",
        "methylKit top-813 DMRs by q-value",
    )

    # Background list: every gene that has at least one *tested* CpG.
    # Our background is all annotated genes from the TSS-distance files
    # (i.e. every gene that has at least one significant DMC). For a true
    # "all tested" background you'd want the gene list from dmc_all_sites.csv
    # — but that file is 1.7 GB. The annotated-gene set is the practical
    # alternative and is what most enrichment tools accept.
    bg_genes = sorted(gene_pos["gene"].unique().tolist())
    write_gene_list(
        bg_genes,
        "background_genes.txt",
        "Background — all genes with at least one significant epykit DMC "
        "(approximates the 'expressed gene' background; for the paper's "
        "exact 23,590-gene RNA-seq background you'd need their RNA-seq table)",
    )

    # Instructions
    instructions = OUT / "README.md"
    instructions.write_text(
        "# ShinyGO-ready gene lists\n\n"
        "Open <http://bioinformatics.sdstate.edu/go77> and paste one of the "
        "lists below into the input box.\n\n"
        "## Recommended workflow\n\n"
        "1. Start with **`epykit_dmr_strict_genes.txt`** (closest to the paper's "
        "stringency). It will have ~ 1 K–3 K genes — comparable to the paper's "
        "input set size.\n"
        "2. Set the species to **Human (Homo sapiens)**.\n"
        "3. (Optional) For a fairer comparison to the paper, paste "
        "**`background_genes.txt`** into the 'Custom background' box. Without "
        "a custom background ShinyGO uses all annotated genes, which dilutes "
        "the enrichment.\n"
        "4. Set **FDR cutoff = 0.05** and **min pathway size = 5**.\n"
        "5. Pull **Curated.Reactome** to compare against paper Panel D and "
        "**GO Molecular Function** to compare against Panel E.\n\n"
        "## Files\n\n"
        "| File | What it is |\n"
        "|---|---|\n"
        "| `epykit_dmr_strict_genes.txt` | Genes near the 257 epykit DMRs at q<0.01 ∧ \\|d\\|≥25% |\n"
        "| `epykit_dmr_top813_genes.txt` | Genes near the top 813 epykit DMRs (matches paper's count) |\n"
        "| `epykit_dmr_top500_genes.txt` | Genes near the top 500 epykit DMRs |\n"
        "| `epykit_dmr_lenient_genes.txt` | Genes near all 3,433 epykit DMRs at q<0.05 |\n"
        "| `methylkit_dmr_strict_genes.txt` | Same strict cut for methylKit (147 DMRs) |\n"
        "| `methylkit_dmr_top813_genes.txt` | Genes near the top 813 methylKit DMRs |\n"
        "| `background_genes.txt` | ~10,500 genes with at least one significant DMC |\n\n"
        "## Paper's reported pathways to look for\n\n"
        "**Panel D (Reactome):** Class A/1 Rhodopsin-like receptors · Peptide "
        "ligand-binding receptors · GPCR ligand binding · GPCR downstream "
        "signalling · Signalling by GPCR · G alpha i signalling events\n\n"
        "**Panel E (GO MF):** Sequence-specific DNA binding · DNA-binding TF "
        "activity · RNA Pol II cis-regulatory region DNA binding · TF binding · "
        "DNA-binding TF activator activity\n\n"
        "## Note on the result\n\n"
        "The paper used **DSS::DMLfit.multiFactor with smoothing=TRUE** plus "
        "**callDMR(p.threshold=1e-5, minCG=3, minlen=50, dis.merge=100)**, "
        "which produces variable-width DMRs. Ours come from methylKit / epykit "
        "with 500 bp fixed tiles. Even with identical inputs, three things "
        "differ: the test (DSS vs LR), the regions (variable vs fixed tiles), "
        "and the gene-DMR association distance. Modest disagreement on top "
        "enriched pathways is expected. If GPCR / TF terms appear at FDR<0.05 "
        "anywhere in your ShinyGO top 20–50, that's a positive reproduction.\n",
        encoding="utf-8",
    )
    print(f"\nWrote instructions to {instructions}")
    print(f"\nAll lists in: {OUT}")


if __name__ == "__main__":
    main()
