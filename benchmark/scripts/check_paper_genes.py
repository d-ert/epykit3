"""Cross-check the original paper's named genes against our DMC/DMR call lists.

Source: AKAP11 paper Figure 6 panels B and E. Genes that the original
authors highlighted as the biologically interesting DMR-associated hits.
"""

from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT.parent / "epykit_vs_methylkit(GSE263850)"
MK_SOURCE = ROOT.parent.parent / "methylkıt_realResults" / "scripts_and_results" / "methylkit_results"

# Genes the original paper named in the figure
HEATMAP_GENES = [
    "KC6", "LOC100506858", "RPLP0P2", "OSBPL8", "TMEM242", "TMEM242.1",
    "PDK3", "LOC100131655", "NAALADL2", "PAX7", "IRX2", "IRX2.1",
    "CCDC177", "CCDC177.1", "GREB1L", "ENPP2", "GNG11", "OTX2",
    "NR2E1", "OTX1",
]
PATHWAY_GENES = [
    "TFAP2B", "FOXC1", "PAX2", "CREM", "OTX1", "SIX3", "DMRTA2", "VAX1",
    "NFIA", "OTX2", "TERF2IP", "SHOX2", "IRX2", "BCL2", "OLIG3",
    "ENPP2",
]
ALL_PAPER_GENES = sorted(set(HEATMAP_GENES + PATHWAY_GENES))


def load_ep_dmc_annotated():
    """Return epykit significant DMCs joined with their gene annotation."""
    sig = pd.read_csv(SOURCE / "epykit_results" / "dmc_significant_qval05.csv")
    sig["pos_1based"] = sig["pos"].astype(int) + 1

    hyper_subset = sig[sig["meth_diff"] > 0].reset_index(drop=True)
    hypo_subset = sig[sig["meth_diff"] < 0].reset_index(drop=True)

    rows = []
    for fn, subset, direction in [
        ("tss_distance_hyper.csv", hyper_subset, "hyper"),
        ("tss_distance_hypo.csv", hypo_subset, "hypo"),
    ]:
        p = SOURCE / "epykit_results" / fn
        t = pd.read_csv(p)
        for _, r in t.iterrows():
            ridx = int(r["target.row"]) - 1
            if 0 <= ridx < len(subset):
                site = subset.iloc[ridx]
                rows.append({
                    "chrom": site["chrom"],
                    "pos_1based": int(site["pos"]) + 1,
                    "meth_diff": float(site["meth_diff"]),
                    "qvalue": float(site["qvalue"]),
                    "gene": str(r["feature.name"]),
                    "dist_to_tss": int(r["dist.to.feature"]),
                    "direction": direction,
                })
    return pd.DataFrame(rows)


def dmrs_overlapping_gene(dmr_df, gene_dmcs, source_label):
    """For each DMR, count how many DMCs near a given gene fall inside."""
    matches = []
    for _, dmr in dmr_df.iterrows():
        chrom = dmr["chrom"]
        start = int(dmr["start"])
        end = int(dmr["end"])
        for _, dmc in gene_dmcs.iterrows():
            if dmc["chrom"] == chrom and start + 1 <= dmc["pos_1based"] <= end:
                matches.append({
                    "dmr_chrom": chrom,
                    "dmr_start": start,
                    "dmr_end": end,
                    "dmr_qvalue": dmr["qvalue"],
                    "dmr_meth_diff": dmr["meth_diff"],
                    "dmr_rank": dmr["rank"],
                    "gene": dmc["gene"],
                    "dist_to_tss": dmc["dist_to_tss"],
                    "dmc_pos": dmc["pos_1based"],
                    "dmc_qvalue": dmc["qvalue"],
                    "source": source_label,
                })
    return pd.DataFrame(matches)


def main():
    print("Loading epykit DMC annotation...")
    ep_ann = load_ep_dmc_annotated()
    print(f"  annotated DMCs: {len(ep_ann):,}")
    print(f"  unique gene names: {ep_ann['gene'].nunique():,}")

    print("\nLoading DMR lists...")
    ep_dmr = pd.read_csv(ROOT / "data" / "study3" / "dmr_significant_lenient.csv")
    mk_dmr = pd.read_csv(MK_SOURCE / "dmr_significant_lenient.csv")
    ep_dmr = ep_dmr.sort_values("qvalue").reset_index(drop=True)
    ep_dmr["rank"] = ep_dmr.index + 1
    mk_dmr = mk_dmr.sort_values("qvalue").reset_index(drop=True)
    mk_dmr["rank"] = mk_dmr.index + 1
    print(f"  epykit DMRs (lenient): {len(ep_dmr):,}")
    print(f"  methylKit DMRs (lenient): {len(mk_dmr):,}")

    lines = []
    lines.append("# Cross-check vs the original AKAP11 paper (Figure 6)\n")
    lines.append(
        "The paper reports **813 total DMRs** (638 hyper + 175 hypo) and highlights "
        "specific genes in panels B (heatmap) and E (pathway enrichment). "
        "We have **3,433 DMRs (epykit)** and **2,661 DMRs (methylKit)** at the "
        "lenient threshold (q < 0.05, |Δ| ≥ 10 %) — more permissive than the paper, "
        "as expected given different parameter choices.\n"
    )
    lines.append("Below we check whether the paper's named genes appear in our calls.\n")

    # ---- For each named gene, find associated DMCs and DMRs ---------------
    lines.append("## Per-gene look-up\n")
    lines.append(
        "For each paper-named gene, we list (a) significant DMCs annotated to that "
        "gene in epykit's call set, and (b) DMRs whose 500 bp tile contains at "
        "least one such DMC. Ranks are by q-value within each tool's lenient set.\n"
    )
    lines.append(
        "| Gene | # epykit DMCs | top DMC \\|d\\| | epykit DMR rank | methylKit DMR rank |"
    )
    lines.append("|---|---:|---:|---:|---:|")

    summary = []
    for gene in ALL_PAPER_GENES:
        gene_dmcs = ep_ann[ep_ann["gene"] == gene].copy()
        if len(gene_dmcs) == 0:
            lines.append(f"| **{gene}** | 0 | — | — | — |")
            summary.append((gene, 0, None, None, None))
            continue
        gene_dmcs = gene_dmcs.sort_values("qvalue")
        top_d = abs(gene_dmcs.iloc[0]["meth_diff"])

        # find DMRs containing these DMCs
        ep_match = dmrs_overlapping_gene(ep_dmr, gene_dmcs, "epykit")
        mk_match = dmrs_overlapping_gene(mk_dmr, gene_dmcs, "methylkit")
        ep_best = ep_match["dmr_rank"].min() if len(ep_match) else None
        mk_best = mk_match["dmr_rank"].min() if len(mk_match) else None

        ep_str = f"#{int(ep_best):,}" if ep_best else "—"
        mk_str = f"#{int(mk_best):,}" if mk_best else "—"
        lines.append(
            f"| **{gene}** | {len(gene_dmcs)} | {top_d:.1f}% | {ep_str} | {mk_str} |"
        )
        summary.append((gene, len(gene_dmcs), top_d, ep_best, mk_best))

    # ---- Summary ---------------------------------------------------------
    n_with_dmc = sum(1 for _, n, _, _, _ in summary if n > 0)
    n_in_ep_dmr = sum(1 for _, _, _, e, _ in summary if e is not None)
    n_in_mk_dmr = sum(1 for _, _, _, _, m in summary if m is not None)
    n_in_top100_ep = sum(1 for _, _, _, e, _ in summary if e is not None and e <= 100)
    n_in_top100_mk = sum(1 for _, _, _, _, m in summary if m is not None and m <= 100)

    lines.append("\n## Summary\n")
    lines.append(
        f"- {n_with_dmc} of {len(ALL_PAPER_GENES)} paper-named genes have **at least one significant DMC** in epykit's call set.\n"
        f"- {n_in_ep_dmr} of {len(ALL_PAPER_GENES)} are inside an **epykit DMR**.\n"
        f"- {n_in_mk_dmr} of {len(ALL_PAPER_GENES)} are inside a **methylKit DMR**.\n"
        f"- {n_in_top100_ep} appear within epykit's **top-100 DMRs** by q-value.\n"
        f"- {n_in_top100_mk} appear within methylKit's **top-100 DMRs** by q-value.\n"
    )

    # ---- Direction agreement ---------------------------------------------
    lines.append("## Direction agreement on paper genes\n")
    lines.append(
        "The paper's panel A reports 638 hyper- and 175 hypo-DMRs. For each "
        "gene we find in our calls, do we agree on the sign?\n"
    )
    direction_agree = 0
    direction_disagree = 0
    direction_lines = ["| Gene | epykit best meth_diff | direction |", "|---|---:|---|"]
    for gene in ALL_PAPER_GENES:
        gene_dmcs = ep_ann[ep_ann["gene"] == gene].copy()
        if len(gene_dmcs) == 0:
            continue
        gene_dmcs = gene_dmcs.sort_values("qvalue")
        best = gene_dmcs.iloc[0]
        d = best["meth_diff"]
        direction_lines.append(
            f"| {gene} | {d:+.1f}% | "
            f"{'**hyper**' if d > 0 else '**hypo**'} |"
        )

    lines.extend(direction_lines)
    lines.append("")

    out_path = ROOT / "paper_gene_check.md"
    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nWrote {out_path}")
    print(f"  {n_with_dmc} / {len(ALL_PAPER_GENES)} paper-named genes have epykit DMCs")
    print(f"  {n_in_ep_dmr} / {len(ALL_PAPER_GENES)} are inside an epykit DMR")
    print(f"  {n_in_mk_dmr} / {len(ALL_PAPER_GENES)} are inside a methylKit DMR")
    print(f"  {n_in_top100_ep} in epykit top-100 DMRs")
    print(f"  {n_in_top100_mk} in methylKit top-100 DMRs")


if __name__ == "__main__":
    main()
