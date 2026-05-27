"""Parse Reactome AnalysisService JSON output and compare to paper Panel D.

For each of the three input gene lists, sort pathways by FDR and surface
the top hits + flag any term related to the paper's GPCR signature.
"""

from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REACTOME = ROOT / "shinygo_lists" / "outputs" / "reactome"

INPUTS = [
    ("methylKit strict (147 DMRs · 639 genes)",
     REACTOME / "methylkit_dmr_strict_genes" / "methylkit_dmr_strict_genes"),
    ("epykit strict (257 DMRs · 1,079 genes)",
     REACTOME / "epykit_dmr_strict_genes" / "epykit_dmr_strict_genes"),
    ("epykit top-500 (1,840 genes)",
     REACTOME / "epykit_dmr_top500_genes" / "epykit_dmr_top500_genes"),
]

# Paper Panel D headline terms
PANEL_D_KEYWORDS = [
    "rhodopsin", "peptide ligand", "gpcr ligand", "gpcr downstream",
    "signaling by gpcr", "signalling by gpcr", "g alpha i",
    "g alpha (i)", "g-protein", "g protein-coupled",
    "class a", "class a/1", "amine ligand-binding", "amine receptor",
]


def parse_file(path: Path):
    data = json.load(path.open(encoding="utf-8"))
    rows = []
    for p in data["pathways"]:
        stats = p["data"]["statistics"]
        # use TOTAL resource (the default)
        s = next((x for x in stats if x.get("resource") == "TOTAL"), stats[0])
        rows.append({
            "name": p["name"],
            "stId": p["stId"],
            "n_found": s["entitiesFound"],
            "n_total": s["entitiesCount"],
            "pvalue": s["entitiesPValue"],
            "fdr": s["entitiesFDR"],
            "ratio": s["entitiesRatio"],
        })
    return sorted(rows, key=lambda r: r["fdr"])


def matches_panel_d(name: str) -> bool:
    n = name.lower()
    return any(kw in n for kw in PANEL_D_KEYWORDS)


def fmt_fdr(f):
    if f < 1e-4:
        return f"{f:.2e}"
    return f"{f:.4f}"


def main():
    lines = []
    lines.append("# Reactome enrichment — comparison to paper Panel D\n")
    lines.append(
        "Reactome AnalysisService output parsed for each of the three gene "
        "lists. Top 25 pathways by FDR; ★ marks terms matching the paper's "
        "GPCR signature.\n"
    )
    lines.append("Paper Panel D: Class A/1 Rhodopsin-like receptors · Peptide ligand-binding receptors · GPCR ligand binding · GPCR downstream signalling · Signalling by GPCR · G alpha (i) signalling events\n")

    summary = []

    for label, path in INPUTS:
        if not path.exists():
            print(f"missing: {path}")
            continue
        rows = parse_file(path)
        n_sig = sum(1 for r in rows if r["fdr"] < 0.05)
        n_panel_d_sig = sum(1 for r in rows if r["fdr"] < 0.05 and matches_panel_d(r["name"]))
        n_panel_d_top50 = sum(1 for r in rows[:50] if matches_panel_d(r["name"]))

        lines.append(f"\n## {label}\n")
        lines.append(
            f"- Total pathways tested: **{len(rows):,}**\n"
            f"- Pathways at FDR < 0.05: **{n_sig:,}**\n"
            f"- GPCR-related at FDR < 0.05: **{n_panel_d_sig}**\n"
            f"- GPCR-related in top 50: **{n_panel_d_top50}**\n"
        )

        lines.append("### Top 25 pathways by FDR\n")
        lines.append("| Rank | Pathway | n / N | FDR | match? |")
        lines.append("|---:|---|---:|---:|:---:|")
        for i, r in enumerate(rows[:25], 1):
            mark = " ★ " if matches_panel_d(r["name"]) else " "
            lines.append(
                f"| {i} | {r['name']} | {r['n_found']}/{r['n_total']} | "
                f"{fmt_fdr(r['fdr'])} | {mark} |"
            )
        lines.append("")

        # Also dump every GPCR-related hit (regardless of rank)
        gpcr_hits = [r for r in rows if matches_panel_d(r["name"])][:20]
        if gpcr_hits:
            lines.append("### All GPCR-related pathways (sorted by FDR)\n")
            lines.append("| Pathway | n / N | FDR | rank in full list |")
            lines.append("|---|---:|---:|---:|")
            for r in gpcr_hits:
                full_rank = next(
                    (i + 1 for i, x in enumerate(rows) if x["stId"] == r["stId"]), -1
                )
                lines.append(
                    f"| {r['name']} | {r['n_found']}/{r['n_total']} | "
                    f"{fmt_fdr(r['fdr'])} | #{full_rank:,} |"
                )
            lines.append("")

        summary.append((label, len(rows), n_sig, n_panel_d_sig, n_panel_d_top50))
        print(f"\n{label}")
        print(f"  total pathways: {len(rows):,}")
        print(f"  FDR<0.05: {n_sig:,}")
        print(f"  GPCR-related at FDR<0.05: {n_panel_d_sig}")
        print(f"  GPCR-related in top 50: {n_panel_d_top50}")

    # Summary table
    lines.append("\n## Summary across the three gene sets\n")
    lines.append("| Gene set | Total pathways | FDR<0.05 | GPCR @ FDR<0.05 | GPCR in top 50 |")
    lines.append("|---|---:|---:|---:|---:|")
    for s in summary:
        lines.append(f"| {s[0]} | {s[1]:,} | {s[2]:,} | {s[3]} | {s[4]} |")

    out = ROOT / "reactome_vs_paper.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
