"""Compare Reactome enrichment results across all our runs against the paper's Panel D.

Inputs:
  - Reactome JSON for each of the smoothed-DMR gene lists
  - Reactome JSON for the original tile-DMR gene lists
Output:
  - Markdown table comparing top GPCR-related terms with the paper's exact
    panel D values
"""

from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REACTOME = ROOT / "shinygo_lists" / "outputs" / "reactome"

# Paper Panel D: Supp Table 6 values
PAPER_PANEL_D = [
    ("Class A/1 (Rhodopsin-like receptors)", 18, 321, 3.297, 5.71e-3),
    ("Peptide ligand-binding receptors",     13, 188, 4.139, 5.71e-3),
    ("GPCR ligand binding",                  21, 454, 2.659, 1.20e-2),
    ("GPCR downstream signalling",           31, 1199, 2.023, 3.14e-2),
    ("Signaling by GPCR",                    32, 1263, 1.933, 3.74e-2),
    ("G alpha (i) signalling events",        19, 470, 2.461, 3.74e-2),
]

INPUTS = [
    ("methylKit strict (tile, 639 genes)",
     REACTOME / "methylkit_dmr_strict_genes" / "methylkit_dmr_strict_genes"),
    ("epykit strict (tile, 1,079 genes)",
     REACTOME / "epykit_dmr_strict_genes" / "epykit_dmr_strict_genes"),
    ("epykit top-500 (tile, 1,840 genes)",
     REACTOME / "epykit_dmr_top500_genes" / "epykit_dmr_top500_genes"),
    ("epykit smoothed alpha=1e-5 nearest (516 genes)",
     REACTOME / "epykit_smoothed_alpha1e-5_nearest_genes" / "epykit_smoothed_alpha1e-5_nearest_genes"),
    ("epykit smoothed alpha=1e-5 100kb (2,448 genes)",
     REACTOME / "epykit_smoothed_alpha1e-5_100kb_genes" / "epykit_smoothed_alpha1e-5_100kb_genes"),
]


def parse(path):
    data = json.load(path.open(encoding="utf-8"))
    rows = []
    for p in data["pathways"]:
        stats = p["data"]["statistics"]
        s = next((x for x in stats if x.get("resource") == "TOTAL"), stats[0])
        rows.append({
            "name": p["name"],
            "n_found": s["entitiesFound"],
            "n_total": s["entitiesCount"],
            "pvalue": s["entitiesPValue"],
            "fdr": s["entitiesFDR"],
        })
    return rows


def find(rows, name_substr):
    for r in rows:
        if name_substr.lower() in r["name"].lower():
            return r
    return None


def fmt_fdr(f):
    if f is None or f >= 0.9999:
        return "—"
    if f < 1e-4:
        return f"{f:.2e}"
    if f < 0.01:
        return f"{f:.4f}"
    return f"{f:.3f}"


def main():
    by_input = {}
    for label, path in INPUTS:
        if not path.exists():
            print(f"missing: {path}")
            continue
        rows = parse(path)
        rows.sort(key=lambda r: r["fdr"])
        by_input[label] = rows
        # diagnostics
        n_sig = sum(1 for r in rows if r["fdr"] < 0.05)
        print(f"{label}: {len(rows):,} pathways, {n_sig} at FDR<0.05")

    lines = []
    lines.append("# Reactome panel-D reproduction — final comparison\n")
    lines.append(
        "For each paper Panel D pathway, look up the same Reactome term in our "
        "five enrichment runs. **★** = our FDR < 0.05 (reproduces); **+** = "
        "FDR < 0.1 (marginal); blank = not reproduced.\n"
    )

    # Build a wide table
    headers = ["Paper Panel D pathway", "Paper FDR", "Paper n/N (Fold)"]
    for label in by_input:
        headers.append(label)
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("|" + "|".join(["---"] * len(headers)) + "|")

    for term, paper_n, paper_N, paper_fold, paper_fdr in PAPER_PANEL_D:
        row = [
            f"**{term}**",
            f"{paper_fdr:.2e}",
            f"{paper_n}/{paper_N} ({paper_fold:.2f}×)",
        ]
        for label, rows in by_input.items():
            r = find(rows, term)
            if r is None or r["fdr"] is None:
                row.append("—")
            else:
                fdr_str = fmt_fdr(r["fdr"])
                mark = ""
                if r["fdr"] < 0.05:
                    mark = " ★"
                elif r["fdr"] < 0.10:
                    mark = " +"
                row.append(f"{r['n_found']}/{r['n_total']} · FDR {fdr_str}{mark}")
        lines.append("| " + " | ".join(row) + " |")

    lines.append("")

    # Also list ALL pathways at FDR<0.05 for the smoothed runs
    for label in by_input:
        if "smoothed" not in label:
            continue
        rows = by_input[label]
        sig = [r for r in rows if r["fdr"] < 0.05]
        lines.append(f"\n## All FDR < 0.05 pathways for: {label}\n")
        if not sig:
            lines.append("_(none)_\n")
            continue
        lines.append("| Pathway | n / N | FDR |")
        lines.append("|---|---:|---:|")
        for r in sig[:40]:
            lines.append(f"| {r['name']} | {r['n_found']}/{r['n_total']} | {fmt_fdr(r['fdr'])} |")
        lines.append("")

    out = ROOT / "panel_d_reproduction.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
