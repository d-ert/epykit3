"""Export the exact data each three-way figure (F1-F10) plots, as TSVs,
into ONE folder for manual inspection.

Output: benchmark/figures/study3_real_GSE263850/three_way/source_data_tsv/
  README.md                         <- figure -> TSV -> original source map
  F1_upset_intersections.tsv
  F1_set_sizes.tsv
  F2_dis_merge_sweep.tsv
  F3_dmr_lengths_long.tsv
  F4_top_named_gene_hits.tsv
  F5_annotation_distribution.tsv
  F6_methylation_heatmap.tsv
  F7_resources.tsv
  F8_enrichment_top20.tsv
  F9_per_dmr_concordance.tsv

Each TSV is the literal data behind the figure's bars / points / cells,
so you can open them in Excel and see exactly what every figure shows.
"""

from __future__ import annotations
import sys
import io
import json
from pathlib import Path
import pandas as pd
import polars as pl

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

BENCH = Path(__file__).resolve().parents[1]
STUDY3 = BENCH / "data" / "study3"
SWEEP = BENCH / "data" / "multi_thread_and_chain_sweep" / "chain_merge_dis_merge_sweep"
COMP = STUDY3 / "comparisons"
THREE_WAY = BENCH / "figures" / "study3_real_GSE263850" / "three_way"
OUT = THREE_WAY / "source_data_tsv"

# External reference inputs (not part of the committed benchmark tree)
PAPER_T5 = Path(r"D:/Coding/Projeler/methyl_lib/epykit2/GSE263850_RAW/Paper resources/DMR_total_list.xlsx")
MK_TILE = Path(r"D:/Coding/Projeler/methyl_lib/methylkıt_realResults/scripts_and_results/methylkit_results/dmr_significant_lenient.csv")
MK_STEP = Path(r"D:/Coding/Projeler/methyl_lib/methylkıt_realResults/scripts_and_results/methylkit_results/benchmark/step_benchmarks.csv")
EK_STEP = Path(r"D:/Coding/Projeler/methyl_lib/benchmarkin_merges/epykit_vs_methylkit(GSE263850)/epykit_results/benchmark/step_benchmarks.csv")


def w(df, name):
    """Write a pandas/polars frame as TSV."""
    if isinstance(df, pl.DataFrame):
        df = df.to_pandas()
    df.to_csv(OUT / name, sep="\t", index=False)
    print(f"  {name}: {len(df)} rows")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Writing TSVs to {OUT}")

    # ---- F1: 3-way DMR overlap (UpSet) ---------------------------------
    f1 = pd.read_csv(THREE_WAY / "F1_upset_data.csv")
    w(f1, "F1_upset_intersections.tsv")
    paper = pd.read_excel(PAPER_T5).rename(columns={"chr": "chrom"})
    ek100 = pl.read_parquet(STUDY3 / "chain_merge" / "dmr_chain_merge.parquet").to_pandas()
    dss = pd.read_csv(STUDY3 / "dss" / "dmr_dss.csv")
    w(pd.DataFrame({
        "set": ["paper (Supp Table 5)", "epykit-chain_merge-100", "DSS-from-scratch"],
        "n_dmrs": [len(paper), len(ek100), len(dss)],
    }), "F1_set_sizes.tsv")

    # ---- F2: dis.merge sweep -------------------------------------------
    w(pd.read_csv(SWEEP / "sweep_summary.csv"), "F2_dis_merge_sweep.tsv")

    # ---- F3: DMR length distributions (long form) ----------------------
    def lengths(df, label):
        d = df.copy()
        d["length_bp"] = d["end"].astype(int) - d["start"].astype(int)
        d["caller"] = label
        return d[["caller", "length_bp"]]
    ek250 = pl.read_parquet(SWEEP / "dis_merge_250" / "dmr.parquet").to_pandas()
    f3 = pd.concat([
        lengths(paper, "paper-DSS (Supp Table 5)"),
        lengths(pd.read_csv(MK_TILE), "methylKit-tile"),
        lengths(pd.read_csv(STUDY3 / "dmr_significant_lenient.csv"), "epykit-tile"),
        lengths(ek100, "epykit-chain_merge-100"),
        lengths(ek250, "epykit-chain_merge-250"),
        lengths(dss, "DSS-from-scratch"),
    ], ignore_index=True)
    w(f3, "F3_dmr_lengths_long.tsv")
    # plus a per-caller length summary for quick reading
    w(f3.groupby("caller")["length_bp"].describe().reset_index(),
      "F3_dmr_lengths_summary.tsv")

    # ---- F4: top named gene hits ---------------------------------------
    w(pd.read_csv(THREE_WAY / "F4_top_named_gene_hits_data.csv"),
      "F4_top_named_gene_hits.tsv")

    # ---- F5 / F5b: annotation distribution -----------------------------
    w(pd.read_csv(COMP / "annotation_distribution.csv"),
      "F5_annotation_distribution.tsv")

    # ---- F6: per-sample methylation heatmap of top-20 genes ------------
    p6 = paper.copy()
    p6["abs_dB"] = p6["diff.meth_mean"].abs()
    p6["gene"] = p6["Gene.Name"].astype(str)
    hyper = p6[p6["diff.meth_mean"] > 0].sort_values("abs_dB", ascending=False).drop_duplicates("gene").head(10)
    hypo = p6[p6["diff.meth_mean"] < 0].sort_values("abs_dB", ascending=False).drop_duplicates("gene").head(10)
    sel = pd.concat([hyper, hypo], ignore_index=True)
    sample_cols = ["SBP009_1_mean", "SBP009_2_mean", "SBP009_3_mean",
                   "Het-AKAP11-KO-Clone16_mean", "Het-AKAP11-KO-Clone20_mean",
                   "Het-AKAP11-KO-Clone21_mean"]
    sel["direction"] = ["hyper"] * len(hyper) + ["hypo"] * len(hypo)
    sel["coord"] = sel["chrom"].astype(str) + ":" + sel["start"].astype(str) + "-" + sel["end"].astype(str)
    cols6 = ["gene", "direction", "coord", "diff.meth_mean"] + [c for c in sample_cols if c in sel.columns]
    w(sel[cols6], "F6_methylation_heatmap.tsv")

    # ---- F7: resource bars ---------------------------------------------
    w(pd.read_csv(THREE_WAY / "F7_resources_data.csv"), "F7_resources.tsv")

    # ---- F8: enrichment top-20 (flattened) -----------------------------
    data = json.loads((COMP / "enrichment_three_way.json").read_text(encoding="utf-8"))
    rows = []
    for caller, entry in data.items():
        for lib, libdata in entry.get("libraries", {}).items():
            for rank, t in enumerate(libdata.get("top20", []), start=1):
                rows.append({
                    "caller": caller, "library": lib, "rank": rank,
                    "term": t["term"], "p_value": t["p_value"], "p_adj": t["p_adj"],
                    "n_overlap": t["n_overlap"],
                    "paper_term_match": t.get("paper_term_match"),
                })
    w(pd.DataFrame(rows), "F8_enrichment_top20.tsv")

    # ---- F9: per-DMR effect-size concordance ---------------------------
    w(pd.read_csv(COMP / "per_dmr_stat_concordance.csv"),
      "F9_per_dmr_concordance.tsv")

    # ---- README --------------------------------------------------------
    readme = """# Source data for the Study 3 three-way figures

Every TSV here is the **exact data** plotted by the corresponding figure
in `../` (the `figures/study3_real_GSE263850/three_way/*.png`). Open them
in Excel / a text editor to inspect what each figure shows.

Regenerate with: `uv run python benchmark/scripts/export_figure_source_tsv.py`

| Figure (png) | What it shows | TSV here | Original source file(s) |
|---|---|---|---|
| **F1a/F1b** upset | 3-way DMR overlap (paper / ek-cm-100 / DSS), any-bp & J>=0.5 | `F1_upset_intersections.tsv`, `F1_set_sizes.tsv` | paper Supp Table 5 xlsx; `data/study3/chain_merge/dmr_chain_merge.parquet`; `data/study3/dss/dmr_dss.csv` |
| **F2** dis.merge sweep | recall/precision/length/panel-E vs dis.merge (vs DSS-922) | `F2_dis_merge_sweep.tsv` | `data/multi_thread_and_chain_sweep/chain_merge_dis_merge_sweep/sweep_summary.csv` |
| **F3** length distributions | DMR length per caller (violin) | `F3_dmr_lengths_long.tsv` (+ `_summary`) | paper xlsx; methylKit & epykit `dmr_significant_lenient.csv`; chain_merge + dis_merge_250 parquets; dss.csv |
| **F4** top named genes | 20 paper Fig-3B genes, best-overlap Jaccard per caller | `F4_top_named_gene_hits.tsv` | paper xlsx; methylKit tile; chain_merge / dis_merge_250 parquets; dss.csv |
| **F5 / F5b** annotation | HOMER genomic-feature distribution of DMRs (pie + stacked bar) | `F5_annotation_distribution.tsv` | `data/study3/comparisons/annotation_distribution.csv` |
| **F6** methylation heatmap | per-sample beta of the 20 named genes (paper-side) | `F6_methylation_heatmap.tsv` | paper Supp Table 5 xlsx |
| **F7** resources | wall / CPU / peak-RSS per caller (Linux pivoine) | `F7_resources.tsv` | methylKit & epykit `step_benchmarks.csv`; `data/study3/dss/resources.json` + `step_timings_resume.tsv` |
| **F8** enrichment | top-20 Reactome/KEGG/GO-MF terms per caller | `F8_enrichment_top20.tsv` | `data/study3/comparisons/enrichment_three_way.json` |
| **F9** per-DMR concordance | per-matched-DMR ek vs DSS effect size + Jaccard | `F9_per_dmr_concordance.tsv` | `data/study3/comparisons/per_dmr_stat_concordance.csv` |
| **F10** summary composite | reuses F2 + F3 + F4 + F5 + F7 data (6 panels) | (see those TSVs) | all of the above |

## Notes
- "paper xlsx" = `epykit2/GSE263850_RAW/Paper resources/DMR_total_list.xlsx`
  (the source paper's Supp Table 5) — an EXTERNAL reference input, not in
  the committed benchmark tree.
- methylKit-tile DMRs (`dmr_significant_lenient.csv`) and its
  `step_benchmarks.csv` live in the external `methylkit_realResults` repo.
- Everything else is in `benchmark/data/` (the committed source tree).
- F7 numbers are the Linux host (pivoine) measurements that match
  paper Table 5b.
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")
    print(f"  README.md")
    print(f"\nDone. {len(list(OUT.glob('*.tsv')))} TSVs + README in {OUT}")


if __name__ == "__main__":
    main()
