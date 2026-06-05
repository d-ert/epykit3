# Source data for the Study 3 three-way figures

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
