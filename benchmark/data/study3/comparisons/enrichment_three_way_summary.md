# Three-way pathway enrichment summary

Gene lists from paper Supp Table 5 (705 paper-named genes), our DMR callers (100 kb linkage for chain_merge/DSS, nearest-TSS for methylKit-tile). Enrichr REST API; Reactome_2022 + KEGG_2021_Human + GO_MF_2023; top 20 by p-value per library.


## n_paper_term_matches in top-20 (per library, per caller)

| Caller | Reactome | KEGG | GO MF | n_genes |
|---|---:|---:|---:|---:|
| paper_Table5 | 2 | 2 | 2 | 705 |
| methylKit_tile | 1 | 1 | 1 | 2111 |
| ek_chain_merge_100 | 2 | 1 | 0 | 1290 |
| ek_chain_merge_250 | 1 | 1 | 0 | 1645 |
| DSS_from_scratch | 2 | 1 | 0 | 1467 |

(See enrichment_three_way.json for the full top-20 + overlap genes per cell.)