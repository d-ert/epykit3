# Comparison against the paper's supplementary tables

Loaded **Supp Table 5** (n = 813 DMRs, 705 unique genes) and **Supp Table 8** (n = 59 DMR-DEG rows, 46 unique genes).

## DMR size distribution

| Source | n DMRs | median len | mean len | min | max |
|---|---:|---:|---:|---:|---:|
| Paper (DSS callDMR) | 813 | 240 | 293 | 52 | 2,352 |
| epykit lenient (500 bp tiles) | 3,433 | 500 | 500 | 500 | 500 |
| epykit strict | 257 | 500 | 500 | 500 | 500 |
| methylKit lenient | 2,661 | 500 | 500 | 500 | 500 |
| methylKit strict | 147 | 500 | 500 | 500 | 500 |

## DMR coordinate overlap

For each paper DMR (Supp Table 5), is there any overlap with our called DMRs (interval intersection)?

| Our call set | n DMRs | recall of paper (paper ⊆ ours) | precision (ours overlap paper) |
|---|---:|---:|---:|
| epykit lenient | 3,433 | 63 / 813 = **7.7%** | 63 / 3,433 = 1.8% |
| epykit strict | 257 | 39 / 813 = **4.8%** | 39 / 257 = 15.2% |
| methylKit lenient | 2,661 | 72 / 813 = **8.9%** | 74 / 2,661 = 2.8% |
| methylKit strict | 147 | 20 / 813 = **2.5%** | 20 / 147 = 13.6% |

## DMR-associated gene-list overlap

Paper Supp Table 5 has **705** unique DMR-associated genes (HOMER nearest gene). How many of these appear in our DMR-associated gene sets (built using 100 kb-TSS-to-DMR-midpoint as in the paper)?

| Our gene set | size | & paper 705 | recall of paper |
|---|---:|---:|---:|
| epykit strict (1,079 genes) | 1,079 | 133 | **18.9%** |
| epykit top-813 (2,493 genes) | 2,493 | 305 | **43.3%** |
| epykit top-500 (1,840 genes) | 1,840 | 226 | **32.1%** |
| methylKit strict (639 genes) | 639 | 83 | **11.8%** |
| methylKit top-813 (2,319 genes) | 2,319 | 275 | **39.0%** |

## Panel E gene-list capture (the 46 critical genes)

Panel E (TF binding enrichment) was computed on **46 DMR-near-DEG genes** (Supp Table 8). Do we have these genes in our calls?

| Our gene set | & Table 8 (46 genes) | recall |
|---|---:|---:|
| epykit strict (1,079 genes) | 14 / 46 | **30%** |
| epykit top-813 (2,493 genes) | 22 / 46 | **48%** |
| epykit top-500 (1,840 genes) | 19 / 46 | **41%** |
| methylKit strict (639 genes) | 9 / 46 | **20%** |
| methylKit top-813 (2,319 genes) | 23 / 46 | **50%** |

### Captured by epykit top-813 (22 / 46)

BCL2, CREM, DCT, DMRTA2, EGFLAM, FOXC1, GABRR1, IRX2, KANK1, OLIG3, OTX1, OTX2, PAX2, RBFOX1, RELL1, RGS3, SLN, STOX2, TFAP2B, TOM1L1, VAX1, WDR1

### Missed by epykit top-813 (24)

ANXA1, CLEC19A, CNR1, COL1A2, COL6A4P2, CTTN, CXCR4, EDNRB, ENPP2, NECTIN3, NFIA, PCDH18, PGAP1, RABGEF1, RBM47, RPE, SCGB2B2, SHOX2, SIX3, SPRY2, ST6GALNAC5, TERF2IP, TNS3, TSPYL5


## Paper heatmap genes — direct DMR coordinate check

For each gene named in the paper's Figure 6B heatmap, look up its DMR coordinate in Supp Table 5 and check whether we called a DMR at that exact region.

| Gene | Paper DMR coord | length | epykit DMR (any q) overlap? | methylKit DMR (any q) overlap? |
|---|---|---:|:---:|:---:|
| **NR2E1** | chr6:108,174,360–108,176,711 | 2,352 bp |   |   |
| **OTX1** | chr2:63,058,510–63,059,731 | 1,222 bp |   | ✓ |
| **IRX2** | chr5:2,746,408–2,748,285 | 1,878 bp |   |   |
| **OTX2** | chr14:56,803,751–56,805,135 | 1,385 bp |   |   |
| **ENPP2** | chr8:119,671,771–119,672,992 | 1,222 bp |   |   |
| **GREB1L** | chr18:21,243,968–21,244,841 | 874 bp |   |   |
| **CCDC177** | chr14:69,572,158–69,572,337 | 180 bp |   |   |
| **PAX7** | chr1:18,638,903–18,639,471 | 569 bp |   |   |
| **NAALADL2** | chr3:174,439,464–174,439,838 | 375 bp |   |   |
| **PDK3** | chrX:24,494,783–24,494,855 | 73 bp |   |   |
| **TMEM242** | chr6:157,097,394–157,097,482 | 89 bp |   |   |
| **OSBPL8** | chr12:76,421,512–76,422,630 | 1,119 bp |   |   |
| **GNG11** | chr7:93,923,359–93,923,518 | 160 bp |   |   |
| **KC6** | chr18:40,939,564–40,939,816 | 253 bp |   |   |
| **RPLP0P2** | chr11:61,603,371–61,603,520 | 150 bp |   |   |
| **LOC100506858** | chr5:2,112,078–2,112,146 | 69 bp |   |   |
| **LOC100131655** | chr18:76,771,220–76,771,894 | 675 bp |   |   |

## Headline numbers

- **epykit lenient recalls 8% of paper DMRs** (interval overlap)
- **methylKit lenient recalls 9% of paper DMRs**
- **epykit top-813 captures 43% of the paper's 705 DMR-associated genes**
- **epykit top-813 captures 48% of the 46 Panel E genes**
