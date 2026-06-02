# epykit chain_merge / smoothed DMRs vs the paper

epykit run with **smoothing enabled** and **alpha = 1e-5** (matching the paper's `DSS::callDMR(p.threshold = 1e-5)` parameter).

- DMRs called: **702** (paper: 813; ratio = 0.86×)
- Median length: **123 bp** (paper: 240 bp — much closer than fixed tiles' 500 bp)
- Hyper / Hypo: 579 / 123 (paper: 638 / 175 — same direction bias)


## DMR coordinate overlap

| Source | DMRs | Recall of paper | Precision (ours ∩ paper) |
|---|---:|---:|---:|
| **epykit smoothed alpha=1e-5** | 702 | **429/813 = 52.8%** | **453/702 = 64.5%** |
| (for reference) epykit lenient 500 bp tiles | 3,433 | 63/813 = 7.7% | 63/3,433 = 1.8% |
| (for reference) methylKit lenient 500 bp tiles | 2,661 | 72/813 = 8.9% | 74/2,661 = 2.8% |


## Gene-list overlap with paper Supp Table 5

| Gene set | size | ∩ paper 705 | recall of paper |
|---|---:|---:|---:|
| epykit smoothed: nearest-gene | 516 | 126 | **17.9%** |
| epykit smoothed: 100 kb assoc | 2,448 | 288 | **40.9%** |


## Panel E gene-list capture (46 critical genes)

| Gene set | ∩ Table 8 (46) | recall |
|---|---:|---:|
| epykit smoothed: nearest-gene | 14 | **30%** |
| epykit smoothed: 100 kb assoc | 25 | **54%** |

### Panel E genes captured (25 / 46)

ANXA1, CNR1, COL6A4P2, CTTN, CXCR4, DCT, DMRTA2, EDNRB, EGFLAM, FOXC1, IRX2, KANK1, NFIA, OTX1, OTX2, RBFOX1, RBM47, RELL1, RGS3, SCGB2B2, SLN, STOX2, TOM1L1, TSPYL5, VAX1

### Panel E genes missed (21)

BCL2, CLEC19A, COL1A2, CREM, ENPP2, GABRR1, NECTIN3, OLIG3, PAX2, PCDH18, PGAP1, RABGEF1, RPE, SHOX2, SIX3, SPRY2, ST6GALNAC5, TERF2IP, TFAP2B, TNS3, WDR1



## Paper heatmap genes — direct coordinate match

| Gene | Paper DMR coord | smoothed-epykit overlap? |
|---|---|:---:|
| NR2E1 | chr6:108,174,360-108,176,711 | ✓ |
| OTX1 | chr2:63,058,510-63,059,731 | ✓ |
| IRX2 | chr5:2,746,408-2,748,285 | ✓ |
| OTX2 | chr14:56,803,751-56,805,135 | ✓ |
| ENPP2 | chr8:119,671,771-119,672,992 | ✓ |
| GREB1L | chr18:21,243,968-21,244,841 | ✓ |
| CCDC177 | chr14:69,572,158-69,572,337 | ✓ |
| PAX7 | chr1:18,638,903-18,639,471 |   |
| NAALADL2 | chr3:174,439,464-174,439,838 |   |
| PDK3 | chrX:24,494,783-24,494,855 |   |
| TMEM242 | chr6:157,097,394-157,097,482 |   |
| OSBPL8 | chr12:76,421,512-76,422,630 |   |
| GNG11 | chr7:93,923,359-93,923,518 |   |


## Summary

With smoothing + alpha=1e-5 (matching the paper's parameters), epykit:

- recovers **53% of paper DMRs** at the coordinate level (vs 7.7% with fixed tiles — a 6.9× improvement)
- captures **41% of paper's 705 genes** (vs 43% with tiles)
- captures **54% of the 46 Panel E genes** (vs 48% with tiles)

The remaining gap is plausibly DSS-vs-LR test differences (DSS smooths dispersion across CpGs in a different way than epykit's combined-pvalue approach). But qualitatively, when the region model matches, the call set converges.
