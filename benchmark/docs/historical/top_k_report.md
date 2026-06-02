# Top-K intersection — epykit vs methylKit on GSE263850

Computed empirically from the real call files. Coordinates aligned (epykit 0-based + 1 → methylKit 1-based).

## DMC intersection by top-K

- epykit significant DMCs: **30,965**
- methylKit significant DMCs: **51,792**


| K  | epykit ∩ methylKit | recall (∩/K) |
|---:|---:|---:|
| 5 | 0 | 0.0% |
| 10 | 2 | 20.0% |
| 25 | 7 | 28.0% |
| 50 | 19 | 38.0% |
| 100 | 38 | 38.0% |
| 250 | 95 | 38.0% |
| 500 | 203 | 40.6% |

## Top-10 DMCs by q-value — methylKit

| # | chrom:pos | gene (Δbp from TSS) | meth_diff | q-value | also in epykit top-10? |
|---:|---|---|---:|---:|:---:|
| 1 | chr5:158,978,234 | KCNQ5-AS1 (-37427 bp) | +93.5% | 3.42e-31 | ✓ |
| 2 | chr6:108,175,374 | GLRA3 (+20301 bp) | +100.0% | 7.79e-28 | ✓ |
| 3 | chr14:48,561,422 | MIR4637 (+19331 bp) | +82.0% | 3.57e-27 |   |
| 4 | chr6:79,335,552 | LOC101927908 (-104887 bp) | -77.3% | 1.33e-26 |   |
| 5 | chr6:55,874,728 | C10orf90 (-1155 bp) | +81.5% | 4.99e-26 |   |
| 6 | chr5:158,978,263 | MIR3914-2 (-2201 bp) | +87.8% | 1.10e-25 |   |
| 7 | chr5:158,978,245 | KRTAP5-7 (-22679 bp) | +85.9% | 2.99e-25 |   |
| 8 | chr1:82,426,116 | AA06 (-28120 bp) | +80.6% | 8.34e-25 |   |
| 9 | chr4:170,990,366 | CTNNA2 (+163292 bp) | +73.0% | 1.17e-24 |   |
| 10 | chr9:109,426,798 | CAMKMT (+15085 bp) | +81.4% | 2.89e-24 |   |

## Top-10 DMCs by q-value — epykit

| # | chrom:pos | gene (Δbp from TSS) | meth_diff | q-value | also in methylKit top-10? |
|---:|---|---|---:|---:|:---:|
| 1 | chr6:98,601,453 | EDIL3 (+270752 bp) | -78.9% | 3.96e-22 |   |
| 2 | chr2:184,945,323 | LINC02740 (+15765 bp) | -72.5% | 1.12e-21 |   |
| 3 | chr4:155,674,883 | NUP42 (+23896 bp) | +74.2% | 1.78e-21 |   |
| 4 | chr3:172,112,504 | LIMCH1 (+17164 bp) | +85.4% | 5.01e-21 |   |
| 5 | chr10:51,607,103 | MIR5584 (-55521 bp) | -76.2% | 1.85e-20 |   |
| 6 | chr6:88,006,571 | MIR151A (-41925 bp) | +81.1% | 1.10e-19 |   |
| 7 | chr8:18,016,860 | LOC93463 (+16321 bp) | -73.5% | 2.44e-19 |   |
| 8 | chr5:158,978,234 | KCNQ5-AS1 (-37427 bp) | +90.7% | 5.30e-19 | ✓ |
| 9 | chr6:55,874,490 | NLGN1-AS1 (-394183 bp) | +86.1% | 5.60e-19 |   |
| 10 | chr6:108,175,374 | GLRA3 (+20301 bp) | +100.0% | 7.95e-19 | ✓ |

## DMR intersection by top-K (lenient threshold)

- epykit significant DMRs: **3,433**
- methylKit significant DMRs: **2,661**


| K  | epykit ∩ methylKit | recall (∩/K) |
|---:|---:|---:|
| 5 | 3 | 60.0% |
| 10 | 6 | 60.0% |
| 25 | 12 | 48.0% |
| 50 | 27 | 54.0% |
| 100 | 53 | 53.0% |
| 250 | 143 | 57.2% |
| 500 | 291 | 58.2% |

## Top-10 DMRs by q-value — methylKit

| # | chrom:start-end | gene(s) in tile | meth_diff | q-value | also in epykit top-10? |
|---:|---|---|---:|---:|:---:|
| 1 | chr6:27,977,501–27,978,000 | FRMD4A (-124579 bp), SNORA20B (-20961 bp) | +56.0% | 2.59e-89 | ✓ |
| 2 | chr3:169,660,001–169,660,500 | SELENOF (+91472 bp), CPZ (+57507 bp) | +24.9% | 4.36e-70 | ✓ |
| 3 | chr10:122,530,001–122,530,500 | FGF1 (+118894 bp), C6orf118 (-85040 bp) | +46.9% | 2.05e-67 | ✓ |
| 4 | chr6:28,001,501–28,002,000 | MBP (-26513 bp), LINC02005 (-79028 bp) | +37.2% | 1.35e-64 | ✓ |
| 5 | chr4:19,193,001–19,193,500 | PRPF18 (+46273 bp), KLF13 (+465437 bp) | +29.6% | 4.33e-64 | ✓ |
| 6 | chr6:75,156,001–75,156,500 | PDE4B (+3917 bp), PPP2R3A (-306073 bp) | +45.9% | 1.42e-63 | ✓ |
| 7 | chr6:155,121,501–155,122,000 | MAPK10 (+63286 bp), OR13A1 (+75695 bp) | +34.8% | 3.69e-56 |   |
| 8 | chr19:50,781,001–50,781,500 | CTAGE11P (+40966 bp), C12orf56 (+61473 bp) | +41.2% | 6.17e-52 |   |
| 9 | chr2:9,961,001–9,961,500 | LOC158435 (+18181 bp), TMEM47 (+6629 bp) | +30.2% | 1.26e-46 |   |
| 10 | chr9:973,001–973,500 | (no DMC-level annotation in tile) | +23.9% | 6.30e-45 |   |

## Top-10 DMRs by q-value — epykit

| # | chrom:start-end | gene(s) in tile | meth_diff | q-value | also in methylKit top-10? |
|---:|---|---|---:|---:|:---:|
| 1 | chr6:27,977,501–27,978,000 | FRMD4A (-124579 bp), SNORA20B (-20961 bp) | +55.3% | 3.10e-60 | ✓ |
| 2 | chr4:9,106,001–9,106,500 | FAM155A (-139878 bp), GCOM1 (-59365 bp) | -37.4% | 3.14e-51 |   |
| 3 | chr1:217,925,001–217,925,500 | LINC00355 (-280990 bp), LINC01747 (-116165 bp) | +19.4% | 9.40e-50 |   |
| 4 | chr10:122,530,001–122,530,500 | FGF1 (+118894 bp), C6orf118 (-85040 bp) | +46.9% | 3.42e-47 | ✓ |
| 5 | chr3:169,660,001–169,660,500 | SELENOF (+91472 bp), CPZ (+57507 bp) | +24.6% | 3.42e-47 | ✓ |
| 6 | chr4:19,193,001–19,193,500 | PRPF18 (+46273 bp), KLF13 (+465437 bp) | +30.0% | 4.51e-46 | ✓ |
| 7 | chr6:28,001,501–28,002,000 | MBP (-26513 bp), LINC02005 (-79028 bp) | +37.1% | 3.47e-45 | ✓ |
| 8 | chr6:75,156,001–75,156,500 | PDE4B (+3917 bp), PPP2R3A (-306073 bp) | +46.1% | 3.47e-45 | ✓ |
| 9 | chr18:463,501–464,000 | P2RX4 (+1117 bp), TENM2 (+2149 bp) | +38.3% | 4.85e-41 |   |
| 10 | chr3:157,631,501–157,632,000 | OTOL1 (+63786 bp), LOC100128079 (+18014 bp) | +36.8% | 1.35e-39 |   |

## Summary

- **Top 5 DMC overlap:** 0 / 5 = 0%
- **Top 10 DMC overlap:** 2 / 10 = 20%
- **Top 5 DMR overlap:** 3 / 5 = 60%
- **Top 10 DMR overlap:** 6 / 10 = 60%
