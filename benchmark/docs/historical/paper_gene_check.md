# Cross-check vs the original AKAP11 paper (Figure 6)

The paper reports **813 total DMRs** (638 hyper + 175 hypo) and highlights specific genes in panels B (heatmap) and E (pathway enrichment). We have **3,433 DMRs (epykit)** and **2,661 DMRs (methylKit)** at the lenient threshold (q < 0.05, |Δ| ≥ 10 %) — more permissive than the paper, as expected given different parameter choices.

Below we check whether the paper's named genes appear in our calls.

## Per-gene look-up

For each paper-named gene, we list (a) significant DMCs annotated to that gene in epykit's call set, and (b) DMRs whose 500 bp tile contains at least one such DMC. Ranks are by q-value within each tool's lenient set.

| Gene | # epykit DMCs | top DMC \|d\| | epykit DMR rank | methylKit DMR rank |
|---|---:|---:|---:|---:|
| **BCL2** | 9 | 32.6% | #1,439 | #904 |
| **CCDC177** | 1 | 44.7% | — | — |
| **CCDC177.1** | 0 | — | — | — |
| **CREM** | 3 | 38.7% | #157 | — |
| **DMRTA2** | 21 | 44.7% | #325 | — |
| **ENPP2** | 5 | 38.2% | — | — |
| **FOXC1** | 13 | 44.6% | #437 | #248 |
| **GNG11** | 2 | 68.0% | — | — |
| **GREB1L** | 21 | 56.8% | #43 | #23 |
| **IRX2** | 18 | 44.4% | #1,136 | #899 |
| **IRX2.1** | 0 | — | — | — |
| **KC6** | 11 | 30.1% | — | — |
| **LOC100131655** | 0 | — | — | — |
| **LOC100506858** | 18 | 34.8% | — | — |
| **NAALADL2** | 3 | 44.7% | — | — |
| **NFIA** | 8 | 39.4% | #1,033 | #590 |
| **NR2E1** | 17 | 57.3% | — | — |
| **OLIG3** | 12 | 58.2% | #398 | — |
| **OSBPL8** | 3 | 46.5% | — | — |
| **OTX1** | 14 | 53.6% | — | — |
| **OTX2** | 24 | 64.8% | #689 | #396 |
| **PAX2** | 4 | 35.5% | — | — |
| **PAX7** | 13 | 37.2% | — | — |
| **PDK3** | 0 | — | — | — |
| **RPLP0P2** | 1 | 14.3% | — | — |
| **SHOX2** | 4 | 55.7% | #877 | — |
| **SIX3** | 0 | — | — | — |
| **TERF2IP** | 0 | — | — | — |
| **TFAP2B** | 20 | 44.9% | — | — |
| **TMEM242** | 3 | 40.1% | — | — |
| **TMEM242.1** | 0 | — | — | — |
| **VAX1** | 9 | 42.1% | #227 | — |

## Summary

- 25 of 32 paper-named genes have **at least one significant DMC** in epykit's call set.
- 11 of 32 are inside an **epykit DMR**.
- 6 of 32 are inside a **methylKit DMR**.
- 1 appear within epykit's **top-100 DMRs** by q-value.
- 1 appear within methylKit's **top-100 DMRs** by q-value.

## Direction agreement on paper genes

The paper's panel A reports 638 hyper- and 175 hypo-DMRs. For each gene we find in our calls, do we agree on the sign?

| Gene | epykit best meth_diff | direction |
|---|---:|---|
| BCL2 | -32.6% | **hypo** |
| CCDC177 | +44.7% | **hyper** |
| CREM | -38.7% | **hypo** |
| DMRTA2 | -44.7% | **hypo** |
| ENPP2 | -38.2% | **hypo** |
| FOXC1 | +44.6% | **hyper** |
| GNG11 | +68.0% | **hyper** |
| GREB1L | +56.8% | **hyper** |
| IRX2 | +44.4% | **hyper** |
| KC6 | -30.1% | **hypo** |
| LOC100506858 | -34.8% | **hypo** |
| NAALADL2 | -44.7% | **hypo** |
| NFIA | -39.4% | **hypo** |
| NR2E1 | +57.3% | **hyper** |
| OLIG3 | +58.2% | **hyper** |
| OSBPL8 | -46.5% | **hypo** |
| OTX1 | +53.6% | **hyper** |
| OTX2 | +64.8% | **hyper** |
| PAX2 | -35.5% | **hypo** |
| PAX7 | -37.2% | **hypo** |
| RPLP0P2 | -14.3% | **hypo** |
| SHOX2 | +55.7% | **hyper** |
| TFAP2B | -44.9% | **hypo** |
| TMEM242 | +40.1% | **hyper** |
| VAX1 | +42.1% | **hyper** |
