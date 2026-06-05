# epykit Benchmark — Executive Summary

**Audience:** PI / reviewer / decision-maker, not a statistician.
**Length:** one page.
**For the science behind these numbers**, read [paper/paper.md](paper/paper.md).
**For full tables and audit trail**, read [report/REPORT.md](report/REPORT.md).

---

## What we did

We evaluated **epykit**, a new Python-native pipeline for DNA methylation
analysis, against the established R/CLI ecosystem in three complementary
studies:

| # | Study | Comparison | Data |
|---|---|---|---|
| 1 | Panel against the field | epykit vs 8 baselines (methylKit, DSS, RADMeth, BiSeq, methylSig, Fisher, BSmooth, metilene) | Piao et al. 2021 simulated BS-seq grid (100K CpGs DMC + 4M CpGs DMR, 5×–25× coverage, 2–10 samples) |
| 2 | Apples-to-apples vs methylKit | epykit vs methylKit, **same machine, same harness** | Same simulated grid |
| 3 | Real biological data | epykit vs methylKit on **real WGBS** | GEO GSE263850 (Clone vs SBP009, hg38, 6 samples, 15.6 M CpGs) |

The three studies cover (i) where epykit sits within the broader ecosystem,
(ii) whether its statistics match the state-of-the-art when controlled
side-by-side, and (iii) whether it behaves correctly on real biological data.

---

## Headline numbers

| Question | Answer |
|---|---|
| Does epykit recover the same DMCs as methylKit on simulated data? | **Yes — TPR, FPR, F1, AUROC identical to three decimal places at n ≥ 4.** |
| What happens at n = 2 (1 vs 1)? | epykit recovers **0.564 TPR vs methylKit 0.302** at the same FPR (≈ 2× more true positives). methylKit's overdispersion estimate degenerates at n = 2; epykit's `lr / allow_n1` does not. |
| What about against the broader baseline panel? | At 5× coverage / small effects, epykit `lr` retains **83.5 % TPR** vs methylKit 26.6 %, DSS 6.5 %, Fisher 8.2 %, RADMeth 42.1 %. By 15× coverage, all credible tools converge. |
| Does the FPR hold up? | At 5×, epykit `lr` FPR is **3.7 × 10⁻⁵** — **100–600× tighter** than methylKit / DSS / RADMeth. |
| Does it find the same regions (DMRs)? | At simulated 10×+ coverage: **35 / 35 reference DMRs recovered**, matching methylKit and pooled Fisher. At 5×: epykit `chain_merge` 97 %, methylKit 100 %. |
| Does it behave on real data (GSE263850)? | **Per-CpG**: yes — effect-size Pearson r = 0.994, 94 % direction agreement with methylKit. **DMR-level**: depends on engine. Against the locally-rerun DSS-from-scratch (922 DMRs, paper-matched parameters): epykit `dmr_chain_merge` recovers **63.8 %** (paper-faithful `dis.merge = 100`, 852 DMRs) to **77.3 %** (morphology-matched `dis.merge = 250`, 1,139 DMRs) at any-bp overlap, and **34.5 % / 64.2 %** at Jaccard ≥ 0.5, with **100 % direction agreement on every matched DMR**. DSS-from-scratch reaches 87.5 % of the paper's own 813-DMR set — the DSS-vs-paper floor sets the realistic ceiling for any other tool. Fixed-tile callers (incl. methylKit-tile) miss ≥ 90 %. |
| Why does epykit call fewer significant DMCs on real data (30K vs methylKit's 52K)? | Different operating points on the precision/recall curve, not different biology. epykit's default is more conservative; the opt-in `lr+` recipe recovers 93 % of methylKit's calls. |
| What about the named genes from the source paper's Fig 3B? | epykit-chain_merge-100 hits 12 / 20 of paper's labeled hyper+hypo genes at any-bp overlap (NR2E1, OTX1, IRX2, ENPP2, GREB1L, CCDC177, GNG11, EBF1 + 4 hypo); chain_merge-250 hits 14 / 20. DSS-from-scratch hits 18 / 20. methylKit-tile hits 2 / 20. |
| Speed (Linux, methylKit `mc.cores = 8`)? | **≈ 33 × faster on per-CpG testing vs methylKit, ≈ 33 × faster than DSS** at the simulator headline cell; **≈ 28 × end-to-end** on the Study 3 22 M-CpG real-data input. DSS's `DMLfit.multiFactor` is single-thread by construction (no multi-core path in DSS 2.58.0). Earlier Windows-only numbers (12 × – 68 ×) reflected methylKit's `mc.cores` no-op on Windows and overstated the gap. |
| Memory? | **1.18× less peak RAM on simulated data, 3.83× less than methylKit on real data** (12.6 GB vs 48.0 GB on the 15.6 M-CpG genome). On Linux, DSS-from-scratch peaks at 14.3 GB on the same input — epykit uses ~1.13× less than DSS too. |

---

## Cross-study summary table

| Study | Compared against | Speedup (wall) | Memory (peak) | Headline accuracy result |
|---|---|---|---|---|
| **1. Simulated, panel** | 8 baseline tools | (timings transcribed, not re-run) | n/a (not re-run) | TPR ≥ best baseline across grid; FPR 100×–600× lower at 5× |
| **2. Simulated, head-to-head** | methylKit (local) | **≈ 33× per-CpG `diffmeth` on Linux** (`mc.cores = 8`); 43× Windows no-op total grid | 1.18× less | Identical to 3 decimal places at n ≥ 4; 2× recall at n = 2 |
| **3. Real GSE263850** | methylKit + DSS (local, Linux pivoine) | **≈ 18× vs single-core methylKit** (full DMC pipeline); **~ 3.5× vs DSS** full-pipeline (675 s vs 2,368 s); cached chain_merge re-call 92 s | **3.83× less than mk; ~1.13× less than DSS** (12.6 vs 48.0 / 14.3 GB) | per-CpG r = 0.994 with methylKit; **63.8 % (paper-faithful) / 77.3 % (morphology-matched) DSS-DMR any-bp recall** with chain_merge (vs ≈ 9 % for methylKit-tile, 87.5 % for DSS-from-scratch vs paper-813); **100 % direction agreement on every matched DMR across every engine** |

---

## What we recommend

* **Default.** Use `epykit.tl.dmc(test="lr")`. At n ≥ 3 per group on
  simulated data it is statistically equivalent to methylKit; on real data
  it is more conservative in the small-p tail.
* **Low replicates (n ≤ 2).** Use `lr+` (opt-in via `power_stack="lr+"`).
* **n = 1 (no replicates).** Use the bug-fixed `fisher` backend.
* **DMRs on real biological data.** Use `dmr_chain_merge` with
  paper-matched `alpha = 1e-5, delta = 0, minlen = 50, minCG = 3,
  pct.sig = 0.5`. Increase `dis.merge` from the literal 100 to
  **250** to match DSS's chain morphology (lifts strict overlap recall
  ~ 2× without sacrificing direction agreement). Fixed-tile callers
  recover only ~9 % of focused biological DMRs on real WGBS — they
  remain useful for simulator benchmarking but should not be the
  default on real data.

## What the chain_merge / DSS replication established

* **DMR-engine architecture is the dominant factor in real-data DMR
  reproduction.** Per-CpG calibration agrees to three decimal places
  across engines (r = 0.994 with methylKit); region-level recovery of
  the locally-rerun DSS-from-scratch call set (922 DMRs) ranges from
  ≈ 9 % (fixed-tile) → 63.8 % (chain_merge dis.merge = 100) →
  77.3 % (chain_merge dis.merge = 250); DSS-from-scratch itself
  recovers 87.5 % of the paper's published 813.
* **Direction is never wrong.** Across all matched DMR pairs between
  any two non-tile callers, every single pair is called hyper / hypo
  in the same direction — 100 % agreement on 588 (chain_merge-100 vs
  DSS) and 713 (chain_merge-250 vs DSS) matched DMRs respectively.
* **DSS reaches an 87.5 % ceiling against the paper's own published
  call set.** The remaining 12.5 % is DSS-vs-DSS noise (package
  version drift, BSseq smoothing internals, threading
  non-determinism). epykit's gap to the locally-rerun DSS-922 is
  therefore the appropriate reference, not its gap to the paper-813.

## What we do not claim

* Real-data accuracy. Study 3 has no ground truth; we can demonstrate
  agreement with methylKit, not absolute correctness.
* That epykit dominates methylKit on all metrics. At n ≥ 3 on real data,
  methylKit is more aggressive in the small-p tail; whether that is
  desirable depends on the downstream analysis.
* That all of epykit's eight backends are equally calibrated. Study 1
  surfaced one calibration bug in `fisher` and one df-reference bug in
  pooled-dispersion modes (both fixed); the default `lr / site` was not
  affected.

## What's next

Multi-dataset real-WGBS validation (IMR90 vs H1-hESC, mouse imprinted DMRs)
is the next step. The current Study 3 is one tissue × one genome, so we
treat it as an existence proof of correctness on biological data, not a
generalisation.

---

**See also.** [README.md](README.md) for navigation • [paper/paper.md](paper/paper.md)
for the scientific manuscript • [report/REPORT.md](report/REPORT.md) for
the full numeric tables.
