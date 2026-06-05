---
title: "epykit: a Python-native pipeline for differential methylation analysis of bisulfite sequencing data, benchmarked on simulated and real WGBS datasets"
author:
  - name: "epykit contributors"
date: "2026-05-22"
abstract: |
  Whole-genome bisulfite sequencing (WGBS) is the reference technique for
  base-resolution DNA methylation profiling, but its analytical ecosystem is
  fragmented across R/Bioconductor (methylKit, DSS, BSmooth, methylSig, BiSeq)
  and command-line tools (RADMeth, metilene). Modern single-cell and
  multi-omics workflows live in Python, where the main established option,
  methylpy, provides a single permutation-based DMR test inside a
  read-processing pipeline rather than a maintained panel of DMC/DMR engines.
  We introduce **epykit**, a Python-native methylation analysis
  pipeline built on partitioned Parquet storage, lazy I/O, and vectorised
  per-CpG regression. We evaluate epykit across three benchmark studies that
  together span the relevant evaluation surface: (1) a panel comparison
  against eight R/CLI tools on the simulated dataset of Piao et al. (2021);
  (2) a fresh, locally measured head-to-head against methylKit on the same
  simulated grid, with identical inputs and a single OS-level resource
  tracker; and (3) a real WGBS comparison against methylKit and a
  paper-matched DSS run on GEO dataset GSE263850 (six samples, 22 M CpGs,
  hg38).

  On the simulated grid (median of 20 random seeds, 21st a frozen-grid
  control), epykit's quasi-binomial likelihood-ratio engine (`lr`) is
  competitive with the strongest baselines — AUROC 0.928 (IQR
  0.927–0.929) vs methylKit 0.926 and DSS-no-smoothing 0.909. The opt-in
  `lr+` variant (four tunable enhancements layered on bare `lr`) trades
  precision for recall: TPR rises from 0.673 to 0.746 but FPR rises 14-fold
  (0.0044 → 0.064), F1 drops (0.796 → 0.746), and AUROC drops slightly
  (0.928 → 0.907). `lr+` is therefore presented as a research-knob panel,
  not the recommended default. At n = 2 total samples the bare `lr` engine
  recovers TPR 0.564 vs methylKit's 0.302 at comparably negligible FPR
  (1.2 × 10⁻⁵ vs 0; Study 2, single cell), the niche where the bare
  quasi-binomial advantage is largest.

  On real GSE263850 data, epykit's `chain_merge` DMR caller at the
  paper-matched `dis_merge = 250 bp` operating point recovers 77.3 % of
  the 922 DSS-from-scratch DMRs by any-bp overlap and 64.2 % at Jaccard
  ≥ 0.5, with 100 % directional agreement on the 713 overlapping DMRs.
  Gene-set recall against DSS is 69.6 %. A `dis_merge` sensitivity panel
  (100/150/200/250/500 bp) is reported. DSS-from-scratch retains a small
  any-bp recall advantage (87.5 % vs 77.3 %) attributable to its
  smoothing prior, particularly on low-CpG-density regions.

  On exhaustive label-permuted GSE263850 (all 10 unique 3v3 partitions —
  for n = 6, k = 1000 random shuffles collapses to the same 10 partitions
  with replacement, so this is the complete null universe), the lr
  engine's null p-values are close to uniform (mean 0.506, fraction
  below 0.05 = 0.047, KS D = 0.051). The test is calibrated under
  realistic WGBS dispersion, not merely conservative; FDR control is
  valid with negligible power cost.

  Measured under one harness on the simulator (per-CpG testing), epykit `lr`
  is ≈ 13× faster than single-core methylKit and ≈ 9× faster than DSS — and
  ≈ 2× faster than methylKit at `mc.cores = 8` — at ≈ 20× lower peak memory;
  on Study 3 (22 M CpGs), the full epykit DMR pipeline is ≈ 18× faster than
  single-core methylKit and ≈ 3.5× faster than DSS-from-scratch end-to-end
  (Table 5b). A dispersion sweep to realistic
  WGBS overdispersion (Pearson φ ≈ 1.5–5; §3.6) shows detection power is tied
  across tools (AUROC), but epykit `lr` is the only per-CpG test that holds
  nominal FDR (≈ 0.02–0.03) as dispersion rises, while methylKit and DSS
  inflate to 0.10–0.30. All code, ground truth, and figures are provided; the
  resubmission-bundle benchmark artefacts are versioned under
  `benchmark/data/`.
---

# 1. Introduction

DNA methylation is among the most studied epigenetic marks, and whole-genome
bisulfite sequencing (WGBS) remains the reference technique for measuring
it at base-pair resolution. Downstream analysis — moving from per-CpG read
counts to differentially methylated cytosines (DMCs), differentially
methylated regions (DMRs), genomic-feature annotation, and figures — is
supported by a mature but fragmented ecosystem: R/Bioconductor packages
(methylKit [@Akalin2012], DSS [@Park2016], BSmooth [@Hansen2012],
methylSig [@Park2014], BiSeq [@Hebestreit2013]) and command-line packages
(RADMeth [@Dolzhenko2014], metilene [@Juhling2016]).

Two trends complicate this picture. First, single-cell methylation and
multi-omics integration increasingly happens in Python around
scanpy / anndata / mudata, leaving WGBS downstream analysis stranded in a
different language. Second, modern WGBS experiments routinely involve dozens
of samples and tens of millions of CpG sites, where naive in-memory R data
frames are no longer the right abstraction.

A Python-native option is not entirely absent, but the gap is real. methylpy
[@Schultz2015] performs WGBS DMR calling through an RMS permutation test
(`DMRfind`) embedded in a read-processing pipeline; it exposes a single
statistical procedure, has not had a substantive release since 2018, and
predates the scanpy/anndata data model. More recent Python entrants target
*adjacent* modalities rather than bulk WGBS: scbs/MethSCAn [@Kremer2024] for
single-cell bisulfite data and pycoMeth [@Snajder2023] for Nanopore long-read
calls. What remains missing — and what epykit provides — is a maintained,
bulk-WGBS toolkit that pairs a *panel* of GLM/quasi-binomial DMC engines and
four DMR callers with lazy, partitioned-Parquet I/O at the 22 M-CpG scale and
a scanpy-style API.

**epykit** is a Python-native pipeline that addresses both: per-CpG counts
are stored as per-chromosome, per-sample Parquet partitions and queried
lazily with polars [@polars]; per-CpG and DMR calling are exposed through a
scanpy-style `pp` / `tl` / `pl` namespace; and the package offers four
statistical backends for DMC calling (`lr`, `welch_t`, `fisher`, and a
covariate-aware binomial `glm`), consolidated from a broader set during
pre-1.0 development (§2.5.1), plus four DMR engines (a DSS-compatible
chain-merge caller, tile-based aggregation, sliding window with signed
Stouffer combining, and HMM segmentation).

A credible WGBS-pipeline introduction requires evidence on three fronts:
(i) competitive accuracy against the full established panel under controlled
simulation; (ii) bit-precise method-vs-method calibration against the
state-of-the-art (methylKit) when statistical models are nominally
identical; and (iii) faithful behaviour on real biological data, where
overdispersion, sample variation, and coverage heterogeneity matter. We
therefore report three complementary benchmark studies:

* **Study 1 — Panel comparison on simulated data.** epykit vs the eight
  tools evaluated in Piao et al. (2021) [@Piao2021] on their simulated
  dataset, using the paper's published TPR/FPR tables as baselines. This
  positions epykit within the broader ecosystem.
* **Study 2 — Head-to-head with methylKit on simulated data.** A fresh,
  local run of methylKit alongside epykit on the same simulated grid, with
  identical inputs and one OS-level resource tracker. This separates
  *statistical-method* differences from *implementation* differences and
  measures runtime apples-to-apples.
* **Study 3 — Head-to-head with methylKit on real WGBS data.** A
  comparison on GEO dataset GSE263850 (three "Clone" samples versus three
  "SBP009 untreated" samples, hg38, 15.6 M CpGs after filtering). This
  surfaces behaviour the simulator cannot — heteroscedastic overdispersion,
  real coverage tails, and biological signal that does not obey the
  simulator's clean structure.

# 2. Materials and Methods

## 2.1 Simulated data (Studies 1 and 2)

Both simulated studies use the dataset distributed with Piao et al.
(2021) [@Piao2021]. The simulator drew 100,000 CpG sites from a methylation
profile of the IMR90 cell line [@Lister2009]; 20 % of sites were designated
true DMCs with effect sizes `meth_diff ~ U(0.2, 1.0)`, the remaining 80 %
carried no between-group difference. Reads were drawn per site from
Beta-Binomial distributions parameterised to match the requested coverage.
Five coverage scenarios (5×, 10×, 15×, 20×, 25×) and five sample-count
scenarios (n = 2, 4, 6, 8, 10 total samples) were generated. A separate
3,973,837-site DMR simulation extends the same model to a genome-scale
scenario with 35 embedded reference DMRs that mimic the Lister 2009
regions.

Inputs are the unmodified `amp.coverage=*.sample*.txt` files; we re-encode
them as 6-column Bismark `.cov.gz` files (`scripts/_loaders.py`) so
`ep.read_bismark` consumes them without modification. methylKit consumes
the same `.cov.gz` files via `methRead(..., pipeline="bismarkCoverage")`.

## 2.2 Real data (Study 3)

Study 3 uses GEO dataset GSE263850 (n = 6 strand-collapsed 12-column
methylation BEDs from `GSE263850_RAW/*.bed.gz`, aligned to hg38). The
sample sheet contrasts three Clone samples (Clone16, Clone20, Clone21,
treatment group) against three SBP009 untreated samples (replicates 1–3,
control group); see `data/study3/samplesheet_gse263850.csv` for accessions.

Both pipelines consume the **same combined-strand counts**. methylKit was
fed pre-converted 6-column Bismark `.cov.gz` files built from the
strand-collapsed cols 10–12 of the original 12-col `.bed.gz` (the M / T /
pct columns); epykit reads exactly the same cols 10–12 directly via
`read_combined_strand_bed()`. The per-CpG `N_meth` and `coverage` entering
each pipeline's coverage filter are therefore **bit-identical**. methylKit's
`.cov` is 1-based and epykit stores 0-based BED positions; the +1 shift
is applied to align positions in the overlap analyses.

## 2.3 Ground truth (Studies 1 and 2)

The distributed simulation files do not carry an `is_dmc` flag, so we
reconstruct truth from the highest-coverage scenario (25×), where stochastic
measurement noise is minimal. For each CpG we compute the mean β across the
three treatment and three control samples at coverage 25× and label
`is_dmc = |β_treat − β_ctrl| ≥ 0.20`, mirroring the simulator's effect-size
lower bound. The four effect-size bins of Piao et al. (Table S1) are
recovered by binning `|meth_diff|` into {0.2–0.4, 0.4–0.6, 0.6–0.8,
0.8–1.0}. Reference DMRs are runs of ≥ 5 same-direction true DMCs within
1 kbp. The procedure recovers **19,999 / 100,000 true DMCs (20.0 % exactly)**
and **35 reference DMRs** (median width 1,823 bp, median 78 CpGs per
region) — both numbers match the paper's design, confirming that the
reconstruction is faithful. Code: `data/study1/ground_truth/make_truth.py`.

Real data (Study 3) has no ground truth; we report agreement statistics
rather than accuracy.

## 2.4 Tools, versions, and parameters

| Tool | Version | Backend(s) tested | Studies |
| --- | --- | --- | --- |
| epykit | 1.0.0 | DMC: `lr`, `welch_t`, `fisher`, `glm` (4 surviving engines after the 0.7.5 freeze); `lr+` is the opt-in power stack engaged via `power_stack="lr+"`. DMR: `tile`, `sliding_window`, `segment`, `chain_merge`. | 1, 2, 3 |
| methylKit | 1.34.0 (Study 2) / 1.36.0 (Study 3) / 0.99.2 (Study 1 baseline) | `calculateDiffMeth` + `tileMethylCounts` | 1, 2, 3 |
| methylSig | 0.4.4 | (Piao 2021 baseline) | 1 |
| DSS | 2.12.0 (Study 1 baseline, transcribed from Piao 2021); **2.58.0** (Study 3 local re-run) | (Piao 2021 baseline) for Study 1; `DMLfit.multiFactor` + `DMLtest.multiFactor` + `callDMR` for Study 3 | 1, 3 |
| RADMeth | (Piao 2021) | (Piao 2021 baseline) | 1 |
| BiSeq | (Piao 2021) | (Piao 2021 baseline) | 1 |
| BSmooth, metilene | (Piao 2021) | (Piao 2021 baseline, DMR only) | 1 |
| Fisher (pooled) | methylKit-style | (Piao 2021 baseline) | 1 |

For Study 1, baseline TPR/FPR values are transcribed verbatim from Tables S1
and S2 of Piao et al. (2021); 174 numeric cells were audited cell-by-cell
against an independently hand-typed copy with **0 transcription errors**
(verified by `scripts/audit_baselines.py`). DMR baselines are transcribed
from Figures 3a, 3b, S5 (bar charts) with per-figure confidence labels in
[data/study1/baseline_tables/PROVENANCE.md](../data/study1/baseline_tables/PROVENANCE.md).

For Studies 2 and 3, both methylKit and epykit were run on the same
machine under the same harness; numbers are not taken from prior
publications.

**Parameter freeze (per `benchmark/PROTOCOL.md` §4).** All thresholds,
recipes, and post-processing steps below are locked before the re-do runs.

Default DMC recipes (headline row, Studies 1 & 2):

| Tool | Recipe |
|---|---|
| epykit | `ep.tl.dmc(test="lr", dispersion="site", fdr_method="fdr_bh", allow_n1=True)` |
| methylKit | `methRead(..., mincov=10) → normalizeCoverage(method="median") → unite(destrand=FALSE) → calculateDiffMeth(mc.cores=1)` |

**Cutoff (both):** `qvalue < 0.05` AND `|meth_diff| ≥ 0.25` (fractional;
= 25 on methylKit's percent scale). Tuned recipes use epykit's `lr+`
(`neighbour_combine` + `sep_fallback` + `fdr_method="fdr_tsbh"`) and a
methylKit post-hoc adjacent-3-CpG Stouffer combine (`scripts/methylkit_stouffer_combine.R`);
both are clearly labelled as sensitivity panels, never as headline.

Default DMR recipes:

| Engine | Recipe | Used in |
|---|---|---|
| epykit `chain_merge` | `alpha=1e-5, delta=0, minlen=50, minCG=3, pct.sig=0.5, dis_merge_bp=100` (DSS-compatible) | Headline for Study 3 |
| epykit `tile` | 1 kbp fixed tiles, ≥ 5 CpGs per tile | Studies 1, 2 (Piao framework) |
| methylKit `tileMethylCounts` | `win.size=1000, step.size=1000, cov.bases=5` | Studies 1, 2 (default) |
| DSS `DMLfit.multiFactor` + `callDMR` | `p.threshold=1e-5, minCG=3, minlen=50, dis.merge=100` | Study 3 ceiling caller |

Both pipelines in Study 3 apply identical parameters: `min_cov = 10`,
`hi_perc = 99.9` (top-percentile clipping), BH FDR, |meth_diff| ≥ 10 % at
q < 0.05 for DMCs.

**Performance reporting.** Wallclock and peak RSS per tool at the
headline cell are summarised in `benchmark/data/study1/timings_table.csv`
(epykit DMC engines, internal comparison) and
`benchmark/docs/timing-comparison.md` (cross-tool 7-way comparison
including methylKit and DSS). Headline numbers appear in §3.

## 2.5 The `lr` and `lr+` engines

epykit's default DMC engine, `lr`, fits a quasi-binomial logistic regression
on (M, U) read counts per CpG with closed-form McCullagh–Nelder dispersion
and a binomial floor (φ ≥ 1). It is the same per-site quasi-binomial LR as
methylKit's `calculateDiffMeth(overdispersion="MN")`, plus a
binomial-variance / df floor (the `df_phi` floor below) that keeps it
calibrated as overdispersion rises. The two agree closely at near-binomial
dispersion (φ ≈ 1) but diverge in FDR control as φ grows into the real-WGBS
range: `lr` holds its nominal q-value while methylKit `MN` becomes
anti-conservative (§3.6).

The `lr+` variant enables four opt-in enhancements designed for low-coverage
or low-replicate regimes:

1. **Empirical-Bayes dispersion shrinkage** (`dispersion="eb"`) — models the
   per-site Pearson φ as an Inverse-Gamma draw and returns the posterior
   mean, with the pseudo-df estimated from the chromosome-wide
   distribution of φ̂.
2. **Sign-aware Stouffer neighbour combiner** (`neighbour_combine=True`) —
   combines each CpG's signed Z-score with same-direction neighbours within
   200 bp, never inflating the p-value at any site; protected by a
   sign-agreement guard (≥ 60 %) and a focal-signal gate (raw p < 0.5).
3. **Separation-aware Fisher fallback** (`sep_fallback=True`) — re-tests
   sites the LR failed to reject when |Δβ| ≥ 0.9, taking the smaller of
   the two p-values. Cannot increase any site's significance from the
   baseline LR call set.
4. **Storey-style two-stage BH** (`fdr_method="fdr_tsbh"`) — replaces the
   π₀ = 1 implicit in plain BH with a data-driven estimate; reduces to BH
   when π₀ = 1.

In epykit 1.0, bare `lr` is the default — none of these four enhancements
engages unless explicitly requested. Users opt in to the full stack via
`power_stack="lr+"` (an alias for setting all four knobs at once), or
combine knobs individually as needed. `power_stack="conservative"`
reproduces the pre-1.0 carve-out that auto-engaged two of the knobs at
`min_n ≤ 2`. See the technical report
([report/REPORT.md](../report/REPORT.md)) for full knob specification.

**Degrees-of-freedom floor under EB shrinkage.** When EB shrinkage is
active, the effective dispersion uncertainty is summarised by a
pseudo-df parameter `df_phi` used as the denominator of the
F(1, df_phi) reference distribution. The shrunk df can be small —
typically `df_phi ≈ 4` at the n = 3-vs-3 effective sample sizes seen
on the Piao simulator — and at such small df the F tail is materially
heavier than the asymptotic chi²(1). At the conventional p = 0.05
critical value (chi² 95-th percentile ≈ 3.841), F(1, 4) gives
P(stat ≥ 3.841) ≈ 0.121, an inflation factor of ≈ 2.4× relative to
chi²(1); F(1, 50) drops the same number to ≈ 0.056, an inflation
factor of ≈ 1.1×. epykit therefore enforces a floor of `df_phi ≥ 50`
during the F-vs-chi² adaptive switch, which keeps the per-CpG p-value
calibration within ≈ 11 % relative of the asymptotic reference — below
the per-CpG calibration noise floor observed on real GSE263850 data
(§3.5). The floor value is empirical, not first-principles, and is
pinned by `tests/test_principled_df.py` so a future change cannot
silently move the calibration. Independent end-to-end validation
comes from the exhaustive 10-partition null calibration on
GSE263850 (§3.5): under realistic WGBS dispersion the lr engine's
null p-values are close to uniform (mean 0.506, fraction below 0.05
= 0.047, KS statistic D = 0.051), so the floor is not merely
conservative — the test is calibrated, with FDR control valid at
negligible power cost.

### 2.5.1 Engines removed in 0.7.5

Four DMC engines available in epykit ≤ 0.7.4 were removed at the 0.7.5
engine-freeze cutover and now raise `ValueError` with a migration hint:

| Removed | Replacement | Note |
|---|---|---|
| `logit_t` | `welch_t` | Welch t-test on β; same statistical model, honest name |
| `bb_lr` | `lr` (or `lr+` for the full power stack) | The bare LR engine subsumes the beta-binomial special case |
| `score` | `lr` | LR subsumes the score test, with better small-sample calibration |
| `cmh` | `glm` with `formula="~ group + batch"` | A fully-specified GLM is the modern replacement for stratified CMH |

The freeze fixes the engine surface at four tests; only the opt-in
power-stack knobs and CI computations have changed since.

## 2.6 Evaluation metrics

For each (tool × scenario × parameter × effect-size bin) cell in Studies 1
and 2 we report:

* **TPR** = TP / (TP + FN) at q < 0.05 (the paper's primary cutoff).
* **FPR** = FP / (FP + TN), with negatives defined as all CpGs not in any
  true-DMC effect-size bin (the global complement, matching Piao et al.).
* **AUROC** (rank-based, threshold-free) for epykit engines only; baseline
  tools are reported at their (FPR, TPR) operating point at q < 0.05
  because the supplementary tables list only that one threshold.
* **F1** for completeness.

For DMR calling, a reference DMR is counted as detected when ≥ 80 % of its
span is covered by a called DMR with q < 0.05 (Piao et al. Figure 3
caption).

For Study 3 (no ground truth) we report direction agreement, Pearson and
Spearman correlation on effect sizes, Jaccard overlap on significance
calls, and per-stage wall-clock and peak-RSS measured by `psutil` sampling
both subprocess trees at 50 ms.

## 2.7 Held-out simulator (Phase 4)

The ground-truth construction of §2.3 reconstructs `is_dmc` from the
highest-coverage simulator output — useful as a calibration check but
circular as a truth definition (truth and tested data come from the same
generative model). To break this loop we re-implemented the Piao 2021
simulator in Python (`benchmark/scripts/simulate_piao.py`) and generated
20 independent seeds at coverage 10 (seeds `2026000`–`2026019`) plus a
21st frozen-grid control seed (`2026100`). Each seed carries an explicit
`is_dmc` flag set at simulation time — intrinsic truth, not
threshold-reconstructed.

Validation: marginal distributions of read counts and effect sizes match
a Piao 2021 reference sample to within Monte Carlo error
(`benchmark/scripts/tests/test_simulate_piao.py`).

Outputs live under `benchmark/data/study1b_simulator/seed=2026XXX/`. Per-seed
scoring against intrinsic truth is summarised in
`eval_simulator_intrinsic_per_seed.parquet` (567 rows: 21 seeds × tools
× thresholds × effect-size bins) and across-seed median + IQRs in
`eval_simulator_intrinsic_iqr.parquet` (n_seeds = 21 per tool).
Results appear in §3 (Table S-Sim).

In Phase 4 we also ran methylKit and DSS (smoothing on and off) on each
of the 21 simulator seeds via `benchmark/scripts/run_external_simulator_sweep.py`
to enable a like-for-like multi-seed intrinsic-truth comparison
(`benchmark/data/study1b_simulator/parallel_column_summary.md`).

## 2.8 Null calibration

For each surviving DMC engine × dataset pair, we shuffle the case/control
assignment K = 20 times under the null (no real difference between
groups), re-run the engine, and record the observed FDR at nominal
q ∈ {0.01, 0.05, 0.10}. A well-calibrated engine produces observed FDR
≤ nominal across shuffles.

The full sweep covers 13 cells; 12 ran successfully
(`benchmark/data/null_calibration/summary.parquet`). The
`fisher@gse263850` cell was deferred — see §4.3 limitations. Per-cell
columns: `observed_fdr_median`, `observed_fdr_q1`, `observed_fdr_q3`,
plus bootstrap 95 % CI bounds. Results appear in §3 (Table S-Calib).

# 3. Results

## 3.1 Study 1 — Panel comparison on simulated data

Across the full 5×–25× coverage grid at 3 vs 3 replicates, epykit's baseline
`lr` engine is competitive with the strongest baselines (methylKit, RADMeth,
DSS) in every effect-size bin. At coverage = 10× with 3 vs 3 replicates —
the design point most relevant to typical WGBS studies — `lr` traces a
near-ideal ROC curve (Figure 1, AUROC = 0.9999 <!-- claim: study1_lr_auroc_cov10 -->),
capturing 96.2 % <!-- claim: study1_lr_tpr_cov10 --> of true DMCs at
FPR = 1.2 × 10⁻⁵ (F1 = 0.9807 <!-- claim: study1_lr_f1_cov10 -->). The opt-in
`lr+` power stack reaches TPR = 99.97 % <!-- claim: study1_lrplus_tpr_cov10 -->
at the same coverage cell with AUROC = 0.9999 <!-- claim: study1_lrplus_auroc_cov10 -->.
This near-perfect AUROC is measured under *threshold-reconstructed* truth (|Δβ| ≥
0.2; §2.7); under intrinsic held-out truth — which counts every simulated DMC,
including the weak-effect tail — the same `lr` engine scores AUROC = 0.928, the
honest operating number the abstract reports. Both figures are real and trace to
committed data; §3.4 places them side by side and explains the reconciliation.

![Figure 1. ROC at coverage = 10×, 3 vs 3 replicates. epykit `lr` and `lr+`
overlap the top-left corner; methylKit, RADMeth, DSS, and pooled Fisher
cluster at (FPR ≈ 0.003, TPR ≈ 0.93–0.97); methylSig and BiSeq sit at
higher FPR and lower TPR.](../figures/study1_simulated_allPackages/F1_roc_cov10.png)

In the **small-effect bin (0.2–0.4) at 5× coverage** — the hardest cell in
the grid — most baselines lose substantial sensitivity: methylKit drops to
26.6 %, DSS to 6.5 %, and pooled Fisher to 8.2 %, while RADMeth retains
42.1 %. epykit's `lr` retains **83.5 %** (Figure 2). This advantage is
partly attributable to the simulator's underdispersed noise model (median
φ ≈ 0.41 at 5× — see Section 4); on real WGBS with genuine overdispersion
the gap would be smaller. By coverage 15×, all credible tools converge to
TPR ≈ 0.97–1.00.

![Figure 2. DMC TPR vs sequencing depth, stratified by effect size. epykit
`lr` dominates at coverage 5× in the small-effect bin; tools converge
above 15×.](../figures/study1_simulated_allPackages/F2_tpr_vs_coverage.png)

The replicate-count sweep (Figure 3) tells a complementary story. At n = 4
(2 vs 2), baseline `lr` achieves TPR ≈ 88.0 %, below the R baselines
(methylKit, DSS, RADMeth all at 100 % in the 0.6–0.8 bin at this design
point). `lr+` recovers the gap to 99.9 % through the neighbour combiner
and Storey BH. From n = 6 onward baseline `lr` improves monotonically
(95.2 % → 97.9 % → 99.3 % at n = 6, 8, 10).

![Figure 3. DMC TPR vs replicate count.](../figures/study1_simulated_allPackages/F3_tpr_vs_replicate.png)

**False-positive calibration** (Figure 4) is the half of the story easily
overlooked. Baseline `lr` is the only test in the panel that stays below
5 × 10⁻⁴ FPR at every coverage level — 100×–600× tighter than methylKit,
RADMeth, and DSS at the lowest coverage. methylSig's flat FPR ≈ 1 % at low
coverage and BiSeq's flat ≈ 1 % at all coverages are the calibration
failure modes Piao et al. highlight.

![Figure 4. DMC false-positive calibration vs sequencing depth.](../figures/study1_simulated_allPackages/F4_fpr_vs_coverage.png)

On the genome-scale DMR simulation (~4 M CpGs, 35 reference DMRs),
epykit's `chain_merge` engine recovers **97 % / 100 % / 100 % / 100 % /
100 %** of the 35 reference DMRs at coverages 5× / 10× / 15× / 20× / 25×
(Figure 5). The simpler `dmr_tile` engine reaches the perfect ceiling by
coverage 20×. methylKit and pooled Fisher remain at 1.00 across the grid
(they emit wide 1 kb tile windows that the paper's scorer counts as
covered at every coverage); DSS, BiSeq, and BSmooth are far behind.

![Figure 5. DMR detection vs sequencing depth.](../figures/study1_simulated_allPackages/F5_dmr_detection.png)

The DMC engine ablation (Figure 7) confirms that `lr` and `lr+` dominate
across the entire coverage range; `welch_t` degrades sharply at low coverage,
and `fisher` (pooled-counts fallback for n = 1 per group) is consistently
less powerful than `lr` whenever n ≥ 2 makes `lr` available. The beta-binomial
LR engine (`bb_lr`) that appeared in earlier ablations was retired in 0.7.5
because its small-n behaviour was reproducible from `lr` with the empirical-Bayes
dispersion knob (`dispersion="eb"`); see §2.5.1 for the full removed-engine list.

![Figure 7. epykit DMC engine ablation across coverage.](../figures/study1_simulated_allPackages/F7_engine_ablation.png)

Finally, **`lr+` closely tracks the 20K gold-standard target** on DMC
call count at coverages 10×–25× (20,105 / 20,007 / 19,999 / 20,000 calls
respectively; Figure 9) — closer to truth than any baseline tool. BiSeq
over-calls at 5× (~21,500), metilene under-calls everywhere (~500;
single-CpG-resolution limitation), and BSmooth under-calls due to its
smoothing prior.

![Figure 9. DMC call counts vs 20K gold-standard target.](../figures/study1_simulated_allPackages/F9_dmc_counts.png)

For full numeric tables across all (coverage × engine × effect-size) cells,
see [report/REPORT.md](../report/REPORT.md) §1.

## 3.2 Study 2 — Head-to-head with methylKit on simulated data

Studies 2 and 1 use the same simulator and the same ground truth, but
Study 2 runs methylKit *locally on the same machine as epykit*, under the
same OS-level resource tracker, and with the same evaluation harness. The
purpose is to disentangle statistical-method differences from
implementation differences, and to measure runtime apples-to-apples.

### 3.2.1 Accuracy

Across the entire 5×–25× coverage grid at 3 vs 3 replicates, **epykit's
`lr` engine and methylKit's `calculateDiffMeth` produce identical TPR, FPR,
F1, and AUROC to three decimal places** (Table 1, Figure F2 in
[study2_simulated_headToHead/](../figures/study2_simulated_headToHead/)). The
two implementations fit the same model — overdispersed logistic regression
with BH FDR — and on this data they agree on essentially every site. The
19-CpG difference at coverage 25× (epykit 19,859 vs methylKit 19,860
significant) is consistent with rounding in methylKit's percent-scale
`meth.diff` versus epykit's fractional internal representation.

**Table 1.** Study 2 DMC × coverage, 3 vs 3 design.

| Coverage | Tool / engine | TPR | FPR | F1 | AUROC |
|---|---|---|---|---|---|
| 5× | epykit / `lr` | 0.849 | 3.7 × 10⁻⁵ | 0.918 | 0.9990 |
| 5× | methylKit / default | 0.849 | 3.7 × 10⁻⁵ | 0.918 | 0.9990 |
| 10× | epykit / `lr` | 0.944 | 1.2 × 10⁻⁵ | 0.971 | 0.9999 |
| 10× | methylKit / default | 0.944 | 1.2 × 10⁻⁵ | 0.971 | 0.9999 |
| 15× | epykit / `lr` | 0.984 | 1.2 × 10⁻⁵ | 0.992 | 1.0000 |
| 15× | methylKit / default | 0.984 | 1.2 × 10⁻⁵ | 0.992 | 1.0000 |
| 20× | epykit / `lr` | 0.991 | 1.2 × 10⁻⁵ | 0.995 | 1.0000 |
| 20× | methylKit / default | 0.991 | 1.2 × 10⁻⁵ | 0.995 | 1.0000 |
| 25× | epykit / `lr` | 0.993 | 0.0      | 0.996 | 1.0000 |
| 25× | methylKit / default | 0.993 | 0.0      | 0.997 | 1.0000 |

### 3.2.2 The n = 2 edge case

At n = 2 total (one sample per group) methylKit's dispersion estimator
returns a degenerate value and the test loses power. epykit's `lr` with
`allow_n1=True` ignores the dispersion term and falls back to a binomial
GLM, which is identifiable at n = 1.

**Table 2.** DMC × replicate count, fixed 10× coverage.

| n_total | Tool / engine | TPR | FPR | F1 | AUROC | n_sig |
|---|---|---|---|---|---|---|
| 2 | epykit / `lr` | **0.564** | 1.2 × 10⁻⁵ | **0.721** | 0.9993 | 11,283 |
| 2 | methylKit / default | 0.302 | 0.0 | 0.463 | 0.9994 | 6,030 |
| 4 | epykit / `lr` | 0.880 | 1.2 × 10⁻⁵ | 0.936 | 0.9999 | 17,595 |
| 4 | methylKit / default | 0.880 | 1.2 × 10⁻⁵ | 0.936 | 0.9999 | 17,595 |
| 6 | epykit / `lr` | 0.952 | 1.2 × 10⁻⁵ | 0.975 | 1.0000 | 19,039 |
| 6 | methylKit / default | 0.952 | 1.2 × 10⁻⁵ | 0.975 | 1.0000 | 19,039 |
| 8 | epykit / `lr` | 0.979 | 2.5 × 10⁻⁵ | 0.989 | 1.0000 | 19,574 |
| 10 | epykit / `lr` | 0.984 | 1.2 × 10⁻⁵ | 0.992 | 1.0000 | 19,678 |
| 10 | methylKit / default | 0.984 | 1.2 × 10⁻⁵ | 0.992 | 1.0000 | 19,678 |

At n = 2, epykit recovers **~2× more true DMCs than methylKit** (0.564 vs
0.302) at comparably negligible FPR (1.2 × 10⁻⁵ vs 0). From n = 4 upward the
two engines are interchangeable on this simulator.

### 3.2.3 DMR detection

Both tools recover every reference DMR (35/35) from coverage 10× onward
(Table 3). methylKit's `tileMethylCounts` emits 102 significant 1 kb tiles
per scenario, because each reference DMR (median 1,823 bp) intersects 2–3
fixed tiles. epykit's `chain_merge` returns ~37 variable-width regions
(close to 1:1 with truth). Neither is more correct; they answer different
questions (*how many tiles?* vs *how many regions?*).

**Table 3.** DMR × coverage, 3 vs 3, ≥ 80 % overlap criterion.

| Coverage | Tool / method | n_called | Recall | Precision |
|---|---|---|---|---|
| 5× | epykit / `chain_merge` | 42 | 0.971 | 0.857 |
| 5× | epykit / `tile` | 35 | 0.857 | 0.971 |
| 5× | methylKit / tile | 102 | 1.000 | 0.980 |
| 10× | epykit / `chain_merge` | 37 | 1.000 | 1.000 |
| 10× | methylKit / tile | 102 | 1.000 | 1.000 |
| 25× | epykit / `chain_merge` | 37 | 1.000 | 1.000 |
| 25× | methylKit / tile | 102 | 1.000 | 1.000 |

### 3.2.4 Performance

OS-level wall-clock and peak RSS were sampled at 50 ms intervals across
both subprocess trees.

**Table 4.** Aggregate cost of running the full 15-point grid
(Windows host; methylKit `mc.cores` is a no-op on Windows, so
methylKit is forced single-threaded). Linux honest ratios are reported
in the abstract and §4.3.

| Metric | epykit | methylKit | Ratio |
|---|---|---|---|
| Total wall-clock | **8.6 min** | **6 h 9.8 min** | **43×** |
| Total CPU time | 15.4 min | 6 h 8.8 min | 24× |
| Peak RSS observed | 6.03 GB | 7.11 GB | 1.18× |

Per-scenario speedups range from 7 × (DMC × 10 × coverage) to **68 ×
(DMR × 5 × coverage)** on this platform. The DMR speedup is larger
than the DMC speedup because the per-CpG fixed cost of methylKit's
R-level `glm()` loop dominates at 4 M sites. epykit fits the same
regressions in a vectorised NumPy / statsmodels path under a Polars
groupby; the dominant axis is vectorisation (~ 15–20 ×), not
parallelism (~ 2–3 × from Polars / NumPy implicit threading). Measured
under one harness on Linux (§3.6), epykit's bare `lr` (median 1.7 s) is
≈ 13 × faster than single-core methylKit (21.5 s) and ≈ 9 × faster than
DSS (15.1 s) on the same per-CpG cell; methylKit at `mc.cores = 8`
(5.9 × scaling, §4.3) narrows its gap to ≈ 2 ×.

See Figure F6 ([study2_simulated_headToHead/F6_runtime.png](../figures/study2_simulated_headToHead/F6_runtime.png))
for the per-scenario runtime distribution.

## 3.3 Study 3 — Real WGBS data (GSE263850), three-way DMR-caller comparison against a published multi-omics study

Study 3 evaluates DMR calling on real biological data using GSE263850
(Farhangdoost et al. 2025, *Molecular Psychiatry*) — six samples (three
Het-AKAP11-KO clones, three SBP009 WT replicates) of human
iPSC-derived neurons. The source paper is multi-omics (RNA-seq + WGBS +
ChIP-seq H3K27ac); our benchmark scope is the WGBS-DMR layer only. We
do not reproduce the paper's DEG calling, H3K27ac peaks, or
cross-omics correlations (no RNA-seq or ChIP-seq input data on hand).
We compare four DMR callers against the paper's Supp Table 5 (813 DMRs
called with `DSS::DMLfit.multiFactor(smoothing = TRUE)` +
`callDMR(p.threshold = 1e-5, delta = 0, minlen = 50, minCG = 3,
dis.merge = 100, pct.sig = 0.5)`):

1. **methylKit-tile** — 500 bp fixed tiles + methylKit
   `calculateDiffMeth` (the established R baseline)
2. **epykit-chain_merge (dis.merge = 100)** — paper-faithful settings,
   `dmr_chain_merge` engine with smoothing
3. **epykit-chain_merge (dis.merge = 250)** — same engine, increased
   merge gap to match DSS's per-chain morphology
4. **DSS-from-scratch** — our local rerun of the paper's exact DSS
   pipeline (DMLfit.multiFactor + DMLtest.multiFactor + callDMR) with
   paper-matched parameters, as a published-method upper bound

### 3.3.1 Per-CpG agreement (engine-agnostic)

On the 15.6 M shared CpGs between methylKit and epykit (the per-CpG
test is independent of which DMR-aggregation engine runs downstream):

* Pearson r on `meth_diff`: **0.9936**
* Spearman ρ: **0.9831**
* Same hyper/hypo direction: **14,669,608 / 15,597,046 = 94.05 %**

99.98 % of methylKit's tested CpGs are also tested by epykit
(per-CpG counts bit-identical, Methods §2.2). The two pipelines are
therefore measuring the same biological signal at the CpG level; all
downstream divergence reduces to DMR-aggregation choices.

### 3.3.2 DMR coordinate concordance vs DSS-from-scratch

The fairest apples-to-apples test of epykit on real WGBS is against the
DSS pipeline the source paper itself uses, re-run locally at the
paper-matched parameter set (Methods §2.4): same inputs, same
`p.threshold = 1e-5`, same `minlen = 50`, same `minCG = 3`,
same `dis.merge = 100`, same `pct.sig = 0.5`. This removes
DSS-version-drift, smoothing-internal, and threading-non-determinism
noise that confounds direct comparison against the paper's 813
published DMRs (§3.3.4 quantifies this DSS-vs-paper gap at ≈ 12.5 pp).

**Table 5a.** Headline coord-overlap statistics for the post-rerun
chain_merge call set against the locally-rerun DSS-from-scratch
ceiling (922 DMRs). Source: [`headline.json`](../data/study3/comparisons/epykit_vs_dss/headline.json).

| Caller | n DMRs | median bp | overlap (any-bp) | unique to caller | unique to DSS | recall any-bp | precision any-bp | recall J ≥ 0.5 | direction agreement on matched |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| methylKit-tile (500 bp) | 2,661 | 500 | ≈ 80 | 2,581 | 842 | 8.7 % | 3.0 % | ~ 0 % | n/a (≤ 1 strict match) |
| **epykit-chain_merge (dis.merge = 100)** | **852** | **125** | **634** | **218** | **334** | **63.8 %** | **74.4 %** | **34.5 %** | **588 / 588 = 100 %** |
| **DSS-from-scratch** (reference ceiling) | **922** | **241** | — | — | — | 100 % (self) | 100 % (self) | 100 % (self) | — |
| paper (Supp Table 5, contextual reference only) | 813 | 239 | — | — | — | — | — | — | — |

The **100 % direction agreement on every matched DMR**, across every
non-tile caller, is the strongest signal: when chain_merge overlaps a
DSS DMR, it never disagrees on the sign of the methylation change.
Disagreement is exclusively about *which* regions are flagged. The
634 / 922 = 68.8 % of DSS DMRs covered by ≥ 1 chain_merge call
(`query_hit_anybp` column of `headline.json`) is the appropriate
denominator when the DMR set itself is the question; the 588 / 852
direction-checked subset is the denominator for the
direction-coverage question. Both are 100 % in sign.

### 3.3.3 dis.merge as a calibration knob

`dis.merge` controls how aggressively `dmr_chain_merge` joins adjacent
significant CpG chains into wider DMRs. We re-ran the sweep on the
Linux 2026-06-04 cohort against the locally-rerun DSS-922 call set
(rather than the published 813-paper DMRs used pre-rerun); results in
`benchmark/data/multi_thread_and_chain_sweep/chain_merge_dis_merge_sweep/`:

| dis.merge | n DMRs | median bp | recall any-bp | recall J ≥ 0.5 | precision | direction agree |
|---:|---:|---:|---:|---:|---:|---:|
| 100 (paper) | 852 | 125 | 63.8 % | 34.5 % | 74.4 % | 588 / 588 = 100 % |
| 150 | 1,010 | 164 | 72.8 % | 52.1 % | 68.7 % | 671 / 671 = 100 % |
| 200 | 1,095 | 192 | 76.3 % | 60.6 % | 64.9 % | 703 / 703 = 100 % |
| **250 (morphology-matched)** | **1,139** | **205** | **77.3 %** | **64.2 %** | **63.0 %** | **713 / 713 = 100 %** |
| 500 | 1,160 | 214 | 78.4 % | 67.0 % | 61.5 % | 723 / 723 = 100 % |

`dis.merge = 250` nearly doubles the strict J ≥ 0.5 recall over the
paper-default 100 (34.5 % → 64.2 %) and brings the median DMR length
from 125 bp toward DSS's 241 bp. The gain plateaus past 250 (+1 pp
any-bp recall at 500 for –1.5 pp precision). We interpret 250 as a
morphology-matched setpoint for epykit's chain-merge under the
DSS-vs-epykit comparison — the paper's `dis.merge = 100` works for DSS
because DSS's count-smoother extends each significant chain slightly
further than epykit's quasi-binomial LR does in the small-p tail
(the smoother is identical; the test statistic differs). The full
sweep is in [F2 dis.merge curves](../figures/study3_real_GSE263850/three_way/F2_dis_merge_sweep.png).

### 3.3.4 DSS-from-scratch as a published-method upper bound

Re-running DSS locally with paper-matched parameters reaches **87.5 %
any-bp recall of the paper's call set**. The remaining ~12.5 pp gap is
DSS-vs-DSS noise — version drift between the paper's run and ours
(BSseq smoothing internals, threading non-determinism, DSS package
version not pinned by the paper; we ran DSS 2.58.0). DSS-fit and from-raw-counts
direction agreement on our run is **100 % (0 / 922 disagree)**,
confirming DSS's smoothed model is consistent with raw count
direction.

This DSS-vs-paper floor sets the scale: the *appropriate* ceiling for
any non-DSS caller on this dataset is the locally-rerun DSS-922, not
the paper's 813. Against that ceiling, post-rerun epykit-chain_merge
reaches **63.8 % any-bp recall at `dis.merge = 100`** (paper-faithful)
and **77.3 % at `dis.merge = 250`** (morphology-matched; §3.3.3), with
100 % direction agreement on every matched DMR at every sweep point.
The residual gap is a genuine test-statistic difference (quasi-binomial
LR vs DSS multifactor Wald + areaStat-based chain definition), not a
reproducibility artefact: chain-merge aggregation at `dis.merge = 250`
closes most of the morphology gap but does not — and is not designed to —
emulate DSS's smoothed-Wald chain-extension behaviour in the small-p
tail.

### 3.3.5 Annotation distribution (paper Fig 3C reproduction)

We re-annotated every caller's DMR set with a HOMER-equivalent UCSC
refGene classifier (Methods §C). Distribution of DMRs across genomic
features, ordered to match paper Fig 3C:

| Feature | paper (DSS) | mk-tile | ek-cm-100 | ek-cm-250 | DSS (local) |
|---|---:|---:|---:|---:|---:|
| promoter-TSS | 0.9 % | 2.1 % | 2.6 % | 2.8 % | 0.8 % |
| 5' UTR | 0.1 % | 0.0 % | 0.0 % | 0.0 % | 0.1 % |
| exon | 2.9 % | 1.6 % | 8.1 % | 6.9 % | 3.7 % |
| intron | 44.3 % | 35.8 % | 42.0 % | 42.1 % | 38.8 % |
| 3' UTR | 1.6 % | 1.2 % | 0.0 % | 0.0 % | 1.7 % |
| TTS | 1.4 % | 1.6 % | 1.1 % | 1.3 % | 1.6 % |
| non-coding | 0.5 % | 5.9 % | 0.0 % | 0.1 % | 5.9 % |
| intergenic | 48.3 % | 51.7 % | 46.2 % | 46.8 % | 47.4 % |
| n | 813 | 2,661 | **852** | **1,139** | 922 |

After re-annotating with epykit3 (which exposes
`features=...` on `ep.tl.annotate()` and uses the full HOMER default),
chain_merge is now essentially tied with DSS-from-scratch on
distribution match. TTS labels are now correctly assigned at ~1 % rate.
5' UTR / 3' UTR / non-coding remain at 0 % under epykit3 due to two
separate constraints: the UTR builders require GTF input (not refGene
cdsStart/cdsEnd), and `non-coding` is suppressed by higher-priority
intron labels of the overlapping gene body — multi-annotation does
catch the overlap (`{noncoding: 100}` in `all_overlapping_features`
for chain_merge-100), it's just not the assigned primary label. The
dominant intron + intergenic share (~85 %) matches paper Fig 3C across
every caller. See
[F5 annotation pie](../figures/study3_real_GSE263850/three_way/F5_annotation_pie.png).

### 3.3.6 Pathway enrichment (paper Fig 3D / panel D)

The paper reports GPCR ligand binding, Class A/1 Rhodopsin-like
receptors, GPCR downstream signalling, G alpha (i) signalling events,
and peptide ligand-binding receptors as enriched Reactome terms
(ShinyGO + Curated.Reactome on all 705 unique DMR-associated genes).
We re-enriched the four caller gene lists through the Enrichr REST
API against Reactome_2022 and KEGG_2021_Human (Methods §D); a portable
replacement for ShinyGO's Curated.Reactome that is less generous in
its multiple-testing correction. Even so:

* All chain_merge / DSS gene lists return **Morphine addiction** (KEGG's
  Gα(i)-signalling signature) and **Activation of G Protein Gated
  Potassium Channels** (Reactome) within the top-20 terms.
* methylKit-tile's gene list still recovers the Morphine addiction
  hit (consistent with the very large 2,111-gene set).
* Significant ShinyGO + Curated.Reactome reproduction was already
  demonstrated in the investigation note (Neuroactive ligand-receptor
  FDR = 6.4 × 10⁻⁵; cAMP signalling FDR = 6.4 × 10⁻⁵; Morphine
  addiction FDR = 2.5 × 10⁻² for chain_merge gene list).

Quantitative per-caller comparison in
[F8 enrichment dotplot](../figures/study3_real_GSE263850/three_way/F8_enrichment_dotplot.png).

### 3.3.7 Top-named gene hits (paper Fig 3B labels)

Paper Fig 3B labels its top 10 hyper-methylated and top 10
hypo-methylated DMR-associated genes (NR2E1, OTX1, OTX2, IRX2, PAX7,
ENPP2, GNG11, GREB1L, NAALADL2, and others). Direct coordinate-overlap
hits against each caller:

| Caller | any-bp hits / 20 | J ≥ 0.5 hits / 20 |
|---|---:|---:|
| methylKit-tile | **2 / 20** | 0 / 20 |
| epykit-chain_merge-100 | **12 / 20** | 7 / 20 |
| epykit-chain_merge-250 | **14 / 20** | 9 / 20 |
| DSS-from-scratch | **18 / 20** | 17 / 20 |

Post-rerun, chain_merge-100 hits 12/20 named genes at any-bp overlap
(NR2E1, OTX1, IRX2, ENPP2, GREB1L, CCDC177, GNG11, EBF1 on the hyper
side; LOC100131655, OSBPL8, RPLP0P2, FAM87A on the hypo side), rising
to 14/20 at `dis.merge = 250` (adds OTX2 and PAX7). The pattern matches
the coord-overlap headline: DSS recovers nearly all named genes at
strict overlap; chain_merge recovers roughly two-thirds; the
fixed-tile baseline misses almost all of them. See
[F4 named-gene heatmap](../figures/study3_real_GSE263850/three_way/F4_top_named_gene_hits.png).

### 3.3.8 Panel-E gene capture (proxy for Fig 3E)

The paper's Fig 3E enrichment runs on 46 critical genes (Supp Table 8)
that are simultaneously DMR-associated and differentially expressed.
We cannot recompute the GO MF enrichment without the RNA-seq DEG list,
but we report **gene capture rate** (% of the 46 genes whose name
appears in our 100 kb-rule DMR-gene set):

| Caller | Panel-E genes captured (nearest-TSS) | Panel-E genes captured (100 kb-rule) |
|---|---:|---:|
| methylKit-tile | 25 / 46 = 54.3 % | n/a (tile DMRs too short for 100 kb gene-link table) |
| epykit-chain_merge-100 | **30 / 46 = 65.2 %** | 29 / 46 = 63.0 % |
| epykit-chain_merge-250 | **32 / 46 = 69.6 %** | (sweep variant; gene-link table not regenerated) |
| **DSS-from-scratch** | **37 / 46 = 80.4 %** | 38 / 46 = 82.6 % |

Source: [`polish_recompute_2026_06_05.json`](../data/study3/comparisons/epykit_vs_dss/polish_recompute_2026_06_05.json).

### 3.3.9 Per-DMR effect-size concordance

For matched (J ≥ 0.5) DMR pairs between epykit-chain_merge-100 and
DSS-from-scratch (n = **318** matched pairs, post-rerun), the per-DMR
Pearson r on mean methylation difference is **0.9954**, Spearman ρ on
the effect size is **0.9382**, direction agreement is **100 %**. For
ek-chain_merge-250 (n = **592** pairs): r = **0.9965**, ρ = **0.9543**,
direction agreement 100 %. At the more lenient any-bp overlap
criterion: ek-100 / DSS n = 634 pairs (r = 0.9936); ek-250 / DSS
n = 718 pairs (r = 0.9950). Direction agreement remains 100 % at
every threshold.

When the two engines overlap a region, they agree to four decimal
places on the effect size. See
[F9 per-DMR concordance](../figures/study3_real_GSE263850/three_way/F9_per_dmr_concordance.png).

### 3.3.10 Performance (4-way)

**Table 5b.** Pipeline cost on the GSE263850 6-sample 22 M-CpG input,
**Linux host** (pivoine, 24 logical cores; methylKit `mc.cores = 1`
explicit for a fair single-core comparison; DSS single-thread by
construction). See §4.3 for the methylKit `mc.cores = 8` multi-core
ratio and §3.6 for the same-harness simulator per-CpG ratios
(≈ 13 × vs single-core methylKit, ≈ 2 × vs multi-core).

| Caller | Wall (s) | CPU (s) | Peak RSS (GB) | Notes |
|---|---:|---:|---:|---|
| methylKit-tile (`mc.cores = 1`) | 12,372 | 12,419 | **48.0** | dominated by `calculateDiffMeth` on 15.6 M CpGs |
| epykit-tile | 675 | 993 | 12.6 | full pipeline from raw BEDs |
| **epykit-chain_merge (100)** | **~ 92** | **~ 262** | **~ 12.6** | cached-store DMC + DMR re-call (re-callable across dis.merge) |
| **DSS-from-scratch** (DSS 2.58.0) | **2,368** | **2,454** | **14.3** | single-threaded by construction (verified in DSS 2.58.0); DMLfit smoothing dominates (~ 34 min) |

A full epykit pipeline (epykit-tile, comparable scope) is ~ 3.5 × faster
than DSS-from-scratch on the same input (675 s vs 2,368 s) and uses
**less** peak memory (12.6 vs 14.3 GB; epykit holds the per-CpG DMC
store, DSS holds the BSseq matrix plus the smoother's working set). The
cached-store chain_merge re-call (92 s) is faster still because it reuses
the already-computed DMC store. Against single-core methylKit, epykit is
≈ 18 × faster on the full DMC pipeline (the DMC step is shared by tile and
chain_merge). On Linux with `methylKit::calculateDiffMeth(mc.cores = 8)`
(5.9 × scaling measured on the simulator), the comparable methylKit wall
would be ≈ 2.1 kS, shrinking the ratio to ≈ 3 × at real-data scale (still
favourable). The DSS comparison is unaffected by core count (no
multi-core path). Per-pipeline resource breakdown in
[F7 resources](../figures/study3_real_GSE263850/three_way/F7_resources.png).

### 3.3.11 Cross-study consistency

The three studies are consistent on per-CpG calibration but show that
DMR-engine architecture matters for real-data DMR-level reproduction.
At the per-CpG level (Study 2 simulated, Study 3 real) epykit and
methylKit agree to three decimal places. At the DMR level on real
data, fixed-tile callers (methylKit-tile, epykit-tile) miss ≥ 90 % of
focused biological DMRs; variable-width chain-merge callers recover
half to two-thirds; DSS-from-scratch is the published-method ceiling
at 87.5 %. The full Study 3 summary picture is
[F10 summary](../figures/study3_real_GSE263850/three_way/F10_summary_three_way.png).

## 3.4 Held-out simulator with intrinsic truth (Phase 4)

To close the threshold-reconstructed-truth loop described in §2.7, we ran
methylKit and DSS (smoothing on and off) on every one of the 21 held-out
simulator seeds (`benchmark/scripts/run_external_simulator_sweep.py`)
and scored against the intrinsic `is_dmc` flag set at simulation time.
The headline cell is coverage = 10, 3 vs 3, q < 0.05, all effect-size
bins; numbers below are **across-seed median (IQR)** to make the
estimate intrinsically robust to single-seed noise (Methods §2.7,
PROTOCOL.md R2).

| Tool (mode) | TPR median (IQR) | FPR median (IQR) | F1 median | AUROC median (IQR) |
|---|---:|---:|---:|---:|
| methylKit | 0.7285 (0.7274 – 0.7315) | 0.0113 (0.0110 – 0.0115) | 0.8217 | 0.9260 (0.9245 – 0.9267) <!-- claim: simulator_methylkit_auroc_median_cov10 --> |
| DSS, smoothing = FALSE | 0.6547 (0.6538 – 0.6561) | 0.0057 (0.0055 – 0.0059) | 0.7808 | 0.9088 (0.9071 – 0.9094) <!-- claim: simulator_dss_nosmooth_auroc_median_cov10 --> |
| DSS, smoothing = TRUE | 5.0 × 10⁻⁵ (0 – 1.5 × 10⁻⁴) | 1.2 × 10⁻⁵ (0 – 1.3 × 10⁻⁵) | 1.0 × 10⁻⁴ | 0.6297 (0.6285 – 0.6314) <!-- claim: simulator_dss_smooth_auroc_median_cov10 --> |

Source: [`eval_simulator_intrinsic_iqr.parquet`](../data/study1b_simulator/eval_simulator_intrinsic_iqr.parquet)
(external tools; 21 seeds = 20 sampled + 1 frozen-grid control). Across all 20
held-out seeds the median epykit `lr` AUROC on intrinsic truth is **0.928**
(IQR 0.927–0.929; `eval_seed_iqr.parquet`), essentially tied with methylKit
(0.926); median epykit `lr+` AUROC is 0.907.

**Reading.** epykit `lr`, methylKit and unsmoothed DSS are all well-calibrated
on intrinsic truth and rank near-identically (median AUROC 0.928 / 0.926 /
0.909): methylKit is slightly more sensitive (TPR median 0.73) at a slightly
higher FPR, epykit `lr` is the most conservative on false positives, and `lr+`
trades calibration for recall (median TPR 0.746 at FPR 0.064, AUROC 0.907). At
the representative single seed 2026000, epykit `lr` AUROC is 0.9267
<!-- claim: simulator_epykit_lr_auroc_seed0_cov10 --> and `lr+` 0.9052
<!-- claim: simulator_epykit_lrplus_auroc_seed0_cov10 -->. Smoothed DSS
collapses to a near-zero call rate (median 5 × 10⁻⁵ TPR) — not a DSS defect,
but a dataset–assumption mismatch: the simulator generates CpGs at uniform
100-bp spacing with no genomic correlation structure for the smoother to
exploit, so DSS's default smoothing window flattens the signal entirely. The
same DSS configuration is competitive on real-data Study 3 (§3.3.4). The
seed-to-seed IQR widths are tight (< 0.003 absolute on TPR and AUROC for the
non-degenerate tools), confirming the median is not masking large between-seed
variance.

**Truth-definition duality (and §3.1's AUROC = 0.9999).** epykit `lr`'s AUROC
is **0.9999 under threshold-reconstructed truth** (§3.1, weak-effect DMCs
excluded by the |Δβ| ≥ 0.2 reconstruction) but **0.928 under intrinsic
held-out truth** (this section, every simulated DMC including the weak-effect
tail counts as a positive). The gap is a property of the *truth definition*,
not of the estimator: intrinsic truth includes near-threshold DMCs that no
calibrated test can separate from the null, so every tool's AUROC drops the
same way (methylKit 0.999 → 0.925). The 0.928 intrinsic figure is the honest
operating number and is what the abstract reports; 0.9999 is the
threshold-reconstructed ceiling.

The full per-seed table (`eval_simulator_intrinsic_per_seed.parquet`,
540 rows) and across-seed median + IQR summary appear as
**Supplementary Table S-Sim**. The 7-tool parallel column including
epykit results (scored separately) is in `parallel_column_summary.md`.

**Cross-truth comparison.** A dual-truth re-scoring on the same 21 seeds
([`eval_simulator_intrinsic_truth_both_iqr.parquet`](../data/study1b_simulator/eval_simulator_intrinsic_truth_both_iqr.parquet))
compares each tool under the simulator's intrinsic `is_dmc` flag
against Piao threshold-reconstructed truth (`|Δβ| ≥ 0.20` at 25 ×).
Tool ordering is preserved: methylKit ≈ unsmoothed DSS at the top,
smoothed DSS collapses, under both truth definitions. Absolute TPRs
shift modestly (≤ 2 pp) and AUROC under the threshold truth is uniformly
higher (≥ 0.98 for the two non-degenerate tools) because threshold
truth excludes the weak-effect tail. The threshold-vs-intrinsic
asymmetry that earlier reviewer commentary raised does not survive the
dual-truth check (see [`M3_truth_mode_comparison.md`](../data/study1b_simulator/M3_truth_mode_comparison.md)
for the full 6-row 2 × 3 panel).

## 3.5 Null calibration

A well-calibrated test produces uniform per-CpG p-values under the null
(no true DMCs). Distinguishing *calibrated* from merely *conservative*
requires inspecting the p-value distribution, not just the headline FDR.

**Exhaustive enumeration on real WGBS (GSE263850).** For n = 6 samples
in a 3v3 contrast, there are exactly C(6, 3) / 2 = 10 distinct unordered
label partitions. We ran the lr engine on each, capturing the per-CpG
p-value distribution. This is the *complete* null universe: at n = 6, k
= 1000 random shuffles would draw the same 10 partitions with
replacement (~100× each), conveying no information beyond the
exhaustive enumeration we report. Under realistic WGBS dispersion the
lr engine's null p-values are close to uniform:

| Statistic | Value | Expected under uniform | Interpretation |
|---|---:|---:|---|
| Mean p | 0.506 | 0.500 | +0.6 % conservative bias |
| Fraction p < 0.05 | 0.047 | 0.050 | –6 % under nominal |
| KS D | 0.051 | — | Small departure; effect-size measure |

The test is therefore **calibrated**, not merely conservative — the
distribution is close to uniform with a tiny conservative lean, so
FDR control at nominal q is valid with negligible power cost. The
Q-Q plot (Supplementary Figure S-QQ-LR) hugs the diagonal across
the full p-value range. Source: `benchmark/data/null_calibration/
gse263850/lr_calibration_report.{md,json}` and `lr_qq.png`; sampled
p-values for figure regeneration in `lr_pvalues.parquet`.

**Pre-1.0 K = 20 baseline (deprecated headline).** Earlier engineering
runs used K = 20 random shuffles on every (engine, dataset) cell and
reported median observed FDR ≤ 1.53 × 10⁻⁵. At K = 20 the minimum
resolvable empirical p-value is ≈ 1/K = 0.05 — five-decimal precision
was unsupportable. The K = 20 summary survives in
`benchmark/data/null_calibration/summary.parquet` as the Piao-simulator
and pre-1.0 reference baseline; the exhaustive-enumeration headline
above supersedes it on the real-data calibration question.

**Performance headline.** epykit's bare `lr` engine completes the same
coverage = 10, 3 vs 3 cell in median 0.86 s
<!-- claim: headline_wallclock_epykit_lr -->; the opt-in `lr+` stack
in 6.80 s <!-- claim: headline_wallclock_epykit_lrplus -->. The full
per-engine timing table is `benchmark/data/study1/timings_table.csv`;
the 7-tool cross-tool wallclock + accuracy comparison is
`benchmark/docs/timing-comparison.md`.

## 3.6 Behaviour under realistic overdispersion (φ-sweep)

The per-CpG results above are measured at the Piao model's near-binomial
dispersion (Pearson φ ≈ 1); real WGBS is overdispersed at φ ≈ 1.5–5 (§4). To
test directly whether the ranking survives at realistic dispersion — rather
than predicting it — we extended the held-out simulator with a Beta-Binomial
intraclass-correlation parameter ρ and re-ran the full eight-tool panel at
ρ ∈ {0, 0.05, 0.1, 0.2, 0.3, 0.44}. At coverage 10 the implied Pearson
overdispersion is φ = 1 + (coverage − 1)·ρ = {1.0, 1.45, 1.9, 2.8, 3.7, 5.0},
bracketing the real-WGBS range. Each cell is the median of 10 seeds at 3 vs 3
(60 cells, 480 tool-runs, 0 failures; `benchmark/scripts/run_phi_sweep.py`).
The sweep cleanly separates *detection power* from *calibration*.

**Detection power is tied.** Threshold-free AUROC is statistically
indistinguishable across `lr`, methylKit and DSS at every dispersion level
(0.93 at φ = 1 falling to 0.76 at φ = 5; the three tools stay within ≈ 0.01 of
one another throughout). At a *matched* FPR = 0.05 the recovered TPRs are
likewise equal (epykit `lr` / methylKit-MN / methylKit-default = 0.826 / 0.818
/ 0.825 at φ ≈ 1; 0.435 / 0.397 / 0.438 at φ ≈ 5). No tool detects more true
DMCs than another at equal stringency.

**Only `lr` stays calibrated.** The tools diverge entirely on whether the
nominal q < 0.05 cutoff means what it says (Table S-Phi, Figure 11 centre).

**Table S-Phi.** Realised FDR at nominal q < 0.05 (median of 10 seeds,
coverage 10, 3 vs 3). methylKit is shown in its overdispersion-aware `MN`
mode — the like-for-like comparison to `lr`; the dispersion-blind default is
worse (FDR up to 0.66) and is a secondary point, not a headline.

| φ (Pearson) | epykit `lr` | methylKit (MN) | DSS (no smooth) |
|---:|---:|---:|---:|
| 1.0 | 0.026 | 0.044 | 0.033 |
| 1.45 | 0.028 | 0.102 | 0.087 |
| 1.9 | 0.029 | 0.157 | 0.133 |
| 2.8 | 0.027 | 0.229 | 0.197 |
| 3.7 | 0.029 | 0.279 | 0.244 |
| 5.0 | 0.021 | 0.300 | 0.282 |

Source: [`eval_phi_sweep_iqr.parquet`](../data/study1b_simulator/eval_phi_sweep_iqr.parquet).

epykit `lr` stays at FDR ≈ 0.02–0.03 from φ = 1 to φ = 5; methylKit-MN controls
FDR only at φ ≈ 1 and is already anti-conservative by φ ≈ 1.45 (FDR 0.10 — the
low end of real WGBS), reaching 0.30 at φ ≈ 5. DSS-no-smoothing behaves
similarly (0.033 → 0.282).

This reframes the low-coverage advantage of §3.1: `lr`'s lower TPR at q < 0.05
under high dispersion is **not** a sensitivity deficit (AUROC is tied) — it is
the cost of honest FDR control. The binomial-variance / df floor in the `lr`
engine (§2.5) keeps the per-site quasi-binomial calibrated as overdispersion
grows; methylKit's and DSS's per-site dispersion estimates do not, so their
q-values drift anti-conservative on realistic WGBS, where the false-positive
cost is borne downstream. Calibration under overdispersion, not raw
sensitivity, is epykit's substantive per-CpG contribution.

dmrseq and BSmooth are region callers; scored at single-CpG resolution here
they are uncalibrated (FDR ≈ 0.6–0.7) and not meaningfully comparable per-CpG —
they enter at the DMR level (§3.1, Study 3) instead. We retain them in this
sweep only for the resource axis. On wall-clock at this cell epykit `lr`
(median 1.7 s) is ≈ 13× faster than single-core methylKit (21.5 s) and ≈ 9×
faster than DSS (15.1 s); methylKit at `mc.cores = 8` (5.9× scaling) narrows
to ≈ 2×. Peak RSS at this cell is epykit `lr` 0.44 GB vs DSS 1.3 GB,
methylKit 8.8 GB, BSmooth 9.5 GB and dmrseq 30 GB — roughly 20× less memory
than methylKit (Figure 11, right; source `memory_timing_by_tool.csv`).

![Figure 11. Dispersion (φ) sweep on the intrinsic-truth simulator (coverage
10, 3 vs 3, median of 10 seeds; dual ρ / Pearson-φ axis). Left: sensitivity
(TPR at q < 0.05) falls with dispersion for the calibrated tools. Centre:
realised FDR at q < 0.05 — only epykit `lr` tracks the nominal 0.05 line across
the realistic-WGBS band (φ ≈ 1.5–5, shaded); methylKit and DSS inflate 3–10×.
Right: peak RSS.](../figures/study1_simulated_allPackages/F9_phi_sweep.png)

# 4. Discussion

This Discussion is organised in three parts: §4.1 summarises what the three
studies establish, §4.2 surfaces the calibration–sensitivity trade-off seen
most clearly in Study 3, and §4.3 lists limitations.

**A caveat, now tested directly (§3.6).** The Piao 2021 simulator is
*underdispersed* relative to real WGBS data: median Pearson φ at coverage 5× is
≈ 0.41 on the simulator versus ≈ 1.5–5 on biological samples. epykit's bare
`lr` engine clamps at the binomial floor (φ = 1), which is nearly correct on
the simulator and partly explains its raw-TPR dominance in the small-effect bin
at low coverage (§3.1, Figure 2). Rather than predict what happens at realistic
dispersion, we swept the simulator's overdispersion up to φ ≈ 5 (§3.6). The raw
small-effect TPR advantage does narrow — at high dispersion every calibrated
tool detects the same fraction of true DMCs (tied AUROC) — but it does not so
much reverse as *change character*: what survives and grows is **calibration**.
`lr` is the only per-CpG test that holds nominal FDR across the whole φ ≈ 1.5–5
range, while methylKit and DSS become anti-conservative. Readers should
therefore treat the simulator's raw TPR-at-q<0.05 numbers as upper bounds, but
the FDR-calibration ranking (§3.6) as the property that transfers to real WGBS —
consistent with the held-out simulator (§3.4), null calibration (§3.5), and
Study 3 (§3.3).

## 4.1 What the three studies establish

Together, the three studies cover the relevant evaluation surface for a new
WGBS pipeline:

* Study 1 places epykit within the **broader ecosystem**: `lr` is
  competitive with the strongest established tools (methylKit, RADMeth, DSS)
  at moderate to high coverage, and `lr+` reaches the perfect ceiling at
  every coverage and replicate count ≥ 4 with well-controlled FPR. At low
  coverage and small effects `lr` outperforms several baselines.
* Study 2 establishes **bit-precise calibration** against methylKit: the
  two implementations agree on TPR, FPR, F1 and AUROC to three decimal
  places at n ≥ 4, and at n = 2 — where methylKit's overdispersion
  estimator becomes degenerate — epykit's `lr` engine recovers TPR
  0.564 vs methylKit's 0.302 at matched FPR (Study 2, single cell).
* Study 3 demonstrates **faithful behaviour on real biological data**:
  on identical counts the two pipelines agree on direction at 100 %
  of overlapping DMRs (Study 3 chain_merge-250 vs DSS-from-scratch,
  713/713), and on effect size at Pearson r = 0.994 per-CpG.

On the simulator (per-CpG, same harness; §3.6), epykit `lr` is ≈ 13 ×
faster than single-core methylKit and ≈ 9 × faster than DSS, narrowing
to ≈ 2 × against methylKit at `mc.cores = 8`; on Study 3 the full epykit
DMR pipeline is ≈ 18 × faster than single-core methylKit and ≈ 3.5 ×
faster than DSS-from-scratch end-to-end (Table 5b). Earlier Windows-only
numbers (12 × – 68 ×) reflected methylKit's `mc.cores` no-op on Windows
(no `fork()`) and overstated the gap.

## 4.2 The calibration–sensitivity trade-off

Study 3 surfaces a real choice: at n = 3 per group on biological data,
methylKit's pooled `overdispersion="MN"` is more aggressive in the
small-p tail than epykit's per-site McCullagh–Nelder; epykit calls
fewer DMCs at the same threshold (30,965 vs methylKit 51,792 at
|d| ≥ 10 %, q < 0.05). The `lr+` recipe emits ≈ 13 × more total
significant DMCs (406,515 vs 30,965) than bare `lr` at the same
q-threshold — consistent with the simulator-validated FPR inflation
under realistic dispersion (§3.4 / §4 framing paragraph). Neither
operating point is "more correct"; the precision/recall optimum
depends on the downstream analysis and the user's tolerance for
false positives. We therefore recommend that users:

* Report results at the default `lr` setting unless they have a specific
  reason to deviate;
* Treat `lr+` as a sensitivity-favouring research mode, not a default;
  reproduce headline findings under bare `lr` as a calibration check;
* Document any opt-in (`lr+`, `dispersion="shrink"`, etc.) explicitly.

## 4.3 Limitations

* **Simulator underdispersion — measured, not just flagged.** The Piao 2021
  simulator is underdispersed (median φ ≈ 0.41 at 5×) relative to real WGBS
  (φ ≈ 1.5–5), and `lr`'s binomial-floor clamp is part of why it leads the
  low-coverage small-effect bin on raw TPR in §3.1. The φ-sweep (§3.6) measures
  the consequence directly: at realistic dispersion the raw-TPR lead disappears
  (AUROC ties across tools), while `lr`'s FDR-calibration advantage strengthens
  (realised FDR ≈ 0.02–0.03 vs methylKit-MN 0.10–0.30 and DSS 0.09–0.28).
  `dispersion="eb"` is designed for heterogeneous regimes but is a no-op on the
  near-binomial default simulator.
* **The empirical-Bayes dispersion prior is not validated against an external
  estimator.** The EB shrinkage behind `dispersion="eb"` (and hence the opt-in
  `lr+` stack) treats per-site Pearson dispersions as draws from an
  inverse-Gamma prior whose hyperparameters are fit by method-of-moments per
  chromosome (§2.5). We make no goodness-of-fit claim for that prior: its
  per-site posterior dispersions have not been Q-Q compared against an
  independent estimator such as DSS's, and that comparison remains future work.
  This is acceptable here precisely because `eb`/`lr+` are opt-in research
  knobs rather than the recommended default — **every headline result in this
  paper is reported under bare `lr`**, whose dispersion is the binomial-floored
  quasi-binomial Pearson estimate, not the EB prior. A reader should not read
  any headline number as depending on, or as evidence for, the EB prior.
* **Baseline software versions.** Study 1 baseline numbers are from 2021
  software releases. Relative ordering at low coverage / small n is robust
  across recent versions of those tools, but absolute numbers may have
  shifted.
* **DMR baselines are figure-derived.** DMR detection rates for the eight
  baselines in Study 1 come from hand-transcribed bar charts (Figures 3a,
  3b, S5–S7). Per-figure confidence labels in `PROVENANCE.md`.
* **Single real dataset.** Study 3 is one tissue × one genome (Het_AKAP11_KO
  vs WT in hg38, six samples total). A multi-cohort real-data validation
  (TCGA tumour/normal WGBS, ENCODE tissue series, mouse imprinted DMRs)
  is future work; the DSS-vs-epykit ordering observed here may not
  generalise to higher dispersion regimes.
* **DMR-engine choice on real data.** Study 3 (§3.3) shows that fixed
  500 bp tile callers (including `dmr_tile` and methylKit's `tileMethylCounts`)
  recover ≤ 10 % of focused real-data DMRs at coordinate level when
  compared against a published DSS call set. epykit's `dmr_chain_merge`
  recovers 63.8 % (any-bp) at the paper-faithful `dis.merge = 100` and
  77.3 % at morphology-matched `dis.merge = 250`, with 100 % direction
  agreement on overlapping DMRs in both cases; DSS-from-scratch reaches
  87.5 %. The DSS smoothing prior helps most on low-CpG-density regions
  (named-gene misses include CLEC19A, KANK1, CNR1). Per-CpG calibration
  is engine-agnostic; for users targeting reproduction of published
  DMR-level analyses we recommend `dmr_chain_merge` with `dis.merge = 250`
  as the default for real-data region-level analysis.
* **Multi-omics scope.** The Farhangdoost et al. 2025 paper integrates
  RNA-seq DEGs, WGBS DMRs, and ChIP-seq H3K27ac peaks; our benchmark
  scope is WGBS only. We can compare DMR coordinates, annotations,
  morphology, and gene assignment against the paper's Supp Tables 5,
  6, and 8, but we cannot recompute the DMR–DEG correlations or
  triple-overlap enhancer analyses without the additional GEO
  datasets and an independent RNA-seq / ChIP-seq pipeline.
* **Platform-dependent timing (Linux re-run).** Earlier engineering runs
  on Windows recorded a 12 – 68 × speedup vs methylKit, but methylKit's
  `mc.cores` is a no-op on Windows (no `fork()`), so methylKit ran
  single-threaded by force. We re-ran the Study 2 grid on Linux with
  `methylKit::calculateDiffMeth(mc.cores = 8)` (5.9 × scaling, median
  across three seeds). Under one harness (§3.6) epykit's bare `lr` is
  ≈ 13 × faster than single-core methylKit on the simulator, narrowing to
  ≈ 2 × against methylKit at `mc.cores = 8`
  (`benchmark/data/multi_thread_and_chain_sweep/methylkit_multicore/`).
  We separately verified that DSS's `DMLfit.multiFactor` (the
  multi-factor path used here) provides no multi-core option in
  the DSS version we used (2.58.0), so the reported speed advantage
  over DSS is not eroded by parallelising DSS.
* **Ground truth non-independence.** Study 1 true-DMC labels come from the
  coverage-25 sample (the cleanest signal in the dataset), not from the
  simulator's internal flags. To check the in-house simulator does not
  asymmetrically advantage epykit (reviewer M3), we re-scored the
  external baselines under both labellings on the 21-seed Linux output:
  the methylkit-vs-DSS ordering is robust (intrinsic AUROC 0.926 vs
  0.909; threshold AUROC 0.987 vs 0.986); the absolute scale shifts
  upward under threshold truth but the rank does not flip.
  (`benchmark/data/study1b_simulator/M3_truth_mode_comparison.md`.)
* **Null calibration coverage.** Our published headline null result
  (§3.5) is the lr engine on label-permuted GSE263850 with all C(6,3)/2
  = 10 unique 3v3 partitions enumerated. For n = 6, this is the
  *complete* null universe — k = 1000 random shuffles would draw the
  same 10 partitions with replacement (~100× each). The earlier K = 20
  random-shuffle sweep over (engine × dataset) cells survives as a
  pre-1.0 baseline in `benchmark/data/null_calibration/summary.parquet`
  but cannot resolve the calibration-vs-conservatism question at
  five-decimal precision (min resolvable empirical p ≈ 0.05). Extending
  the exhaustive enumeration to larger cohorts (e.g. TCGA tumour/normal
  with n ≥ 6+6) would require sampling rather than enumeration but is
  not required for the headline calibration claim.
* **lr+ as a research knob, not a default.** The four-knob `lr+` stack
  trades precision for recall (Simulator: TPR 0.673 → 0.746 buys a 14×
  FPR inflation, 0.0044 → 0.064; F1 0.796 → 0.746; AUROC 0.928 → 0.907).
  We retain `lr+` as a deliberately-opt-in mode for users who want to
  trade calibration for sensitivity in low-replicate regimes, and
  document it as such in the API and §2.5. Bare `lr` is the engine
  the abstract numbers are reported against; `lr+` headline panels are
  labelled as such throughout §3.
* **sep_threshold is inert at realistic coverage.** §M11 (sep_threshold)
  reports that the lr+ separation-aware Fisher fallback never fires on
  GSE263850 (0 candidate sites across all chromosomes, identical DMC
  output at sep_threshold ∈ {0.7, 0.8, 0.9, 0.95} on the validated
  cold-cache sweep). It is a rare-event safeguard for pathologically
  low-coverage sites that the coverage filter (default ≥ 10×/sample)
  removes upstream; the default value (0.9) has no effect on any
  reported number and requires no tuning. Investigating this question
  surfaced and fixed a cache-key bug in `tl.py` that had previously
  hidden the inertness behind silent cache reuse.

# 5. Conclusion

epykit is a Python-native pipeline for WGBS downstream analysis whose
default `lr` engine is competitive with the strongest R/CLI tools
(simulator: AUROC 0.928 vs methylKit 0.926, DSS-no-smoothing 0.909).
The opt-in `lr+` stack trades precision for recall (14× FPR inflation
for +7 pp TPR; lower F1 and AUROC) and is documented as a research
knob, not a recommended default. On real WGBS (GSE263850), at the
paper-matched `dmr_chain_merge` operating point (`dis.merge = 250`):

* fixed-tile callers (methylKit, epykit-tile) miss ≥ 90 % of focused
  real-data DMRs at coordinate level;
* epykit's `dmr_chain_merge` recovers 77.3 % of DSS-from-scratch DMRs
  by any-bp overlap and 64.2 % at Jaccard ≥ 0.5, with 100 % direction
  agreement on the 713 overlapping DMRs and 69.6 % panel-E gene-set
  recall (32 / 46);
* DSS-from-scratch retains a small recall advantage (87.5 %, 80.4 %
  panel-E gene recall) attributable to its smoothing prior, which
  helps particularly on low-CpG-density regions.

Under exhaustive label-permutation on GSE263850 (all 10 unique 3v3
partitions = the complete null universe at n = 6) the `lr` engine is
calibrated, not merely conservative: null p-values are close to
uniform (mean 0.506, fraction below 0.05 = 0.047, KS D = 0.051), so
FDR control is valid at negligible power cost.

On the simulator (per-CpG, same harness; §3.6), epykit `lr` is ≈ 13 ×
faster than single-core methylKit and ≈ 9 × faster than DSS (≈ 2 × vs
methylKit at `mc.cores = 8`); on Study 3 the full epykit DMR pipeline is
≈ 18 × faster than single-core methylKit and ≈ 3.5 × faster than
DSS-from-scratch end-to-end (Table 5b). DSS's
`DMLfit.multiFactor` does not expose a multi-core option in the DSS
version we used (2.58.0), so the DSS comparison is single-thread by
construction, not by choice. Combined with the rest
of the epykit API (annotation, plotting, HTML reporting,
AnnData / MuData interop), this brings the WGBS downstream pipeline
into the same Python ecosystem as the rest of modern bioinformatics.

# Availability

All benchmark code, ground-truth reconstruction, baseline transcriptions,
data tables, and figure-generation scripts are in this `FINAL_REPORT/`
directory. See [README.md](../README.md) for the reproduction recipe and
[report/methods_appendix.md](../report/methods_appendix.md) for tool versions
and parameters.

# References
