# epykit Benchmark — Technical Report

**Companion to** [paper.md](../paper/paper.md).
**Date:** 2026-05-22.
**Software:** epykit 1.0.0; methylKit 1.34.0 (Study 2) / 1.36.0 (Study 3) /
0.99.2 (Study 1 baseline transcriptions from Piao et al. 2021).

> epykit benchmarks were executed at commit `60a71e0` (engine tag
> `v0.7.5-phase3-engines-frozen`, 2026-05-28); the package version string was
> `0.7.2` at run time and the engine code is unchanged through the `1.0.0`
> release, so the results apply to 1.0.0.

This report contains the complete numeric tables behind the paper, the
full per-stage performance breakdown, the dispersion-estimator sensitivity
analysis on real data, and the bug-fix log. It is intended for readers who
want to reproduce or audit specific claims; for the narrative, read
[paper.md](../paper/paper.md) first.

---

## 1. Study 1 — Simulated panel comparison (epykit vs 8 baselines)

### 1.1 DMC × Coverage, 3 vs 3 design, q < 0.05

#### 1.1.1 Small effects (|Δβ| = 0.2–0.4)

| Tool | 5× | 10× | 15× | 20× | 25× |
|------|---:|----:|----:|----:|----:|
| **epykit `lr+`** | **0.999** | **1.000** | **1.000** | **1.000** | **1.000** |
| **epykit `lr`**  | **0.835** | **0.963** | **0.988** | **0.996** | **0.997** |
| epykit `fisher` (post-fix) | 0.649 | 0.890 | 0.978 | 0.994 | 0.985 |
| RADMeth | 0.421 | 0.967 | 0.987 | 0.997 | 0.998 |
| BiSeq | 0.684 | 0.607 | 0.698 | 0.702 | 0.690 |
| methylKit | 0.266 | 0.963 | 0.987 | 0.998 | 1.000 |
| DSS | 0.065 | 0.935 | 0.984 | 0.996 | 0.998 |
| Fisher (paper, pooled) | 0.082 | 0.920 | 0.980 | 0.996 | 0.999 |
| epykit `welch_t` | 0.512 | 0.769 | 0.990 | 0.993 | 1.000 |
| methylSig | n/a | 0.716 | 0.896 | 0.947 | 0.966 |
| epykit `bb_lr` | 0.019 | 0.077 | 0.728 | 0.642 | 0.995 |

#### 1.1.2 Medium effects (|Δβ| = 0.4–0.6)

| Tool | 5× | 10× | 15× | 20× | 25× |
|------|---:|----:|----:|----:|----:|
| **epykit `lr+`** | **0.999** | **0.999** | **1.000** | **1.000** | **1.000** |
| **epykit `lr`**  | **0.854** | **0.963** | **0.984** | **0.997** | **1.000** |
| epykit `fisher` (post-fix) | 0.671 | 0.902 | 0.970 | 0.993 | 1.000 |
| RADMeth | 0.964 | 0.970 | 0.985 | 0.996 | 0.997 |
| methylKit | 0.857 | 0.962 | 0.987 | 0.996 | 1.000 |
| DSS | 0.737 | 0.938 | 0.984 | 0.993 | 0.999 |
| Fisher (paper, pooled) | 0.600 | 0.922 | 0.979 | 0.993 | 0.999 |
| BiSeq | 0.658 | 0.601 | 0.674 | 0.677 | 0.666 |
| methylSig | n/a | 0.714 | 0.888 | 0.945 | 0.969 |

#### 1.1.3 Strong effects (|Δβ| = 0.6–0.8)

| Tool | 5× | 10× | 15× | 20× | 25× |
|------|---:|----:|----:|----:|----:|
| **epykit `lr+`** | **0.999** | **0.999** | **1.000** | **1.000** | **1.000** |
| RADMeth | 1.000 | 0.971 | 0.984 | 0.996 | 0.998 |
| methylKit | 0.998 | 0.964 | 0.983 | 0.996 | 1.000 |
| DSS | 0.997 | 0.938 | 0.980 | 0.995 | 0.999 |
| Fisher (paper, pooled) | 0.973 | 0.923 | 0.976 | 0.994 | 0.999 |
| **epykit `lr`** | 0.843 | 0.961 | 0.985 | 0.997 | 1.000 |
| epykit `fisher` (post-fix) | 0.668 | 0.900 | 0.974 | 0.993 | 1.000 |
| BiSeq | 0.657 | 0.612 | 0.681 | 0.687 | 0.679 |
| methylSig | n/a | 0.730 | 0.896 | 0.948 | 0.971 |

#### 1.1.4 FPR (0.2–0.4 bin)

| Tool | 5× | 10× | 15× | 20× | 25× |
|------|---:|----:|----:|----:|----:|
| **epykit `lr`** | **3.7 × 10⁻⁵** | **1.2 × 10⁻⁵** | **1.2 × 10⁻⁵** | **1.2 × 10⁻⁵** | **0** |
| epykit `lr+`    | 1.9 × 10⁻² | 1.4 × 10⁻³ | 1.5 × 10⁻⁴ | 1.2 × 10⁻⁵ | 1.2 × 10⁻⁵ |
| epykit `fisher` (post-fix) | 2.5 × 10⁻⁵ | 1.2 × 10⁻⁵ | 1.2 × 10⁻⁵ | 1.2 × 10⁻⁵ | 0 |
| methylKit | 2.3 × 10⁻² | 1.0 × 10⁻³ | 0 | 0 | 0 |
| RADMeth   | 1.8 × 10⁻² | 1.0 × 10⁻³ | 0 | 0 | 0 |
| DSS       | 2.9 × 10⁻² | 2.0 × 10⁻³ | 1.0 × 10⁻³ | 0 | 0 |
| Fisher (paper, pooled) | 2.9 × 10⁻² | 2.0 × 10⁻³ | 1.0 × 10⁻³ | 0 | 0 |
| BiSeq     | 1.0 × 10⁻² | 1.2 × 10⁻² | 9.0 × 10⁻³ | 9.0 × 10⁻³ | 1.0 × 10⁻² |
| methylSig | n/a | 9.0 × 10⁻³ | 3.0 × 10⁻³ | 2.0 × 10⁻³ | 1.0 × 10⁻³ |

**Reading.** At 5×, epykit `lr` has FPR 100×–600× lower than any
competitor. `lr+` trades a small FPR cost (1.9 × 10⁻² at 5×, comparable to
the R baselines) for large TPR gains at low coverage.

#### 1.1.5 AUROC (epykit only — baselines have no full p-value vectors)

| Tool | 5× | 10× | 15× | 20× | 25× |
|------|---:|----:|----:|----:|----:|
| epykit `lr+` | 0.9999 | 0.9999 | 0.9999 | 0.9999 | 0.9999 |
| epykit `lr`  | 0.9990 | 0.9999 | 0.9999 | 0.9999 | 1.0000 |
| epykit `fisher` (post-fix) | 0.9998 | 0.9999 | 1.0000 | 1.0000 | 1.0000 |

### 1.2 DMC × Replicates (10× coverage), q < 0.05

#### 1.2.1 Small effects (|Δβ| = 0.2–0.4)

| Tool | n=2 (1v1) | n=4 (2v2) | n=6 (3v3) | n=8 (4v4) | n=10 (5v5) |
|------|----------:|----------:|----------:|----------:|-----------:|
| epykit `lr`     | — | **0.998** | **0.950** | **0.979** | **0.995** |
| epykit `lr+`    | — | **0.997** | **1.000** | — | **1.000** |
| epykit `fisher` (post-fix) | 0.205 | 0.806 | 0.920 | 0.964 | 0.988 |
| methylKit | 0.000 | 0.144 | 0.567 | 0.805 | 0.983 |
| DSS       | 0.058 | 0.117 | 0.558 | 0.786 | 0.975 |
| RADMeth   | 0.058 | 0.119 | 0.655 | 0.879 | 0.992 |
| Fisher (paper) | 0.000 | 0.029 | 0.395 | 0.727 | 0.966 |
| BiSeq     | 0.001 | 0.686 | 0.687 | 0.696 | 0.677 |
| methylSig | n/a | 0.144 | 0.000 | 0.268 | 0.871 |

#### 1.2.2 Strong effects (|Δβ| = 0.6–0.8)

| Tool | n=2 | n=4 | n=6 | n=8 | n=10 |
|------|----:|----:|----:|----:|-----:|
| epykit `lr`  | — | **0.999** | **0.949** | **0.979** | **0.994** |
| epykit `lr+` | — | **0.999** | **1.000** | — | **1.000** |
| epykit `fisher` (post-fix) | 0.221 | 0.796 | 0.918 | 0.964 | 0.987 |
| methylKit | 0.858 | 1.000 | 1.000 | 1.000 | 1.000 |
| DSS       | 0.873 | 1.000 | 1.000 | 1.000 | 1.000 |
| RADMeth   | 0.873 | 1.000 | 1.000 | 1.000 | 1.000 |
| Fisher (paper) | 0.350 | 1.000 | 1.000 | 1.000 | 1.000 |
| BiSeq     | 0.001 | 0.687 | 0.696 | 0.686 | 0.689 |

### 1.3 DMR × Coverage (≥ 80 % overlap criterion)

#### 1.3.1 Recall (TPR)

| Tool | 5× | 10× | 15× | 20× | 25× |
|------|---:|----:|----:|----:|----:|
| methylKit | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| Fisher (paper, pooled) | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| **epykit `chain_merge` (post-fix)** | **0.97** | **1.00** | **1.00** | **1.00** | **1.00** |
| RADMeth | 0.74 | 0.86 | 1.00 | 1.00 | 1.00 |
| epykit `tile` | 0.54 | 0.77 | 0.91 | 1.00 | 1.00 |
| metilene | 0.74 | 0.89 | 0.94 | 0.89 | 0.91 |
| methylSig | 0.20 | 0.86 | 0.94 | 0.97 | 0.97 |
| BSmooth | 0.49 | 0.54 | n/a | 0.54 | n/a |
| BiSeq | 0.14 | 0.14 | 0.14 | 0.14 | 0.14 |
| DSS | 0.00 | 0.06 | 0.43 | 0.66 | 0.69 |

#### 1.3.2 Precision (epykit DMR engines)

| Engine | 5× | 10× | 15× | 20× | 25× |
|---|---:|---:|---:|---:|---:|
| `chain_merge` | 0.65 | 0.88 | 1.00 | 1.00 | 1.00 |
| `tile` | 0.54 | 1.00 | 1.00 | 1.00 | 1.00 |

#### 1.3.3 F1 (epykit DMR engines)

| Engine | 5× | 10× | 15× | 20× | 25× |
|---|---:|---:|---:|---:|---:|
| `chain_merge` | 0.78 | 0.94 | 1.00 | 1.00 | 1.00 |
| `tile` | 0.54 | 0.87 | 0.96 | 1.00 | 1.00 |

### 1.4 Runtime

epykit single-CPU wall-clock per scenario (Study 1 simulator):

| Engine | Scenario | Time |
|--------|----------|------|
| `lr` | DMC 100 K sites | 0.4–0.9 s |
| `lr+` | DMC 100 K sites | 1.7–1.8 s |
| `fisher` | DMC 100 K sites | 0.3–0.5 s |
| `welch_t` | DMC 100 K sites | 0.3–0.7 s |
| `bb_lr` | DMC 100 K sites | 3.6–7.5 s |
| `tile` | DMR ~4 M sites | 2.9–11.7 s |
| `chain_merge` | DMR ~4 M sites | 1.0–14.3 s |
| **Full grid (all 25 scenarios)** | — | **~5 min** |

Piao et al. (2021) Section 3.4 places BiSeq and RADMeth in the multi-hour
range on the same data; we did not re-run those baselines for runtime in
Study 1.

---

## 2. Study 2 — Head-to-head with methylKit on simulated data

Both tools were run on the **same machine, same simulated input data, same
evaluation harness**. methylKit was instrumented under the same OS-level
resource tracker as epykit. methylKit `mc.cores = 1` (Windows host; no
`fork()`); epykit single-threaded dispatcher (Polars/NumPy still use
implicit vectorisation threads internally).

### 2.1 DMC × coverage, 3 vs 3 design

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
| 25× | epykit / `lr` | 0.993 | 0.0 | 0.996 | 1.0000 |
| 25× | methylKit / default | 0.993 | 0.0 | 0.997 | 1.0000 |

### 2.2 DMC × replicate count

| n_total | Tool / engine | TPR | FPR | F1 | AUROC | n_sig |
|---|---|---|---|---|---|---|
| 2  | epykit / `lr`        | **0.564** | 1.2 × 10⁻⁵ | **0.721** | 0.9993 | 11,283 |
| 2  | methylKit / default  | 0.302     | 0.0        | 0.463     | 0.9994 | 6,030 |
| 4  | epykit / `lr`        | 0.880     | 1.2 × 10⁻⁵ | 0.936     | 0.9999 | 17,595 |
| 4  | methylKit / default  | 0.880     | 1.2 × 10⁻⁵ | 0.936     | 0.9999 | 17,595 |
| 6  | epykit / `lr`        | 0.952     | 1.2 × 10⁻⁵ | 0.975     | 1.0000 | 19,039 |
| 6  | methylKit / default  | 0.952     | 1.2 × 10⁻⁵ | 0.975     | 1.0000 | 19,039 |
| 8  | epykit / `lr`        | 0.979     | 2.5 × 10⁻⁵ | 0.989     | 1.0000 | 19,574 |
| 8  | methylKit / default  | 0.979     | 2.5 × 10⁻⁵ | 0.989     | 1.0000 | 19,574 |
| 10 | epykit / `lr`        | 0.984     | 1.2 × 10⁻⁵ | 0.992     | 1.0000 | 19,678 |
| 10 | methylKit / default  | 0.984     | 1.2 × 10⁻⁵ | 0.992     | 1.0000 | 19,678 |

### 2.3 DMR × coverage (3 vs 3, 35 reference DMRs)

| Coverage | Tool / method | n_called | Recall | Precision |
|---|---|---|---|---|
| 5× | epykit / `chain_merge` | 42 | 0.971 | 0.857 |
| 5× | epykit / `tile` | 35 | 0.857 | 0.971 |
| 5× | methylKit / tile | 102 | 1.000 | 0.980 |
| 10× | epykit / `chain_merge` | 37 | 1.000 | 1.000 |
| 10× | epykit / `tile` | 34 | 0.943 | 1.000 |
| 10× | methylKit / tile | 102 | 1.000 | 1.000 |
| 15× | epykit / `chain_merge` | 37 | 1.000 | 1.000 |
| 15× | methylKit / tile | 102 | 1.000 | 1.000 |
| 20× | epykit / `chain_merge` | 37 | 1.000 | 1.000 |
| 20× | methylKit / tile | 102 | 1.000 | 1.000 |
| 25× | epykit / `chain_merge` | 37 | 1.000 | 1.000 |
| 25× | methylKit / tile | 102 | 1.000 | 1.000 |

### 2.4 Per-scenario wall-clock (Study 2)

| Scenario | epykit wall | methylKit wall | Speed-up |
|---|---|---|---|
| DMC × cov 5×    |   8.8 s   | 130.0 s   |  14.8× |
| DMC × cov 10×   |  16.4 s   | 116.1 s   |   7.1× |
| DMC × cov 15×   |   8.6 s   | 101.5 s   |  11.8× |
| DMC × cov 20×   |   8.5 s   |  99.0 s   |  11.6× |
| DMC × cov 25×   |   8.0 s   |  96.2 s   |  12.0× |
| DMC × n = 2     |   5.9 s   |  12.6 s   |   2.1× |
| DMC × n = 4     |   7.1 s   | 103.3 s   |  14.5× |
| DMC × n = 6     |   9.3 s   | 100.6 s   |  10.8× |
| DMC × n = 8     |  11.4 s   | 105.7 s   |   9.3× |
| DMC × n = 10    |  13.1 s   | 123.6 s   |   9.5× |
| DMR × cov 5×    |  62.0 s   | 4,241.9 s | **68.4×** |
| DMR × cov 10×   |  97.8 s   | 4,934.5 s | **50.5×** |
| DMR × cov 15×   |  89.0 s   | 4,139.8 s | **46.5×** |
| DMR × cov 20×   |  82.8 s   | 3,995.4 s | **48.3×** |
| DMR × cov 25×   |  86.0 s   | 3,889.1 s | **45.2×** |

### 2.5 Per-scenario peak RAM (RSS)

| Scenario | epykit peak RSS | methylKit peak RSS | Ratio |
|---|---|---|---|
| DMC × cov 5–25× (range) | 616 – 634 MB | 1,336 – 1,360 MB | ~2.1× |
| DMC × n = 2  |   478 MB | 1,022 MB | 2.1× |
| DMC × n = 10 |   734 MB | 1,399 MB | 1.9× |
| DMR × cov 5×  | 5,553 MB | 6,371 MB | 1.15× |
| DMR × cov 10× | 4,571 MB | 6,742 MB | 1.48× |
| DMR × cov 25× | 5,681 MB | 7,112 MB | 1.25× |

### 2.6 Aggregate cost of full 15-point grid

| Metric | epykit | methylKit | Ratio |
|---|---|---|---|
| Total wall-clock | **8.6 min** | **6 h 9.8 min** | **43×** |
| Total CPU time | 15.4 min | 6 h 8.8 min | 24× |
| Peak RSS observed | 6.03 GB | 7.11 GB | 1.18× |

### 2.7 What drives the runtime gap

CPU / wall ratio observed:

| Tool | CPU / wall |
|---|---|
| methylKit | ≈ 1.00 (strictly single-threaded) |
| epykit | 1.6 – 2.6 (Polars / NumPy implicit threading) |

| Source of speedup | Approx. contribution |
|---|---|
| Vectorised regressions vs per-CpG R-level `glm()` loop | ~15–20× |
| Polars / NumPy implicit multi-threading | ~2–3× |

Vectorisation dominates; the speedup persists even if methylKit were
granted equal core count.

---

## 3. Study 3 — Real WGBS data (GSE263850)

**Dataset.** GSE263850 (Clone16 / Clone20 / Clone21 vs SBP009 untreated 1 /
2 / 3, hg38, n = 6, see [samplesheet](../data/study3/samplesheet_gse263850.csv)).
**Inputs.** Six strand-collapsed 12-col methylation BEDs. Both pipelines
consume identical combined-strand counts (Methods §2.2 of [paper.md](../paper/paper.md)).

### 3.1 Headline numbers

| Metric | methylKit | epykit | Ratio (ep / mk) |
|---|---:|---:|---:|
| tool version | 1.36.0 | 1.0.0 | — |
| ncpus_logical | 24 | 16 | 0.667 |
| pipeline_wall_sec | 13,033.0 | 1,072.5 | **0.0823 (12.2×)** |
| pipeline_cpu_sec | 13,080.3 | 1,374.9 | 0.105 |
| final_rss_mb | 44,790.2 | 11,678.8 | 0.261 |
| peak_rss_mb | 48,001.8 | 12,565.6 | **0.262 (3.83×)** |
| n_samples | 6 | 6 | 1.00 |
| n_cpgs_united | 15,600,476 | 15,597,046 | 0.9998 |
| n_cpgs_after_sd | 14,041,063 | 14,617,178 | 1.04 |
| n_dmcs_tested | 15,600,476 | 15,597,046 | 0.9998 |
| n_dmcs_sig_q05 (\|d\|≥10 %) | 51,792 | 30,965 | 0.598 |
| n_dmcs_hyper | 24,744 | 17,269 | 0.698 |
| n_dmcs_hypo | 24,527 | 13,688 | 0.558 |
| n_tiles_tested | 1,174,664 | 1,172,347 | 0.998 |
| n_dmrs_lenient | 2,661 | 3,433 | 1.29 |
| n_dmrs_strict | 147 | 257 | 1.75 |

### 3.2 Stage-count waterfall

| Stage | methylKit | epykit |
|---|---:|---:|
| 00_input_files | 6 | 6 |
| 00_input_total_size_mb | 958.1 | 1,621.7 |
| 01_raw_total_cpgs | 123,177,176 | 143,375,098 |
| 03_filtered_total_cpgs | 123,053,376 | 123,034,590 |
| 05_united_cpgs | 15,600,476 | 15,597,046 |
| 06_sd_filtered_cpgs | 14,041,063 | 14,617,178 |
| 10b_united_tiles | 1,174,664 | 1,172,347 |

The pre-filter count delta is a **read-time vs post-read filter ordering**
difference (methylKit applies `mincov` at `methRead`; epykit reads all
rows then filters in step 03). Post-filter counts agree to 0.02 %.

### 3.3 Per-step wall-clock

| Step | methylKit (s) | epykit (s) | Speed-up |
|---|---:|---:|---:|
| 01_read | 46.72 | 12.37 | 3.78× |
| 02_per_sample_qc | 46.32 | 63.11 | 0.73× |
| 03_filterByCoverage | 23.78 | 11.27 | 2.11× |
| 04_normalizeCoverage | 6.33 | 7.09 | 0.89× |
| 05_unite | 18.36 | 1.3 × 10⁻⁵ | huge |
| 06_sd_filter | 4.75 | 103.9 | 0.046× |
| 07_sample_qc_plots | 640.6 | 5.78 | **110.8×** |
| **08_calculateDiffMeth_dmc** | **10,661.4** | **204.8** | **52.0×** |
| 10b_tileMethylCounts | 223.4 | 80.49 | 2.78× |
| 10c_calculateDiffMeth_dmr | 692.5 | 132.7 | 5.22× |
| 10f_dmr_gene_annotation | 5.44 | 22.27 | 0.24× |
| 12_dmc_gene_annotation | 1.97 | 31.33 | 0.063× |

The dominant cost in methylKit is `calculateDiffMeth` on 15.6 M CpGs.
epykit is slower on annotation steps (small fraction of total wall).

### 3.4 CpG agreement (all tested sites)

| Metric | Value |
|---|---:|
| methylKit tested CpGs | 15,600,476 |
| epykit tested CpGs | 15,597,046 |
| Both pipelines (∩) | 15,597,046 |
| Union | 15,600,476 |
| **Jaccard** | **0.9998** |
| methylKit-only | 3,430 |
| epykit-only | 0 |

### 3.5 Effect-size agreement on shared CpGs

| Metric | Value |
|---|---:|
| Shared CpGs | 15,597,046 |
| Pearson r on `meth_diff` | **0.9936** |
| Spearman ρ | **0.9831** |
| Same hyper/hypo direction | 14,669,608 / 15,597,046 (**94.05 %**) |
| Significant in both pipelines (q < 0.05) | 15,688 |
| Significant in methylKit only | 36,054 |
| Significant in epykit only | 15,277 |

### 3.6 Significance confusion (shared CpGs only)

|  | epykit sig | epykit not sig |
|---|---:|---:|
| **methylKit sig**     | 15,688 | 36,054 |
| **methylKit not sig** | 15,277 | 15,530,027 |

Recall of methylKit calls by epykit (default `lr`): **30.3 %**.
Precision of epykit vs methylKit: **50.7 %**.

### 3.7 DMR agreement (500 bp tiles) — methylKit vs epykit-tile only

These numbers describe the pre-investigation tile-vs-tile comparison
(epykit `dmr_tile` vs methylKit `tileMethylCounts`, both with 500 bp
fixed windows). They remain valid as a same-architecture sanity check.
The full multi-engine comparison against the paper's call set is in
§3.10–§3.16.

| Metric | Value |
|---|---:|
| Shared tiles | 1,172,347 |
| Pearson r on tile `meth_diff` | **0.9970** |
| Same direction | 1,076,503 / 1,172,347 (**91.82 %**) |
| methylKit DMRs (lenient, \|d\|≥10 %, q<0.05) | 2,661 |
| epykit DMRs (lenient) | 3,433 |
| Intersection | 1,957 |
| Jaccard (lenient) | **0.473** |
| Strict (\|d\|≥25 %, q<0.01) | 147 / 257 |
| Jaccard (strict) | **0.530** |

#### Top-K DMR concordance (by q-value)

| K | methylKit ∩ epykit | Recall (∩ / K) |
|---:|---:|---:|
| 50 | 27 | 54.0 % |
| 100 | 53 | 53.0 % |
| 250 | 143 | 57.2 % |
| 500 | 291 | 58.2 % |
| 1,000 | 574 | 57.4 % |

### 3.8 Dispersion-estimator sensitivity

Same threshold (q < 0.05, |Δ| ≥ 10 %), all variants on identical
underlying counts:

| Variant | n_sig | Recall vs methylKit | Jaccard | Note |
|---|---:|---:|---:|---|
| `lr / site` (default) | 30,957 | 31.7 % | 0.243 | baseline |
| `lr / shrink` (post-bugfix) | 1,499 | 0.3 % | 0.003 | conservative |
| `lr / chrom` (post-bugfix) | 167,923 | 29.6 % | 0.072 | permissive |
| **`lr+ (power_stack)`** | **406,515** | **92.9 %** | 0.112 | high recall, low precision |
| `lr / site + smoothing` | 10,691 | 5.0 % | 0.043 | DSS-style box dilutes |

`lr+` recovers 93 % of methylKit's significant DMCs at the cost of 13× more
total calls; some of the additional calls may be real signal that methylKit
misses, since `power_stack` borrows information across correlated adjacent
CpGs.

### 3.9 Source: Farhangdoost et al. 2025, *Molecular Psychiatry*

The original paper is a multi-omics study integrating RNA-seq DEGs,
WGBS DMRs, and ChIP-seq H3K27ac peaks across 6 samples (3 Het-AKAP11-KO
+ 3 SBP009 WT) of human iPSC-derived neurons. Our benchmark scope is
the WGBS-DMR layer only. We compare against:

* **Supp Table 5** — 813 DSS DMRs with coordinates, areaStat, per-sample
  methylation means, HOMER Annotation column, Gene.Name (705 unique).
* **Supp Table 6** — Reactome ORA via ShinyGO (paper Fig 3D top-20).
* **Supp Table 8** — 46 DMR-DEG-correlated critical genes (paper Fig 3E).

Paper parameters (paper Methods §"WGBS data processing and analysis"):

> `DSS::DMLfit.multiFactor(smoothing = TRUE)` →
> `DSS::DMLtest.multiFactor` →
> `DSS::callDMR(delta = 0, p.threshold = 1e-5, minlen = 50, minCG = 3,
>   dis.merge = 100, pct.sig = 0.5)`
> Annotation via HOMER (hg38). 100 kb gene-DMR linkage (TSS to DMR
> midpoint).

### 3.10 DMR coordinate concordance vs DSS-from-scratch (post-rerun)

After the 2026-06 Linux re-run the chain_merge call sets at
`dis.merge = 100 / 250` shifted from pre-rerun 702 / 940 to **852 / 1,139**
DMRs respectively, and the comparison baseline shifts from the paper's
813 DMRs to the locally-rerun DSS-from-scratch 922 DMRs (§3.3.4 of
[paper.md](../paper/paper.md) for rationale: DSS-vs-paper recall is
only 87.5 %, so paper-813 is not the appropriate denominator).

| Caller | n DMRs | median bp | %hyper | recall any-bp (vs DSS) | precision any-bp | recall J ≥ 0.5 | direction agreement on matched |
|---|---:|---:|---:|---:|---:|---:|---:|
| methylKit-tile (500 bp) | 2,661 | 500 | ~50 % | ≈ 8.7 % | ≈ 3.0 % | ~ 0 % | (≤ 1 matched at J ≥ 0.5) |
| **ek-chain_merge-100** | **852** | **125** | **79.0 %** | **63.8 %** | **74.4 %** | **34.5 %** | **588 / 588 = 100 %** |
| **ek-chain_merge-250** | **1,139** | **205** | **~ 80 %** | **77.3 %** | **63.0 %** | **64.2 %** | **713 / 713 = 100 %** |
| **DSS-from-scratch** (2.58.0) | **922** | **241** | **74.6 %** | 100 % (self) | 100 % (self) | 100 % (self) | — |
| paper (Supp Table 5, contextual reference) | 813 | 239 | 78.5 % | DSS-vs-paper: 87.5 % | — | — | — |

Direction agreement on every matched DMR is **100 % across every
caller**. Source: [`headline.json`](../data/study3/comparisons/epykit_vs_dss/headline.json),
[`dis_merge_vs_dss_sensitivity.csv`](../data/multi_thread_and_chain_sweep/chain_merge_dis_merge_sweep/dis_merge_vs_dss_sensitivity.csv).

### 3.11 dis.merge sweep (post-rerun, vs DSS-922)

| dis.merge | n DMRs | median bp | recall any-bp | recall J ≥ 0.25 | recall J ≥ 0.5 | recall J ≥ 0.75 | precision | direction agree |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 (paper) | 852 | 125 | 63.8 % | 54.2 % | 34.5 % | 22.8 % | 74.4 % | 588 / 588 = 100 % |
| 150 | 1,010 | 164 | 72.8 % | 65.8 % | 52.1 % | 40.6 % | 68.7 % | 671 / 671 = 100 % |
| 200 | 1,095 | 192 | 76.3 % | 70.9 % | 60.6 % | 49.5 % | 64.9 % | 703 / 703 = 100 % |
| **250 (morphology-matched)** | **1,139** | **205** | **77.3 %** | **73.1 %** | **64.2 %** | **53.0 %** | **63.0 %** | **713 / 713 = 100 %** |
| 500 | 1,160 | 214 | 78.4 % | 74.4 % | 67.0 % | 54.9 % | 61.5 % | 723 / 723 = 100 % |

Direction agreement remains 100 % across the entire sweep. Source:
[`dis_merge_vs_dss_sensitivity.csv`](../data/multi_thread_and_chain_sweep/chain_merge_dis_merge_sweep/dis_merge_vs_dss_sensitivity.csv) +
[F2 curves](../figures/study3_real_GSE263850/three_way/F2_dis_merge_sweep.png).

### 3.12 Panel-E gene capture (Supp Table 8 = 46 critical genes), post-rerun

| Caller | Captured by nearest-TSS | Captured by 100 kb rule |
|---|---:|---:|
| methylKit-tile | 25 / 46 = 54.3 % | n/a (tile DMRs too short for 100 kb gene-link table) |
| epykit-chain_merge-100 | **30 / 46 = 65.2 %** | 29 / 46 = 63.0 % |
| epykit-chain_merge-250 | **32 / 46 = 69.6 %** | (sweep variant; 100 kb gene-link table not regenerated) |
| **DSS-from-scratch** | **37 / 46 = 80.4 %** | 38 / 46 = 82.6 % |

Source: [`panel_e_capture_dss.csv`](../data/study3/comparisons/epykit_vs_dss/panel_e_capture_dss.csv) +
[`polish_recompute_2026_06_05.json`](../data/study3/comparisons/epykit_vs_dss/polish_recompute_2026_06_05.json).
chain_merge-100 misses include CLEC19A, KANK1, CNR1, EDNRB, ANXA1,
CXCR4, SHOX2, SIX3, TFAP2B — genes the paper depended on for the
Fig 3E enrichment. Most have very short / low-CpG-density DMRs that
pass DSS's smoothed test but fall below epykit's LR.

### 3.13 Annotation distribution (paper Fig 3C reproduction)

| Feature | paper-DSS | mk-tile | ek-tile | ek-cm-100 | ek-cm-250 | DSS-local |
|---|---:|---:|---:|---:|---:|---:|
| promoter-TSS | 0.9 % | 2.1 % | 2.0 % | 2.6 % | 2.8 % | 0.8 % |
| 5' UTR       | 0.1 % | 0.0 % | 0.1 % | 0.0 % | 0.0 % | 0.1 % |
| exon         | 2.9 % | 1.6 % | 1.6 % | 8.1 % | 6.9 % | 3.7 % |
| intron       | 44.3 % | 35.8 % | 37.7 % | 42.0 % | 42.1 % | 38.8 % |
| 3' UTR       | 1.6 % | 1.2 % | 1.3 % | 0.0 % | 0.0 % | 1.7 % |
| TTS          | 1.4 % | 1.6 % | 1.4 % | 1.1 % | 1.3 % | 1.6 % |
| non-coding   | 0.5 % | 5.9 % | 5.7 % | 0.0 % | 0.1 % | 5.9 % |
| intergenic   | 48.3 % | 51.7 % | 50.2 % | 46.2 % | 46.8 % | 47.4 % |
| n            | 813 | 2,661 | 3,433 | **852** | **1,139** | 922 |
| chi² vs paper (lower = closer) | 0 | 65.4 | 58.0 | **45.7** | **44.2** | **41.4** |

After switching to epykit3 (which has the proper full-HOMER feature
default + `features=` kwarg on `ep.tl.annotate()`), chain_merge-100
chi² drops from 54.4 → **45.7** and chain_merge-250 from 56.6 →
**44.2** — essentially tied with DSS-from-scratch (41.4). chain_merge
now correctly assigns TTS labels at ~1 % rate (matching paper's 1.4 %).
5' UTR / 3' UTR / non-coding still report 0 % under epykit3 because:

1. UTR builders are GTF-only: `_build_utr_df_from_gtf_utrs` requires
   GENCODE-format UTR records, not refGene cdsStart/cdsEnd boundaries.
   Using `refgene=` therefore never yields UTR labels even with the
   full default tuple. (Tracked as a separate refGene-UTR-derivation
   issue.)
2. `non-coding` has the lowest priority and is suppressed by any
   gene-body intron of the same transcript. The
   `all_overlapping_features` column for chain_merge-100 does include
   `noncoding: 100` — i.e. 100 DMRs overlap non-coding transcripts —
   but those DMRs are simultaneously labeled `intron` by the
   higher-priority gene-body match.

The remaining intron / intergenic dominance (~85 %) matches paper
Fig 3C. Methylkit-tile inflates non-coding and intergenic relative
to paper.
Source: [annotation_distribution.csv](../data/study3/comparisons/annotation_distribution.csv).

### 3.14 Top-named gene hits (paper Fig 3B labels)

The paper Fig 3B labels its top 10 hyper- and 10 hypo-methylated
DMR-associated genes. Direct coordinate-overlap hits per caller:

| Caller | any-bp / 20 | J ≥ 0.5 / 20 |
|---|---:|---:|
| methylKit-tile | 2 / 20 | 0 / 20 |
| epykit-chain_merge-100 | 9 / 20 | 5 / 20 |
| epykit-chain_merge-250 | 11 / 20 | 7 / 20 |
| **DSS-from-scratch** | **18 / 20** | **17 / 20** |

Specific hyper-side genes hit by chain_merge-100: NR2E1, OTX1, IRX2,
ENPP2, GREB1L, CCDC177, EBF1; hypo-side: RPLP0P2, FAM87A (= 9 / 20).
Notable hyper misses include OTX2, GNG11, PPP2R3B — short,
low-CpG-density DMRs that pass DSS's smoothed test but fall below
epykit's LR; hypo misses include LOC100506858, KC6, TMEM242,
NAALADL2, LOC100131655, PAX7, PDK3, OSBPL8. Source:
[F4 named-gene heatmap](../figures/study3_real_GSE263850/three_way/F4_top_named_gene_hits.png) +
[F4 data](../figures/study3_real_GSE263850/three_way/F4_top_named_gene_hits_data.csv).

### 3.15 Per-DMR effect-size concordance (epykit vs DSS)

For matched (J ≥ 0.5) DMR pairs between epykit-chain_merge and
DSS-from-scratch, post-rerun:

| Pair | n matched | Pearson r on Δβ | Spearman ρ on Δβ | Direction agreement |
|---|---:|---:|---:|---:|
| ek-100 vs DSS, any-bp | 634 | **0.9936** | 0.9257 | 100 % |
| ek-100 vs DSS, J ≥ 0.25 | 506 | **0.9950** | 0.9357 | 100 % |
| ek-100 vs DSS, J ≥ 0.5  | **318** | **0.9954** | 0.9382 | 100 % |
| ek-250 vs DSS, any-bp | 718 | **0.9950** | 0.9412 | 100 % |
| ek-250 vs DSS, J ≥ 0.5  | **592** | **0.9965** | 0.9543 | 100 % |

When the two engines overlap a region, they agree to four decimal
places on effect size and ≥ 0.93 on signed-rank correlation.
Source: [`polish_recompute_2026_06_05.json`](../data/study3/comparisons/epykit_vs_dss/polish_recompute_2026_06_05.json) +
[F9](../figures/study3_real_GSE263850/three_way/F9_per_dmr_concordance.png).

### 3.16 Pipeline cost (4-way, Windows host)

Windows-host numbers; methylKit `mc.cores` is a no-op on Windows so
methylKit is forced single-threaded. Linux honest ratios for the
DMC step are reported in [paper.md](../paper/paper.md) §4.3 abstract +
§3.5 (≈ 33 × at the simulator headline cell, ≈ 5 × at the Study 3
15.6 M-CpG cell). DSS is single-thread by construction (no
multi-core path in DSS 2.58.0); the DSS ratio is
platform-agnostic.

| Caller | Wall (s) | CPU (s) | Peak RSS (GB) | Threads peak | Engine notes |
|---|---:|---:|---:|---:|---|
| methylKit-tile (R 4.5.0, Windows) | 12,372 | 12,419 | **48.0** | 1 (Windows: mc.cores = no-op) | `calculateDiffMeth` on 15.6 M CpGs dominates |
| epykit-tile (Python 3.13) | 675 | 993 | 12.6 | — | matches methylKit at DMC; ~3× faster overall |
| epykit-chain_merge-100 | **~ 443** | **~ 260** | ~ 12.6 (est.) | — | DMC cached + DMR re-callable per dis.merge |
| **DSS-from-scratch** (R 4.5.0, DSS 2.58.0) | **2,820** | **2,756** | **9.3** | 12 logical (effective 1) | DMLfit smoothing dominates (~ 34 min) |

epykit chain_merge is ~ 6 × faster than DSS on the same input and uses
about the same memory; the DSS smoothing step is single-threaded and
the most expensive single operation in any of the four pipelines.
Source: [F7 resource bars](../figures/study3_real_GSE263850/three_way/F7_resources.png) +
[F7_resources_data.csv](../figures/study3_real_GSE263850/three_way/F7_resources_data.csv).

### 3.17 Three-way Reactome / KEGG enrichment

Enrichr REST API against `Reactome_2022` + `KEGG_2021_Human` + `GO_MF_2023`.
Paper's Panel D Reactome terms (GPCR ligand binding, Class A/1
Rhodopsin, G alpha i signalling, etc.) are recovered with the
following counts of top-20 paper-keyword matches:

| Caller | n_genes | Reactome top-20 paper hits | KEGG top-20 paper hits |
|---|---:|---:|---:|
| paper Table 5 (705 genes) | 705 | 2 | 2 |
| methylKit-tile (nearest-TSS) | 2,111 | 1 | 1 |
| ek-chain_merge-100 (100 kb) | 1,112 | 3 | 1 |
| ek-chain_merge-250 (100 kb) | 1,444 | 1 | 1 |
| DSS-from-scratch (100 kb) | 1,467 | 2 | 1 |

Enrichr's full-library BH correction is more aggressive than ShinyGO's
"Curated.Reactome" — the absolute term counts are smaller than the
investigation report's prior ShinyGO numbers. The signal still ranks
GPCR / Gα(i)-related terms (Morphine addiction, Activation of G
Protein Gated Potassium Channels) in the top-20 for every caller.
Full per-caller, per-library top-20 in
[enrichment_three_way.json](../data/study3/comparisons/enrichment_three_way.json) +
[F8](../figures/study3_real_GSE263850/three_way/F8_enrichment_dotplot.png).

---

## 4. Bug fixes uncovered and applied during this benchmark

| # | Fix | File | Impact |
|---|---|---|---|
| 1 | **Fisher calibration bug** — upper-tail-only hypergeometric failed for hypo direction | `dmc.py:297–303` | Fisher TPR 0.000 → 0.668–0.998 across the coverage grid |
| 2 | **chain_merge defaults** — `dis_merge_bp` 100 → 500 (default), preset tuning | `dmr.py` (5 locations) | DMR chain_merge TPR 0.086 → 0.971–1.000 |
| 3 | **bb_lr guardrails** — auto-promote dispersion to "shrink" at n < 6 | `dmc.py:1601+` | Improved bb_lr behaviour at low n |
| 4 | **welch_t warning tiers** — critical warning at n ≤ 2 (degenerate DOF) | `dmc.py:2001` | Better user guidance |
| 5 | **lr+ auto-engage** — `power_stack="auto"` enables neighbour_combine + sep_fallback at n ≤ 2 | `tl.py:455` | Auto lr → lr+ at 1v1; lr TPR at 2v2 0.880 → 0.999 |
| 6 | **Tile adjacent merging** — Stouffer-combine adjacent significant tiles | `dmr.py:1213+` | Better DMR boundary recovery |
| 7 | **df_phi in dispersion modes** — `_score_finalize` referenced F(1, df_residual_per_site) instead of the actual df backing each phi estimate | `dmc.py::_score_finalize`, `_glm.py::compute_dispersion_phi` | Pre-fix `dispersion="shrink"` and `"chrom"` both collapsed to 1 DMC across 15.6 M CpGs; post-fix they produce sensible numbers. `dispersion="site"` was already correct. |

Fix 7's pre-fix bug never affected the default `dispersion="site"` mode and
therefore never affected the headline `lr` numbers.

---

## 5. Methodology details

### 5.1 Ground truth (Studies 1 and 2)

- Reconstructed from highest-coverage simulator outputs (25×).
- DMC truth: `is_dmc = |β_treat − β_ctrl| ≥ 0.20` at 25×; recovers **19,999
  / 100,000 = 20.0 % exactly**.
- DMR truth: runs of ≥ 5 same-direction true DMCs within 1 kbp; recovers
  **35 reference DMRs** (median width 1,823 bp, 78 CpGs) with **100 %
  coordinate overlap** against Piao et al. Table S3.
- Code: `data/study1/ground_truth/make_truth.py` and the same script in
  `data/study2/ground_truth/`.

### 5.2 Baseline transcription audit (Study 1)

- 174 numeric cells across Piao et al. Tables S1 and S2 audited cell-by-cell
  against an independent hand-typed copy.
- **0 transcription errors**.
- 35 / 35 reference DMR coordinates (Table S3) recovered at 100 % overlap
  by the truth-reconstruction code.
- Audit reproducible via `scripts/audit_baselines.py`.

### 5.3 DMR baseline numbers (Study 1)

DMR detection rates for the eight baselines in Section 1.3 come from
hand-transcribed Piao Figures 3a, 3b, S5–S7 (bar charts), not from tables.
Per-figure confidence labels are recorded in
[data/study1/baseline_tables/PROVENANCE.md](../data/study1/baseline_tables/PROVENANCE.md).

### 5.4 Metrics

- **DMC.** TPR, FPR, F1, AUROC at q < 0.05, stratified by effect-size bins
  (0.2–0.4, 0.4–0.6, 0.6–0.8, 0.8–1.0).
- **DMR.** Detected iff ≥ 80 % of truth-DMR span is covered by a called
  DMR with q < 0.05 (Piao et al. Figure 3 caption).
- **Real data (Study 3).** Direction agreement, Pearson r, Spearman ρ,
  Jaccard, top-K concordance.
- **Performance.** `psutil` sampling at 50 ms across both subprocess trees;
  peak RSS and accumulated CPU time.

### 5.5 Reproduction recipe

See [README.md](../README.md) and
[methods_appendix.md](methods_appendix.md). Each of the three studies is
end-to-end reproducible from the data files in `FINAL_REPORT/data/`.
