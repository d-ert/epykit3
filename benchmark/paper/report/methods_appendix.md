# Methods Appendix — Reproducibility details

Companion to [paper.md](../paper/paper.md) and [REPORT.md](REPORT.md). This
file collects exact tool versions, parameters, command lines, file
schemas, and known limitations needed to reproduce the three benchmark
studies.

---

## 1. Software versions

| Tool | Version | Where used |
|---|---|---|
| Python | 3.12 | All studies |
| epykit | 1.0.0 (all studies) | All studies |
| polars | 1.x | All studies (lazy I/O) |
| NumPy, SciPy, statsmodels | latest at time of run | All studies |
| psutil | ≥ 5.9 | Study 2 perf sampling |
| R | 4.5.0 | Studies 2 and 3 |
| methylKit | 1.34.0 (Study 2); 1.36.0 (Study 3); 0.99.2 (Study 1 baseline tables, transcribed) | All studies |
| **DSS** | **2.58.0 (Study 3 local re-run; per `data/study3/dss/dss_session_info.txt`)**; 2.12.0 (Study 1 baseline, transcribed from Piao 2021) | Studies 1 (transcribed), 3 (local) |
| Other Study 1 baselines | methylSig 0.4.4, RADMeth, BiSeq, Fisher (pooled), BSmooth, metilene — versions from Piao et al. 2021 | Study 1 (transcribed) |

All three studies were executed at commit `60a71e0` (engine tag
`v0.7.5-phase3-engines-frozen`, 2026-05-28); the package version string was
`0.7.2` at run time and the engine code is unchanged through the `1.0.0`
release, so the results apply to 1.0.0. (An earlier Study-3 pass on epykit
0.6.0 predated the `df_phi` dispersion-estimator fix; it was superseded by
the frozen-engine re-run reported here.)

---

## 2. Hardware

* **Studies 1 & 2.** Single Windows workstation. psutil samples both
  subprocesses at 50 ms; numbers are OS-level (not the tools'
  self-reports). methylKit's `mc.cores` is a no-op on Windows (no
  `fork()`), so methylKit ran single-threaded by force.
* **Study 3.** Two hosts: methylKit on `pivoine` (24 logical CPUs); epykit
  on `DESKTOP-0GUMIA4` (16 logical CPUs). Same input files; same parameter
  choices. Comparison is per-pipeline pipeline_wall_sec, not per-core.

---

## 3. Parameters — simulated benchmark (Studies 1 and 2)

### 3.1 DMC pipeline

```python
ep.tl.dmc(
    md,
    test="lr",                  # or "lr+", "bb_lr", "welch_t", "fisher"
    dispersion="site",          # default; see §4 for "shrink"/"chrom"/"eb"
    fdr_method="fdr_bh",        # or "fdr_tsbh" inside lr+
    allow_n1=True,              # enables n=1 fallback for lr
)
# Significance cut:
sig = (df["qvalue"] < 0.05) & (df["meth_diff"].abs() >= 0.25)
```

### 3.2 DMR pipeline

* **`dmr_tile`**: 1 kbp fixed tiles, read-pooled binomial test per tile,
  ≥ 5 CpGs per tile.
* **`dmr_chain_merge`**: DSS-callDMR semantics; default after fix
  `dis_merge_bp = 500`, `min_cpgs = 5`, q < 0.05 per merged region.

### 3.3 Ground truth construction (Studies 1 and 2)

```python
# In data/study{1,2}/ground_truth/make_truth.py
TRUTH_THRESHOLD = 0.20         # effect-size threshold
DMR_GAP_BP      = 1000         # max gap between same-direction true DMCs
DMR_MIN_CPGS    = 5            # min CpGs per reference DMR
```

Recovered: 19,999 DMCs (20.0 %); 35 reference DMRs.

### 3.4 methylKit driver (Study 2)

```r
mk <- methRead(samplesheet, pipeline = "bismarkCoverage",
               treatment = c(rep(1, n_trt), rep(0, n_ctrl)),
               mincov = 10, sample.id = sample_ids)
mk <- normalizeCoverage(mk, method = "median")
mk <- unite(mk, destrand = FALSE)
diff <- calculateDiffMeth(mk, mc.cores = 1)
dmc <- getMethylDiff(diff, difference = 25, qvalue = 0.05)
```

For DMRs:

```r
tiles <- tileMethylCounts(mk, win.size = 1000, step.size = 1000,
                          cov.bases = 5)
tiles_united <- unite(tiles, destrand = FALSE)
diff_t <- calculateDiffMeth(tiles_united, mc.cores = 1)
dmr <- getMethylDiff(diff_t, difference = 25, qvalue = 0.05)
```

Same `mincov = 10`, same `q < 0.05`, same `difference = 25` (= 0.25 on
fractional scale) for both pipelines.

---

## 4. Parameters — real data (Study 3, GSE263850)

### 4.1 Input

Six 12-column strand-collapsed BEDs (cols 10–12 = M / T / pct after
strand merging). methylKit was fed pre-converted 6-col Bismark `.cov.gz`
derived from these BEDs; epykit reads cols 10–12 directly via
`read_combined_strand_bed()`. Per-CpG counts entering each pipeline's
coverage filter are **bit-identical**.

### 4.2 Both pipelines

| Parameter | Value |
|---|---|
| min coverage | 10 |
| top-percentile clipping | 99.9 |
| FDR method | Benjamini–Hochberg |
| DMC threshold | q < 0.05, \|meth_diff\| ≥ 0.10 |
| DMR tile size | 500 bp |
| DMR min CpGs per tile | 5 |
| DMR threshold (lenient) | q < 0.05, \|d\| ≥ 0.10 |
| DMR threshold (strict) | q < 0.01, \|d\| ≥ 0.25 |
| Reference genome | hg38 |

### 4.3 Coordinate convention

methylKit `.cov` is 1-based (`pos = bed.end`); epykit stores 0-based BED
positions (`pos = bed.start`). The two differ by exactly 1 bp at every CpG.
Overlap analyses apply a `+1` shift on epykit's `pos`/`start` to align.

### 4.4 Sample sheet

[data/study3/samplesheet_gse263850.csv](../data/study3/samplesheet_gse263850.csv).
Groups:

* **sbp009** (control, n = 3): SBP009 untreated 1, 2, 3
  (GSM8200109, GSM8200110, GSM8200111).
* **clone** (treatment, n = 3): Clone16, Clone20, Clone21.

---

## 5. Evaluation metrics — implementation

### 5.1 TPR / FPR / F1

* TPR = TP / (TP + FN), FN = true DMC not called at q < 0.05.
* FPR = FP / (FP + TN), TN = CpG not in any true-DMC effect-size bin and
  not called.
* F1 = 2 · precision · recall / (precision + recall), with
  precision = TP / (TP + FP).

### 5.2 AUROC

Mann–Whitney U statistic with `1 − pvalue` as the score, normalised to
[0, 1]. Reported only for epykit engines (baselines have only the single
q < 0.05 operating point).

### 5.3 DMR detection (Study 1)

A reference DMR is detected when the **union** of called DMRs overlapping
it covers ≥ 80 % of its span. This matches Piao et al. Figure 3 caption.
Precision is (truth-DMRs covered) / (total DMRs called).

### 5.4 Real-data agreement (Study 3)

* **Direction agreement.** sign(meth_diff_mk) == sign(meth_diff_ep) on
  shared CpGs / tiles.
* **Pearson r, Spearman ρ.** Computed on shared CpG / tile keys after
  the 1-bp coordinate alignment.
* **Jaccard.** |A ∩ B| / |A ∪ B| on the set of significant call IDs.
* **Top-K concordance.** Take each tool's top K calls by q-value, count
  the intersection.

---

## 6. File schemas

### 6.1 `eval_summary.parquet` (Study 1, [data/study1/eval_summary.parquet](../data/study1/eval_summary.parquet))

Tidy long-format: one row per (tool, scenario, scenario_value,
effect_size_bin, metric). Columns:

| Column | Type | Description |
|---|---|---|
| tool | str | `epykit_lr`, `methylKit`, `DSS`, `RADMeth`, … |
| scenario | str | `coverage` or `replicate` |
| scenario_value | int | 5, 10, 15, 20, 25 (coverage) or 2, 4, 6, 8, 10 (replicate) |
| effect_size_bin | str | `0.2-0.4`, `0.4-0.6`, `0.6-0.8`, `0.8-1.0`, or `all` |
| metric | str | `tpr`, `fpr`, `f1`, `auroc` |
| value | float | metric value (NaN where not reported) |
| source | str | `our_run` (epykit) or `piao2021_table_s1` (baseline) |

### 6.2 `timings.parquet` (Study 1, [data/study1/timings.parquet](../data/study1/timings.parquet))

Per (tool, scenario, scenario_value) wall-clock and CPU time.

### 6.3 methylKit results (Study 2, [data/study2/methylkit_results/](../data/study2/methylkit_results/))

One `.tsv` per scenario:

| Column | Type | Description |
|---|---|---|
| chr, start, end | str, int, int | CpG / tile coordinates (1-based) |
| pvalue | float | raw p-value |
| qvalue | float | BH-adjusted q-value |
| meth.diff | float | difference in percent methylation (−100 … +100) |

### 6.4 Study 3 benchmark CSVs (`data/study3/benchmark/`)

| File | Content |
|---|---|
| `run_summary.csv` | Headline numbers (Table 1 in REPORT §3.1) |
| `step_benchmarks.csv` | Per-step wall-clock and peak RSS |
| `stage_counts.csv` | Stage-count waterfall |
| `dmc_threshold_sweep.csv`, `dmr_threshold_sweep.csv` | Sweep over q and Δβ |
| `dmc_pvalue_histogram.csv`, `dmr_pvalue_histogram.csv` | p-value distributions |
| `package_versions.csv`, `session_info.txt` | Reproducibility metadata |

---

## 7. Reproduction recipe (Study 1)

```bash
# Repository assumed: epykit benchmarks at .../benchmark/
cd .../epykit_vs_allPackages\(simulated_approxData\)

# 1. Ground truth
uv run python ground_truth/make_truth.py

# 2. Run all scenarios
uv run python scripts/run_epykit_dmc_coverage.py
uv run python scripts/run_epykit_dmc_replicate.py
uv run python scripts/run_epykit_dmr_coverage.py

# 3. Evaluate + join baselines
uv run python scripts/load_figure_estimates.py
uv run python scripts/evaluate.py --scenario all

# 4. Figures
uv run python scripts/make_figures.py

# 5. Audit baseline transcriptions (should print "0 mismatches")
uv run python scripts/audit_baselines.py
```

Total runtime: ~5 minutes on a single CPU core.

## 8. Reproduction recipe (Study 2)

```bash
cd .../epykit_vs_methylkit\(simulated_realRun\)

# 0. Ground truth (cached after first run)
uv run python scripts/_make_truth.py

# 1. Full grid through the perf harness
uv run python scripts/run_bench_all.py \
    --tools epykit methylkit \
    --scenarios dmc_coverage dmc_replicate dmr_coverage

# 2. Score + render
uv run python scripts/compare_tools.py
```

Selective re-runs: `--only K` (single grid value), `--tools epykit` or
`--tools methylkit` (single tool).

## 9. Reproduction recipe (Study 3)

```bash
cd .../epykit_vs_methylkit\(GSE263850\)

# epykit pipeline
uv run python epykit_analysis.py

# methylKit pipeline (separately, R)
Rscript methylkit_pipeline.R   # produced by Study 3 driver

# Compare
uv run python comparison.py
```

epykit analysis runtime: ~18 minutes. methylKit: ~3 h 37 min on this
hardware.

---

## 10. Cross-study summary figures

A small set of cross-benchmark figures was generated as part of
consolidating the three studies into this FINAL_REPORT folder:

* **S1.** Wall-clock runtime, epykit vs methylKit, across all three
  studies. Log-y bar chart.
* **S2.** Best-engine TPR vs methylKit TPR across coverage levels for
  Studies 1 and 2 — a sanity check that both simulated studies tell the
  same story.
* **S3.** Effect-size scatter on Study 3 alongside an agreement summary
  (Jaccard, Pearson, recall) heatmap.

Generated by [scripts/make_summary_figures.py](../scripts/make_summary_figures.py);
PNGs in [figures/summary/](../figures/summary/).

---

## A. epykit `dmr_chain_merge` paper-matched parameter set

The Study 3 chain_merge replication maps DSS::callDMR parameters to
epykit's `ep.tl.dmr(method='chain_merge', ...)` one-for-one:

| DSS parameter | epykit parameter | Value |
|---|---|---|
| `p.threshold = 1e-5` | `alpha=1e-5` | per-CpG p-value threshold |
| `delta = 0` | `min_abs_meth_diff=0.0` | no Δβ cutoff |
| `minlen = 50` | `minlen_bp=50` | min DMR length |
| `minCG = 3` | `min_cpgs=3` | min CpGs per DMR |
| `dis.merge = 100` | `dis_merge_bp=100` | max gap to merge adjacent chains |
| `pct.sig = 0.5` | `pct_sig=0.5` | min fraction of CpGs in DMR significant |
| (post-hoc) | `min_mean_qvalue=0.05` | BH-q filter on the DMR table |

Per-CpG test (matches DSS::DMLfit.multiFactor(smoothing=TRUE)):

```python
ep.tl.dmc(md,
    test='lr', dispersion='site',
    smoothing=True, smoothing_span_bp=500,
)
```

The smoother is a uniform-box average over CpGs within
`±smoothing_span_bp//2` (250 bp on each side), implemented in
[`_smooth_sample_counts_box`](../../epykit2/epykit2/src/epykit/dmc.py).
This matches DSS's `smooth.chr(..., method="avg")` exactly with the
default `smoothing.span = 500`.

Script: [run_chain_merge_replication.py](../scripts/run_chain_merge_replication.py).
Output: [data/study3/chain_merge/](../data/study3/chain_merge/).
Parameter sweep: [sweep_dis_merge.py](../scripts/sweep_dis_merge.py) →
[chain_merge_dis_merge_sweep/](../data/study3/chain_merge_dis_merge_sweep/).

## B. DSS-from-scratch replication

Script: [run_dss_replication.R](../scripts/run_dss_replication.R)
(driven by [run_dss_replication.py](../scripts/run_dss_replication.py)
for psutil-based resource sampling, then
[resume_dss_from_dmltest.R](../scripts/resume_dss_from_dmltest.R)
+ [resume_dss_from_dmltest.py](../scripts/resume_dss_from_dmltest.py)
to finish the post-callDMR steps after a column-name bug in the
initial run).

Pipeline:

```r
samplesheet  # chrom, sample_id, group, path; n=6 (3 KO + 3 WT)
# 1. Read each 12-col BED (cols 1-3 chrom/start, 10-11 M/T combined-strand)
#    keeping CpGs with coverage >= 5
# 2. Intersect on (chrom, pos) across all 6 samples
# 3. Build BSseq(chr, pos+1, M, Cov)  -- BED is 0-based, BSseq is 1-based
# 4. DMLfit.multiFactor(BSobj, design=data.frame(group=...),
#                       formula=~group, smoothing=TRUE)
#    -- this is the slow step (2,044 s wall / 2,029 s CPU on the
#       15.4 GB Windows host, 22.0 M CpGs; single-threaded)
# 5. DMLtest.multiFactor(fit, coef='groupWT')
# 6. callDMR(test_res, delta=0, p.threshold=1e-5, minlen=50, minCG=3,
#            dis.merge=100, pct.sig=0.5)
```

Note: `callDMR` on multifactor test output returns
`[chr, start, end, length, nCG, areaStat]` — it does NOT include
`diff.Methy` / `meanMethy1` / `meanMethy2` (those require the
two-group `DMLtest` path). Per-DMR per-group mean methylation
(matching the paper's Supp Table 5 `diff.meth_mean` column) is derived
in the resume script directly from the 6 per-sample BEDs via
`data.table::foverlaps`.

Versions used (Windows host):

| Component | Version |
|---|---|
| R | 4.5.0 (2025-04-11) |
| DSS | **2.58.0** (Bioconductor; per `data/study3/dss/dss_session_info.txt`) |
| bsseq | **1.46.0** (Bioconductor) |
| data.table | 1.x |
| optparse | 1.x |

Full session info: [dss_session_info.txt](../data/study3/dss/dss_session_info.txt).
Resource trace: [resources.csv](../data/study3/dss/resources.csv) +
[resources.json](../data/study3/dss/resources.json).
Cached per-CpG DMLtest table: `dmr_dss/dmltest_per_cpg.tsv.gz` (461 MB,
22 M rows × 5 cols); allows resume of just the callDMR + annotation
steps in ~ 8 minutes without redoing the 34-minute DMLfit.

## C. HOMER-equivalent annotation (UCSC refGene)

HOMER is not available on the Windows host; we re-implemented HOMER's
`annotatePeaks.pl` annotation rules in pure Python (epykit) and R
(DSS pipeline) using HOMER's bundled UCSC refGene catalog
(`refGene.txt.gz`, hg38).

Feature priority (HOMER):

```
promoter-TSS (0) > TTS (1) > 5'UTR (2) > 3'UTR (3)
> exon (4) > intron (5) > non-coding (6) > intergenic (7)
```

Promoter window: TSS −1000 / +100 bp (strand-aware).
TTS window: TTS −100 / +1000 bp (strand-aware).
Coding vs non-coding distinction: refGene `acc` starts with `NM_`.

For each DMR midpoint, the feature with the lowest priority code
that contains the midpoint becomes the assigned `feature_type`.
Nearest TSS is computed strand-aware: `distance = midpoint − tss`
on `+`-strand, negated on `−`-strand.

100 kb DMR-gene linkage (paper Methods): for each DMR midpoint,
all genes whose canonical TSS (most upstream per gene, per strand)
is within ±100 kb of the midpoint are linked. Outputs in
[dmr_gene_links_100kb.csv](../data/study3/chain_merge/dmr_gene_links_100kb.csv)
(and parallel files under [dss/](../data/study3/dss/) and
[chain_merge_dis_merge_sweep/dis_merge_*/](../data/study3/chain_merge_dis_merge_sweep/)).

Annotation reimplementation:
[compare_homer_refseq.py](../../epykit2/GSE263850_RAW/compare_homer_refseq.py)
(epykit's `ep.tl.annotate(refgene=...)` also uses this catalog); R-side
in [run_dss_replication.R](../scripts/run_dss_replication.R) and
[resume_dss_from_dmltest.R](../scripts/resume_dss_from_dmltest.R).

## D. Comparison metric definitions

Used throughout §3.10–§3.17 of the technical report:

| Metric | Definition |
|---|---|
| any-bp overlap | Two DMRs from different callers overlap if `min(end_a, end_b) > max(start_a, start_b)`. |
| Jaccard | `intersection_bp / union_bp` for a single matched pair. Range [0, 1]. |
| Reciprocal-overlap fraction | `intersection_bp / min(length_a, length_b)`. Range [0, 1]. |
| Recall (any-bp) of target by query | Fraction of target-set DMRs with ≥ 1 any-bp-overlapping query DMR. |
| Recall (J ≥ t) | Fraction of target-set DMRs with at least one query DMR at Jaccard ≥ t. Stricter than any-bp. |
| Precision (any-bp) | Fraction of query-set DMRs with ≥ 1 any-bp-overlapping target DMR. |
| Direction agreement on matched | For each matched (query, target) pair, do their `hyper`/`hypo` directions agree? Counted on matched pairs only. |
| Nearest-TSS gene capture | Fraction of target-set gene names (`nearest_tss_gene` column) recovered in the query-set gene list. |
| 100 kb gene capture | Fraction of target-set gene names recovered in the query-set's 100 kb-rule gene-link table. |

Best-match rule: for each row in the query set, the best target match
is the target row with the highest Jaccard (or 0 if no any-bp
overlap). Implemented as a sorted-sweep over per-chromosome interval
indices in [compare_chain_merge_to_paper.py](../scripts/compare_chain_merge_to_paper.py)
and [compare_epykit_to_dss.py](../scripts/compare_epykit_to_dss.py).

## E. Three-way enrichment methodology

Gene lists (Methods §B above for 100 kb-rule construction):

| List | Source | n genes |
|---|---|---:|
| paper Table 5 (705 unique) | Supp Table 5 `Gene.Name` column | 705 |
| methylKit-tile | HOMER-equivalent nearest TSS of every methylKit DMR | 2,111 |
| ek-chain_merge-100 (100 kb rule) | dmr_gene_links_100kb.csv (post-rerun, 852 DMRs) | 1,290 |
| ek-chain_merge-250 (100 kb rule) | dmr_gene_links_100kb.csv (post-rerun sweep, 1,139 DMRs) | 1,645 |
| DSS-from-scratch (100 kb rule) | DSS dmr_gene_links_100kb.csv | 1,467 |

Backend: Enrichr REST API (`maayanlab.cloud/Enrichr`), no auth.
Libraries:

* `Reactome_2022` (1,818 pathways; closest to paper's
  ShinyGO Curated.Reactome but with full-library BH correction)
* `KEGG_2021_Human` (KEGG-side reproduction of GPCR / Gα(i) signals)
* `GO_Molecular_Function_2023` (paper Panel E equivalent)

Each gene list is submitted as a fresh user list; top-20 terms by raw
p-value per library are retained with overlap-gene lists. Paper Panel D
+ KEGG-equivalent keyword matches are flagged in the output JSON.

Note: the paper used ShinyGO V0.77 + Curated.Reactome with the paper's
own 23,590-gene "expressed gene" background. Enrichr does not allow
the same background substitution via API; this limits direct FDR
comparison. Term-rank comparisons remain interpretable; absolute
significance does not transfer.

Script: [run_enrichment_three_way.py](../scripts/run_enrichment_three_way.py).
Output: [enrichment_three_way.json](../data/study3/comparisons/enrichment_three_way.json) +
[summary.md](../data/study3/comparisons/enrichment_three_way_summary.md).
