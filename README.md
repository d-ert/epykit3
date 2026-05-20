# epykit

A Python-native WGBS methylation analysis pipeline built on Parquet partitioning and lazy I/O — from Bismark `.cov` (or MethylDackel `.bedGraph`) files to differentially methylated cytosines (DMCs), differentially methylated regions (DMRs), gene-feature annotation, and shareable HTML reports.

epykit ingests Bismark / MethylDackel coverage output into a partitioned Parquet **methylstore** and runs the whole downstream analysis (QC → filtering → DMC → DMR → annotation → plotting → report) over that store with [polars](https://pola.rs) and lazy I/O. The Python API is organised in a scanpy-style `pp` / `tl` / `pl` namespace; a CLI mirrors the same operations for scripting.

> **Status:** version 0.7.0 (unreleased), pre-1.0. API may change. MIT licensed.

[Documentation](https://d-ert.github.io/epykit2/) | [Changelog](CHANGELOG.md)


---

## Highlights

- **Partitioned Parquet methylstore.** Per-chromosome, per-sample columnar storage — never load a whole genome into RAM. DMC results follow the same convention: `process_chromosomes_dmc(..., return_store=True)` returns a `DMCStore` handle backed by per-chromosome parquet files under `<methylstore>/.cache/dmc/<test>/`, and `apply_multiple_testing_correction` / `call_dmr_sliding_window` stream from it so peak memory stays at O(largest chromosome) on whole-genome inputs (~22 M CpGs).
- **Statistical engines.** Per-CpG DMC tests cover the quasi-binomial likelihood-ratio (`lr`, the default at n ≥ 2; closed-form with McCullagh-Nelder dispersion), Pearson score, full IRLS binomial GLM with covariates, Welch t on logit(β), Welch t on raw β (`welch_t`), a true quasi-binomial LRT (`bb_lr`), Cochran-Mantel-Haenszel, and pooled Fisher exact. Every test surfaces 95 % Wald CIs on `meth_diff`. Permutation empirical FDR is available end-to-end: `tl.dmc(..., empirical_fdr=True)` and `tl.dmr(..., empirical_fdr=True)` shuffle labels, re-run the engine, and add `empirical_pvalue` / `empirical_qvalue` columns.
- **`lr+` power stack (since 0.7.1).** Four opt-in enhancements to the `lr` engine that close the asymptotic-quasi-binomial gap to methylKit / RADMeth / DSS at low coverage and small cohorts: empirical-Bayes dispersion shrinkage (`dispersion="eb"`), sign-aware RADMeth-style Stouffer combiner over neighbouring CpGs (`neighbour_combine=True`, default ±200 bp), separation-aware Fisher fallback for near-perfect-separation 2 × 2 tables (`sep_fallback=True`), and Storey/two-stage BH q-values (`fdr_method="fdr_tsbh"` or `"fdr_storey"`). All four enabled together push TPR ≥ 0.999 on the Piao et al. 2021 simulated benchmark at every coverage ≥ 10× and every replicate count ≥ 4 while keeping FPR strictly tighter than every R baseline. See `comparison_test/benchmark/REPORT.md` for the head-to-head.
- **Multi-group & covariate contrasts.** `tl.dmc(formula="~ group + age", contrast="group")` runs a joint F-test across factor levels; `contrast="age"` runs a Wald test on a continuous covariate as the primary effect.
- **Four DMR engines plus permutation FDR.** Tile-based (read-pooled, default), per-CpG sliding-window with signed Stouffer's combining, HMM segmentation over `meth_diff`, and a DSS-compatible **chain-merge** caller (`tl.dmr(method="chain_merge", preset="strict" | "default" | "permissive")`) that mirrors DSS `callDMR` semantics. `tl.dmr(..., empirical_fdr=True, n_perm=100)` re-runs the engine on shuffled labels and reports empirical p- and q-values. `tl.diagnose_dmr_calling(md, reference_dmrs)` buckets unrecovered reference DMRs into actionable categories (coverage loss vs. weak test vs. structural filter) for triage.
- **Differential variability.** `tl.dvc(md)` finds CpGs whose between-replicate variance differs between groups even when the means don't — the iEVORA signal that mean-based DMC misses.
- **Clinical / cohort QC.** Opt-in `qc.sex_check` (chrX mean β), `qc.contamination_estimate` (β-distribution bimodality), `qc.sample_correlation` (sample-swap detection), and `qc.power` (sample-size calculator). Bisulfite conversion rate is reported (CHH context, dashboard + MultiQC) but **not applied** to per-CpG counts — matching `bsseq` / `methylKit` defaults. A poorly converted library should be re-prepped, not papered over with a multiplicative count adjustment.
- **Replicate-aware throughout.** Per-site `min_samples_treatment` / `min_samples_control` guards, per-site or chromosome-level McCullagh-Nelder dispersion, optional covariate design matrices via Wilkinson formulas.
- **Annotation.** Gene features (promoter / 5'UTR / exon / intron / 3'UTR) from GENCODE / Ensembl **GTF** or UCSC **`refGene.txt`** (HOMER's default catalog), plus CpG-island context (island / shore / shelf / open-sea). Opt-in `gene_type_filter="protein_coding"` drops lincRNAs / pseudogenes. `multi_annotation=True` (default) adds annotatr-style `nearest_tss_gene` / `nearest_tss_distance` and one-to-many `all_overlapping_genes` / `all_overlapping_features` columns so a site that's intronic for one gene AND in another gene's promoter window is faithfully represented.
- **Visualisation pack.** matplotlib volcano, MA, Manhattan, coverage histogram, methylation heatmap, PCA, UMAP, sample-correlation heatmap, QC dashboard, DMR boxplot, genomic-context bar, CpG-island pie, TSS metaplot — plus Plotly twins for the HTML report.
- **Interop.** Self-contained HTML report (`md.report(out.html)`), AnnData (`md.to_anndata()`), MuData (`md.to_mudata()`), methylKit-compatible tabix tables (`md.to_methylkit_tabix(dir)`), MultiQC custom-content JSON (`ep.report_multiqc(md, dir)`), nf-core/methylseq QC ingestion (`ep.read_nfcore_methylseq_qc(...)`).
- **CLI.** `epykit convert | filter | dmc | dmr | annotate | qc-report | smooth | report | aggregate-regions | export` — every stage scriptable from the shell.

---

## Installation

Requires Python ≥ 3.9.

```bash
# from the repo checkout
pip install -e .

# or with uv
uv pip install -e .

# dev install
pip install -e ".[dev]"

# full feature install (report + export + anndata + viz)
pip install -e ".[all]"
```

Core dependencies: `polars`, `pyarrow`, `numpy`, `scipy`, `numba`, `bioframe`, `pyfaidx`, `statsmodels`, `patsy`, `psutil`, `scikit-learn`, `matplotlib`, `seaborn`. Optional extras: `report` (Jinja2 + Plotly), `export` (pyBigWig), `anndata` (anndata + mudata), `viz` (umap-learn), `methylkit` (pysam, for tabix indexing on Linux/macOS). The CLI is installed as the `epykit` console script.

---

## Quickstart

### 1. Samplesheet

epykit reads a CSV with three required columns. Any extra columns are kept on `md.obs` and are available as GLM covariates.

```csv
sample_id,group,path
ctrl_1,control,raw_data/bismark/ctrl_1.bismark.cov.gz
ctrl_2,control,raw_data/bismark/ctrl_2.bismark.cov.gz
cd55_1,cd55,raw_data/bismark/cd55_1.bismark.cov.gz
cd55_2,cd55,raw_data/bismark/cd55_2.bismark.cov.gz
```

### 2. End-to-end analysis (Python API)

```python
import epykit as ep
import polars as pl

# Ingest: converts each .cov to per-chromosome Parquet under
# methyl_store/.cache/raw/ and returns a MethylData object.
md = ep.read_bismark(
    "samplesheet.csv",
    treatment_group="cd55",
    control_group="control",
    assembly="hg38",
    store_dir="methyl_store",
)
print(md)

# Preprocessing (pp.*) — each step repoints md.store at a cached store.
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)
ep.pp.normalize_coverage(md, method="median")
ep.pp.unite(md, type="intersect")        # or "union" + min_samples_* guards

# Tools (tl.*) — populate md.obs / md.varm / md.uns.
ep.tl.qc(md)                              # populates md.obs with QC metrics
ep.tl.dmc(md, test="auto")                # md.varm["dmc_lr"] (n≥2) or "dmc_fisher" (n=1)
ep.tl.dmr(md, tile_size_bp=500, min_cpgs_per_tile=5)   # tile-based, default

# Inspect.
total = len(md.dmc)
sig   = md.dmc.filter(pl.col("qvalue") < 0.05).height
print(f"DMCs: {sig:,} / {total:,} ({100 * sig / total:.2f}%)")
print(md.uns["dmr"].filter(pl.col("qvalue") < 0.05))

# Annotation.
ep.tl.annotate(
    md,
    gtf="raw_data/gencode.v49.annotation.gtf",
    cpg_islands="raw_data/hg38_cpg_islands.bed",
)

# Persist the analysis and emit a shareable HTML report.
md.save("cd55_analysis")
md.report("cd55_report.html")             # interactive Plotly + Jinja2

# Plotting (pl.*) — works on a freshly loaded MethylData.
md = ep.load("methyl_store/results/cd55_analysis")
ep.pl.volcano(md,              save="volcano")
ep.pl.ma_plot(md,              save="ma")
ep.pl.manhattan(md,            save="manhattan")
ep.pl.coverage_histogram(md,   save="coverage_hist")
ep.pl.pca(md,                  save="pca")
ep.pl.umap(md,                 save="umap")
ep.pl.qc_dashboard(md,         save="qc_dashboard")
ep.pl.dmr_boxplot(md, top_n=6, save="dmr_boxplot")
```

`scratch_plan2.py` at the repo root exercises every Plan 2 feature on real Bismark data and is the canonical real-data smoke test.

### 3. Covariate-adjusted analysis

When `md.obs` has additional columns (sex, batch, age, donor, …), pass them through a formula. The engine uses a binomial GLM internally:

```python
# Covariate-adjusted binary contrast
ep.tl.dmc(md, formula="~ group + donor", contrast="group")

# Multi-group joint F-test (3+ levels)
ep.tl.dmc(md, formula="~ group", contrast="group")

# Continuous covariate as the primary effect
ep.tl.dmc(md, formula="~ age", contrast="age")
```

### 4. Permutation-based empirical FDR for DMRs

Asymptotic q-values can be miscalibrated on small-n WGBS. For trustworthy DMR-level inference:

```python
ep.tl.dmr(md, method="tile", empirical_fdr=True, n_perm=100, perm_seed=42)
# md.uns["dmr"] now carries empirical_pvalue / empirical_qvalue columns
```

### 5. The `lr+` power stack — closing the gap to methylKit / RADMeth / DSS

The default `lr` test is a *conservative* quasi-binomial likelihood ratio: it controls FPR aggressively, sometimes at the cost of 5-16 pp of TPR at low coverage on strong-effect sites or at 2 v 2 replicate designs (see `comparison_test/benchmark/REPORT.md` Section 4). Four targeted opt-in enhancements close that gap without inflating FPR.

```python
# Drop-in "lr+" recipe: enable all four enhancements at once
ep.tl.dmc(
    md,
    test="lr",
    fdr_method="fdr_tsbh",       # Storey/two-stage BH; estimates pi0 from data
    neighbour_combine=True,       # signed-Stouffer combiner across nearby CpGs
    neighbour_bp=200,             # half-window for the combiner
    sep_fallback=True,            # Fisher fallback for near-separated 2x2 tables
    sep_threshold=0.9,            # |meth_diff| threshold for triggering Fisher
    dispersion="eb",              # empirical-Bayes shrinkage of per-site phi
)
```

Each enhancement is independently controllable so you can adopt the ones that suit your data:

| Option | What it does | When it pays off |
|---|---|---|
| `fdr_method="fdr_tsbh"` (or `"fdr_storey"`) | Two-stage BH / Storey q-values; estimates the proportion of true nulls (π₀) from the p-value histogram and uses it to scale BH. Reduces to plain BH when π₀ = 1. | Whenever a meaningful fraction of the genome carries real signal (typical WGBS). Cheap. |
| `neighbour_combine=True`, `neighbour_bp=200` | Signed-Stouffer Z combiner over CpGs within ±`neighbour_bp` bp. Gates on ≥ 60 % sign agreement and raw p < 0.5 at the focal site so spatially isolated false positives are not amplified. Output `pvalue` becomes the combined p; the raw p is preserved as `pvalue_raw`. | Real WGBS data, where adjacent CpGs in true DMRs share biological signal. Biggest TPR contributor of the four. |
| `sep_fallback=True`, `sep_threshold=0.9` | When the LR p-value fails to reject (p > 0.05) but the observed `|meth_diff|` exceeds `sep_threshold`, re-test with `scipy.stats.fisher_exact` on pooled counts and take the more powerful of the two. Never inflates p. | Very-low-coverage data where strong-effect sites generate near-perfect-separation 2 × 2 tables. |
| `dispersion="eb"` | Empirical-Bayes shrinkage of the per-site quasi-binomial dispersion toward a chromosome-wide inverse-Gamma prior whose pseudo-df is estimated from the data (vs. the fixed pseudo-df = 4 of `"shrink"` mode). | Real cohorts with genomic-context-dependent dispersion (CpG islands vs. gene bodies vs. repeats). No-op when there is no true overdispersion. |

**Result columns when `neighbour_combine=True`.** The `pvalue` / `qvalue` columns now refer to the *combined* p-value; the raw per-CpG values are preserved as `pvalue_raw` / `qvalue_raw`. Two additional columns are added: `pvalue_combined` (identity-equal to `pvalue`, kept for audit clarity) and `pvalue_combined_n_neighbours` (how many CpGs contributed to the combination at each site).

**Footgun.** Downstream code that expects `pvalue` to be the raw per-CpG p-value should switch to `pvalue_raw` after enabling `neighbour_combine`, or skip the combiner.

**Recipe matrix.**

| Your situation | Recommended setup |
|---|---|
| Quick exploratory run | `ep.tl.dmc(md)` — old default, ~0.4 s on 100 K sites |
| Whole-genome WGBS, best power | `lr+` recipe above |
| Strict-FDR validation cohort | `fdr_method="fdr_bh"`, `neighbour_combine=False` |
| Very-low-coverage data (< 10×) | add `sep_fallback=True` and `neighbour_combine=True` |
| Cohort with strong context-dependent dispersion | add `dispersion="eb"` |
| Reproducing a published `lr` baseline | leave all four off |

### 6. Clinical / cohort QC

```python
ep.tl.qc(
    md,
    run_sex_check=True,           # infers sex from chrX β; flags swaps
    run_contamination=True,        # β-distribution bimodality score
    run_sample_correlation=True,   # all-vs-all sample correlation
)
ep.qc.power(meth_diff=0.10, coverage=20, power=0.80)   # minimum n per group
```

---

## CLI

The `epykit` script mirrors the Python pipeline. Every subcommand takes `--methylstore` (the partitioned Parquet directory) and writes Parquet output unless otherwise noted.

| Subcommand          | Purpose |
|---------------------|---------|
| `convert`           | Bismark `.cov[.gz]` → partitioned Parquet |
| `filter`            | Coverage / blacklist filtering |
| `summary`           | Per-sample summary statistics |
| `dmc`               | Per-CpG differential methylation. `--test {lr,score,glm,logit_t,welch_t,bb_lr,cmh,fisher}`, plus `--formula` / `--contrast` / `--covariates` for covariate-adjusted and multi-group designs. The `lr+` power-stack options (`fdr_method`, `neighbour_combine`, `sep_fallback`, `dispersion="eb"`) are currently Python-API-only; CLI flags are pending. |
| `dmr`               | DMR calling — `--method tile` (default) or `--method sliding_window`. Supports `--empirical-fdr --n-perm N`. |
| `annotate`          | Add gene-feature (`--gtf`) and CpG-island (`--cpg-islands`) annotation. |
| `qc-report`         | QC + coverage uniformity report. |
| `smooth`            | Gaussian-kernel β smoothing along the genome. |
| `report`            | Render a self-contained interactive HTML report from a saved analysis. |
| `aggregate-regions` | Aggregate per-CpG counts to user-supplied BED regions. |
| `export`            | Sub-commands: `bedgraph`, `bigwig`, `dmcs-bed`, `dmrs-bed`, `mudata`, `methylkit-tabix`, `multiqc`. |

Run `epykit <subcommand> --help` for the full flag list.

---

## Input formats

- **Bismark `.cov` / `.cov.gz`** — 6-column 0-based BED-like:
  `chrom`, `start`, `end`, `methylation_percent`, `count_methylated`, `count_unmethylated`. Read with `ep.read_bismark(...)` or `epykit convert --format bismark`.
- **MethylDackel `.bedGraph` / `.bedGraph.gz`** — same 6 columns as Bismark with a single `track type="bedGraph" ...` header line that is skipped automatically. Read with `ep.read_methyldackel(...)` or `epykit convert --format methyldackel`.
- **Samplesheet** (CSV) — required columns `sample_id`, `group`, `path`. Any extra column is preserved on `md.obs` and can be referenced as a GLM covariate.
- **GTF** — Ensembl / GENCODE / UCSC; gene features are extracted via [bioframe](https://github.com/open2c/bioframe). `gene_type` (GENCODE) and `gene_biotype` (Ensembl) are both honoured.
- **UCSC `refGene.txt[.gz]`** — HOMER's default gene catalog. Pass `ep.tl.annotate(md, refgene=...)` (Python API; not yet wired into `epykit annotate`). Schema-compatible with the GTF path.
- **CpG-island BED** — UCSC `cpgIslandExt` 4-column BED.

---

## Output layout

`read_bismark(..., store_dir="methyl_store")` produces:

```
methyl_store/
├── .cache/
│   ├── raw/                      # converted .cov → Parquet
│   │   ├── sample=ctrl_1/
│   │   │   └── chrom=chr1/part-0.parquet
│   │   └── sample=cd55_1/
│   │       └── chrom=chr1/part-0.parquet
│   ├── filtered/                 # after pp.filter_coverage
│   ├── normalized/               # after pp.normalize_coverage
│   └── dmc/
│       └── lr/                   # after tl.dmc(test="lr")
│           ├── .epykit_dmc_manifest.json
│           ├── chrom=chr1.parquet
│           └── chrom=chr2.parquet
└── results/
    └── cd55_analysis/            # md.save() target
        ├── obs.parquet
        ├── varm_dmc_lr_annotated.parquet
        ├── uns_dmr.parquet
        └── methyldata.json
```

DMC frames carry: `chrom`, `pos`, `strand`, `n_case`, `n_control`, `mean_beta_case`, `mean_beta_control`, `meth_diff`, `meth_diff_ci_lo`, `meth_diff_ci_hi`, `pvalue`, `qvalue`, `log2_odds_ratio`, plus per-test extras (`coef_treatment` / `coef_se` for GLM and `bb_lr`; `f_stat` / `df1` / `df2` / per-level `mean_beta_<level>` / `meth_diff_max` for multi-group contrasts) and, after `tl.annotate`, `feature_type` / `gene_id` / `cpg_context`. Tile-DMR frames add `start`, `end`, `n_cpgs`, `dmr_type ∈ {hyper, hypo, mixed}`; permutation FDR adds `empirical_pvalue` / `empirical_qvalue`.

---

## Module map

| Module             | Role |
|--------------------|------|
| `methyldata.py`    | `MethylData` dataclass — `obs`, `store`, `varm`, `uns`; `.dmc` / `.treatment_ids` / `.control_ids` properties; `save()` / `load()`; `region_beta()` per-region query |
| `io.py`            | `read_bismark`, `read_nfcore_methylseq`, `load` |
| `convert.py`       | `.cov` → partitioned Parquet |
| `filter.py`        | Coverage filter, coverage normalisation, blacklist intersect |
| `pp.py`            | Preprocessing wrappers (`filter_coverage`, `normalize_coverage`, `unite`, `smooth`, `aggregate_regions`) |
| `dmc.py`           | Streaming per-CpG accumulators + statistical engines (`lr`, `score`, `glm`, `logit_t`, `welch_t`, `bb_lr`, `cmh`, `fisher`), BH correction |
| `_dmc_store.py`    | `DMCStore` handle — persistent per-chromosome DMC parquet directory + manifest; lets BH and sliding-window DMR stream from disk so peak memory is O(largest chrom), not O(genome) |
| `dmr.py`           | `call_dmr_tile_based`, `call_dmr_sliding_window`, `empirical_fdr_for_dmr`, `smooth_methylation_gaussian` |
| `dvc.py`           | Differentially Variable CpG calling (iEVORA-style) |
| `annotate.py`      | `annotate_features` (GTF), `annotate_cpg_islands` (island / shore / shelf / open-sea) |
| `qc.py`            | `bisulfite_conversion_rate`, `global_methylation_report`, `coverage_uniformity`, `sex_check`, `contamination_estimate`, `sample_correlation`, `power` |
| `tl.py`            | High-level orchestrators: `tl.qc`, `tl.dmc`, `tl.dmr`, `tl.dvc`, `tl.annotate` |
| `pl/`              | Plotting — `qc`, `differential`, `genomic`, `clustering`, `metaplot`, `embedding`, `correlation`, `dashboard`, `dmr_boxplot`, plus Plotly twins |
| `report.py`        | Self-contained interactive HTML report (Jinja2 + Plotly) |
| `export.py`        | BedGraph / BigWig / DMC-BED / DMR-BED export |
| `anndata_io.py`    | AnnData export |
| `mudata_io.py`     | MuData export (multi-omics bundling) |
| `methylkit_io.py`  | methylKit-compatible tabix tables |
| `multiqc_export.py`| MultiQC custom-content JSON emitter |
| `nfcore_qc.py`     | nf-core/methylseq run-dir QC ingestion |
| `cli.py`           | `epykit` CLI entry point |
| `_glm.py`          | Wilkinson formula → design matrix, batched IRLS binomial GLM, Wald test on contrasts |
| `_style.py`        | Shared matplotlib palette / theme |

---

## Tests

```bash
pip install -e ".[dev]"
pytest tests/
```

---

## License

MIT.
