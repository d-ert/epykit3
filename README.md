<p align="center">
  <img src="docs/assets/logo.png" alt="epykit logo" width="340">
</p>

<p align="center">
  <strong>A Python-native WGBS methylation analysis pipeline.</strong><br>
  Bismark / MethylDackel coverage files → DMC → DMR → annotation → shareable HTML reports,<br>
  over a partitioned-Parquet methylstore with lazy I/O.
</p>

<p align="center">
  <a href="https://github.com/d-ert/epykit3/actions/workflows/test.yml"><img src="https://github.com/d-ert/epykit3/actions/workflows/test.yml/badge.svg" alt="CI"></a>
  <img src="https://img.shields.io/badge/version-1.1.0-success" alt="version 1.1.0">
  <img src="https://img.shields.io/badge/python-3.10%E2%80%933.13-blue" alt="Python 3.10–3.13">
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-yellow.svg" alt="License: MIT"></a>
  <a href="https://d-ert.github.io/epykit3/"><img src="https://img.shields.io/badge/docs-d--ert.github.io-informational" alt="Documentation"></a>
</p>

<p align="center">
  <a href="https://d-ert.github.io/epykit3/">Documentation</a> ·
  <a href="CHANGELOG.md">Changelog</a> ·
  <a href="benchmark/paper/report/REPORT.md">Benchmark</a>
</p>

---

epykit ingests Bismark / MethylDackel coverage output into a partitioned Parquet **methylstore** and runs the whole downstream analysis (QC → filtering → DMC → DMR → annotation → plotting → report) over that store with [polars](https://pola.rs) and lazy I/O. Whole-genome data (~22 M CpGs) is never held in RAM as a single frame — every step streams over the partition tree. The Python API is organised in a scanpy-style `pp` / `tl` / `pl` namespace; a CLI mirrors the same operations for scripting.

> **Status:** version 1.1.0 — stable API. MIT licensed. Linux / macOS / Windows.

---

## Why epykit

The WGBS analysis ecosystem is fragmented across R/Bioconductor (methylKit, DSS, BSmooth, methylSig, BiSeq) and command-line tools (RADMeth, metilene), while most modern single-cell and multi-omics workflows live in Python. epykit brings a maintained panel of DMC/DMR engines to Python with three design choices that set it apart:

- **Parquet-first, RAM-light.** Per-chromosome, per-sample columnar storage means whole-genome cohorts run on a laptop — DMC/DMR *computation* keeps peak memory at O(largest chromosome), not O(genome). The default `tl.dmc(md)` then assembles the full per-CpG result onto `md.varm` (~700 MB at 22 M CpGs) for plotting/report/export; pass `tl.dmc(md, materialize=False)` to stay O(largest chromosome) end-to-end via the streaming `DMCStore`.
- **A panel, not a single test.** Four per-CpG engines and four DMR callers, all emitting one canonical schema, so you can switch methods without rewriting downstream code.
- **Familiar, scanpy-style API.** `ep.pp.*` / `ep.tl.*` / `ep.pl.*` mirrors the conventions Python bioinformaticians already know, with a CLI that exposes the same pipeline for scripting.

---

## Highlights

- **Partitioned Parquet methylstore.** Per-chromosome, per-sample columnar storage — never load a whole genome into RAM. DMC results follow the same convention: `tl.dmc` writes per-chromosome parquet files under `<methylstore>/.cache/dmc/<test>/` and exposes a streaming `DMCStore` handle (`md.dmc_store`), so DMC computation and the downstream sliding-window DMR caller stay at O(largest chromosome) on whole-genome inputs (~22 M CpGs). By default the per-CpG table is then materialised onto `md.varm` for plotting/report/export; pass `tl.dmc(md, materialize=False)` to keep peak memory O(largest chromosome) end-to-end (`md.dmc` then materialises on demand). Advanced users can drive the streaming engines directly via `from epykit.dmc import process_chromosomes_dmc, apply_multiple_testing_correction`.
- **Statistical engines.** Four per-CpG DMC backends: `lr` (quasi-binomial likelihood-ratio, the default at n ≥ 2; closed-form with McCullagh-Nelder dispersion), `welch_t` (Welch t on raw β), `fisher` (pooled Fisher exact, n = 1 fallback), and `glm` (full IRLS binomial GLM with covariates). `auto` resolves to `fisher` at n < 2 and `lr` at n ≥ 2. Every test surfaces 95 % Wald CIs on `meth_diff`. Permutation empirical FDR shuffles labels, re-runs the engine, and adds `empirical_pvalue` / `empirical_qvalue` columns: `tl.dmc(..., empirical_fdr=True)` on the binary treatment/control path, and `tl.dmr(method="tile", empirical_fdr=True)` for regions (the tile caller is the only DMR method with the permutation harness wired in; the others raise `NotImplementedError`).
- **Multi-group & covariate contrasts.** `tl.dmc(formula="~ group + age", contrast="group")` runs a joint F-test across factor levels; `contrast="age"` runs a Wald test on a continuous covariate as the primary effect.
- **Four DMR engines plus permutation FDR.** A DSS-compatible **chain-merge** caller (`tl.dmr(method="chain_merge", preset="strict" | "default" | "permissive")`, the default) that mirrors DSS `callDMR` semantics, plus tile-based (read-pooled) aggregation, per-CpG sliding-window with signed Stouffer's combining, and rule-based 3-state segmentation over `meth_diff` (`method="segment"`). `tl.dmr(method="tile", empirical_fdr=True, n_perm=100)` re-runs the tile engine on shuffled labels and reports empirical p- and q-values (tile only). `tl.diagnose_dmr_calling(md, reference_dmrs)` buckets unrecovered reference DMRs into actionable categories (coverage loss vs. weak test vs. structural filter) for triage.
- **Differential variability.** `tl.dvc(md)` finds CpGs whose between-replicate variance differs between groups even when the means don't — the iEVORA signal that mean-based DMC misses.
- **Clinical / cohort QC.** Opt-in `qc.sex_check` (chrX mean β), `qc.contamination_estimate` (β-distribution bimodality), `qc.sample_correlation` (sample-swap detection), and `qc.power` (sample-size calculator). Bisulfite conversion rate is reported (CHH context, dashboard + MultiQC) but **not applied** to per-CpG counts — matching `bsseq` / `methylKit` defaults. A poorly converted library should be re-prepped, not papered over with a multiplicative count adjustment.
- **Replicate-aware throughout.** Per-site `min_samples_treatment` / `min_samples_control` guards, per-site or chromosome-level McCullagh-Nelder dispersion, optional covariate design matrices via Wilkinson formulas.
- **Annotation.** Gene features (promoter / 5'UTR / exon / intron / 3'UTR) from GENCODE / Ensembl **GTF** or UCSC **`refGene.txt`** (HOMER's default catalog), plus CpG-island context (island / shore / shelf / open-sea). Opt-in `gene_type_filter="protein_coding"` drops lincRNAs / pseudogenes. `multi_annotation=True` (default) adds annotatr-style `nearest_tss_gene` / `nearest_tss_distance` and one-to-many `all_overlapping_genes` / `all_overlapping_features` columns so a site that's intronic for one gene AND in another gene's promoter window is faithfully represented.
- **Visualisation pack.** matplotlib volcano, MA, Manhattan, coverage histogram, methylation heatmap, PCA, UMAP, sample-correlation heatmap, QC dashboard, DMR boxplot, genomic-context bar, CpG-island pie, TSS metaplot — plus Plotly twins for the HTML report.
- **Result tables as TSV/CSV.** The main analyses write a human-readable TSV by default: `ep.tl.dmc(md)`, `ep.tl.dmr(md)`, and `ep.tl.annotate(md)` drop `dmc.significant.tsv` / `dmr.tsv` / `dmc_annotated.tsv` into `<methylstore>/results/` (override with `tsv="path"`, disable with `tsv=False` or `EPYKIT_NO_AUTO_TSV=1`). Or dump **all** result tables (DMC / DMR / DVC / QC) at once with `md.export_tables("results/")`. Tab-delimited by default (the genomics convention); a `.csv` suffix switches to commas. The individual writers (`ep.export.dmc_to_tsv`, `dmr_to_tsv`, `dvc_to_tsv`, `qc_to_tsv`) are there for fine control, and the CLI auto-emits a sibling `.tsv` next to every parquet output.
- **Interop.** Self-contained HTML report (`md.report("out.html")`), AnnData (`md.to_anndata()`), MuData (`md.to_mudata()`), methylKit-compatible tabix tables (`md.to_methylkit_tabix(dir)`), MultiQC custom-content JSON (`ep.report_multiqc(md, dir)`), nf-core/methylseq QC ingestion (`ep.read_nfcore_methylseq_qc(...)`).
- **CLI.** `epykit convert | filter | dmc | dmr | annotate | qc-report | smooth | report | aggregate-regions | export` — every stage scriptable from the shell.

---

## Installation

Requires Python ≥ 3.10.

```bash
# from the repo checkout
pip install -e .

# or with uv
uv pip install -e .

# dev install (test + lint tools live in the `dev` dependency group)
uv sync --group dev --extra all
pip install -e . --group dev          # pip >= 25.1 equivalent

# full feature install (report + export + anndata + viz)
pip install -e ".[all]"
```

Core dependencies: `polars`, `pyarrow`, `numpy`, `scipy`, `numba`, `bioframe`, `pyfaidx`, `statsmodels`, `patsy`, `psutil`, `scikit-learn`, `matplotlib`, `seaborn`. Optional extras: `report` (Jinja2 + Plotly), `export` (pyBigWig), `anndata` (anndata + mudata), `viz` (umap-learn), `methylkit` (pysam, for tabix indexing on Linux/macOS). The CLI is installed as the `epykit` console script.

> **Platform note.** The pure-Python core and the `report` / `anndata` / `viz` extras run on Linux, macOS, and Windows (CI covers `{ubuntu, windows} × {py3.10, py3.12, py3.13}`). A few extras are Linux/macOS only because they have no Windows wheel: `export` (pyBigWig) and the `pysam`-backed `methylkit` / `bam` extras.

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
ep.pp.set_unite_type(md, type="intersect")  # or "union" + min_samples_* guards

# Tools (tl.*) — populate md.obs / md.varm / md.uns.
ep.tl.qc(md)                                # populates md.obs with QC metrics
ep.tl.dmc(md, test="auto")                  # md.varm["dmc_lr"] (n≥2) or "dmc_fisher" (n=1, allow_n1=True)
ep.tl.dmr(md)                               # chain_merge default (DSS callDMR semantics)

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
# tl.dmc / tl.dmr / tl.annotate already auto-wrote human-readable TSVs to
# methyl_store/results/ by default (dmc.significant.tsv, dmr.tsv, dmc_annotated.tsv;
# pass tsv=False or set EPYKIT_NO_AUTO_TSV=1 to opt out, or tsv="path" to redirect).
md.export_tables("results/tables")        # or grab everything (incl. DVC/QC) in one call

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

### 5. Clinical / cohort QC

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
| `dmc`               | Per-CpG differential methylation. `--test {auto,lr,glm,welch_t,fisher}`, plus `--formula` / `--contrast` / `--covariates` for covariate-adjusted and multi-group designs. `--fdr-method` and `--dispersion` are exposed. The `lr+` power stack (`power_stack`, `neighbour_combine`, `sep_fallback`) is Python API only; no CLI flags are planned. |
| `dmr`               | DMR calling — `--method chain_merge` (default), `tile`, `sliding_window` or `segment`. `--empirical-fdr --n-perm N` is supported with `--method tile`. |
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

DMC frames carry: `chrom`, `pos`, `strand`, `n_case`, `n_control`, `mean_beta_case`, `mean_beta_control`, `meth_diff`, `meth_diff_ci_lo`, `meth_diff_ci_hi`, `pvalue`, `qvalue`, `log2_odds_ratio_pooled` (pooled-count tests: `lr`, `fisher`) or `coef_treatment_log2` (GLM backend; logit coefficient in log2 units, not log2 of an odds ratio), plus per-test extras (`coef_treatment` / `coef_se` for GLM; `f_stat` / `df1` / `df2` / per-level `mean_beta_<level>` / `meth_diff_max` for multi-group contrasts) and, after `tl.annotate`, `feature_type` / `gene_id` / `cpg_context`. Tile-DMR frames add `start`, `end`, `n_cpgs`, `dmr_type ∈ {hyper, hypo, mixed}`; permutation FDR adds `empirical_pvalue` / `empirical_qvalue`.

---

## Module map

The canonical architecture reference is [`docs/advanced/architecture.md`](docs/advanced/architecture.md); this table is a summary.

| Module             | Role |
|--------------------|------|
| `methyldata.py`    | `MethylData` dataclass — `obs`, `store`, `varm`, `uns`; `.dmc` / `.treatment_ids` / `.control_ids` properties; `save()` / `load()`; `region_beta()` per-region query |
| `io.py`            | `read_bismark`, `read_nfcore_methylseq`, `load` |
| `convert.py`       | `.cov` → partitioned Parquet |
| `filter.py`        | Coverage filter, coverage normalisation, blacklist intersect |
| `pp.py`            | Preprocessing wrappers (`filter_coverage`, `normalize_coverage`, `unite`, `smooth`, `aggregate_regions`) |
| `dmc.py`           | Streaming per-CpG accumulators + statistical engines (`lr`, `glm`, `welch_t`, `fisher`), BH correction |
| `_dmc_store.py`    | `DMCStore` handle — persistent per-chromosome DMC parquet directory + manifest; lets BH and sliding-window DMR stream from disk so peak memory is O(largest chrom), not O(genome) |
| `dmr.py`           | `call_dmr_chain_merge` (default, DSS-style), `call_dmr_tile_based`, `call_dmr_sliding_window`, `empirical_fdr_for_dmr`, `smooth_methylation_gaussian`, `smooth_methylation_bsmooth` |
| `dvc.py`           | Differentially Variable CpG calling (iEVORA-style) |
| `annotate.py`      | `annotate_features` (GTF), `annotate_cpg_islands` (island / shore / shelf / open-sea) |
| `qc.py`            | `bisulfite_conversion_rate`, `global_methylation_report`, `coverage_uniformity`, `sex_check`, `contamination_estimate`, `sample_correlation`, `power` |
| `tl.py`            | High-level orchestrators: `tl.qc`, `tl.dmc`, `tl.dmr`, `tl.dvc`, `tl.annotate` |
| `pl/`              | Plotting — `qc`, `differential`, `genomic`, `clustering`, `metaplot`, `embedding`, `correlation`, `dashboard`, `dmr_boxplot`, plus Plotly twins |
| `report.py`        | Self-contained interactive HTML report (Jinja2 + Plotly) |
| `export.py`        | TSV/CSV result tables (`dmc_to_tsv`, `dmr_to_tsv`, `dvc_to_tsv`, `qc_to_tsv`, `export_tables`) + BedGraph / BigWig / DMC-BED / DMR-BED |
| `anndata_io.py`    | AnnData export |
| `mudata_io.py`     | MuData export (multi-omics bundling) |
| `methylkit_io.py`  | methylKit-compatible tabix tables |
| `multiqc_export.py`| MultiQC custom-content JSON emitter |
| `nfcore_qc.py`     | nf-core/methylseq run-dir QC ingestion |
| `cli.py`           | `epykit` CLI entry point |
| `_glm.py`          | Wilkinson formula → design matrix, batched IRLS binomial GLM, Wald test on contrasts |
| `_style.py`        | Shared matplotlib palette / theme |

---

## Benchmark

`benchmark/` reproduces a head-to-head against eight published DMC/DMR tools on the Piao et al. 2021 simulated dataset, plus a real-data cohort (GSE263850). The canonical TPR / FPR / F1 record is [`benchmark/paper/report/REPORT.md`](benchmark/paper/report/REPORT.md); the manuscript lives in [`benchmark/paper/paper.md`](benchmark/paper/paper.md). Raw simulated data and run caches are not bundled — see [`benchmark/README.md`](benchmark/README.md) for the bootstrap.

The headline claims are made around the bare `lr` engine (the 1.0 default). The `lr+` power stack is positioned as an exploratory research knob — see [`docs/analysis/lr-plus.md`](docs/analysis/lr-plus.md) and the [engine architecture map](docs/advanced/architecture.md).

---

## Tests

```bash
uv sync --group dev --extra all   # or: pip install -e . --group dev  (pip >= 25.1)
pytest -m "not slow"        # fast tier (the CI invocation)
pytest -m slow              # opt-in slow tier (>~5s tests)
```

---

## Citing epykit

If epykit contributes to your work, please cite it. A benchmark manuscript —
*"epykit: a Python-native pipeline for differential methylation analysis of
bisulfite sequencing data, benchmarked on simulated and real WGBS datasets"* —
is in preparation ([`benchmark/paper/paper.md`](benchmark/paper/paper.md)). In
the meantime you can cite the software:

```bibtex
@software{epykit,
  author  = {Ertuğrul, Deniz},
  title   = {epykit: a Python-native WGBS methylation analysis pipeline},
  year    = {2026},
  version = {1.1.0},
  url     = {https://github.com/d-ert/epykit3}
}
```

---

## License

[MIT](LICENSE).
