# epykit

**Parquet-backed WGBS methylation analysis pipeline.**

epykit is a Python-native whole-genome bisulfite sequencing (WGBS) analysis toolkit
that takes you from Bismark or MethylDackel output to annotated differentially
methylated regions in a single, scriptable pipeline. It is built on partitioned
Parquet storage and never loads the whole genome into RAM.

## Highlights

- **Partitioned Parquet store** -- hive-partitioned by `sample=*/chrom=*`, enabling random-access queries and chromosome-streaming analysis on whole-genome data without loading everything into memory.
- **4 DMC backends** -- `lr` (quasi-binomial likelihood-ratio, default at n >= 2), `glm` (IRLS binomial with covariates), `welch_t` (Welch t on raw beta), and `fisher` (pooled Fisher exact, n = 1 fallback). `auto` resolves to `lr` at n >= 2 and `fisher` at n = 1. Multi-group and continuous-covariate contrasts via patsy formulas.
- **4 DMR methods** -- `chain_merge` (DSS-style, default), `tile`, `sliding_window`, and `segment`. Preset bundles (`strict`, `default`, `permissive`) for chain-merge.
- **`lr+` power stack (opt-in)** -- four enhancements on top of `lr` (empirical-Bayes dispersion, neighbour combining, separation-aware fallback, two-stage BH q-values) engaged via `tl.dmc(..., power_stack="lr+")`. Off by default; bare `lr` is what you get out of the box.
- **25+ plots** -- volcano, MA, Manhattan, PCA, UMAP, QC dashboard, karyogram, TSS/gene-body metaplots, DMR boxplots/violins/heatmaps, annotation stacked bars, and more.
- **Clinical QC** -- bisulfite conversion rate, coverage uniformity, sex check, contamination estimation, sample correlation, and power analysis.
- **Interop** -- export to AnnData, MuData, methylKit tabix, BedGraph, BigWig, BED, and self-contained interactive HTML reports.
- **Scanpy-style API** -- `ep.pp.*` (preprocessing), `ep.tl.*` (analysis tools), `ep.pl.*` (plotting), `ep.query.*` (random-access queries). Functions modify the `MethylData` object in place.

## Documentation

| Section | Description |
|---------|-------------|
| [Installation](getting-started/installation.md) | Install epykit and optional extras |
| [Quickstart](getting-started/quickstart.md) | End-to-end pipeline walkthrough |
| [Core Concepts](getting-started/concepts.md) | MethylData, methylstore, API namespaces, caching |
| [Samplesheet Format](getting-started/samplesheet.md) | How to prepare your input CSV |

## Minimal example

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor", control_group="normal",
                     assembly="hg38", store_dir="my_store")
ep.pp.filter_coverage(md)
ep.pp.normalize_coverage(md)
ep.pp.set_unite_type(md, type="union")
ep.tl.dmc(md)                  # default test="auto" -> "lr" at n >= 2
ep.tl.dmr(md)                  # default method="chain_merge"
ep.tl.annotate(md, gtf="gencode.v44.gtf.gz")
md.save("results/my_analysis")
```
