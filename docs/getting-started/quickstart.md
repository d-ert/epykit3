# Quickstart

This page walks through a complete WGBS analysis, from raw Bismark coverage
files to annotated DMRs and an HTML report.

## 1. Prepare a samplesheet

Create a CSV file (`samples.csv`) with at least three columns: `sample_id`,
`group`, and `path`. Each row points to one Bismark `.cov` or `.cov.gz` file.

```csv
sample_id,group,path
tumor_1,tumor,/data/bismark/tumor_1.cov.gz
tumor_2,tumor,/data/bismark/tumor_2.cov.gz
tumor_3,tumor,/data/bismark/tumor_3.cov.gz
normal_1,normal,/data/bismark/normal_1.cov.gz
normal_2,normal,/data/bismark/normal_2.cov.gz
normal_3,normal,/data/bismark/normal_3.cov.gz
```

See [Samplesheet Format](samplesheet.md) for details on optional columns.

## 2. Read input

```python
import epykit as ep

md = ep.read_bismark(
    "samples.csv",
    treatment_group="tumor",
    control_group="normal",
    assembly="hg38",
    store_dir="methyl_store",
)
```

This converts each Bismark coverage file into a partitioned Parquet
methylstore under `methyl_store/` and returns a `MethylData` object.

For MethylDackel bedGraph input, use `ep.read_methyldackel()` with the same
signature.

## 3. Preprocess

Apply coverage filtering, normalisation, and site alignment:

```python
# Remove sites with coverage < 10 or above the 99.9th percentile
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)

# Equalise per-sample coverage depth (median normalisation)
ep.pp.normalize_coverage(md)

# Record site-alignment strategy (union keeps all sites; intersect
# restricts to sites covered in every sample)
ep.pp.unite(md, type="union")
```

The recommended order is **filter -> normalize -> unite**. Each step updates
`md.store` to point at the processed Parquet partition.

## 4. Quality control

```python
ep.tl.qc(md)
```

This populates `md.obs` with per-sample metrics (global methylation, mean
coverage, coverage breadth at 1x/5x/10x, low-coverage flags) and stores
detailed QC tables in `md.uns`.

Optional clinical QC checks:

```python
ep.tl.qc(
    md,
    run_sex_check=True,
    run_contamination=True,
    run_sample_correlation=True,
)
```

## 5. Differential methylation -- CpG level (DMC)

```python
ep.tl.dmc(md, test="auto")
```

`test="auto"` selects the quasi-binomial likelihood-ratio test (`lr`) when
you have 2+ replicates per group, or Fisher exact when n=1 (with a warning).
The result is stored in `md.varm["dmc_lr"]` and accessible via `md.dmc`.

Available tests: `lr`, `score`, `glm`, `logit_t`, `welch_t`, `bb_lr`, `cmh`,
`fisher`.

For covariate-adjusted analysis:

```python
ep.tl.dmc(md, formula="~ group + age + sex", contrast="group")
```

## 6. Differential methylation -- region level (DMR)

```python
ep.tl.dmr(md)
```

The default method is `chain_merge` (DSS-style), which merges significant
neighbouring CpGs into contiguous regions. The result is stored in
`md.uns["dmr"]` as a Polars DataFrame.

Other methods:

```python
# Tile-based DMR (aggregates counts within fixed-width tiles)
ep.tl.dmr(md, method="tile", tile_size_bp=1000)

# Sliding-window DMR (combines per-CpG p-values with signed Stouffer's Z)
ep.tl.dmr(md, method="sliding_window", window_bp=500, step_bp=250)
```

Use the `preset` parameter with chain-merge for curated parameter bundles:

```python
ep.tl.dmr(md, preset="strict")    # validation-ready
ep.tl.dmr(md, preset="default")   # balanced (recommended)
ep.tl.dmr(md, preset="permissive") # recall-oriented
```

## 7. Annotate

```python
ep.tl.annotate(
    md,
    gtf="gencode.v44.annotation.gtf.gz",
    cpg_islands="cpg_islands_hg38.bed",
)
```

This adds gene-feature context (`promoter`, `exon`, `intron`, `intergenic`)
and CpG island context (`island`, `shore`, `shelf`, `open_sea`) to both DMC
and DMR tables. Annotated DMC tables are stored as
`md.varm["dmc_lr_annotated"]`.

## 8. Save and generate report

```python
md.save("results/my_analysis")
```

This writes `obs.parquet`, per-varm parquet files, and `methyldata.json` to
the output directory.

Generate a self-contained interactive HTML report:

```python
md.report("results/report.html")
```

## 9. Plotting

```python
import matplotlib.pyplot as plt

# Volcano plot of DMC results
ep.pl.volcano(md)

# MA plot
ep.pl.ma_plot(md)

# Manhattan plot across chromosomes
ep.pl.manhattan(md)

# PCA of sample-level methylation
ep.pl.pca(md)

# QC dashboard (coverage, global methylation, conversion rate)
ep.pl.qc_dashboard(md)

plt.show()
```

Additional plots include `tss_metaplot`, `gene_body_metaplot`, `umap`,
`karyogram`, `genomic_context_bar`, `cpg_island_pie`, `dmr_boxplot`,
`dmr_violin`, `dmr_heatmap`, `sample_correlation`, `methylation_heatmap`,
`coverage_histogram`, and `figure_grid` (a layout composer).

## 10. Load a saved analysis

```python
md = ep.load("results/my_analysis")
print(md)
```

The loaded `MethylData` object has all `obs`, `varm`, and `uns` restored,
including DMC/DMR results and QC metrics. Preprocessing state (`filtered`,
`united`, etc.) is derived automatically from the loaded metadata.
