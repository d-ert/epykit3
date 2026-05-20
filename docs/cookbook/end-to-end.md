# End-to-end Pipeline Walkthrough

This tutorial walks through a complete WGBS analysis from samplesheet
preparation to exported results. Each step builds on the previous one. By the
end you will have annotated DMRs, publication-ready plots, and an interactive
HTML report.

## 1. Prepare the Samplesheet

Create a CSV with at least `sample_id`, `group`, and `path`. Add any covariates
(age, sex, batch) as extra columns -- they are carried through the entire
pipeline and available for covariate-adjusted analysis.

```csv
sample_id,group,path,age,sex,batch
tumor_1,tumor,/data/bismark/tumor_1.cov.gz,55,M,batch1
tumor_2,tumor,/data/bismark/tumor_2.cov.gz,62,F,batch1
tumor_3,tumor,/data/bismark/tumor_3.cov.gz,48,M,batch2
normal_1,normal,/data/bismark/normal_1.cov.gz,57,F,batch1
normal_2,normal,/data/bismark/normal_2.cov.gz,71,M,batch2
normal_3,normal,/data/bismark/normal_3.cov.gz,44,F,batch2
```

See [Samplesheet Format](../getting-started/samplesheet.md) for the full
specification.

## 2. Ingest

```python
import epykit as ep

md = ep.read_bismark(
    "samples.csv",
    treatment_group="tumor",
    control_group="normal",
    assembly="hg38",
    store_dir="methylstore",
)
```

This converts each Bismark coverage file into a partitioned Parquet methylstore
and returns a `MethylData` object. For MethylDackel output, use
`ep.read_methyldackel()` with the same interface.

## 3. Preprocessing

Apply the standard preprocessing pipeline in order: filter, normalize, unite.

```python
# Remove low- and high-coverage sites
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)

# Equalise per-sample coverage depth
ep.pp.normalize_coverage(md)

# Align site sets across samples
ep.pp.unite(md, type="union")
```

Each step updates the methylstore in place. The recommended order ensures that
normalisation operates on filtered data, and that the united site set reflects
normalised coverage.

## 4. Quality Control

```python
ep.tl.qc(md)
```

This populates `md.obs` with per-sample QC metrics: global methylation level,
mean coverage, coverage breadth at 1x/5x/10x, and low-coverage flags.

For clinical-grade checks:

```python
ep.tl.qc(
    md,
    run_sex_check=True,
    run_contamination=True,
    run_sample_correlation=True,
    run_power_analysis=True,
)
```

Review the results:

```python
print(md.obs[["sample_id", "global_meth", "mean_coverage", "sex_inferred"]])
```

## 5. DMC Calling with lr+ Power Stack

Run per-CpG differential methylation calling using the `lr` test with the full
lr+ power stack for maximum sensitivity:

```python
ep.tl.dmc(
    md,
    power_stack=True,
    fdr_method="fdr_tsbh",
)
```

This enables empirical-Bayes dispersion, neighbour-aware p-value combining,
and separation-aware Fisher fallback alongside two-stage BH correction.

```python
dmc = md.dmc  # shorthand for md.varm["dmc_lr"]
sig = dmc.filter(dmc["qvalue"] < 0.05)
print(f"Tested {len(dmc)} CpGs, {len(sig)} significant at q < 0.05")
```

See [lr+ Power Stack](../analysis/lr-plus.md) for details on each enhancement.

## 6. DMR Calling with chain_merge

Aggregate significant CpGs into differentially methylated regions:

```python
ep.tl.dmr(md, method="chain_merge", preset="default")
```

```python
dmrs = md.uns["dmr"]
print(f"Found {len(dmrs)} DMRs")
print(dmrs.head())
```

For validation experiments, use `preset="strict"`. For exploratory analyses,
use `preset="permissive"`. See [DMR Calling](../analysis/dmr.md) for tuning
guidance.

## 7. Annotation

Add gene-feature and CpG island context to the DMC and DMR results:

```python
ep.tl.annotate(
    md,
    gtf="gencode.v44.annotation.gtf.gz",
    cpg_islands="cpg_islands_hg38.bed",
)
```

This adds `gene_feature` (promoter, exon, intron, intergenic) and
`cpg_context` (island, shore, shelf, open_sea) columns to both the DMC
DataFrame and the DMR DataFrame.

```python
# Check annotation distribution
annotated_dmc = md.varm["dmc_lr_annotated"]
print(annotated_dmc["gene_feature"].value_counts())
```

## 8. Visualization

```python
import matplotlib.pyplot as plt

# Volcano plot
ep.pl.volcano(md)

# Manhattan plot
ep.pl.manhattan(md)

# PCA of sample-level methylation
ep.pl.pca(md)

# DMR boxplot for the top DMR
ep.pl.dmr_boxplot(md, dmr_index=0)

# QC dashboard
ep.pl.qc_dashboard(md)

plt.show()
```

Additional plots available: `ma_plot`, `umap`, `karyogram`, `tss_metaplot`,
`gene_body_metaplot`, `dmr_violin`, `dmr_heatmap`, `genomic_context_bar`,
`cpg_island_pie`, `sample_correlation`, `coverage_histogram`, and
`figure_grid` for multi-panel layouts. See [Plotting](../plotting/index.md).

## 9. Export

### HTML Report

```python
md.report("results/report.html")
```

Generates a self-contained interactive HTML report with summary statistics,
plots, and downloadable tables.

### AnnData

```python
adata = ep.to_anndata(md, layer="beta")
adata.write_h5ad("results/methylation.h5ad")
```

### BigWig

```python
ep.export.to_bigwig(md, output_dir="results/bigwig/")
```

Writes one BigWig file per sample for genome browser visualisation.

### MultiQC

```python
ep.report_multiqc(md, output_dir="results/multiqc/")
```

## 10. Save, Load, and Continue

Save the full analysis state:

```python
md.save("results/my_analysis")
```

This writes `obs.parquet`, per-varm parquet files, `uns` data, and a
`methyldata.json` metadata file.

Load it back in a new session:

```python
md = ep.load("results/my_analysis")

# Everything is restored -- continue from where you left off
print(md.state)
print(md.completed_stages)
print(f"DMC results: {len(md.dmc)} sites")
print(f"DMR results: {len(md.uns['dmr'])} regions")
```

You can re-run any downstream step (e.g., with different DMR parameters) without
repeating preprocessing:

```python
# Try stricter DMR calling on the same DMC results
ep.tl.dmr(md, method="chain_merge", preset="strict")
```
