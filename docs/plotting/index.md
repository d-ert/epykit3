# Plotting Overview

epykit provides 25 plotting functions through the `ep.pl` namespace for visualizing
DNA methylation data. All plots accept a `MethylData` object and produce
publication-ready matplotlib figures.

## Common Parameters

Every plotting function shares these parameters:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | The methylation data object |
| `ax` | `matplotlib.axes.Axes` | `None` | Axes to draw on. If `None`, a new figure is created |
| `figsize` | `tuple[float, float]` | varies | Figure dimensions in inches `(width, height)` |
| `save` | `str` | `None` | Filename to save the figure. Format inferred from extension |

When `ax` is provided, the plot is drawn onto that axes, which is useful for
composing multi-panel figures. When `save` is provided, the figure is written to
disk and the axes object is still returned.

```python
import epykit as ep

md = ep.read("methylation_data/")

# Draw onto a new figure
ep.pl.volcano(md, save="volcano.png")

# Draw onto an existing axes
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ep.pl.volcano(md, ax=ax)
plt.show()
```

## Styling

### Themes

`apply_theme()` sets a consistent visual style across all subsequent plots:

```python
ep.pl.apply_theme()          # Apply the default epykit theme
ep.pl.apply_theme("minimal") # Use the minimal theme variant
```

### Color Palettes

`set_palette()` configures the color scheme used for group distinctions:

```python
ep.pl.set_palette("Set2")           # Use a named matplotlib colormap
ep.pl.set_palette(["#e41a1c", "#377eb8", "#4daf4a"])  # Custom colors
```

## Plot Catalog

All 25 functions organized by category:

### QC (4 plots)

| Function | Description | Prerequisite |
|----------|-------------|--------------|
| [`ep.pl.coverage_histogram(md)`](qc-plots.md) | Per-sample coverage distribution | -- |
| [`ep.pl.methylation_heatmap(md)`](qc-plots.md) | Heatmap of methylation levels | -- |
| [`ep.pl.mbias_plot(mbias_data)`](qc-plots.md) | M-bias across read positions | -- |
| [`ep.pl.qc_dashboard(md)`](qc-plots.md) | Composite QC figure | `tl.qc()` |

### Differential (3 plots)

| Function | Description | Prerequisite |
|----------|-------------|--------------|
| [`ep.pl.volcano(md)`](differential.md) | Volcano plot of DMCs | `tl.dmc()` |
| [`ep.pl.ma_plot(md)`](differential.md) | MA plot of DMCs | `tl.dmc()` |
| [`ep.pl.manhattan(md)`](differential.md) | Manhattan plot across chromosomes | `tl.dmc()` |

### Clustering and Embedding (3 plots)

| Function | Description | Prerequisite |
|----------|-------------|--------------|
| [`ep.pl.pca(md)`](clustering.md) | PCA of beta matrix | -- |
| [`ep.pl.umap(md)`](clustering.md) | UMAP embedding | `pip install 'epykit[viz]'` |
| [`ep.pl.sample_correlation(md)`](clustering.md) | Clustered correlation heatmap | `tl.qc(run_sample_correlation=True)` |

### Metaplots (2 plots)

| Function | Description | Prerequisite |
|----------|-------------|--------------|
| [`ep.pl.tss_metaplot(md, gtf=...)`](metaplots.md) | Methylation around TSS | GTF file |
| [`ep.pl.gene_body_metaplot(md, gtf=...)`](metaplots.md) | Methylation across gene bodies | GTF file |

### Genomic Context (3 plots)

| Function | Description | Prerequisite |
|----------|-------------|--------------|
| [`ep.pl.genomic_context_bar(md)`](genomic.md) | Feature category bar chart | `tl.annotate()` |
| [`ep.pl.cpg_island_pie(md)`](genomic.md) | CpG island context pie chart | `tl.annotate(cpg_islands=...)` |
| [`ep.pl.karyogram(md)`](genomic.md) | Chromosomal ideogram overlay | -- |

### DMR (4 plots)

| Function | Description | Prerequisite |
|----------|-------------|--------------|
| [`ep.pl.dmr_boxplot(md)`](dmr-plots.md) | Per-sample DMR beta distributions | `tl.dmr()` |
| [`ep.pl.dmr_violin(md)`](dmr-plots.md) | Violin plots of DMR methylation | `tl.dmr()` |
| [`ep.pl.dmr_heatmap(md)`](dmr-plots.md) | DMR methylation heatmap | `tl.dmr()` |
| [`ep.pl.dmr_overlap(md)`](dmr-plots.md) | Overlap of DMR sets | `tl.dmr()` |

### Annotation (4 plots)

| Function | Description | Prerequisite |
|----------|-------------|--------------|
| [`ep.pl.plot_annotation_counts(md)`](annotation-plots.md) | Annotation category bar chart | `tl.annotate()` |
| [`ep.pl.plot_numerical_by_annotation(md)`](annotation-plots.md) | Numeric metric by feature type | `tl.annotate()` |
| [`ep.pl.plot_coannotations(md)`](annotation-plots.md) | Annotation co-occurrence heatmap | `tl.annotate()` |
| [`ep.pl.plot_categorical(md)`](annotation-plots.md) | Categorical annotation distribution | `tl.annotate()` |

### Utility (1 plot)

| Function | Description | Prerequisite |
|----------|-------------|--------------|
| [`ep.pl.figure_grid(plots)`](figure-grid.md) | Compose plots into a multi-panel figure | -- |

## Interactive HTML Reports

The `md.report()` method generates a self-contained HTML report with interactive
Plotly versions of many of the plots listed above. This is useful for exploratory
analysis and sharing results with collaborators.

```python
md.report(output="methylation_report.html")
```
