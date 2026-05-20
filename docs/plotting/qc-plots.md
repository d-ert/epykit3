# QC Plots

Four plotting functions for quality control of methylation data. These help
identify problematic samples, coverage issues, and sequencing biases before
downstream analysis.

```python
import epykit as ep

md = ep.read("methylation_data/")
```

## Coverage Histogram

`ep.pl.coverage_histogram()` displays the distribution of read coverage across
CpG sites for each sample. This helps identify samples with low or uneven
coverage that may need filtering.

```python
ep.pl.coverage_histogram(md, save="coverage.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Methylation data object |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(10, 6)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example

```python
# Compare coverage across a large cohort
ep.pl.coverage_histogram(md, figsize=(14, 6), save="coverage_all_samples.png")
```

---

## Methylation Heatmap

`ep.pl.methylation_heatmap()` produces a heatmap of beta values across samples
and CpG sites. Rows represent sites and columns represent samples. Hierarchical
clustering is applied to both axes by default.

```python
ep.pl.methylation_heatmap(md, save="meth_heatmap.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Methylation data object |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(10, 8)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example

```python
ep.pl.methylation_heatmap(md, figsize=(12, 10), save="heatmap_large.png")
```

---

## M-bias Plot

`ep.pl.mbias_plot()` shows the relationship between methylation percentage and
position along the sequencing read. This reveals position-dependent biases
introduced during library preparation (e.g., end-repair artifacts). Separate
curves are drawn for read 1 (R1) and read 2 (R2) across three cytosine
contexts: CpG, CHG, and CHH.

This function accepts either a pre-parsed DataFrame or a file path to
Bismark-style M-bias output.

```python
# From a file path
ep.pl.mbias_plot("QC/mbias_report.txt", save="mbias.png")

# From a pre-parsed DataFrame
ep.pl.mbias_plot(mbias_df, save="mbias.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mbias_data` | `DataFrame` or `str` | *required* | M-bias data as a DataFrame or path to an M-bias file |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(12, 4)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Cytosine contexts

The plot displays three panels, one for each context:

- **CpG** -- the primary context of interest for most methylation studies
- **CHG** -- where H is A, C, or T; typically unmethylated in mammals
- **CHH** -- typically unmethylated in mammals; useful for detecting incomplete bisulfite conversion

### Example

```python
# Typical usage from Bismark output
ep.pl.mbias_plot(
    "bismark_output/sample1.M-bias.txt",
    figsize=(14, 5),
    save="mbias_sample1.png",
)
```

---

## QC Dashboard

`ep.pl.qc_dashboard()` generates a composite figure that combines multiple QC
metrics into a single overview panel. This requires that `tl.qc()` has been
run first.

```python
ep.tl.qc(md)
ep.pl.qc_dashboard(md, save="qc_dashboard.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with QC results from `tl.qc()` |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(16, 12)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Prerequisite

Run `tl.qc()` before calling this function:

```python
ep.tl.qc(md)
ep.pl.qc_dashboard(md, save="qc_dashboard.png")
```

### Example: full QC workflow

```python
import epykit as ep

md = ep.read("methylation_data/")

# Run QC analysis
ep.tl.qc(md, run_sample_correlation=True)

# Generate individual plots
ep.pl.coverage_histogram(md, save="coverage.png")
ep.pl.methylation_heatmap(md, save="heatmap.png")

# Generate the combined dashboard
ep.pl.qc_dashboard(md, save="qc_dashboard.png")
```
