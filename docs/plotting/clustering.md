# Clustering and Embedding Plots

Three functions for exploring sample-level relationships using dimensionality
reduction and correlation analysis. These are useful for identifying batch
effects, outliers, and verifying that samples cluster by biological group.

```python
import epykit as ep

md = ep.read("methylation_data/")
```

## PCA

`ep.pl.pca()` performs principal component analysis on the beta-value matrix
and displays the first two components as a scatter plot. Samples are colored
by their group label.

```python
ep.pl.pca(md, save="pca.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Methylation data object |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(8, 6)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example

```python
# PCA colored by the default group column
ep.pl.pca(md, save="pca_groups.png")
```

---

## UMAP

`ep.pl.umap()` computes a UMAP embedding of the beta-value matrix and displays
the result as a 2D scatter plot. This is helpful when PCA does not separate
groups clearly, as UMAP can capture nonlinear structure.

**Requirement:** UMAP support requires the optional visualization dependencies:

```bash
pip install 'epykit[viz]'
```

This installs `umap-learn` and its dependencies.

```python
ep.pl.umap(md, save="umap.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Methylation data object |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(8, 6)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example

```python
ep.pl.umap(md, figsize=(10, 8), save="umap_large.png")
```

---

## Sample Correlation

`ep.pl.sample_correlation()` draws a clustered heatmap of pairwise Pearson
correlation coefficients between samples, computed from the beta-value matrix.
Hierarchical clustering is applied to both rows and columns to group similar
samples together.

**Prerequisite:** The correlation matrix must be pre-computed by running `tl.qc()`
with `run_sample_correlation=True`:

```python
ep.tl.qc(md, run_sample_correlation=True)
ep.pl.sample_correlation(md, save="sample_corr.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with correlation matrix from `tl.qc()` |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(10, 10)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example: complete clustering workflow

```python
import epykit as ep

md = ep.read("methylation_data/")

# Compute QC metrics including sample correlation
ep.tl.qc(md, run_sample_correlation=True)

# Generate all three clustering/embedding plots
ep.pl.pca(md, save="pca.png")
ep.pl.umap(md, save="umap.png")
ep.pl.sample_correlation(md, save="correlation_heatmap.png")
```

### Interpreting the heatmap

- Diagonal values are always 1.0 (self-correlation).
- Off-diagonal blocks of high correlation suggest sample grouping consistent
  with biological or batch variables.
- Isolated samples with uniformly low correlation may indicate outliers or
  technical failures.
