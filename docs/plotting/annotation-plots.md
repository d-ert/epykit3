# Annotation Plots

Four annotatr-style visualization functions for exploring the relationship
between methylation data and genomic annotations. All four require that
`tl.annotate()` has been run first to assign feature annotations to CpG sites.

```python
import epykit as ep

md = ep.read("methylation_data/")
ep.tl.dmc(md, group_column="condition", comparison=("tumor", "normal"))
ep.tl.annotate(md, gtf="reference/gencode.gtf")
```

## Annotation Counts

`ep.pl.plot_annotation_counts()` displays a bar chart of the number of CpG
sites (or DMCs) falling into each annotation category. Categories include
promoter, 5' UTR, exon, intron, 3' UTR, intergenic, and others depending on
the annotation source.

```python
ep.pl.plot_annotation_counts(md, save="annotation_counts.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with annotations from `tl.annotate()` |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(10, 6)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example

```python
ep.pl.plot_annotation_counts(md, figsize=(12, 6), save="annot_counts.png")
```

---

## Numerical by Annotation

`ep.pl.plot_numerical_by_annotation()` shows the distribution of a numeric
metric (e.g., methylation difference, mean beta) grouped by annotation feature
type. This is useful for understanding how effect sizes or methylation levels
vary across genomic contexts.

```python
ep.pl.plot_numerical_by_annotation(md, column="meth_diff", save="meth_diff_by_annot.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with annotations from `tl.annotate()` |
| `column` | `str` | `"meth_diff"` | Name of the numeric column to plot |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(10, 6)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example: different numeric columns

```python
# Methylation difference by feature type
ep.pl.plot_numerical_by_annotation(
    md, column="meth_diff", save="meth_diff_by_feature.png"
)

# Mean beta value by feature type
ep.pl.plot_numerical_by_annotation(
    md, column="mean_beta", save="mean_beta_by_feature.png"
)
```

---

## Co-annotations

`ep.pl.plot_coannotations()` generates a heatmap showing how often pairs of
annotations co-occur on the same CpG sites. For example, a site might be
annotated as both "promoter" and "CpG island." The heatmap values represent
the number (or proportion) of sites carrying both annotations simultaneously.

```python
ep.pl.plot_coannotations(md, save="coannotations.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with annotations from `tl.annotate()` |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(10, 10)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example

```python
ep.pl.plot_coannotations(md, figsize=(12, 12), save="coannotation_heatmap.png")
```

---

## Categorical Annotations

`ep.pl.plot_categorical()` visualizes the distribution of categorical
annotation values. This provides a breakdown of how sites are distributed
across discrete categories such as chromosome, strand, or custom labels.

```python
ep.pl.plot_categorical(md, save="categorical.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with annotations from `tl.annotate()` |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(10, 6)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example: complete annotation workflow

```python
import epykit as ep

md = ep.read("methylation_data/")
ep.tl.dmc(md, group_column="condition", comparison=("tumor", "normal"))
ep.tl.annotate(
    md,
    gtf="reference/gencode.gtf",
    cpg_islands="reference/cpg_islands.bed",
)

# All four annotation plots
ep.pl.plot_annotation_counts(md, save="annot_counts.png")
ep.pl.plot_numerical_by_annotation(md, column="meth_diff", save="meth_diff_annot.png")
ep.pl.plot_coannotations(md, save="coannotations.png")
ep.pl.plot_categorical(md, save="categorical.png")
```

### Example: composing annotation plots

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(20, 14))

ep.pl.plot_annotation_counts(md, ax=axes[0, 0])
axes[0, 0].set_title("Annotation Counts")

ep.pl.plot_numerical_by_annotation(md, column="meth_diff", ax=axes[0, 1])
axes[0, 1].set_title("Methylation Difference by Feature")

ep.pl.plot_coannotations(md, ax=axes[1, 0])
axes[1, 0].set_title("Co-annotations")

ep.pl.plot_categorical(md, ax=axes[1, 1])
axes[1, 1].set_title("Categorical Distribution")

fig.tight_layout()
fig.savefig("annotation_panel.png", dpi=300)
```
