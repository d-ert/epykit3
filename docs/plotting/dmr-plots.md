# DMR Plots

Four functions for visualizing differentially methylated regions (DMRs). All
require that `tl.dmr()` has been run to identify regions of coordinated
methylation change.

```python
import epykit as ep

md = ep.read("methylation_data/")
ep.tl.dmc(md, group_column="condition", comparison=("tumor", "normal"))
ep.tl.dmr(md)
```

## DMR Boxplot

`ep.pl.dmr_boxplot()` shows per-sample beta-value distributions for the top
DMRs as box plots. Each panel displays one DMR, with samples grouped by
condition. Beta values within each region are computed using `md.region_beta()`.

```python
ep.pl.dmr_boxplot(md, top_n=6, save="dmr_boxplot.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with DMR results from `tl.dmr()` |
| `top_n` | `int` | `6` | Number of top DMRs to display |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(12, 8)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example: showing more regions

```python
# Show the top 12 DMRs in a larger figure
ep.pl.dmr_boxplot(md, top_n=12, figsize=(16, 12), save="dmr_box_top12.png")
```

---

## DMR Violin Plot

`ep.pl.dmr_violin()` displays violin plots of methylation distributions within
DMRs, providing a richer view of the distribution shape compared to box plots.
Each violin shows the density of beta values for a given group within a DMR.

```python
ep.pl.dmr_violin(md, save="dmr_violin.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with DMR results from `tl.dmr()` |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(12, 8)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example

```python
ep.pl.dmr_violin(md, figsize=(14, 8), save="dmr_violin_wide.png")
```

---

## DMR Heatmap

`ep.pl.dmr_heatmap()` produces a heatmap of mean methylation levels across
DMRs (rows) and samples (columns). Hierarchical clustering is applied to
reveal patterns of coordinated methylation change across regions and samples.

```python
ep.pl.dmr_heatmap(md, save="dmr_heatmap.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with DMR results from `tl.dmr()` |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(10, 12)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example

```python
ep.pl.dmr_heatmap(md, figsize=(12, 14), save="dmr_heatmap_large.png")
```

---

## DMR Overlap

`ep.pl.dmr_overlap()` creates Venn diagram or UpSet-style plots comparing
multiple DMR sets. This is useful when you have DMRs from different comparisons,
methods, or parameter settings and want to assess their overlap.

```python
ep.pl.dmr_overlap(md, save="dmr_overlap.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with DMR results |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(10, 8)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example: complete DMR visualization workflow

```python
import epykit as ep

md = ep.read("methylation_data/")
ep.tl.dmc(md, group_column="condition", comparison=("tumor", "normal"))
ep.tl.dmr(md)

# Generate all DMR plots
ep.pl.dmr_boxplot(md, top_n=6, save="dmr_boxplot.png")
ep.pl.dmr_violin(md, save="dmr_violin.png")
ep.pl.dmr_heatmap(md, save="dmr_heatmap.png")
ep.pl.dmr_overlap(md, save="dmr_overlap.png")
```

### Example: composing DMR plots in a grid

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(20, 16))

ep.pl.dmr_boxplot(md, top_n=4, ax=axes[0, 0])
ep.pl.dmr_violin(md, ax=axes[0, 1])
ep.pl.dmr_heatmap(md, ax=axes[1, 0])
ep.pl.dmr_overlap(md, ax=axes[1, 1])

fig.suptitle("DMR Analysis Summary", fontsize=16)
fig.tight_layout()
fig.savefig("dmr_summary_panel.png", dpi=300)
```
