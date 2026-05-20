# Figure Grid

`ep.pl.figure_grid()` composes multiple plots into a single multi-panel figure.
This is the recommended way to build publication-quality composite figures from
individual epykit plots.

## Basic Usage

Pass a list of plotting functions and specify the grid layout:

```python
import epykit as ep

md = ep.read("methylation_data/")
ep.tl.dmc(md, group_column="condition", comparison=("tumor", "normal"))

ep.pl.figure_grid(
    [
        lambda ax: ep.pl.volcano(md, ax=ax),
        lambda ax: ep.pl.ma_plot(md, ax=ax),
        lambda ax: ep.pl.manhattan(md, ax=ax),
    ],
    ncols=3,
    figsize=(18, 5),
    save="differential_grid.png",
)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `plots` | `list` | *required* | List of callables. Each receives an `Axes` argument and draws onto it |
| `ncols` | `int` | `3` | Number of columns in the grid |
| `figsize` | `tuple` | `(15, 10)` | Overall figure size in inches |
| `save` | `str` | `None` | Filename to save the figure |

## How It Works

`figure_grid` creates a `matplotlib.figure.Figure` with a grid of subplots
determined by the number of plots and `ncols`. It then calls each function in
the `plots` list, passing the corresponding `Axes` object. Any unused subplot
cells (when the number of plots is not a multiple of `ncols`) are hidden
automatically.

## Examples

### Two-column QC summary

```python
ep.tl.qc(md, run_sample_correlation=True)

ep.pl.figure_grid(
    [
        lambda ax: ep.pl.coverage_histogram(md, ax=ax),
        lambda ax: ep.pl.methylation_heatmap(md, ax=ax),
        lambda ax: ep.pl.pca(md, ax=ax),
        lambda ax: ep.pl.sample_correlation(md, ax=ax),
    ],
    ncols=2,
    figsize=(16, 12),
    save="qc_summary_grid.png",
)
```

### Mixed analysis panels

Combine plots from different categories into a single overview figure:

```python
ep.tl.dmc(md, group_column="condition", comparison=("tumor", "normal"))
ep.tl.dmr(md)
ep.tl.annotate(md, gtf="reference/gencode.gtf")

ep.pl.figure_grid(
    [
        lambda ax: ep.pl.volcano(md, ax=ax),
        lambda ax: ep.pl.manhattan(md, ax=ax),
        lambda ax: ep.pl.genomic_context_bar(md, ax=ax),
        lambda ax: ep.pl.dmr_heatmap(md, ax=ax),
        lambda ax: ep.pl.pca(md, ax=ax),
        lambda ax: ep.pl.plot_annotation_counts(md, ax=ax),
    ],
    ncols=3,
    figsize=(20, 12),
    save="full_overview.png",
)
```

### Single-column layout

Use `ncols=1` for a vertical stack of plots:

```python
ep.pl.figure_grid(
    [
        lambda ax: ep.pl.tss_metaplot(md, gtf="reference/gencode.gtf", ax=ax),
        lambda ax: ep.pl.gene_body_metaplot(md, gtf="reference/gencode.gtf", ax=ax),
    ],
    ncols=1,
    figsize=(10, 10),
    save="metaplots_stacked.png",
)
```

### Using pre-made axes

If you need finer control over the layout (e.g., different subplot sizes),
use matplotlib directly and pass axes to individual plot functions instead:

```python
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize=(16, 10))
gs = GridSpec(2, 3, figure=fig)

ax1 = fig.add_subplot(gs[0, :2])   # Wide top-left
ax2 = fig.add_subplot(gs[0, 2])    # Narrow top-right
ax3 = fig.add_subplot(gs[1, :])    # Full-width bottom

ep.pl.manhattan(md, ax=ax1)
ep.pl.volcano(md, ax=ax2)
ep.pl.karyogram(md, ax=ax3)

fig.tight_layout()
fig.savefig("custom_layout.png", dpi=300)
```

This manual approach gives full control over subplot proportions and is
complementary to `figure_grid`, which handles the common case of equal-sized
panels.
