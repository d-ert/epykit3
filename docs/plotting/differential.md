# Differential Plots

Three plotting functions for visualizing differentially methylated cytosine (DMC)
results. All three require that `tl.dmc()` has been run on the `MethylData` object
beforehand, since they read from the DMC results table stored in `md`.

```python
import epykit as ep

md = ep.read("methylation_data/")
ep.tl.dmc(md, group_column="condition", comparison=("tumor", "normal"))
```

## Volcano Plot

`ep.pl.volcano()` plots methylation difference on the x-axis against
statistical significance (-log10 p-value) on the y-axis. Points are colored
by direction and significance: hypermethylated sites in red, hypomethylated in
blue, and non-significant sites in grey.

```python
ep.pl.volcano(md, save="volcano.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with DMC results |
| `alpha` | `float` | `0.05` | Significance threshold for adjusted p-value |
| `meth_diff_threshold` | `float` | `0.1` | Minimum absolute methylation difference to color a site |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(8, 6)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Customizing thresholds

Adjust `alpha` and `meth_diff_threshold` to change which sites are highlighted:

```python
# Stricter thresholds
ep.pl.volcano(
    md,
    alpha=0.01,
    meth_diff_threshold=0.2,
    save="volcano_strict.png",
)
```

### Color scheme

- **Red**: hypermethylated sites (meth_diff > threshold and padj < alpha)
- **Blue**: hypomethylated sites (meth_diff < -threshold and padj < alpha)
- **Grey**: non-significant sites

---

## MA Plot

`ep.pl.ma_plot()` shows the relationship between average methylation level
(mean beta, x-axis) and the methylation difference between groups (y-axis).
This is useful for identifying bias related to overall methylation level.

```python
ep.pl.ma_plot(md, save="ma_plot.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with DMC results |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(8, 6)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example with custom figure size

```python
ep.pl.ma_plot(md, figsize=(10, 5), save="ma_wide.png")
```

---

## Manhattan Plot

`ep.pl.manhattan()` arranges sites along the genome (x-axis, grouped by
chromosome) with -log10(p-value) on the y-axis. A horizontal line indicates
genome-wide significance. Chromosomes are displayed in alternating colors for
readability.

```python
ep.pl.manhattan(md, save="manhattan.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with DMC results |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(14, 5)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example: composing differential plots

Use `figure_grid` or manual subplot layout to combine all three differential
plots into a single figure:

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(22, 6))

ep.pl.volcano(md, ax=axes[0])
ep.pl.ma_plot(md, ax=axes[1])
ep.pl.manhattan(md, ax=axes[2])

fig.tight_layout()
fig.savefig("differential_panel.png", dpi=300)
```
