# Metaplots

Two functions for visualizing methylation patterns relative to gene structure.
Both require a GTF annotation file (e.g., from GENCODE) to define gene
coordinates.

```python
import epykit as ep

md = ep.read("methylation_data/")
gtf_path = "reference/gencode.v44.annotation.gtf"
```

## TSS Metaplot

`ep.pl.tss_metaplot()` computes and displays the average methylation level in
a window around transcription start sites (TSS). CpG sites are binned by their
distance from the TSS, and the mean beta value is plotted for each bin. Separate
lines are drawn for each group.

This plot typically reveals a characteristic dip in methylation at active
promoters.

```python
ep.pl.tss_metaplot(md, gtf=gtf_path, save="tss_metaplot.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Methylation data object |
| `gtf` | `str` | *required* | Path to a GTF annotation file |
| `window_bp` | `int` | `5000` | Distance (in base pairs) upstream and downstream of the TSS to include |
| `n_bins` | `int` | `100` | Number of bins across the window |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(10, 5)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example: adjusting the window

```python
# Wider window to capture distal regulatory elements
ep.pl.tss_metaplot(
    md,
    gtf=gtf_path,
    window_bp=10000,
    n_bins=200,
    save="tss_wide.png",
)
```

### Example: narrower focus on promoter

```python
# Focus tightly on the proximal promoter
ep.pl.tss_metaplot(
    md,
    gtf=gtf_path,
    window_bp=2000,
    n_bins=50,
    save="tss_proximal.png",
)
```

---

## Gene Body Metaplot

`ep.pl.gene_body_metaplot()` shows the methylation profile across scaled gene
bodies. Each gene is divided into three regions -- upstream flank / TSS, gene
body, and downstream flank / TES -- and the methylation level is averaged
within bins across these regions. Separate lines are drawn for each group.

This plot is useful for observing the well-known pattern where gene body
methylation correlates with transcriptional activity.

```python
ep.pl.gene_body_metaplot(md, gtf=gtf_path, save="gene_body.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Methylation data object |
| `gtf` | `str` | *required* | Path to a GTF annotation file |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(10, 5)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example: combining both metaplots

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(18, 5))

ep.pl.tss_metaplot(md, gtf=gtf_path, ax=axes[0])
axes[0].set_title("TSS Metaplot")

ep.pl.gene_body_metaplot(md, gtf=gtf_path, ax=axes[1])
axes[1].set_title("Gene Body Metaplot")

fig.tight_layout()
fig.savefig("metaplots_combined.png", dpi=300)
```

### Notes on GTF files

- Use a comprehensive GTF such as GENCODE or Ensembl for best coverage.
- The GTF must contain `gene` or `transcript` feature entries with `gene_id`
  attributes.
- Compressed GTF files (`.gtf.gz`) are supported.
