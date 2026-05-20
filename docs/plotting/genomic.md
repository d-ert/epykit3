# Genomic Context Plots

Three functions for visualizing the genomic context of CpG sites, DMCs, or DMRs.
These plots show where methylation changes occur relative to gene features and
CpG island annotations.

```python
import epykit as ep

md = ep.read("methylation_data/")
```

## Genomic Context Bar Chart

`ep.pl.genomic_context_bar()` displays a bar chart of CpG sites broken down
by genomic feature category: promoter, exon, intron, intergenic, 5' UTR,
3' UTR, and others. This provides an overview of where the measured or
differentially methylated sites fall within the genome.

**Prerequisite:** Run `tl.annotate()` first to assign feature annotations:

```python
ep.tl.annotate(md, gtf="reference/gencode.gtf")
ep.pl.genomic_context_bar(md, save="genomic_context.png")
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
ep.tl.annotate(md, gtf="reference/gencode.gtf")
ep.pl.genomic_context_bar(md, figsize=(12, 6), save="context_bar.png")
```

---

## CpG Island Pie Chart

`ep.pl.cpg_island_pie()` shows the proportion of CpG sites falling into each
CpG island context category:

- **Island** -- within a CpG island
- **Shore** -- within 2 kb flanking a CpG island
- **Shelf** -- within 2-4 kb flanking a CpG island
- **Open sea** -- more than 4 kb from any CpG island

**Prerequisite:** Run `tl.annotate()` with a CpG island annotation file:

```python
ep.tl.annotate(md, gtf="reference/gencode.gtf", cpg_islands="reference/cpg_islands.bed")
ep.pl.cpg_island_pie(md, save="cpg_island_context.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with CpG island annotations from `tl.annotate()` |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(8, 8)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example

```python
ep.tl.annotate(
    md,
    gtf="reference/gencode.gtf",
    cpg_islands="reference/cpg_islands.bed",
)
ep.pl.cpg_island_pie(md, save="cpg_pie.png")
```

---

## Karyogram

`ep.pl.karyogram()` draws a chromosomal ideogram with an overlay showing the
density or magnitude of DMCs, DMRs, or mean beta values along each chromosome.
This provides a genome-wide spatial view of methylation changes.

Unlike the other two plots in this section, the karyogram does not require
prior annotation. It can work directly with positional data from the
`MethylData` object.

```python
ep.pl.karyogram(md, save="karyogram.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Methylation data object |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(14, 10)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example: genomic context workflow

A typical workflow combining all three genomic context plots:

```python
import epykit as ep

md = ep.read("methylation_data/")
ep.tl.dmc(md, group_column="condition", comparison=("tumor", "normal"))
ep.tl.annotate(
    md,
    gtf="reference/gencode.gtf",
    cpg_islands="reference/cpg_islands.bed",
)

# Individual genomic context plots
ep.pl.genomic_context_bar(md, save="feature_bar.png")
ep.pl.cpg_island_pie(md, save="cpg_pie.png")
ep.pl.karyogram(md, save="karyogram.png")
```
