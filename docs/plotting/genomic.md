# Genomic Context Plots

Three functions for visualizing the genomic context of CpG sites, DMCs, or DMRs.
These plots show where methylation changes occur relative to gene features and
CpG island annotations.

```python
import epykit as ep

md = ep.read("methylation_data/")
```

## Genomic Context Bar Chart

`ep.pl.genomic_context_bar()` displays a bar chart broken down by genomic
feature category: promoter, exon, intron, intergenic, 5' UTR, 3' UTR, and
others. This provides an overview of where the differentially methylated
sites — or regions — fall within the genome.

### DMC vs DMR (`level`)

The `level` argument selects which annotated table the chart counts:

- `level="dmc"` (default) reads the per-CpG annotated table on `md.dmc` —
  "where do differential **cytosines** fall?". This view is weighted by CpG
  density, so CpG-rich features (promoters, islands) dominate.
- `level="dmr"` reads the per-region table on `md.uns["dmr"]` — the
  field-standard "fraction of **DMRs** per feature" (the ChIPseeker / DSS /
  dmrseq-style distribution). One count per region, regardless of how many
  CpGs it spans.

For a level-aware chart with bar/pie choice and co-annotation support, see
[`ep.pl.plot_annotation_counts()`](annotation-plots.md) — `genomic_context_bar`
is the lightweight twin.

**Prerequisite:** Run `tl.annotate()` first to assign feature annotations. For
`level="dmr"` you must also have called `tl.dmr()` so `md.uns["dmr"]` exists
before annotating:

```python
# DMC-level (default)
ep.tl.annotate(md, gtf="reference/gencode.gtf")
ep.pl.genomic_context_bar(md, save="genomic_context.png")

# DMR-level — annotate runs after dmr() so the region table is annotated too
ep.tl.dmr(md)
ep.tl.annotate(md, gtf="reference/gencode.gtf")
ep.pl.genomic_context_bar(md, level="dmr", save="genomic_context_dmr.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with annotations from `tl.annotate()` |
| `level` | `str` | `"dmc"` | `"dmc"` (per-CpG, `md.dmc`) or `"dmr"` (per-region, `md.uns["dmr"]`) |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(7, 4)` | Figure size |
| `save` | `str` | `None` | Filename to save the figure |

### Example

```python
ep.tl.annotate(md, gtf="reference/gencode.gtf")
ep.pl.genomic_context_bar(md, figsize=(12, 6), save="context_bar.png")
```

---

## CpG Island Pie Chart

`ep.pl.cpg_island_pie()` shows the proportion of sites (or regions) falling
into each CpG island context category:

- **Island** -- within a CpG island
- **Shore** -- within 2 kb flanking a CpG island
- **Shelf** -- within 2-4 kb flanking a CpG island
- **Open sea** -- more than 4 kb from any CpG island

Like `genomic_context_bar`, it accepts `level="dmc"` (default, per-CpG on
`md.dmc`) or `level="dmr"` (per-region on `md.uns["dmr"]`). As of the DMR
annotation update, `tl.annotate(cpg_islands=...)` assigns CpG-island context
to the DMR region table as well as the per-CpG table, so the `level="dmr"`
pie has data.

**Prerequisite:** Run `tl.annotate()` with a CpG island annotation file:

```python
# DMC-level (default)
ep.tl.annotate(md, gtf="reference/gencode.gtf", cpg_islands="reference/cpg_islands.bed")
ep.pl.cpg_island_pie(md, save="cpg_island_context.png")

# DMR-level — call dmr() first so md.uns["dmr"] is annotated by annotate()
ep.tl.dmr(md)
ep.tl.annotate(md, gtf="reference/gencode.gtf", cpg_islands="reference/cpg_islands.bed")
ep.pl.cpg_island_pie(md, level="dmr", save="cpg_island_context_dmr.png")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | *required* | Data with CpG island annotations from `tl.annotate()` |
| `level` | `str` | `"dmc"` | `"dmc"` (per-CpG, `md.dmc`) or `"dmr"` (per-region, `md.uns["dmr"]`) |
| `ax` | `Axes` | `None` | Axes to draw on |
| `figsize` | `tuple` | `(5, 5)` | Figure size |
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
