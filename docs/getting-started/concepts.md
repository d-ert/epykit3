# Core Concepts

## MethylData

`MethylData` is the central dataclass that holds everything about an analysis
session. It is the single object you pass to every `ep.pp.*`, `ep.tl.*`, and
`ep.pl.*` function.

```python
@dataclass
class MethylData:
    obs: pl.DataFrame          # sample metadata (one row per sample)
    store: str                 # path to the partitioned Parquet methylstore
    assembly: str              # genome assembly, e.g. "hg38"
    context: str               # methylation context, e.g. "CpG"
    varm: dict[str, pl.DataFrame]  # per-site result tables (DMC, DVC, etc.)
    uns: dict                  # unstructured results (DMR, QC, config, etc.)
```

### Key slots

| Slot | Contents |
|------|----------|
| `obs` | One row per sample. Columns always include `sample_id`, `group`, and `treatment` (0/1). QC metrics, covariates, and inferred annotations are added here by `ep.tl.qc()`. |
| `store` | Filesystem path to the current methylstore. Updated in place by preprocessing steps (`filter_coverage` -> `normalize_coverage`). |
| `varm` | Dictionary of per-site (variable-annotation) DataFrames. DMC results land here as `"dmc_lr"`, `"dmc_glm"`, etc. Annotated variants are stored as `"dmc_lr_annotated"`. |
| `uns` | Dictionary of unstructured results. DMR tables, QC summaries, preprocessing config, and internal bookkeeping live here. |

### Useful properties

```python
md.dmc                  # most-recently-written DMC table (annotated if available)
md.significant_dmcs     # md.dmc filtered to qvalue < 0.05
md.treatment_ids        # list of sample_id values where treatment == 1
md.control_ids          # list of sample_id values where treatment == 0
md.n_samples            # number of samples
md.state                # ordered list of preprocessing steps applied, e.g. ["filtered", "normalized", "united"]
```

### Save and load

```python
md.save("results/my_analysis")   # persist to disk
md = ep.load("results/my_analysis")  # restore from disk
```

Saving writes `obs.parquet`, one parquet per `varm` entry, serialised `uns`
(with DataFrame values saved as separate parquets), and a `methyldata.json`
manifest. Large DMC tables backed by a `DMCStore` are saved via hard-link
from the per-chromosome parquets -- no re-encoding, constant memory.

### Exports

`MethylData` provides convenience methods for common export formats:

```python
md.report("report.html")              # interactive HTML report
md.to_anndata()                       # -> AnnData (requires ep.pp.unite first)
md.to_mudata()                        # -> MuData
md.to_methylkit_tabix("output_dir/")  # per-sample methylKit tabix tables
md.to_bedgraph("sample_1", "out.bedgraph")
md.to_bigwig("sample_1", "out.bw")
md.dmcs_to_bed("dmcs.bed")
md.dmrs_to_bed("dmrs.bed")
```

---

## Methylstore

The methylstore is a Hive-partitioned Parquet directory with the layout:

```
methyl_store/
  sample=tumor_1/
    chrom=chr1/
      part-0.parquet
    chrom=chr2/
      part-0.parquet
    ...
  sample=tumor_2/
    ...
```

Each `part-0.parquet` file contains per-CpG rows with columns: `pos`,
`strand`, `N_meth`, `N_unmeth` (or `coverage`), and `context`.

### Design principles

- **Never loads the whole genome into RAM.** Chromosome-level streaming is
  the default for DMC, DMR, and QC operations.
- **Predicate pushdown.** Polars lazy scans read parquet row-group statistics
  (`min(pos)` / `max(pos)`) and prune row groups that do not overlap a query
  before reading any data. This gives tabix-equivalent random access without
  a separate index file.
- **Partition pruning.** Selecting a sample and chromosome is O(1) via
  direct file-path construction (`sample=X/chrom=Y/part-0.parquet`).

Preprocessing steps create new partitioned stores under a `.cache/` directory
(e.g. `.cache/filtered`, `.cache/normalized`) and update `md.store` to point
at the latest version.

---

## Scanpy-style API

epykit organises its public functions into four namespaces that mirror the
scanpy convention:

### `ep.pp.*` -- Preprocessing

Functions that transform the methylstore and prepare data for analysis.
They modify `md` in place (updating `md.store`, `md.uns`, etc.).

| Function | Description |
|----------|-------------|
| `ep.pp.filter_coverage(md)` | Remove sites below a minimum coverage or above a percentile ceiling. Optionally exclude blacklisted regions. |
| `ep.pp.normalize_coverage(md)` | Per-sample median (or mean) coverage normalisation. |
| `ep.pp.unite(md)` | Record the site-alignment strategy (`"intersect"` or `"union"`) for downstream DMC processing. Does not materialise the full join. |
| `ep.pp.smooth(md)` | Smooth per-sample beta values along the genome (Gaussian kernel or BSmooth). |
| `ep.pp.aggregate_regions(md, bed)` | Aggregate per-CpG counts within user-supplied BED regions, producing a region-level store compatible with `ep.tl.dmc()`. |

### `ep.tl.*` -- Tools / Analysis

Functions that run statistical analyses and deposit results in `md.varm`
or `md.uns`.

| Function | Description |
|----------|-------------|
| `ep.tl.qc(md)` | Per-sample QC metrics (coverage, global methylation, conversion rate). Optional clinical checks (sex, contamination, correlation). |
| `ep.tl.dmc(md)` | Per-CpG differential methylation calling. 8 test backends; supports covariates via patsy formulas. |
| `ep.tl.dmr(md)` | Region-level DMR calling. 4 methods: `chain_merge`, `tile`, `sliding_window`, `hmm`. |
| `ep.tl.annotate(md)` | Gene-feature and CpG-island annotation of DMC and DMR results. |

### `ep.pl.*` -- Plotting

Static matplotlib/seaborn plots. Each function takes `md` (and optionally
extra parameters) and returns a matplotlib `Figure` or `Axes`.

| Function | Description |
|----------|-------------|
| `ep.pl.volcano(md)` | Volcano plot (meth_diff vs. -log10 p-value) |
| `ep.pl.ma_plot(md)` | MA plot (mean beta vs. meth_diff) |
| `ep.pl.manhattan(md)` | Genome-wide Manhattan plot |
| `ep.pl.pca(md)` | PCA of sample-level methylation profiles |
| `ep.pl.umap(md)` | UMAP embedding (requires `[viz]` extra) |
| `ep.pl.qc_dashboard(md)` | Multi-panel QC overview |
| `ep.pl.karyogram(md)` | Chromosome ideogram with DMC density |
| `ep.pl.tss_metaplot(md)` | Methylation profile around transcription start sites |
| `ep.pl.gene_body_metaplot(md)` | Methylation profile across scaled gene bodies |
| `ep.pl.coverage_histogram(md)` | Per-sample coverage distribution |
| `ep.pl.methylation_heatmap(md)` | Heatmap of per-sample methylation |
| `ep.pl.sample_correlation(md)` | Pairwise sample correlation heatmap |
| `ep.pl.dmr_boxplot(md)` | Per-region boxplots of methylation by group |
| `ep.pl.dmr_violin(md)` | Violin plots of DMR-level effect sizes |
| `ep.pl.dmr_heatmap(md)` | Heatmap of DMR methylation across samples |
| `ep.pl.genomic_context_bar(md)` | Stacked bar chart of genomic feature distribution |
| `ep.pl.cpg_island_pie(md)` | Pie chart of CpG island context |
| `ep.pl.figure_grid(...)` | Layout composer for multi-panel figures |

### `ep.query.*` -- Random-access queries

Functions for fetching methylation data at specific genomic loci without
scanning the entire store.

| Function | Description |
|----------|-------------|
| `ep.query.query_region(store, chrom, start, end)` | Fetch all CpG rows within a single region |
| `ep.query.query_regions(store, regions_df)` | Batch query across multiple regions |
| `ep.query.query_sites(store, sites_df)` | Fetch data at exact CpG positions |

```python
# Example: fetch methylation in the BRCA1 promoter region
df = ep.query.query_region(md.store, "chr17", 43_044_295, 43_170_245)
```

---

## DMCStore

`DMCStore` is a streaming, per-chromosome DMC result container. When
`ep.tl.dmc()` runs, it writes one parquet file per chromosome to a temporary
directory and wraps it in a `DMCStore`. This design keeps peak memory at
O(largest chromosome) rather than O(genome).

The `DMCStore` supports:

- **Streaming BH correction** -- multiple-testing correction is applied
  per-chromosome and then globally, without materialising the full table.
- **Lazy materialisation** -- `dmc_store.to_dataframe()` collects all
  chromosomes into a single Polars DataFrame only when needed (e.g. for
  `md.varm` back-compatibility).
- **Zero-copy save** -- `md.save()` hard-links the per-chromosome parquet
  files from the DMCStore directory rather than re-encoding a multi-gigabyte
  DataFrame.

---

## Caching

epykit caches intermediate results under a `.cache/` directory relative to
the analysis root (or adjacent to the methylstore when no analysis root is
set). Cached stages include:

| Cache path | Contents |
|------------|----------|
| `.cache/filtered/` | Coverage-filtered methylstore |
| `.cache/normalized/` | Coverage-normalised methylstore |
| `.cache/dmc/` | Per-chromosome DMC parquets |
| `.cache/regions/` | Region-aggregated methylstore (from `pp.aggregate_regions`) |

The `resumable=True` parameter on `ep.tl.dmc()` enables a checkpoint/resume
mechanism: when a prior run with matching inputs and parameters is found in
the pipeline manifest (`.epykit_manifest.json`), the cached result is loaded
directly, skipping the computation. This is useful for iterative workflows
where you re-run a notebook but only want to recompute stages whose inputs
have changed.
