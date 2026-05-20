# MethylData API

`MethylData` is the central data container in epykit. It is a lightweight
dataclass that holds sample metadata, a reference to the on-disk methylstore,
and all analysis results. Every `ep.pp.*`, `ep.tl.*`, and `ep.pl.*` function
takes a `MethylData` object and modifies it in place.

## Slots

| Slot | Type | Description |
|------|------|-------------|
| `obs` | `pl.DataFrame` | Per-sample metadata. One row per sample. Columns include `sample_id`, `group`, and any samplesheet covariates. QC metrics are added here by `ep.tl.qc()`. |
| `store` | `str` | Filesystem path to the partitioned Parquet methylstore. Updated by preprocessing steps as the active cache stage changes. |
| `varm` | `dict[str, pl.DataFrame]` | Per-variable (per-CpG) result tables, keyed by name. DMC results are stored here (e.g., `"dmc_lr"`, `"dmc_glm"`, `"dmc_lr_annotated"`). |
| `uns` | `dict[str, Any]` | Unstructured storage for analysis results that are not per-sample or per-CpG. DMR results, DVC results, annotation metadata, and other outputs are stored here. |

## Properties

| Property | Type | Description |
|----------|------|-------------|
| `dmc` | `pl.DataFrame` | Shortcut to the most recent DMC result in `varm`. Equivalent to `md.varm["dmc_<test>"]` for the last test that was run. |
| `treatment_ids` | `list[str]` | Sample IDs in the treatment group. |
| `control_ids` | `list[str]` | Sample IDs in the control group. |
| `n_samples` | `int` | Total number of samples (treatment + control). |
| `state` | `str` | Current preprocessing state (e.g., `"united"`, `"filtered"`). |
| `completed_stages` | `list[str]` | Ordered list of completed pipeline stages (e.g., `["ingest", "filter", "normalize", "unite"]`). |

## Methods

### save

```python
md.save(path: str)
```

Save the full analysis state to a directory. Writes `obs.parquet`, per-varm
parquet files, serialized `uns` data, and a `methyldata.json` metadata file.
The methylstore itself is not copied -- the saved metadata references the
original store path.

### get_dmc

```python
md.get_dmc(test: str | None = None) -> pl.DataFrame
```

Retrieve a DMC result by test name. If `test` is None, returns the most recent
DMC result (same as the `dmc` property). If `test` is specified, returns
`md.varm["dmc_<test>"]`.

```python
# Get the most recent DMC result
dmc = md.get_dmc()

# Get a specific test's result
dmc_glm = md.get_dmc(test="glm")
```

### region_beta

```python
md.region_beta(chrom: str, start: int, end: int) -> pl.DataFrame
```

Query per-sample methylation beta values within a genomic region. Returns a
DataFrame with columns `sample_id`, `pos`, and `beta` for all CpGs in the
specified region across all samples.

```python
region = md.region_beta("chr1", 1_000_000, 1_010_000)
print(region.head())
```

### report

```python
md.report(output: str)
```

Generate a self-contained interactive HTML report. The report includes summary
statistics, QC metrics, top DMC/DMR tables, and embedded plots.

```python
md.report("results/report.html")
```

### to_anndata

```python
md.to_anndata(layer: str = "beta") -> anndata.AnnData
```

Convert to an AnnData object with a samples-by-sites matrix. Requires a united
site set (`ep.pp.unite()` must have been called). See
[AnnData / MuData](../export/anndata.md) for details.

```python
adata = md.to_anndata(layer="beta")
```

### to_mudata

```python
md.to_mudata(other_modalities: dict | None = None) -> mudata.MuData
```

Convert to a MuData multi-omics container with methylation as the primary
modality. Optionally include additional modalities (e.g., RNA-seq) for
multi-omics integration.

```python
mdata = md.to_mudata()
```

## Result Storage Conventions

Each analysis step stores its output in a specific slot:

| Analysis | Storage Location | Example Access |
|----------|-----------------|----------------|
| QC metrics | `md.obs` | `md.obs["global_meth"]` |
| DMC results | `md.varm["dmc_<test>"]` | `md.varm["dmc_lr"]` or `md.dmc` |
| Annotated DMC | `md.varm["dmc_<test>_annotated"]` | `md.varm["dmc_lr_annotated"]` |
| DMR results | `md.uns["dmr"]` | `md.uns["dmr"]` |
| DVC results | `md.uns["dvc"]` | `md.uns["dvc"]` |
| Annotation metadata | `md.uns["annotation_*"]` | `md.uns["annotation_gtf"]` |

Annotation (`ep.tl.annotate()`) adds `gene_feature` and `cpg_context` columns
directly to the DMC DataFrame and stores the annotated version as
`md.varm["dmc_<test>_annotated"]`. The original un-annotated DMC result is
preserved.

## Creating a MethylData Object

You do not construct `MethylData` directly. It is returned by the reader
functions:

```python
import epykit as ep

# From Bismark
md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")

# From MethylDackel
md = ep.read_methyldackel("samples.csv", treatment_group="tumor",
                          control_group="normal", assembly="hg38")

# From nf-core/methylseq
md = ep.read_nfcore_methylseq(run_dir="results/", treatment_group="tumor",
                              control_group="normal",
                              treatment_samples=["t1"], control_samples=["c1"],
                              assembly="hg38")

# Load a saved analysis
md = ep.load("results/my_analysis")
```

## Printing

Printing a `MethylData` object displays a summary of its contents:

```python
print(md)
# MethylData object
#   Samples: 6 (3 treatment, 3 control)
#   Store: methylstore/ (united)
#   Stages: ingest -> filter -> normalize -> unite
#   varm keys: dmc_lr, dmc_lr_annotated
#   uns keys: dmr, annotation_gtf
```
