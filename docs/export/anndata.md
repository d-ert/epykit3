# AnnData / MuData

Export methylation data to scverse-compatible formats for integration with single-cell
and multi-omics analysis tools such as scanpy and scvi-tools.

## Installation

```bash
pip install 'epykit[anndata]'
```

This installs `anndata` and `mudata` as optional dependencies.

---

## ep.to_anndata

Convert a MethylData object into an AnnData object with a samples-by-sites matrix.

```python
ep.to_anndata(md, layer="beta")
```

**Parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `md` | `MethylData` | MethylData object with a united site set |
| `layer` | `str` | Which value to fill into `X`. `"beta"` (default) for methylation fraction, `"coverage"` for read depth, or `"M"` for M-values |

**Returns**

An `anndata.AnnData` object where:

- `.obs` contains sample metadata from the samplesheet
- `.var` contains site coordinates (`chrom`, `pos`, `strand`)
- `.X` contains the requested value matrix (samples x sites)

**Prerequisites**

You must run `pp.unite()` before calling `to_anndata` so that all samples share
a common set of CpG sites:

```python
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")
ep.pp.filter_coverage(md, min_cov=10)
ep.pp.unite(md)

adata = ep.to_anndata(md, layer="beta")
print(adata)
# AnnData object with n_obs x n_vars = 12 x 4521307
#     obs: 'sample_id', 'treatment', 'age', 'sex'
#     var: 'chrom', 'pos', 'strand'
```

**Memory**

The matrix is filled using a streamed per-sample, per-chromosome strategy. Each
sample-chromosome partition is read from the Parquet store independently, which keeps
peak memory usage well below what a full in-memory join would require. For very large
datasets, the total AnnData object may still be large -- consider subsetting sites
(e.g., variable CpGs only) before conversion.

---

## ep.to_mudata

Create a MuData multi-omics container with methylation as the primary modality.

```python
ep.to_mudata(md, other_modalities=None)
```

**Parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `md` | `MethylData` | MethylData object with a united site set |
| `other_modalities` | `dict[str, AnnData]` or `None` | Optional dictionary mapping modality names to AnnData objects for multi-omics integration |

**Returns**

A `mudata.MuData` object with:

- `"meth"` modality containing the methylation AnnData (equivalent to `ep.to_anndata(md)`)
- Any additional modalities passed via `other_modalities`

**Example**

```python
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")
ep.pp.unite(md)

# Methylation only
mdata = ep.to_mudata(md)
print(mdata)
# MuData object with 1 modality
#   meth: 12 x 4521307

# Multi-omics: combine methylation with RNA-seq
import scanpy as sc
rna = sc.read_h5ad("rna_counts.h5ad")

mdata = ep.to_mudata(md, other_modalities={"rna": rna})
print(mdata)
# MuData object with 2 modalities
#   meth: 12 x 4521307
#   rna:  12 x 22000
```

---

## Use Cases

### Integration with scanpy

```python
import scanpy as sc

adata = ep.to_anndata(md, layer="beta")

# PCA on most variable CpGs
sc.pp.highly_variable_genes(adata, n_top_genes=5000)
adata = adata[:, adata.var.highly_variable]
sc.tl.pca(adata)
sc.pl.pca(adata, color="treatment")
```

### Integration with scvi-tools

```python
import scvi

adata = ep.to_anndata(md, layer="beta")
scvi.model.SCVI.setup_anndata(adata)
model = scvi.model.SCVI(adata)
model.train()
```

### Multi-omics with MuData

```python
import muon as mu

mdata = ep.to_mudata(md, other_modalities={"rna": rna})
mu.pp.intersect_obs(mdata)
mu.tl.mofa(mdata)
mu.pl.mofa(mdata)
```
