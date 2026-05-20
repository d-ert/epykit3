# k-NN Imputation

!!! warning "Experimental"
    This feature is implemented but has not been validated on real biological data.
    The API may change in future releases.

epykit provides two k-nearest-neighbor imputation functions for filling
missing methylation values: one operating directly on the methylstore, and
one operating on an AnnData object.

## Functions

### `ep.impute_knn_beta` -- Methylstore Imputation

Per-chromosome inverse-distance-weighted kNN imputation on the methylstore.
For each missing value, the `k` nearest CpG sites (within
`max_distance_bp`) that have observed values in the same sample are
identified, and the missing value is filled using an inverse-distance
weighted average.

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)

# Impute missing beta values in the methylstore
ep.impute_knn_beta(md.store, k=5, max_distance_bp=1000)
```

### `ep.impute_knn_anndata` -- AnnData Imputation

kNN imputation on the sample-by-site AnnData matrix produced by
`ep.export.to_anndata()`. Neighbors are identified in feature space (not
genomic distance), making this suitable for array-like or post-unite data
where each column is a CpG site.

```python
adata = ep.export.to_anndata(md)

# Impute missing values in the AnnData object
ep.impute_knn_anndata(adata, k=5)
```

## Parameters

### `ep.impute_knn_beta`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `store` | MethylStore | required | Methylstore backing the MethylData object (`md.store`) |
| `k` | int | `5` | Number of nearest neighbors |
| `max_distance_bp` | int | `1000` | Maximum genomic distance (bp) to search for neighbors |
| `chromosomes` | list | None | Restrict to specific chromosomes |

### `ep.impute_knn_anndata`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | AnnData object with missing values (NaN) |
| `k` | int | `5` | Number of nearest neighbors |

## When to Use Imputation

Imputation is useful when downstream methods require a complete matrix but
the raw data is sparse:

- **Single-cell bisulfite sequencing** (scBS-seq): coverage per cell is
  typically 1-5%, leaving the vast majority of CpG sites unmeasured.
- **Low-coverage WGBS** (< 5x): many sites have no reads in some samples.
- **Post-unite matrices**: `ep.pp.unite(md)` retains only sites with
  coverage in at least one sample, but individual samples may still have
  gaps.

### Which function to choose

| Scenario | Function | Rationale |
|----------|----------|-----------|
| Pre-unite, filling sparse methylstore | `impute_knn_beta` | Operates on raw per-chromosome data using genomic distance |
| Post-unite AnnData for PCA, clustering | `impute_knn_anndata` | Operates in feature space on the dense matrix |

## When NOT to Use Imputation

Imputation introduces synthetic values and can distort downstream results.
Avoid imputation in these scenarios:

- **DMC/DMR calling**: The statistical tests in `ep.tl.dmc()` handle
  missing data natively by testing only sites with observed counts. Filling
  missing values would inflate sample sizes and produce anti-conservative
  p-values.
- **Well-covered bulk WGBS** (> 10x): Most sites will have coverage in all
  samples after `ep.pp.unite()`. Imputation adds no benefit and only
  introduces bias.
- **Quantitative effect-size estimation**: Imputed values regress toward
  local averages, attenuating true methylation differences.

!!! warning "Imputation and hypothesis testing"
    Never impute before running `ep.tl.dmc()` or `ep.tl.dmr()`. The
    imputed values are not real observations and will inflate statistical
    power, leading to false positives.

## Implementation Details

- `impute_knn_beta` is a pure-NumPy implementation. It processes each
  chromosome independently, loading the beta matrix into memory one
  chromosome at a time. Memory usage scales with the number of CpG sites
  on the largest chromosome.
- `impute_knn_anndata` uses scikit-learn's `KNNImputer` under the hood
  when available, falling back to a pure-NumPy implementation otherwise.
- Both functions modify the data in place. There is no undo. If you need
  to preserve the original data, copy the store or AnnData object before
  imputing.

## Example: Single-cell Workflow

```python
import epykit as ep

md = ep.read_bismark("sc_samples.csv", treatment_group="celltype_A",
                     control_group="celltype_B", assembly="hg38")
ep.pp.filter_coverage(md, min_coverage=1)
ep.pp.unite(md)

# Impute the methylstore for sparse single-cell data
ep.impute_knn_beta(md.store, k=5, max_distance_bp=500)

# Export to AnnData for clustering
adata = ep.export.to_anndata(md)

# Or: impute directly on the AnnData
# ep.impute_knn_anndata(adata, k=5)
```
