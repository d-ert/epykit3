# Preprocessing Overview

The `ep.pp` module provides functions that prepare raw methylation data for differential
analysis. Each preprocessing step reads from the current methylstore, applies a
transformation, and writes a new cached output that `md.store` is repointed to. This
means preprocessing is non-destructive -- the original data remains on disk, and each step
produces a new layer that downstream functions read from.

## Recommended Call Order

```python
import epykit as ep

md = ep.read_bismark("samplesheet.csv", treatment_group="tumor", control_group="normal", assembly="hg38")

# 1. Remove low-coverage and extreme-coverage sites
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)

# 2. Normalize coverage across samples
ep.pp.normalize_coverage(md, method="median")

# 3. Define the site universe (intersect or union)
ep.pp.unite(md, type="intersect")

# Optional steps (order-independent after unite)
ep.pp.smooth(md, method="gaussian", bandwidth=1000)
ep.pp.aggregate_regions(md, regions_bed="promoters.bed")
```

## Functions at a Glance

| Function | Purpose | When to Use |
|----------|---------|-------------|
| [`ep.pp.filter_coverage()`](filter-coverage.md) | Remove sites with too few or too many reads | Always -- first step after reading |
| [`ep.pp.normalize_coverage()`](normalize-coverage.md) | Equalize coverage depth across samples | Always -- prevents sequencing-depth bias |
| [`ep.pp.unite()`](unite.md) | Define which sites are included in analysis | Always -- required before DMC/DMR |
| [`ep.pp.smooth()`](smooth.md) | Apply spatial smoothing to methylation values | Optional -- for visualization or exploratory analysis |
| [`ep.pp.aggregate_regions()`](aggregate-regions.md) | Collapse per-CpG data into region-level summaries | Optional -- for promoter/enhancer-level analysis |

## How Caching Works

Each preprocessing function writes its output to a new subdirectory within the store and
updates `md.store` to point at the latest result. For example, after running all three
core steps:

```
results/methylstore/
├── raw/                    # Original data from reader
├── filtered/               # After filter_coverage
├── normalized/             # After normalize_coverage
└── united/                 # After unite
```

This layered approach means you can:

- Re-run a step with different parameters without re-reading the original files.
- Inspect intermediate results by loading a specific layer.
- Resume preprocessing from any checkpoint.

## Minimal Working Example

```python
import epykit as ep

# Read data
md = ep.read_bismark(
    samplesheet="samplesheet.csv",
    treatment_group="tumor",
    control_group="normal",
    assembly="hg38",
)

# Preprocess
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)

# Ready for differential analysis
ep.tl.dmc(md)
```

## Next Steps

See the individual preprocessing pages for detailed parameter documentation and usage
examples:

- [Filter Coverage](filter-coverage.md)
- [Normalize Coverage](normalize-coverage.md)
- [Unite Sites](unite.md)
- [Smoothing](smooth.md)
- [Aggregate Regions](aggregate-regions.md)
