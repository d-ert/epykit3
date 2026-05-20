# Aggregate Regions

`ep.pp.aggregate_regions()` collapses per-CpG methylation data into region-level summaries
using a user-supplied BED file of genomic regions. After aggregation, downstream analysis
(including `ep.tl.dmc()`) operates on region-level counts rather than individual CpG
sites.

## Function Signature

```python
ep.pp.aggregate_regions(md, regions_bed="promoters.bed", min_cpgs_per_region=1)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | required | The MethylData object to aggregate |
| `regions_bed` | `str` | required | Path to a BED file defining the regions to aggregate over |
| `min_cpgs_per_region` | `int` | `1` | Minimum number of CpG sites a region must contain to be retained |

## How It Works

1. **Overlap** -- Each CpG site in the methylstore is intersected with the regions defined
   in `regions_bed`.
2. **Sum counts** -- Within each region, `count_methylated` and `count_unmethylated` are
   summed across all overlapping CpG sites for each sample.
3. **Filter** -- Regions with fewer than `min_cpgs_per_region` CpGs are dropped.
4. **Output** -- The aggregated data replaces the per-CpG store. After this step,
   `md.store` points to region-level data, and all downstream functions (DMC, DMR,
   plotting) work at the region level.

## Usage

### Promoter-level analysis

```python
import epykit as ep

md = ep.read_bismark("samplesheet.csv", treatment_group="tumor", control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)

# Aggregate to promoter regions
ep.pp.aggregate_regions(
    md,
    regions_bed="/references/hg38_promoters.bed",
    min_cpgs_per_region=3,
)

# DMC now tests promoter-level methylation differences
ep.tl.dmc(md)
```

### Enhancer-level analysis

```python
ep.pp.aggregate_regions(
    md,
    regions_bed="/references/hg38_enhancers.bed",
    min_cpgs_per_region=5,
)
```

### CpG island analysis

```python
ep.pp.aggregate_regions(
    md,
    regions_bed="/references/hg38_cpg_islands.bed",
    min_cpgs_per_region=10,
)
```

## BED File Format

The `regions_bed` file should be a standard 3-column (or more) BED file with tab-separated
fields:

```
chr1    9873    16357
chr1    28735   29810
chr1    29320   30081
```

At minimum, the first three columns (`chrom`, `start`, `end`) are required. Additional
columns (name, score, strand, etc.) are ignored during aggregation but a `name` column
(column 4) will be used as the region identifier if present.

## Choosing min_cpgs_per_region

| Value | Effect |
|-------|--------|
| `1` | Keep all regions that contain at least one CpG. Maximizes the number of testable regions. |
| `3` | A moderate filter. Removes regions with very sparse CpG coverage where aggregated counts are unreliable. |
| `5-10` | Stricter. Ensures each region has enough CpG sites for a meaningful methylation estimate. Recommended for CpG islands. |

Higher values improve the reliability of per-region methylation estimates but reduce the
number of regions tested.

## Use Cases

- **Promoter methylation** -- Test whether promoter-level methylation differs between
  conditions, which can be more biologically interpretable than individual CpG changes.
- **Enhancer methylation** -- Assess regulatory region methylation changes.
- **CpG island / CpG shore analysis** -- Aggregate across defined CpG island boundaries.
- **Custom regions** -- Any BED file of genomic intervals works, including gene bodies,
  exons, or regions from ChIP-seq peaks.

## Call Order

Aggregation is optional and should follow the core preprocessing steps:

```python
ep.pp.filter_coverage(md)       # Required
ep.pp.normalize_coverage(md)    # Required
ep.pp.unite(md)                 # Required
ep.pp.aggregate_regions(md)     # Optional -- replaces per-CpG data with region-level data
```

After aggregation, `ep.tl.dmc()` will test regions instead of individual CpG sites. If
you want both per-CpG and region-level results, run the per-CpG analysis first, save the
results, then aggregate and run a second analysis:

```python
# Per-CpG analysis
ep.pp.filter_coverage(md)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)
ep.tl.dmc(md)
md.save("results/per_cpg_analysis")

# Region-level analysis (re-read or reload to start fresh)
md = ep.load("results/per_cpg_analysis")
ep.pp.aggregate_regions(md, regions_bed="promoters.bed", min_cpgs_per_region=3)
ep.tl.dmc(md)
md.save("results/promoter_analysis")
```
