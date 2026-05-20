# Filter Coverage

`ep.pp.filter_coverage()` removes methylation sites that have too few reads (unreliable
estimates) or too many reads (likely mapping artifacts or repeat regions). An optional
blacklist BED file can be supplied to mask known problematic genomic regions.

## Function Signature

```python
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9, blacklist_bed=None)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | required | The MethylData object to filter |
| `lo_count` | `int` | `10` | Minimum total read count. Sites with fewer reads in any sample are removed. |
| `hi_perc` | `float` | `99.9` | Upper percentile cutoff for coverage. Sites above this percentile (computed per sample) are removed. |
| `blacklist_bed` | `str` or `None` | `None` | Path to a BED file of regions to exclude (e.g. ENCODE blacklist regions). |

## How It Works

1. **Low-coverage filter** -- For each sample, any site with total coverage (methylated +
   unmethylated reads) below `lo_count` is marked for removal. A site is removed if it
   fails the threshold in **any** sample.

2. **High-coverage filter** -- For each sample, the coverage distribution is computed and
   sites above the `hi_perc` percentile are marked for removal. This catches PCR
   duplicates and multi-mapping artifacts that inflate coverage at specific loci.

3. **Blacklist filter** -- If `blacklist_bed` is provided, any site overlapping a
   blacklisted region is removed regardless of its coverage.

4. **Output** -- The filtered data is written to a new cached partition and `md.store` is
   repointed to it.

## Usage

### Basic filtering

```python
import epykit as ep

md = ep.read_bismark("samplesheet.csv", treatment_group="tumor", control_group="normal", assembly="hg38")

# Remove sites with < 10 reads or in the top 0.1% of coverage
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)
```

### With a blacklist

```python
# Use the ENCODE blacklist for hg38
ep.pp.filter_coverage(
    md,
    lo_count=10,
    hi_perc=99.9,
    blacklist_bed="/references/hg38-blacklist.v2.bed",
)
```

### Stricter thresholds for low-coverage experiments

```python
# For shallow sequencing, a lower minimum may be appropriate
ep.pp.filter_coverage(md, lo_count=5, hi_perc=99.5)
```

## Choosing Parameters

- **`lo_count=10`** is a widely used default in WGBS literature. It ensures that
  methylation percentages are estimated from at least 10 observations. For targeted
  bisulfite sequencing with high depth, you may increase this.

- **`hi_perc=99.9`** removes the top 0.1% of sites by coverage. If your data has
  substantial PCR duplication, consider lowering this to `99.5` or `99.0`.

- **Blacklists** are especially useful for whole-genome data. The ENCODE blacklist regions
  are available for hg19, hg38, mm9, and mm10.

## Effect on the MethylData Object

After filtering:

- `md.store` points to the filtered data layer.
- The number of sites is reduced. The exact count depends on your thresholds and data.
- No sample-level metadata in `md.obs` is changed.

## Next Steps

After filtering, proceed to [Normalize Coverage](normalize-coverage.md):

```python
ep.pp.normalize_coverage(md, method="median")
```
