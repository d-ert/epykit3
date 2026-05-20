# Normalize Coverage

`ep.pp.normalize_coverage()` adjusts read counts across samples so that differences in
sequencing depth do not bias downstream statistical tests. Without normalization,
deeper-sequenced samples contribute more evidence to pooled-count tests (such as the
logistic regression DMC engine), which can inflate false positives.

## Function Signature

```python
ep.pp.normalize_coverage(md, method="median")
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | required | The MethylData object to normalize |
| `method` | `str` | `"median"` | Normalization method: `"median"` or `"mean"` |

## Methods

### Median normalization (default)

```python
ep.pp.normalize_coverage(md, method="median")
```

Scales each sample's counts so that all samples share the same median total coverage. This
is the recommended method because the median is robust to outlier sites with extreme
coverage.

### Mean normalization

```python
ep.pp.normalize_coverage(md, method="mean")
```

Scales each sample's counts so that all samples share the same mean total coverage. This
method is more sensitive to outliers than median normalization, but may be preferred when
the coverage distribution is relatively uniform.

## How It Works

1. For each sample, compute the summary statistic (median or mean) of per-site total
   coverage.
2. Determine a target value (the median of all per-sample summaries).
3. Scale each sample's `count_methylated` and `count_unmethylated` by the ratio of the
   target to that sample's summary statistic.
4. Round scaled counts to the nearest integer (counts must remain whole numbers).
5. Write the normalized data to a new cached partition and repoint `md.store`.

## Usage

```python
import epykit as ep

md = ep.read_bismark("samplesheet.csv", treatment_group="tumor", control_group="normal", assembly="hg38")

# Step 1: Filter first
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)

# Step 2: Normalize
ep.pp.normalize_coverage(md, method="median")
```

## When to Use Each Method

| Scenario | Recommended Method |
|----------|-------------------|
| General WGBS with variable depth | `"median"` |
| Uniform-depth targeted sequencing | `"mean"` or `"median"` |
| Samples with many outlier sites | `"median"` |

In most cases, `"median"` is the safe default. The difference between the two methods is
small when coverage distributions are well-behaved (i.e., after proper filtering).

## Call Order

Normalization should be called **after** `filter_coverage` and **before** `unite`:

```python
ep.pp.filter_coverage(md)       # Remove unreliable sites first
ep.pp.normalize_coverage(md)    # Then normalize remaining sites
ep.pp.unite(md)                 # Then define the site universe
```

Filtering before normalization ensures that extreme-coverage outliers (which are removed by
`filter_coverage`) do not distort the normalization scaling factors.

## Next Steps

After normalization, proceed to [Unite Sites](unite.md):

```python
ep.pp.unite(md, type="intersect")
```
