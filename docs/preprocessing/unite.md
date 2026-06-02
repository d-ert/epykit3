# Unite Sites

`ep.pp.set_unite_type()` defines which CpG sites are included in downstream
analysis by aligning the site universes across all samples. Different samples
may cover different sets of sites, and this step records how to handle that
mismatch.

It is a **state-marker**: it writes the chosen strategy into
`md.uns["unite"]` and does not materialise any data. The actual site
filtering is applied on-the-fly when downstream functions (such as
`ep.tl.dmc()`) query the store. This keeps preprocessing fast and avoids
redundant disk writes.

## Function Signature

```python
ep.pp.set_unite_type(md, type="union")
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | required | The MethylData object to mark |
| `type` | `str` | `"union"` | Strategy for site alignment: `"intersect"` or `"union"` |

## Strategies

### Union (default)

```python
ep.pp.set_unite_type(md, type="union")
```

Retains all sites that appear in **any** sample. Sites that are missing in
some samples are handled downstream by the DMC engine using
`min_samples_treatment` and `min_samples_control` guards, which require a
site to be present in at least a minimum number of samples per group to be
tested.

**Advantages:**

- Maximizes the number of testable sites.
- Better suited for experiments with variable coverage across samples.

**Disadvantages:**

- Some sites will have data in only a subset of samples.
- Statistical power varies across sites depending on how many samples cover them.

### Intersect

```python
ep.pp.set_unite_type(md, type="intersect")
```

Retains only sites that are present in **every** sample. This is the most
conservative approach -- it guarantees complete data at every site with no
missing values.

**Advantages:**

- No missing data to handle downstream.
- Statistical tests operate on the same observations across all samples.

**Disadvantages:**

- Can lose a substantial fraction of sites if samples have uneven coverage,
  especially with low-depth sequencing.

## Usage

```python
import epykit as ep

md = ep.read_bismark(
    "samplesheet.csv",
    treatment_group="tumor",
    control_group="normal",
    assembly="hg38",
)

ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)
ep.pp.normalize_coverage(md, method="median")

# Default: union, paired with min_samples_* guards in tl.dmc
ep.pp.set_unite_type(md)
ep.tl.dmc(md, min_samples_treatment=2, min_samples_control=2)

# Or intersect for a complete site matrix
ep.pp.set_unite_type(md, type="intersect")
ep.tl.dmc(md)
```

## Choosing a Strategy

| Scenario | Recommended Strategy |
|----------|---------------------|
| High-depth WGBS, similar coverage across samples | `"intersect"` |
| Variable depth or many samples | `"union"` |
| Small sample size (n < 4 per group) | `"intersect"` (missing data is costly) |
| Large cohort with heterogeneous sequencing | `"union"` |

## Call Order

`set_unite_type` should be called **after** `filter_coverage` and
`normalize_coverage`:

```python
ep.pp.filter_coverage(md)       # Step 1
ep.pp.normalize_coverage(md)    # Step 2
ep.pp.set_unite_type(md)        # Step 3
```

## Next Steps

After marking the unite type, the data is ready for differential methylation
analysis:

```python
ep.tl.dmc(md)
ep.tl.dmr(md)
```

For optional additional preprocessing, see [Smoothing](smooth.md) and
[Aggregate Regions](aggregate-regions.md).

## Deprecated alias: `pp.unite`

The function was originally called `ep.pp.unite()`, which suggested a verb
that performs the union but the function only writes `md.uns["unite"]`. As
of epykit 1.0, the canonical name is `ep.pp.set_unite_type()`; the old name
continues to work as a deprecation wrapper through the 1.x series and is
scheduled for removal in 2.0.

```python
# Still works in 1.x -- emits DeprecationWarning, scheduled to be removed in 2.0.
ep.pp.unite(md, type="union")
```
