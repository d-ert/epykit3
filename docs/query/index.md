# Query

Random-access queries on the partitioned Parquet store. No extra index is needed --
epykit uses Polars predicate pushdown on Parquet row-group statistics to skip
irrelevant partitions and row groups, making region lookups fast even on whole-genome
stores with billions of rows.

All query functions operate on `md.store` (the store path as a string or Path), not on
the MethylData object itself. This means you can query a store without loading the full
MethylData object into memory.

---

## ep.query.query_region

Retrieve all CpG measurements within a single genomic region.

```python
ep.query.query_region(store, chrom, start, end, samples=None)
```

**Parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `store` | `str` or `Path` | Path to the methylation Parquet store |
| `chrom` | `str` | Chromosome name (e.g., `"chr7"`) |
| `start` | `int` | Region start position (0-based, inclusive) |
| `end` | `int` | Region end position (0-based, exclusive) |
| `samples` | `list[str]` or `None` | Optional subset of sample IDs. `None` returns all samples |

**Returns**

A Polars DataFrame with columns:

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | `str` | Sample identifier |
| `chrom` | `str` | Chromosome |
| `pos` | `int` | CpG position |
| `strand` | `str` | `"+"` or `"-"` |
| `N_meth` | `int` | Methylated read count |
| `coverage` | `int` | Total read coverage |
| `beta` | `float` | Methylation fraction (`N_meth / coverage`) |

**Example**

```python
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")

# Query a region around BRAF
df = ep.query.query_region(md.store, "chr7", 140_453_000, 140_500_000)
print(df)
# shape: (1284, 7)
# +-----------+------+-----------+--------+--------+----------+-------+
# | sample_id | chrom| pos       | strand | N_meth | coverage | beta  |
# +-----------+------+-----------+--------+--------+----------+-------+
# | tumor_1   | chr7 | 140453012 | +      | 45     | 52       | 0.865 |
# | tumor_1   | chr7 | 140453089 | +      | 3      | 48       | 0.063 |
# | ...       | ...  | ...       | ...    | ...    | ...      | ...   |
# +-----------+------+-----------+--------+--------+----------+-------+

# Query for specific samples only
df = ep.query.query_region(
    md.store, "chr7", 140_453_000, 140_500_000,
    samples=["tumor_1", "normal_1"]
)
```

---

## ep.query.query_regions

Batch query for multiple genomic regions. More efficient than calling `query_region`
in a loop because it groups regions by chromosome partition.

```python
ep.query.query_regions(store, regions_df, samples=None)
```

**Parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `store` | `str` or `Path` | Path to the methylation Parquet store |
| `regions_df` | `polars.DataFrame` or `pandas.DataFrame` | Regions with columns `chrom`, `start`, `end`, and optionally `region_id` |
| `samples` | `list[str]` or `None` | Optional subset of sample IDs |

**Returns**

A Polars DataFrame with the same columns as `query_region`, plus:

| Column | Type | Description |
|--------|------|-------------|
| `region_id` | `str` or `int` | Identifier linking each row back to its source region. Uses the `region_id` column from the input if present, otherwise the row index |

**Example**

```python
import polars as pl
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")

# Define regions of interest
regions = pl.DataFrame({
    "chrom": ["chr7", "chr17", "chr1"],
    "start": [140_453_000, 7_571_000, 11_166_000],
    "end":   [140_500_000, 7_590_000, 11_195_000],
    "region_id": ["BRAF", "TP53", "MTOR"],
})

df = ep.query.query_regions(md.store, regions, samples=["tumor_1", "normal_1"])
print(df.group_by("region_id").len())
# shape: (3, 2)
# +-----------+-----+
# | region_id | len |
# +-----------+-----+
# | BRAF      | 214 |
# | TP53      | 89  |
# | MTOR      | 156 |
# +-----------+-----+
```

---

## ep.query.query_sites

Query exact CpG positions. Useful for targeted lookups such as epigenetic clock CpGs,
validation panel sites, or specific CpGs from a publication.

```python
ep.query.query_sites(store, sites_df, samples=None)
```

**Parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `store` | `str` or `Path` | Path to the methylation Parquet store |
| `sites_df` | `polars.DataFrame` or `pandas.DataFrame` | Sites with columns `chrom` and `pos` |
| `samples` | `list[str]` or `None` | Optional subset of sample IDs |

**Returns**

A Polars DataFrame with the same schema as `query_region`.

**Example**

```python
import polars as pl
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")

# Horvath clock CpG sites (example subset)
clock_sites = pl.DataFrame({
    "chrom": ["chr1", "chr2", "chr6", "chr11"],
    "pos":   [15_865_264, 107_050_538, 11_044_894, 67_418_096],
})

df = ep.query.query_sites(md.store, clock_sites)
print(df)
```

---

## Use Cases

### Targeted validation

Look up specific CpGs reported in a publication to verify concordance:

```python
published_sites = pl.read_csv("published_dmcs.csv")  # chrom, pos columns
df = ep.query.query_sites(md.store, published_sites)
```

### Region aggregation

Combine `query_regions` with downstream aggregation:

```python
df = ep.query.query_regions(md.store, promoter_regions, samples=["tumor_1", "normal_1"])
means = df.group_by(["region_id", "sample_id"]).agg(pl.col("beta").mean())
```

### Lightweight access without MethylData

Since queries operate on the store path, you can skip creating a MethylData object
entirely for quick lookups:

```python
df = ep.query.query_region("methylstore/", "chr7", 140_453_000, 140_500_000)
```

---

## Performance Notes

- Queries exploit the `chrom=` partition structure of the Parquet store, so only the
  relevant chromosome directory is scanned.
- Within each partition, Polars predicate pushdown uses row-group min/max statistics
  on the `pos` column to skip row groups that fall outside the query range.
- For batch queries (`query_regions`), regions on the same chromosome are grouped
  together to minimize the number of file opens.
- No additional index files are needed beyond the standard Parquet metadata.
