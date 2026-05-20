# read_bismark

`ep.read_bismark()` reads standard Bismark coverage files and writes them into a
partitioned Parquet methylstore.

## Input Format

Bismark `.cov` (or `.cov.gz`) files contain six tab-separated columns with **0-based**
coordinates:

```
chrom    start    end    methylation_percent    count_methylated    count_unmethylated
```

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | str | Chromosome name (e.g. `chr1`) |
| `start` | int | 0-based start position |
| `end` | int | 0-based end position |
| `methylation_percent` | float | Percent methylation at the site (0--100) |
| `count_methylated` | int | Number of methylated reads |
| `count_unmethylated` | int | Number of unmethylated reads |

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `samplesheet` | `str` | required | Path to a CSV file mapping sample names, file paths, and group labels |
| `treatment_group` | `str` | required | Group label in the samplesheet that identifies treatment samples |
| `control_group` | `str` | required | Group label in the samplesheet that identifies control samples |
| `assembly` | `str` | required | Genome assembly name (e.g. `"hg38"`, `"mm10"`) |
| `store_dir` | `str` | `"methylstore"` | Directory where the partitioned Parquet store will be written |
| `context` | `str` | `"CpG"` | Cytosine context to retain (`"CpG"`, `"CHG"`, or `"CHH"`) |

## Usage

### 1. Prepare a samplesheet

Create a CSV file (`samplesheet.csv`) with your sample layout:

```csv
sample_name,file_path,group
tumor_1,/data/bismark/tumor_1.cov.gz,tumor
tumor_2,/data/bismark/tumor_2.cov.gz,tumor
normal_1,/data/bismark/normal_1.cov.gz,normal
normal_2,/data/bismark/normal_2.cov.gz,normal
```

### 2. Read the data

```python
import epykit as ep

md = ep.read_bismark(
    samplesheet="samplesheet.csv",
    treatment_group="tumor",
    control_group="normal",
    assembly="hg38",
    store_dir="results/methylstore",
    context="CpG",
)
```

## What It Produces

After reading, the function creates a partitioned Parquet directory under `store_dir`:

```
results/methylstore/
├── sample=tumor_1/
│   ├── chrom=chr1.parquet
│   ├── chrom=chr2.parquet
│   └── ...
├── sample=tumor_2/
│   └── ...
├── sample=normal_1/
│   └── ...
└── sample=normal_2/
    └── ...
```

Each Parquet partition contains columns for `start`, `end`, `count_methylated`,
`count_unmethylated`, and `coverage` (the sum of the two counts).

The returned `MethylData` object (`md`) holds:

- `md.obs` -- A DataFrame of sample-level metadata (name, group, file path).
- `md.store` -- A reference to the on-disk Parquet store for lazy, memory-efficient queries.
- `md.uns` -- A dictionary for unstructured metadata (assembly, context, parameters).

## Next Steps

After reading, proceed to preprocessing:

```python
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)
```

See the [Preprocessing overview](../preprocessing/index.md) for the recommended workflow.
