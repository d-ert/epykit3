# read_combined_strand_bed

`ep.read_combined_strand_bed()` reads 12-column strand-collapsed methylation BED files
and writes them into a partitioned Parquet methylstore.

## Input Format

Combined-strand BED files contain 12 tab-separated columns that merge forward- and
reverse-strand methylation calls into a single record per CpG site. This format is
produced by tools that collapse strand-specific calls before output.

The expected columns are:

| # | Column | Type | Description |
|---|--------|------|-------------|
| 1 | `chrom` | str | Chromosome |
| 2 | `start` | int | 0-based start position |
| 3 | `end` | int | End position |
| 4 | `name` | str | Site identifier or `.` |
| 5 | `score` | int | Score (often 0 or 1000) |
| 6 | `strand` | str | Strand (`.` for combined) |
| 7 | `thick_start` | int | Display start (typically equals `start`) |
| 8 | `thick_end` | int | Display end (typically equals `end`) |
| 9 | `item_rgb` | str | RGB color string |
| 10 | `coverage` | int | Total read coverage at the site |
| 11 | `methylation_percent` | float | Percent methylation (0--100) |
| 12 | Additional field | varies | Tool-specific extra column |

The reader extracts `chrom`, `start`, `end`, `coverage`, and `methylation_percent`, then
derives `count_methylated` and `count_unmethylated` from the coverage and percentage.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `samplesheet` | `str` | required | Path to a CSV mapping sample names, file paths, and group labels |
| `treatment_group` | `str` | required | Group label identifying treatment samples |
| `control_group` | `str` | required | Group label identifying control samples |
| `assembly` | `str` | required | Genome assembly name (e.g. `"hg38"`, `"mm10"`) |
| `store_dir` | `str` | `"methylstore"` | Directory for the partitioned Parquet store |
| `context` | `str` | `"CpG"` | Cytosine context to retain |

## Usage

### 1. Prepare a samplesheet

```csv
sample_name,file_path,group
case_1,/data/beds/case_1.combined.bed.gz,case
case_2,/data/beds/case_2.combined.bed.gz,case
ctrl_1,/data/beds/ctrl_1.combined.bed.gz,control
ctrl_2,/data/beds/ctrl_2.combined.bed.gz,control
```

### 2. Read the data

```python
import epykit as ep

md = ep.read_combined_strand_bed(
    samplesheet="samplesheet.csv",
    treatment_group="case",
    control_group="control",
    assembly="hg38",
    store_dir="results/methylstore",
)
```

## Output

The output is the same partitioned Parquet methylstore and `MethylData` object produced by
all other readers. See [read_bismark](read-bismark.md) for details on the store layout.

## When to Use This Reader

Use `ep.read_combined_strand_bed()` when your methylation caller outputs 12-column
strand-collapsed BED files. If your files follow the simpler 6-column Bismark or
MethylDackel format, use the corresponding reader instead.

## Next Steps

After reading, proceed to [Preprocessing](../preprocessing/index.md):

```python
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)
```
