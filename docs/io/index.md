# I/O Overview

epykit provides several reader functions for ingesting whole-genome bisulfite sequencing
(WGBS) data from common output formats. Each reader parses raw methylation calls, validates
the samplesheet, and writes the data into a partitioned Parquet methylstore for fast
downstream analysis.

## Supported Readers

| Function | Input Format | Description |
|----------|-------------|-------------|
| `ep.read_bismark()` | Bismark `.cov[.gz]` | Standard Bismark coverage files (6-column, 0-based) |
| `ep.read_methyldackel()` | MethylDackel `.bedGraph[.gz]` | MethylDackel bedGraph output (auto-skips track header) |
| `ep.read_combined_strand_bed()` | 12-column BED | Strand-collapsed methylation BED files |
| `ep.read_nfcore_methylseq()` | nf-core run directory | Direct ingestion from an nf-core/methylseq pipeline run |
| `ep.load()` | Saved analysis directory | Reconstruct a previously saved `MethylData` object |

## How Readers Work

All readers follow the same general workflow:

1. **Parse the samplesheet** -- A CSV file mapping sample names to file paths and group
   assignments (treatment vs. control).
2. **Read and validate** -- Each coverage file is read, checked for the expected column
   layout, and filtered for the requested cytosine context (default: CpG).
3. **Write to the methylstore** -- Data is written as partitioned Parquet files under the
   specified `store_dir`, organized by sample and chromosome for efficient random access.
4. **Return a `MethylData` object** -- The returned object holds sample metadata (`obs`),
   a pointer to the on-disk store, and slots for downstream results.

## Samplesheet Format

All readers except `ep.load()` require a samplesheet CSV with at least three columns:

```
sample_name,file_path,group
tumor_1,/data/bismark/tumor_1.cov.gz,tumor
tumor_2,/data/bismark/tumor_2.cov.gz,tumor
normal_1,/data/bismark/normal_1.cov.gz,normal
normal_2,/data/bismark/normal_2.cov.gz,normal
```

- `sample_name` -- Unique identifier for each sample.
- `file_path` -- Absolute or relative path to the methylation calls file.
- `group` -- Group label used to define treatment and control via the `treatment_group`
  and `control_group` parameters.

## Choosing a Reader

- **Bismark pipeline** -- Use `ep.read_bismark()`. This is the most common entry point.
- **MethylDackel / PileOMeth** -- Use `ep.read_methyldackel()`. Handles the bedGraph
  track header automatically.
- **Strand-collapsed BED** -- Use `ep.read_combined_strand_bed()` when your upstream tool
  produces 12-column combined-strand BED output.
- **nf-core/methylseq** -- Use `ep.read_nfcore_methylseq()` to read directly from the
  pipeline output directory without manually building a samplesheet.
- **Resuming work** -- Use `ep.load()` to reload a previously saved analysis.

## Next Steps

After reading data, proceed to [Preprocessing](../preprocessing/index.md) to filter,
normalize, and unite your methylation sites before differential analysis.
