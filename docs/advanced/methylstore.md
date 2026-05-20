# Methylstore Internals

The methylstore is the on-disk Parquet structure that backs every `MethylData`
object. It is designed around two principles: never load the whole genome into
memory, and make chromosome-streaming analysis the default path.

## Partitioned Parquet Layout

The store uses Hive-style partitioning with two levels:

```
methylstore/
├── sample=tumor_1/
│   ├── chrom=chr1/
│   │   └── part-0.parquet
│   ├── chrom=chr2/
│   │   └── part-0.parquet
│   └── ...
├── sample=tumor_2/
│   └── ...
└── sample=normal_1/
    └── ...
```

Each `part-0.parquet` file contains a single sample's data for a single
chromosome. This layout enables two key optimisations:

- **Partition pruning** -- querying a specific sample and chromosome reads only
  the relevant parquet file, skipping all others.
- **Predicate pushdown** -- Parquet column statistics and row-group metadata
  allow the engine to skip row groups that cannot match a filter predicate
  (e.g., `pos > 10_000_000`).

## Per-CpG Schema

Each parquet partition stores one row per CpG site with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `pos` | Int32 | Genomic position (0-based) |
| `strand` | str | Strand (`+` or `-`) |
| `context` | str | Cytosine context (`CpG`, `CHG`, `CHH`) |
| `N_meth` | Int32 | Methylated read count |
| `N_unmeth` | Int32 | Unmethylated read count |
| `coverage` | Int32 | Total read count (`N_meth + N_unmeth`) |

Int32 is used instead of Int64 to halve storage and I/O for columns that never
exceed 32-bit range.

## Cache Hierarchy

Preprocessing steps write their results into a `.cache/` directory inside the
store, organised by processing stage:

```
methylstore/
├── .cache/
│   ├── raw/                    # Original ingested data (immutable)
│   ├── filtered/               # After ep.pp.filter_coverage()
│   ├── normalized/             # After ep.pp.normalize_coverage()
│   └── dmc/
│       ├── lr/                 # DMC results for test="lr"
│       │   ├── chrom=chr1.parquet
│       │   ├── chrom=chr2.parquet
│       │   └── ...
│       └── glm/                # DMC results for test="glm"
│           └── ...
├── sample=tumor_1/
│   └── ...
└── .epykit_manifest.json
```

Each stage writes a complete set of partitions. The `MethylData.store` path
always points to the current active stage. Moving backward (e.g., re-running
filtering with different parameters) overwrites the downstream cache.

## DMCStore

`ep.DMCStore` is a persistent per-chromosome DMC result store. Instead of
holding all chromosomes in a single DataFrame, it writes one parquet file per
chromosome and tracks them through a manifest.

```
.cache/dmc/lr/
├── chrom=chr1.parquet
├── chrom=chr2.parquet
├── ...
├── chrom=chrX.parquet
└── .epykit_dmc_manifest.json
```

The manifest (`.epykit_dmc_manifest.json`) records:

- Test name and parameters
- List of completed chromosomes
- Row counts per chromosome
- Timestamp of each chromosome's completion

### Key Methods and Properties

| Member | Description |
|--------|-------------|
| `path` | Root directory of the DMCStore |
| `manifest` | Parsed manifest dict |
| `to_dataframe()` | Concatenate all chromosomes into a single Polars DataFrame |
| `iter_chromosomes()` | Yield one DataFrame per chromosome (streaming) |

### Why DMCStore Exists

The per-chromosome layout enables two operations at O(largest chromosome) memory
instead of O(genome):

1. **Streaming BH correction** -- `ep.apply_multiple_testing_correction()` can
   scan all chromosome files in a first pass to count total tests, then apply
   BH in a second streaming pass without holding the full result in memory.

2. **DMR calling** -- `ep.call_dmr_chain_merge()` processes one chromosome at a
   time. It never needs more than one chromosome's DMC results in memory.

## Manifest

The top-level `.epykit_manifest.json` file tracks the overall state of the
methylstore:

```json
{
  "version": "0.7.1",
  "assembly": "hg38",
  "samples": ["tumor_1", "tumor_2", "normal_1", "normal_2"],
  "chromosomes": ["chr1", "chr2", "...", "chrX"],
  "stages_completed": ["ingest", "filter", "normalize", "unite"],
  "current_stage": "unite",
  "parameters": {
    "filter": {"lo_count": 10, "hi_perc": 99.9},
    "normalize": {"method": "median"},
    "unite": {"type": "union"}
  }
}
```

This manifest enables **checkpoint/resume**: if a long-running pipeline is
interrupted, the next run reads the manifest, identifies which stages and
chromosomes have already been processed, and continues from where it left off.
The `resumable=True` parameter on `ep.tl.dmc()` uses this mechanism.

## Design Principles

1. **Never load the whole genome.** Every function in epykit processes data one
   chromosome at a time. The methylstore layout makes this the natural access
   pattern.

2. **Predicate pushdown for queries.** `ep.query.*` functions push filters
   (chromosome, position range, strand) into the Parquet reader so only
   matching row groups are read from disk.

3. **Partition pruning.** Hive partitioning on `sample=` and `chrom=` means
   that reading one sample's chromosome never touches another sample's files.

4. **Immutable raw data.** The `raw/` cache preserves the original ingested
   data. Preprocessing steps write new partitions rather than modifying
   existing ones.

5. **Deterministic fingerprinting.** Cache keys are derived from the store
   path, sample lists, parameters, and chromosomes. Changing any input
   invalidates downstream caches automatically.
