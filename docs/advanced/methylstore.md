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

- `test`: the engine name
- `chroms`: one entry per chromosome with `name` and `n_sites`, in submission order
- `total_sites`
- `input_sig`: SHA-256 of the inputs that affect the result (store path, sample lists, test, parameters). A rerun with the same signature reuses the store instead of recomputing.
- `bh_qvalues_applied`, `bh_qvalue_col`, `bh_method`: set once `apply_multiple_testing_correction` has written q-values back into the per-chromosome files
- `epykit_version`, `completed_at`

### Key Methods and Properties

| Member | Description |
|--------|-------------|
| `path` | Root directory of the DMCStore |
| `manifest` | Parsed manifest dict |
| `chroms()` | Chromosome names in manifest order |
| `total_sites`, `bh_applied` | Row count across all chromosomes; whether q-values have been written back |
| `iter_chroms(columns=None)` | Yield `(chrom, DataFrame)` pairs one chromosome at a time (streaming) |
| `read_chrom(chrom)`, `scan_chrom(chrom)` | Eager read or lazy scan of one chromosome |
| `update_chrom(chrom, df)` | Atomically rewrite one chromosome file (used by BH correction) |
| `to_dataframe()` | Concatenate all chromosomes into a single Polars DataFrame (whole-genome eager load; use sparingly) |

### Why DMCStore Exists

The per-chromosome layout enables two operations at O(largest chromosome) memory
instead of O(genome):

1. **Streaming BH correction** -- `epykit.dmc.apply_multiple_testing_correction()` can
   scan all chromosome files in a first pass to count total tests, then apply
   BH in a second streaming pass without holding the full result in memory.

2. **DMR calling** -- `ep.call_dmr_chain_merge()` processes one chromosome at a
   time. It never needs more than one chromosome's DMC results in memory.

## Manifests

Two manifest layers live next to the data. Both are written through the
helpers in `_cache.py`.

**Per-sample, per-step manifests** (`.epykit_raw_manifest.json`,
`.epykit_filter_manifest.json`, `.epykit_normalize_manifest.json` inside each
`sample=<id>/` directory) fingerprint the upstream input and the step
parameters. `read_bismark`, `pp.filter_coverage` and `pp.normalize_coverage`
skip a sample whose manifest already matches.

**The pipeline manifest** (`.epykit_manifest.json` at the analysis root)
records completed stages for the checkpoint/resume API:

```json
{
  "epykit_version": "1.0.0",
  "stages": [
    {
      "name": "dmc_lr",
      "params": {"test": "lr", "unite": true, "dispersion": "eb"},
      "input_sig": "<sha256 of the store fingerprint, sample lists and parameters>",
      "output_path": ".epykit_results/dmc_lr.parquet",
      "completed_at": "2026-06-09T12:00:00Z",
      "extra": {"n_sites": 21873452}
    }
  ]
}
```

`ep.tl.dmc(md, resumable=True)` hashes its inputs, looks for a stage with the
same name and `input_sig`, and loads the sidecar Parquet instead of
recomputing. `MethylData.completed_stages` lists the recorded stages and
`MethylData.resume_from(stage)` re-hydrates one. Without `resumable=True` the
pipeline manifest is neither read nor written.

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
