# Load and Save

epykit can persist a `MethylData` object to disk and reload it later, preserving sample
metadata, analysis results, and the connection to the underlying methylstore.

## Saving

### md.save()

```python
md.save("results/my_analysis")
```

This writes the current state of the `MethylData` object to the specified directory. If
the directory does not exist, it is created automatically.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `path` | `str` | required | Directory path where the analysis will be saved |

### What Gets Saved

The save operation writes multiple files capturing different parts of the analysis state:

```
results/my_analysis/
├── obs.parquet                        # Sample-level metadata (names, groups, QC stats)
├── varm_dmc_lr_annotated.parquet      # Annotated DMC results (if computed)
├── uns_dmr.parquet                    # DMR results (if computed)
└── methyldata.json                    # Configuration: assembly, context, store path, parameters
```

- **`obs.parquet`** -- The sample observation table containing sample names, group
  assignments, and any QC or summary statistics added during preprocessing.
- **`varm_dmc_lr_annotated.parquet`** -- Per-site differential methylation results with
  genomic annotations. Present only after running a DMC analysis.
- **`uns_dmr.parquet`** -- Differentially methylated region (DMR) results. Present only
  after running a DMR analysis.
- **`methyldata.json`** -- A JSON file recording the assembly, cytosine context, store
  directory path, and other configuration needed to reconstruct the object.

The underlying Parquet methylstore (the per-sample, per-chromosome files) is **not**
duplicated. The saved JSON records the path to the existing store, so the store directory
must remain accessible when reloading.

## Loading

### ep.load()

```python
import epykit as ep

md = ep.load("results/my_analysis")
```

This reconstructs the `MethylData` object from the saved files. The methylstore is
reconnected using the path recorded in `methyldata.json`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `path` | `str` | required | Path to the directory created by `md.save()` |

### What Gets Restored

- `md.obs` -- Sample metadata from `obs.parquet`.
- `md.store` -- Reconnected to the original Parquet methylstore.
- `md.varm` -- DMC results from `varm_dmc_lr_annotated.parquet` (if present).
- `md.uns` -- Configuration and DMR results (if present).

## Example Workflow

```python
import epykit as ep

# Initial analysis session
md = ep.read_bismark(
    samplesheet="samplesheet.csv",
    treatment_group="tumor",
    control_group="normal",
    assembly="hg38",
)

ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)

# Save progress before running differential analysis
md.save("results/preprocessed")

# Later, in a new session
md = ep.load("results/preprocessed")

# Continue from where you left off
ep.tl.dmc(md)
ep.tl.dmr(md)

# Save the complete analysis
md.save("results/full_analysis")
```

## Notes

- **Store path stability** -- The methylstore directory path is stored as-is. If you move
  the store after saving, update the path in `methyldata.json` or re-read the data.
- **Incremental saves** -- Calling `md.save()` again overwrites the previous save at that
  path. Use different directory names to keep snapshots of intermediate stages.
- **Portability** -- To share an analysis, copy both the save directory and the methylstore
  directory. The recipient may need to update the store path in `methyldata.json` if the
  directory structure differs.
