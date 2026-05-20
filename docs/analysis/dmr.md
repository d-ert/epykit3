# DMR Calling

`ep.tl.dmr(md)` aggregates per-CpG test results into differentially methylated
regions (DMRs). Four methods are available, each suited to different analysis
scenarios. Results are stored in `md.uns["dmr"]`.

## Basic Usage

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)
ep.tl.dmc(md)

# Default: chain_merge (DSS-style)
ep.tl.dmr(md)

dmrs = md.uns["dmr"]
print(f"Found {len(dmrs)} DMRs")
```

## Methods

### chain_merge (default)

DSS `callDMR`-style contiguous CpG chaining. Walks through the DMC results,
identifies CpGs that pass the significance and effect-size thresholds, and
chains consecutive significant CpGs whose gap is within `dis_merge_bp` into
candidate regions. Regions are then filtered by minimum CpG count, minimum
span length, and minimum fraction of significant CpGs.

Requires a prior `ep.tl.dmc(md)` call.

```python
ep.tl.dmr(md, method="chain_merge")

# With a preset
ep.tl.dmr(md, method="chain_merge", preset="strict")
```

#### Presets

Three parameter bundles are available for `chain_merge`. Any explicit kwarg
passed alongside `preset` overrides the bundled value.

| Preset | alpha | min_abs_meth_diff | dis_merge_bp | min_cpgs | pct_sig | minlen_bp | Use case |
|--------|-------|-------------------|--------------|----------|---------|-----------|----------|
| `"strict"` | 1e-6 | 0.20 | 250 | 5 | 0.5 | 100 | Validation-ready DMRs, low FP tolerance |
| `"default"` | 1e-4 | 0.10 | 500 | 3 | 0.5 | 50 | Balanced starting point for general WGBS |
| `"permissive"` | 1e-4 | 0.05 | 1000 | 3 | 0.5 | 50 | Recall-oriented, exploratory analyses |

```python
# Strict preset for validation experiments
ep.tl.dmr(md, method="chain_merge", preset="strict")

# Permissive preset with custom merge distance
ep.tl.dmr(md, method="chain_merge", preset="permissive", dis_merge_bp=500)
```

!!! tip "Tuning chain_merge"
    If recall is too low, loosen `dis_merge_bp` first (highest Pareto
    leverage), then `min_cpgs`. If precision is too low, tighten `alpha`
    or increase `min_abs_meth_diff`.

### tile

Read-pooled tile aggregation. Divides the genome into fixed-width tiles,
pools read counts within each tile, and runs a single test per tile. Does
**not** require a prior DMC call -- it goes directly to the methylstore.

```python
ep.tl.dmr(md, method="tile", tile_size_bp=500, min_cpgs_per_tile=5)
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `tile_size_bp` | 1000 | Tile width in base pairs |
| `min_cpgs_per_tile` | 5 | Minimum CpGs per tile |
| `merge_adjacent` | True | Merge adjacent significant tiles |

### sliding_window

Per-CpG signed Stouffer's Z combining within overlapping windows. A legacy
method that operates on the in-memory DMC table. Faster (no extra I/O) but
lower-power than tile-based or chain-merge aggregation at typical WGBS coverage.

Requires a prior `ep.tl.dmc(md)` call.

```python
ep.tl.dmr(md, method="sliding_window", window_bp=500, step_bp=250)
```

### hmm

HMM segmentation over `meth_diff`. Fits a two-state (hyper/hypo) hidden Markov
model to the per-CpG methylation differences and segments the genome into
contiguous DMRs. Useful when the spatial structure of methylation changes is
complex.

Requires a prior `ep.tl.dmc(md)` call.

```python
ep.tl.dmr(md, method="hmm")
```

## Empirical FDR

Permutation-based FDR is available for the `tile` method:

```python
ep.tl.dmr(md, method="tile", empirical_fdr=True, n_perm=100)
```

!!! note
    Empirical FDR is not supported with covariate designs because label
    shuffling invalidates the stratified design.

## Output Columns

The DMR result DataFrame (`md.uns["dmr"]`) contains:

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | str | Chromosome |
| `start` | int | Region start (0-based) |
| `end` | int | Region end |
| `n_cpgs` | int | Number of CpGs in the region |
| `mean_meth_diff` | float | Mean per-CpG methylation difference |
| `mean_qvalue` | float | Mean per-CpG q-value within the region |
| `dmr_type` | str | Direction: `"hyper"`, `"hypo"`, or `"mixed"` |
| `pvalue` | float | Region-level p-value |
| `qvalue` | float | Region-level q-value |

## DMR Calling Diagnostics

`ep.tl.diagnose_dmr_calling()` compares your DMR results against a reference
set (e.g., from a published study or another pipeline) and classifies each
missed reference DMR into an actionable category.

```python
import polars as pl

reference = pl.read_csv("published_dmrs.bed", separator="\t",
                        has_header=False,
                        new_columns=["chrom", "start", "end"])

diag = ep.tl.diagnose_dmr_calling(md, reference_dmrs=reference)
print(diag["summary"])
```

The five diagnostic buckets:

| Bucket | Meaning | Fix |
|--------|---------|-----|
| `SUCCESS_OVERLAP` | Our DMR overlaps the reference DMR | None needed |
| `H1_NO_CPGS` | Zero united CpGs inside the reference region | Relax coverage filter or use `unite(type="union")` |
| `H2_NO_SIG_CPGS` | CpGs present but none reach q < 0.05 | Need a more powerful test statistic |
| `H3a_WEAK_ALPHA` | Sig CpGs at q < 0.05 but not at chain-merge alpha | Loosen `alpha` (e.g., `preset="permissive"`) |
| `H3b_STRUCTURE` | Sig CpGs exist but chain-merge dropped the candidate | Loosen `dis_merge_bp` first, then `min_cpgs` |

```python
# Example output:
# DMR-calling diagnostic on 150 reference DMRs:
#   alpha_threshold = 1e-05  (matches chain-merge alpha)
#
#   SUCCESS_OVERLAP       98 (65.3%)  -- already recovered (no fix needed)
#   H1_NO_CPGS             5 ( 3.3%)  -- no CpGs present -> coverage/unite issue
#   H2_NO_SIG_CPGS        12 ( 8.0%)  -- no sig CpGs at q<0.05 -> need better test stat
#   H3a_WEAK_ALPHA        20 (13.3%)  -- sig CpGs but not at alpha -> loosen alpha
#   H3b_STRUCTURE          15 (10.0%)  -- sig CpGs exist but chain-merge dropped
```

## Parameters

### chain_merge Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `preset` | str | None | Parameter bundle: `"strict"`, `"default"`, `"permissive"` |
| `alpha` | float | 0.05 | Per-CpG significance threshold |
| `min_abs_meth_diff` | float | 0.1 | Minimum absolute methylation difference |
| `dis_merge_bp` | int | 500 | Maximum gap (bp) between consecutive sig CpGs |
| `min_cpgs` | int | 3 | Minimum CpGs in a region |
| `pct_sig` | float | 0.5 | Minimum fraction of significant CpGs |
| `minlen_bp` | int | 50 | Minimum region span (bp) |
| `use_q_for_sig` | bool | False | Use q-value instead of p-value for significance |
| `min_mean_qvalue` | float | 0.05 | Post-hoc DMR-level q-value filter |

### Shared Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | MethylData | required | Analysis object |
| `method` | str | `"chain_merge"` | DMR method: `"chain_merge"`, `"tile"`, `"sliding_window"`, `"hmm"` |
| `chromosomes` | list | None | Restrict to specific chromosomes |
| `backend` | str | `"sequential"` | Execution backend (tile method only) |
| `empirical_fdr` | bool | False | Permutation FDR (tile method only) |
| `n_perm` | int | 100 | Number of permutations |
