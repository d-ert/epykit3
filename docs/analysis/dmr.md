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
ep.pp.set_unite_type(md)
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

`chain_merge` is also available from the CLI (added in 1.0), so you don't
need to drop into Python just to run the DSS-style caller:

```bash
epykit dmr \
    --methylstore methylstore/ \
    --method chain_merge \
    --dmc-results md/varm/dmc_lr.parquet \
    --preset strict
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
ep.tl.dmr(md, method="chain_merge", preset="permissive", dis_merge_bp=2000)
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
| `canonical_only` | False | Keep only the fixed human-style chromosome set (`1`-`22`, `X`, `Y`, `M`/`MT`, with or without `chr`) of the auto-detected partitions; an explicit `chromosomes=` list wins |

`canonical_only` is a tile-only option because the tile caller enumerates
the store partitions itself. It filters the auto-detected chromosome list
before the tile test and the BH correction, uses the same resolved list for
the observed tiles and every `empirical_fdr` permutation, logs one INFO line
naming the dropped contigs, and is recorded in `md.uns["dmr_params"]`. An
explicit `chromosomes=` list, including an empty one, is used verbatim.
`chain_merge`, `sliding_window` and `segment` work on the DMC table and
inherit the chromosome universe of the upstream `ep.tl.dmc()` run, so they
raise `ValueError` on `canonical_only=True`. Filter upstream instead: ingest
with `canonical_only=True` (see [read_bismark](../io/read-bismark.md#canonical-chromosomes-only))
or restrict `ep.tl.dmc(md, chromosomes=...)`.

### sliding_window

Per-CpG signed Stouffer's Z combining within overlapping windows. A legacy
method that operates on the in-memory DMC table. Faster (no extra I/O) but
lower-power than tile-based or chain-merge aggregation at typical WGBS coverage.

Requires a prior `ep.tl.dmc(md)` call.

```python
ep.tl.dmr(md, method="sliding_window", window_bp=500, step_bp=250)
```

### segment

Rule-based 3-state segmentation over `meth_diff`. Walks the per-CpG result
and segments the genome into contiguous hyper / hypo / null runs, then
filters by minimum CpG count and minimum span. Useful when the spatial
structure of methylation changes is more complex than chain-merge's
"contiguous significant CpGs" assumption.

Requires a prior `ep.tl.dmc(md)` call.

```python
ep.tl.dmr(md, method="segment")
```

!!! note "Renamed from `hmm` in 1.0"
    Previous releases shipped this engine as `method="hmm"`; the name was
    misleading because the implementation is rule-based, not an HMM. The
    underlying segmenter is unchanged. `method="hmm"` was deprecated in
    0.7.5 with `FutureWarning` and now raises `ValueError` -- use
    `method="segment"`.

## Empirical FDR

The region-level `qvalue` (tile) and `combined_qvalue` (chain_merge,
sliding_window, segment) are BH corrections of asymptotic region p-values.
Adjacent WGBS CpGs are positively correlated and real coverage is
overdispersed, so these q-values rank regions well but are not a calibrated
region FDR. Permutation FDR re-runs the caller on shuffled treatment and
control labels and compares the observed regions with the shuffled (decoy)
regions. It is available for `tile` and `chain_merge`:

```python
ep.tl.dmr(md, method="tile", empirical_fdr=True, n_perm=100, perm_seed=42)
ep.tl.dmr(md, method="chain_merge", empirical_fdr=True, n_perm=100, perm_seed=42)
# md.uns["dmr"] gains empirical_pvalue, empirical_qvalue and empirical_fdr_set
```

Threshold `empirical_qvalue`, not `qvalue` or `combined_qvalue`, for FDR
control.

### `fdr_method`: how the empirical FDR is computed

Both callers accept `fdr_method`. The default reproduces the numbers earlier
releases computed; `"region"` is opt-in.

| `fdr_method` | Construction | `empirical_pvalue` | `empirical_qvalue` | `empirical_fdr_set` |
|---|---|---|---|---|
| `"max_t"` (default) | Westfall-Young min-P. Each permutation contributes its genome-wide minimum null p-value. Controls the family-wise error rate and is very conservative at genome scale under realistic dispersion. | Fraction of permutations whose minimum null p-value is at or below the observed p-value, with a pseudo-count. | BH transform of `empirical_pvalue`. | NaN. |
| `"region"` | Count-ratio target-decoy FDR (BSmooth, SAM). At each threshold `t`, the mean number of decoy survivors with `p <= t` divided by the number of observed survivors with `p <= t`, made monotone with a suffix minimum and clipped to `[0, 1]`. The overdispersion that inflates both counts cancels in the ratio. | Fraction of pooled decoy survivors with `p` at or below the observed p-value. A diagnostic, not a calibrated per-region p-value. | The count-ratio estimate at the region's own threshold. | `min(mean decoy survivors / observed survivors, 1)`, the set-level estimate. Also stored in `md.uns["dmr_params"]["empirical_fdr_set"]`. |

The two modes treat permutations differently:

- `max_t` counts every permutation that produced at least one region,
  including assignments equal to the observed split or its mirror. Failed
  permutations and permutations with no region leave the denominator.
- `region` excludes assignments equal to the observed split or its mirror
  (their statistics are the observed ones, not null draws) and failed
  permutations. A permutation that ran cleanly and produced no region counts
  as zero decoys. When no usable assignment remains, all three columns are
  NaN and a `UserWarning` is emitted.
- In both modes an observed region with a non-finite p-value gets NaN in
  `empirical_pvalue` and `empirical_qvalue`.
- `n_perm` must be positive; the method and count are validated before the
  first permutation runs.

```python
ep.tl.dmr(md, method="tile", empirical_fdr=True, n_perm=100, fdr_method="region")
sig = md.uns["dmr"].filter(pl.col("empirical_qvalue") < 0.05)
md.uns["dmr_params"]["empirical_fdr_set"]  # e.g. 0.13
```

`md.uns["dmr_params"]` records `empirical_fdr`, `n_perm`, `perm_seed`,
`fdr_method` and `empirical_fdr_set` (`None` when the estimate is NaN).

!!! warning "Small cohorts"
    Permutation FDR needs enough distinct label assignments. At fewer than
    four samples per group, `fdr_method="region"` emits a `UserWarning`: few
    shuffles exist, and draws adjacent to the true split leak signal into
    the null, so the estimate is conservative (biased high). Read it as a
    floor.

!!! note "Calibration evidence"
    The engine hash gate under `benchmark/` covers selected per-CpG `lr`
    output only. It does not establish that either permutation construction
    is calibrated on your data. The
    [design note](../review/2026-06-08-region-empirical-fdr-design.md)
    records the real-data comparison that motivated region mode.

### chain_merge permutations

`method="chain_merge"` replays the observed analysis for every permutation.
Each permutation recomputes the per-CpG DMC on the shuffled labels with the
knobs recorded in `md.uns["dmc"]` (test, `unite`, minimum sample counts,
dispersion, reference, smoothing, separation fallback), over the same
chromosome universe as the observed DMC, applies the observed
multiple-testing method, chain-merges with the same preset and knobs, and
applies the same `min_mean_qvalue` filter. The surviving `combined_pvalue`
values are the decoys.

Requirements and limits:

- The DMC must be a two-group `lr`, `welch_t` or `fisher` run from
  `ep.tl.dmc`. GLM, formula and contrast designs, and `use_smoothed=True`
  are rejected before any permutation.
- An explicit `chromosomes=` must equal the observed DMC universe. To
  restrict the scan, rerun `ep.tl.dmc(md, chromosomes=...)` first.
- `empirical_strata` must name an `md.obs` column that covers every
  treatment and control sample. A missing or partial column raises instead
  of falling back to an unrestricted shuffle. The same rule applies to the
  tile harness.
- Raw `pvalue` and `qvalue` drive the chain_merge gate, as in the observed
  run. Neighbour-combined columns from `neighbour_combine=True` are never
  substituted.
- Each permutation recomputes a genome-wide DMC, so a whole-genome run with
  `n_perm=100` takes hours. The cost is logged once at INFO. `perm_n_jobs`
  parallelises permutations. Each permutation writes its DMC to a private
  temporary directory and never touches the observed store.
- The CLI (`epykit dmr --empirical-fdr`) stays tile-only. `sliding_window`
  and `segment` raise `NotImplementedError`.
- The per-CpG `empirical_fdr_for_dmc` keeps the min-P construction. Region
  mode exists only in the DMR API.

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
| `empirical_pvalue` | float | With `empirical_fdr=True`: permutation p-value (`max_t`) or pooled-null tail fraction (`region`) |
| `empirical_qvalue` | float | With `empirical_fdr=True`: the permutation FDR to threshold |
| `empirical_fdr_set` | float | With `empirical_fdr=True`: constant set-level estimate in `region` mode, NaN in `max_t` mode |

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
| `min_cpgs` | int or None | None | Minimum CpGs in a region. `None` takes the active preset's value, otherwise 5 |
| `pct_sig` | float | 0.5 | Minimum fraction of significant CpGs |
| `minlen_bp` | int | 50 | Minimum region span (bp) |
| `use_q_for_sig` | bool | False | Use q-value instead of p-value for significance |
| `min_mean_qvalue` | float | 0.05 | Post-hoc DMR-level q-value filter |

### Shared Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | MethylData | required | Analysis object |
| `method` | str | `"chain_merge"` | DMR method: `"chain_merge"`, `"tile"`, `"sliding_window"`, `"segment"` |
| `csv` | str | None | Write the DMR table to this TSV path; auto-derived next to the DMR output when unset. Pass `csv=False` to disable. |
| `chromosomes` | list | None | Restrict to specific chromosomes (tile method); an explicit list wins over `canonical_only` |
| `canonical_only` | bool | False | Tile method only: keep the fixed human-style chromosome set of the auto-detected partitions; the other methods raise |
| `backend` | str | `"sequential"` | Execution backend (tile method only) |
| `empirical_fdr` | bool | False | Permutation FDR (tile and chain_merge) |
| `n_perm` | int | 100 | Number of permutations (must be positive) |
| `perm_seed` | int | 42 | Seed for the label shuffles |
| `perm_n_jobs` | int | 1 | joblib workers for the permutations |
| `empirical_strata` | str | None | `md.obs` column that defines shuffle strata; must cover every sample |
| `fdr_method` | str | `"max_t"` | Permutation FDR construction: `"max_t"` or `"region"` |
