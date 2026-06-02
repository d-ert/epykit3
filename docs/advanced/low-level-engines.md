# Low-level Engines

The `ep.tl.*` wrappers (e.g., `ep.tl.dmc()`, `ep.tl.dmr()`) are the
recommended entry points for most users. They handle parameter validation,
result storage, caching, and MethylData bookkeeping automatically.

The functions documented on this page are the underlying engines that `ep.tl.*`
calls internally. Use them when you need direct control: custom pipelines,
non-standard workflows, benchmarking, or building your own higher-level
abstractions.

!!! warning "1.0 import path"
    Several DMC-level engines (`process_chromosomes_dmc`,
    `apply_multiple_testing_correction`, `empirical_fdr_for_dmc`,
    `fisher_exact_vectorized`, `shrink_meth_diff`) were removed from the
    top-level `epykit.*` namespace in 1.0. They remain importable from
    their owning submodule:

    ```python
    from epykit.dmc import (
        process_chromosomes_dmc,
        apply_multiple_testing_correction,
        empirical_fdr_for_dmc,
        fisher_exact_vectorized,
        shrink_meth_diff,
    )
    ```

    The bare `ep.<name>` access still resolves via a `__getattr__` shim
    that emits `DeprecationWarning`; the shim is scheduled for removal in
    1.2. The samples below use the new explicit import path; older `ep.X`
    snippets will continue to run for the 1.x series.

## DMC Engines

### process_chromosomes_dmc

Per-chromosome streaming DMC accumulator. This is the core function that
iterates over chromosomes, runs the selected statistical test on each, and
writes results to a `DMCStore`.

```python
from epykit.dmc import process_chromosomes_dmc

dmc_store = process_chromosomes_dmc(
    store="methylstore/",
    samples={"treatment": ["t1", "t2"], "control": ["c1", "c2"]},
    test="lr",
    dispersion="eb",
    chromosomes=["chr1", "chr2"],
    output_dir=".cache/dmc/lr/",
)
```

Returns a `DMCStore` with one parquet file per chromosome.

### apply_multiple_testing_correction

BH or Storey correction on a `DMCStore` or a Polars/pandas DataFrame. When
given a `DMCStore`, it streams through chromosomes in two passes (count, then
correct) to avoid loading the full genome into memory.

```python
from epykit.dmc import apply_multiple_testing_correction

# On a DMCStore (streaming, in-place: writes corrected per-chrom parquets)
apply_multiple_testing_correction(dmc_store, method="fdr_bh")

# On a DataFrame (in-memory)
corrected_df = apply_multiple_testing_correction(dmc_df, method="fdr_tsbh")
```

Supported methods: `"fdr_bh"` (Benjamini-Hochberg), `"fdr_tsbh"` (two-stage BH).

### empirical_fdr_for_dmc

Permutation-based empirical FDR. Shuffles treatment/control labels `n_perm`
times, re-runs the test, and estimates an empirical null distribution.

```python
from epykit.dmc import empirical_fdr_for_dmc
empirical_fdr_for_dmc(md, n_perm=100, test="lr")
```

Adds `empirical_pvalue` and `empirical_qvalue` columns to the DMC result.

### fisher_exact_vectorized

Vectorized Fisher exact test for 2x2 contingency tables. Operates on arrays
of (a, b, c, d) values and returns p-values and odds ratios.

```python
import numpy as np
from epykit.dmc import fisher_exact_vectorized

a = np.array([10, 20, 5])
b = np.array([5, 10, 15])
c = np.array([3, 8, 12])
d = np.array([12, 2, 8])

pvalues, odds_ratios = fisher_exact_vectorized(a, b, c, d)
```

Used internally by the `fisher` test backend and the separation fallback in
lr+.

## DMR Engines

### ep.call_dmr_chain_merge

DSS-style chain-merge DMR caller. Takes a DMC DataFrame (or iterates over a
DMCStore) and chains consecutive significant CpGs into regions.

```python
dmr_df = ep.call_dmr_chain_merge(
    dmc_df,
    alpha=1e-4,
    min_abs_meth_diff=0.10,
    dis_merge_bp=500,
    min_cpgs=3,
    pct_sig=0.5,
    minlen_bp=50,
)
```

Returns a Polars DataFrame of DMR regions. See [DMR Calling](../analysis/dmr.md)
for the full parameter reference.

### ep.call_dmr_sliding_window

Sliding-window DMR engine. Combines per-CpG p-values within overlapping windows
using signed Stouffer's Z.

```python
dmr_df = ep.call_dmr_sliding_window(
    dmc_store,
    window_bp=500,
    step_bp=250,
)
```

### ep.DMR_PRESETS

Dictionary of parameter bundles for `call_dmr_chain_merge`. Each preset is a
dict of keyword arguments.

```python
print(ep.DMR_PRESETS.keys())
# dict_keys(['strict', 'default', 'permissive'])

# Use a preset directly
dmr_df = ep.call_dmr_chain_merge(dmc_df, **ep.DMR_PRESETS["strict"])
```

## Smoothing Engines

### ep.smooth_methylation_gaussian

Gaussian-kernel smoothing of methylation values in the methylstore.

```python
ep.smooth_methylation_gaussian(store, bandwidth_bp=200)
```

### ep.smooth_methylation_bsmooth

BSmooth-style local-likelihood smoothing. Fits a local weighted regression to
methylation fractions, accounting for read-count uncertainty.

```python
ep.smooth_methylation_bsmooth(store, bandwidth_bp=1000, min_coverage=1)
```

Both smoothing engines write results to the `.cache/` hierarchy and operate
one chromosome at a time.

## DVC / DVR Engines

### ep.process_chromosomes_dvc

Streaming engine for differential variability calling (DVC). Analogous to
`process_chromosomes_dmc` but tests for differences in methylation variance
rather than mean.

```python
dvc_store = ep.process_chromosomes_dvc(
    store="methylstore/",
    samples={"treatment": ["t1", "t2"], "control": ["c1", "c2"]},
    chromosomes=["chr1", "chr2"],
)
```

### ep.call_dvr_density

DVR density aggregation. Groups significant DVCs into differentially variable
regions using a kernel-density approach.

```python
dvr_df = ep.call_dvr_density(dvc_df, bandwidth_bp=1000, min_dvcs=3)
```

## Utility Engines

### ep.build_design

Build a design matrix from a Wilkinson-style formula string and an observation
DataFrame.

```python
import polars as pl

obs = pl.DataFrame({
    "sample_id": ["t1", "t2", "c1", "c2"],
    "group": ["tumor", "tumor", "normal", "normal"],
    "age": [55, 62, 48, 71],
})

design, contrast_vector = ep.build_design("~ group + age", obs)
```

Used internally by the `glm` test backend when `formula=` is specified in
`ep.tl.dmc()`.

### ep.DMCStore

The persistent per-chromosome DMC result store. See
[Methylstore Internals](methylstore.md) for the full description of its layout
and design rationale.

```python
dmc_store = ep.DMCStore(".cache/dmc/lr/")
df = dmc_store.to_dataframe()            # concatenate all chromosomes
for chrom_df in dmc_store.iter_chromosomes():  # stream one at a time
    print(chrom_df.shape)
```

## When to Use These Engines

Use the low-level engines when:

- You are building a **custom pipeline** that does not follow the standard
  preprocess-test-annotate workflow.
- You need to **benchmark** individual components (e.g., compare the runtime of
  `lr` vs `welch_t` on a single chromosome).
- You want to run a test on **data not managed by MethylData** (e.g., a
  DataFrame you constructed manually).
- You are implementing a **new DMR method** that consumes DMCStore output.
- You need to apply **multiple testing correction separately** from the DMC
  calling step.

For standard analyses, prefer `ep.tl.dmc()`, `ep.tl.dmr()`, and the other
`ep.tl.*` wrappers. They handle caching, result storage, parameter validation,
and MethylData state management automatically.
