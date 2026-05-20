# Effect-size Shrinkage

!!! warning "Experimental"
    This feature is implemented but has not been validated on real biological data.
    The API may change in future releases.

`ep.shrink_meth_diff(dmc_df)` applies apeglm-style empirical-Bayes shrinkage
to the `meth_diff` column of a DMC result DataFrame. This stabilizes effect
estimates at low-coverage sites and reduces noise in downstream
visualizations.

## Basic Usage

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)
ep.tl.dmc(md)

# Shrink the meth_diff values
dmc = md.dmc
dmc_shrunk = ep.shrink_meth_diff(dmc)

# Compare original vs. shrunk volcano plots
ep.pl.volcano(md, save="volcano_original.png")
md.dmc = dmc_shrunk
ep.pl.volcano(md, save="volcano_shrunk.png")
```

## How It Works

1. An empirical-Bayes prior is estimated from the observed distribution
   of `meth_diff` values and their associated standard errors.
2. Each site's effect estimate is shrunk toward zero, with the degree of
   shrinkage proportional to its uncertainty. High-coverage sites with
   precise estimates are barely affected; low-coverage sites with noisy
   estimates are pulled strongly toward zero.
3. The shrunk `meth_diff` replaces the original column. A new column
   `meth_diff_unshrunk` preserves the original values.

This is analogous to the `apeglm` shrinkage estimator used in DESeq2 for
RNA-seq log-fold-change estimation, adapted for methylation differences.

## Input and Output

**Input**: A Polars DataFrame with DMC results, as produced by
`ep.tl.dmc()`. Required columns: `meth_diff`, `pvalue`, and either
`meth_diff_se` (standard error) or `meth_diff_ci_lo` / `meth_diff_ci_hi`
(confidence interval bounds, from which the SE is derived).

**Output**: A Polars DataFrame with the same schema, plus:

| Column | Type | Description |
|--------|------|-------------|
| `meth_diff` | float | Shrunk methylation difference |
| `meth_diff_unshrunk` | float | Original methylation difference (preserved) |

All other columns are passed through unchanged.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dmc_df` | DataFrame | required | DMC result DataFrame from `ep.tl.dmc()` |
| `prior_scale` | float | None | Scale of the Cauchy prior. If None, estimated from the data. |

## When to Use Shrinkage

Shrinkage is most useful when:

- **Low sample sizes** (n = 2--3 per group): effect estimates are noisy,
  and many sites have inflated `meth_diff` purely due to sampling
  variability.
- **Volcano or MA plots**: unshrunk plots show a characteristic "fan"
  shape at low coverage, with extreme meth_diff values that are not
  biologically meaningful. Shrinkage compresses the fan into a more
  interpretable shape.
- **Ranking sites by effect size**: shrunk estimates provide a more
  reliable ranking than raw `meth_diff`, because noisy sites are
  penalized rather than promoted.

## Visualizing the Effect

The primary benefit is visible in volcano and MA plots:

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Original
md.dmc = dmc
ep.pl.volcano(md, ax=axes[0], title="Original")

# Shrunk
md.dmc = dmc_shrunk
ep.pl.volcano(md, ax=axes[1], title="Shrunk")

plt.tight_layout()
plt.savefig("shrinkage_comparison.png", dpi=150)
```

Before shrinkage, low-coverage sites produce extreme effect estimates that
dominate the wings of the volcano plot. After shrinkage, only sites with
strong evidence retain large effect sizes.

## Notes

- Shrinkage does not change p-values or q-values. It only modifies the
  effect-size estimate. Sites that were significant before shrinkage
  remain significant afterward.
- The function operates on a standalone DataFrame and does not modify
  the MethylData object in place. Assign the result back to `md.dmc` if
  you want subsequent plotting functions to use the shrunk values.
- Shrinkage is a post-processing step. Run it after `ep.tl.dmc()`, not
  before.
- If `meth_diff_se` is not present in the DataFrame, the standard error
  is approximated from the confidence interval bounds as
  `(ci_hi - ci_lo) / (2 * 1.96)`.
