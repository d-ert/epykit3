# Choosing a DMC Test

epykit ships 4 statistical backends for per-CpG differential methylation
calling (plus an `auto` dispatcher). This guide helps you pick the right one
for your analysis scenario.

## Decision Guide

### Default: auto-selection

```python
ep.tl.dmc(md, test="auto")
```

The `test="auto"` default inspects your data and selects the appropriate test:

- **n >= 2 per group** -- selects `lr` (quasi-binomial likelihood-ratio).
- **n = 1 per group** -- selects `fisher` (pooled Fisher exact). Requires
  `allow_n1=True`.
- **`formula=` specified** -- selects `glm` regardless of sample size.

For most users, `test="auto"` is the right choice.

### Best power: lr with power stack

```python
ep.tl.dmc(
    md,
    fdr_method="fdr_tsbh",
    neighbour_combine=True,
    sep_fallback=True,
    dispersion="eb",
)
# Or equivalently:
ep.tl.dmc(md, power_stack=True, fdr_method="fdr_tsbh")
```

The lr+ power stack enables four enhancements that close the sensitivity gap
to methylKit, RADMeth, and DSS. This is the recommended configuration when
maximising detection power is the priority. See
[lr+ Power Stack](../analysis/lr-plus.md) for details.

### Covariates needed: GLM

```python
ep.tl.dmc(md, formula="~ group + age + sex", contrast="group")
```

When `formula=` is specified, epykit auto-selects the `glm` test. This fits a
full IRLS binomial GLM per CpG, adjusting for all covariates in the formula.
Use this when confounders (age, sex, batch) must be accounted for.

### Conservative baseline: lr with defaults

```python
ep.tl.dmc(md, test="lr")
```

The `lr` test with default settings (no power stack, standard BH) provides a
conservative baseline with well-controlled FPR. Use this when strict FDR
control matters more than sensitivity.

### Quick exploration: welch_t

```python
ep.tl.dmc(md, test="welch_t")
```

This transformation-based test is fast because it skips the binomial model
entirely. It compares mean beta values using a Welch t-test. Useful for rapid
exploration but less powerful than `lr` for formal analysis.

### Single replicate: fisher

```python
ep.tl.dmc(md, test="auto", allow_n1=True)
# auto-selects fisher at n=1
```

When you have only one sample per group, the `fisher` test pools reads across
the single replicate and runs a Fisher exact test per CpG. This is the only
option at n=1, but it cannot account for between-replicate variance and
produces anti-conservative p-values.

!!! warning "n=1 limitations"
    Fisher exact at n=1 ignores biological variance. P-values should not be
    interpreted as evidence of differential methylation. Use only as an
    exploratory tool.

### Cross-validation with R packages

If you need to reproduce or compare results with R-based tools:

- `ep.tl.dmc(md, test="lr")` produces results comparable to **methylKit**'s
  logistic regression test.

!!! note "Engines removed in 0.7.5"
    `cmh` (use `formula='~ group + batch'`), `bb_lr` (use `lr`),
    `score` (use `lr`), `logit_t` (use `welch_t`). All raise `ValueError`
    with a one-line migration hint.

## Comparison Table

| Test | Speed | Power | FPR Control | Covariate Support | Min Replicates |
|------|-------|-------|-------------|-------------------|----------------|
| `lr` | Fast | High | Good | No (use `glm`) | 2 |
| `glm` | Moderate | High | Good | Yes (full formula) | 2 |
| `welch_t` | Very fast | Moderate | Good | No | 2 |
| `fisher` | Fast | Low (no bio variance) | Anti-conservative | No | 1 |

**Speed**: relative wall time per chromosome. "Very fast" = no iterative model
fitting. "Moderate" = iterative (IRLS for GLM).

**Power**: ability to detect true positives at a given FDR threshold. Assessed
on simulated WGBS data with known ground truth.

**FPR Control**: how well the test controls false positive rate at the nominal
level. "Good" = well-calibrated. "Anti-conservative" = systematically inflated
(fisher at n=1).

## Recipes by Scenario

| Scenario | Recommended Test | Key Parameters |
|----------|-----------------|----------------|
| Standard two-group WGBS (n >= 3) | `lr` with power stack | `power_stack=True, fdr_method="fdr_tsbh"` |
| Two-group with confounders | `glm` | `formula="~ group + age + sex"` |
| Small sample size (n = 2) | `lr` with power stack | `power_stack=True` (EB dispersion stabilises small-n estimates) |
| Single replicate (n = 1) | `fisher` | `allow_n1=True` |
| Quick exploratory look | `welch_t` | No special parameters |
| Batch-stratified / covariate design | `glm` | `formula="~ group + batch"` |
| Reproduce methylKit results | `lr` | Default parameters |
| Publication with strict FDR | `lr` | `fdr_method="fdr_bh"`, no power stack |
