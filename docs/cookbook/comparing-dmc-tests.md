# Choosing a DMC Test

epykit ships 8 statistical backends for per-CpG differential methylation
calling. This guide helps you pick the right one for your analysis scenario.

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

### Quick exploration: welch_t or logit_t

```python
ep.tl.dmc(md, test="welch_t")
# or
ep.tl.dmc(md, test="logit_t")
```

These transformation-based tests are fast because they skip the binomial model
entirely. They compare mean beta values (or logit-transformed betas) using a
Welch t-test. Useful for rapid exploration but less powerful than `lr` for
formal analysis.

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

### Stratified design: cmh

```python
ep.tl.dmc(md, test="cmh")
```

The Cochran-Mantel-Haenszel test is designed for stratified 2x2 tables. Use it
when your samples are grouped into strata (e.g., matched pairs, batches) and
you want to test for a common odds ratio across strata.

### Cross-validation with R packages

If you need to reproduce or compare results with R-based tools:

- `ep.tl.dmc(md, test="lr")` produces results comparable to **methylKit**'s
  logistic regression test.
- `ep.tl.dmc(md, test="bb_lr")` produces results comparable to **DSS**'s
  beta-binomial likelihood-ratio test.

## Comparison Table

| Test | Speed | Power | FPR Control | Covariate Support | Min Replicates |
|------|-------|-------|-------------|-------------------|----------------|
| `lr` | Fast | High | Good | No (use `glm`) | 2 |
| `score` | Fast | Slightly higher than `lr` | Mildly anti-conservative | No | 2 |
| `glm` | Moderate | High | Good | Yes (full formula) | 2 |
| `logit_t` | Very fast | Moderate | Good | No | 2 |
| `welch_t` | Very fast | Moderate | Good | No | 2 |
| `bb_lr` | Moderate | High | Good | No | 2 |
| `cmh` | Fast | Moderate | Good | No (stratified only) | 2 per stratum |
| `fisher` | Fast | Low (no bio variance) | Anti-conservative | No | 1 |

**Speed**: relative wall time per chromosome. "Very fast" = no iterative model
fitting. "Moderate" = iterative (IRLS for GLM, numerical optimisation for
bb_lr).

**Power**: ability to detect true positives at a given FDR threshold. Assessed
on simulated WGBS data with known ground truth.

**FPR Control**: how well the test controls false positive rate at the nominal
level. "Good" = well-calibrated. "Mildly anti-conservative" = slightly
inflated at small n. "Anti-conservative" = systematically inflated (fisher at
n=1).

## Recipes by Scenario

| Scenario | Recommended Test | Key Parameters |
|----------|-----------------|----------------|
| Standard two-group WGBS (n >= 3) | `lr` with power stack | `power_stack=True, fdr_method="fdr_tsbh"` |
| Two-group with confounders | `glm` | `formula="~ group + age + sex"` |
| Small sample size (n = 2) | `lr` with power stack | `power_stack=True` (EB dispersion stabilises small-n estimates) |
| Single replicate (n = 1) | `fisher` | `allow_n1=True` |
| Quick exploratory look | `welch_t` | No special parameters |
| Matched-pair or batch-stratified | `cmh` | Strata defined via samplesheet |
| Reproduce methylKit results | `lr` | Default parameters |
| Reproduce DSS results | `bb_lr` | Default parameters |
| Publication with strict FDR | `lr` | `fdr_method="fdr_bh"`, no power stack |
