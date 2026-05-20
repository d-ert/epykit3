# Covariate & Multi-group Contrasts

`ep.tl.dmc()` supports covariate-adjusted designs, multi-group joint tests,
and continuous-covariate primary effects via Wilkinson-style formulas. When
`formula=` is specified, epykit auto-selects the GLM backend and fits a full
binomial GLM at every CpG site.

## Covariate-Adjusted Binary Contrast

Adjust for nuisance variables (batch, donor, age) while testing a binary
treatment effect:

```python
ep.tl.dmc(md, formula="~ group + donor", contrast="group")
```

This fits `logit(beta) ~ group + donor` at every CpG and tests the `group`
coefficient via a Wald z-test. The `donor` variable absorbs inter-individual
variance that would otherwise inflate the residual, increasing power.

!!! tip "When to adjust for covariates"
    Covariate adjustment is most impactful when the nuisance variable
    explains substantial variance (e.g., batch effects across sequencing
    runs, matched tumor/normal pairs from the same donor). If all samples
    come from a single batch and independent donors, the unadjusted `lr`
    test is usually sufficient.

## Multi-group Joint F-test

Test whether methylation differs across 3 or more groups simultaneously:

```python
ep.tl.dmc(md, formula="~ group", contrast="group")
```

When `group` is a categorical column with 3+ levels (e.g., `"WT"`, `"KO"`,
`"HET"`), epykit constructs the dummy-coded design matrix and runs a joint
F-test across all group coefficients. This identifies CpGs where at least one
group differs from the others.

## Continuous Covariate as Primary Effect

Test whether methylation changes with a continuous variable (age, dose):

```python
ep.tl.dmc(md, formula="~ age", contrast="age")
```

This fits a linear effect of `age` on logit(beta) and tests the slope
coefficient. Useful for identifying age-associated CpGs across a cohort.

## Auto-Selection of GLM Backend

When `formula=` is specified, the test is forced to the GLM path regardless of
the `test=` parameter. This is because the closed-form `lr` and `score` tests
do not support arbitrary design matrices.

```python
# These are equivalent when formula is provided:
ep.tl.dmc(md, formula="~ group + batch", contrast="group")
ep.tl.dmc(md, formula="~ group + batch", contrast="group", test="glm")
```

## Additional Output Columns

### Binary Contrast with Covariates

When `formula=` is used with a binary contrast, the output includes additional
GLM-specific columns alongside the standard DMC columns:

| Column | Description |
|--------|-------------|
| `coef_treatment` | Estimated coefficient for the treatment term |
| `coef_se` | Standard error of the treatment coefficient |

### Multi-group F-test

When the contrast is a categorical factor with 3+ levels, the output adds:

| Column | Description |
|--------|-------------|
| `f_stat` | Joint F-statistic |
| `df1` | Numerator degrees of freedom |
| `df2` | Denominator degrees of freedom |
| `mean_beta_<level>` | Mean beta for each factor level (e.g., `mean_beta_WT`, `mean_beta_KO`) |
| `meth_diff_max` | Maximum pairwise methylation difference across levels |

## Example: Adjusting for Batch Effects

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")

# Ensure batch info is in md.obs
# (typically added via the samplesheet or manually)
print(md.obs.select(["sample_id", "group", "batch"]))

ep.pp.filter_coverage(md)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)

# Adjust for batch while testing group
ep.tl.dmc(md, formula="~ group + batch", contrast="group")

# Results stored at md.varm["dmc_glm_contrast"]
results = md.varm["dmc_glm_contrast"]
sig = results.filter(results["qvalue"] < 0.05)
print(f"{len(sig)} significant DMCs after batch adjustment")
```

## Example: Three-group Comparison

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="KO",
                     control_group="WT", assembly="hg38")

# md.obs has a 'genotype' column with levels: WT, KO, HET
ep.pp.filter_coverage(md)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)

ep.tl.dmc(md, formula="~ genotype", contrast="genotype")

results = md.varm["dmc_glm_contrast"]
print(results.select(["chrom", "pos", "f_stat", "pvalue", "qvalue",
                       "mean_beta_WT", "mean_beta_KO", "mean_beta_HET"]))
```

## Formula Syntax

Formulas follow the Wilkinson notation as implemented by
[patsy](https://patsy.readthedocs.io/). The left-hand side is always the
methylation response (implicit); only the right-hand side is specified.

| Formula | Meaning |
|---------|---------|
| `"~ group"` | Single factor (binary or multi-group) |
| `"~ group + batch"` | Factor with nuisance covariate |
| `"~ group + donor"` | Paired/matched design |
| `"~ age"` | Continuous primary effect |
| `"~ age + sex"` | Continuous effect adjusted for a categorical covariate |
| `"~ group + age + sex"` | Factor adjusted for mixed covariates |

## Contrast Specification

The `contrast` parameter accepts several forms:

| Form | Example | Effect |
|------|---------|--------|
| Column name (continuous) | `contrast="age"` | Single-coefficient Wald z-test |
| Factor name (categorical) | `contrast="group"` | Joint F-test over all dummies |
| Patsy linear combination | `contrast="group[T.KO] - group[T.WT]"` | Single pairwise contrast |
| Raw numpy matrix | `contrast=np.array(...)` | Custom `(k, p)` contrast matrix |

!!! note "Empirical FDR not supported"
    Empirical FDR (`empirical_fdr=True`) is not compatible with the
    covariate / multi-group path. Label shuffling invalidates the stratified
    design encoded by the formula.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `formula` | str | None | Patsy formula on `md.obs` columns |
| `contrast` | str/array | None | Contrast specification (see table above) |
| `covariates` | list[str] | None | Convenience list of nuisance column names |
| `treatment_col` | str | `"treatment"` | Binary 0/1 column for the legacy path |
