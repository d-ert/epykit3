# Quality Control

`ep.tl.qc(md)` populates `md.obs` with per-sample QC metrics and caches
detailed QC tables in `md.uns`. The default call is fast and non-destructive;
clinical checks are opt-in.

## Basic Usage

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.unite(md)

# Default QC: global methylation + coverage uniformity
ep.tl.qc(md)

print(md.obs)
```

The default call computes:

- **Global methylation** -- Mean CpG beta across the genome for each sample.
- **Coverage uniformity** -- Per-sample breadth-of-coverage statistics
  (`mean_coverage`, `frac_ge_1x`, `frac_ge_5x`, `frac_ge_10x`,
  `low_coverage_flag`).
- **Bisulfite conversion rate** -- Requires a separate CHH-context store
  (see below).

## Bisulfite Conversion Rate

Conversion rate estimation requires CHH-context methylation data stored in a
separate Parquet store. Pass the path via `chh_context_store`:

```python
ep.tl.qc(md, chh_context_store="path/to/chh_store")
```

!!! note "Reported, not applied"
    The conversion rate is **reported** in `md.obs`, the QC dashboard, and
    the HTML report, but it is **not** used to rescale read counts before
    testing. This matches the default behaviour of methylKit and the bsseq
    family (`read.bismark`, etc.), which leave count-level correction to
    the user. For a well-converted library (>=99.5%) the correction is
    statistically negligible; for a poorly converted one the right action
    is usually to re-prep the library.

## Clinical Checks

Clinical checks are opt-in via `run_*` flags to keep the default call fast.

### Sex Check

```python
ep.tl.qc(md, run_sex_check=True, expected_sex_col="reported_sex")
```

Computes mean CpG beta on chrX for each sample and infers sex from the
bimodal distribution (females have higher chrX methylation due to X
inactivation). When `expected_sex_col` names a column in `md.obs`, the
inferred sex is compared to the reported value and mismatches are flagged.

Results:

- `md.obs` gains columns `inferred_sex` and `sex_mismatch`.
- `md.uns["qc_sex_check"]` holds the full per-sample sex-check DataFrame.

### Contamination Estimation

```python
ep.tl.qc(md, run_contamination=True)
```

Computes a beta-distribution bimodality score per sample. Well-separated
bimodal peaks at 0 and 1 indicate clean WGBS data; a flat or unimodal
distribution may indicate contamination or poor bisulfite conversion.

Results:

- `md.obs` gains a `contamination_score` column.

### Sample Correlation

```python
ep.tl.qc(md, run_sample_correlation=True, correlation_method="spearman")
```

Computes all-vs-all pairwise correlations across samples for sample-swap
detection. Samples from the same group should correlate highly; an anomalously
low minimum pairwise correlation suggests a mislabelled or contaminated sample.

Results:

- `md.obs` gains a `min_pairwise_corr` column.
- `md.uns["qc_sample_correlation"]` holds the full pairwise correlation
  DataFrame.

### All Clinical Checks at Once

```python
ep.tl.qc(
    md,
    chh_context_store="path/to/chh_store",
    run_sex_check=True,
    expected_sex_col="reported_sex",
    run_contamination=True,
    run_sample_correlation=True,
    correlation_method="spearman",
)
```

## Direct-Access QC Functions

Each QC metric is also available as a standalone function under the `ep.qc.*`
namespace for use outside the standard pipeline:

| Function | Description |
|----------|-------------|
| `ep.qc.bisulfite_conversion_rate(store, sample, chh_store)` | Conversion rate from CHH context |
| `ep.qc.global_methylation_report(store, samples)` | Per-sample, per-context global methylation |
| `ep.qc.coverage_uniformity(store, sample)` | Coverage breadth statistics |
| `ep.qc.sex_check(store, samples)` | chrX-based sex inference |
| `ep.qc.contamination_estimate(store, sample)` | Beta bimodality score |
| `ep.qc.sample_correlation_qc(store, samples)` | All-vs-all correlation matrix |

## Power Calculator

`ep.qc.power_calc()` estimates the minimum number of samples per group needed
to detect a given methylation difference at a specified power level, or
computes the achievable power for a fixed sample size.

```python
# How many samples per group to detect a 10% meth diff at 80% power?
n = ep.qc.power_calc(meth_diff=0.10, coverage=20, power=0.80)
print(f"Need {n} samples per group")

# What power do I have with 5 samples per group?
pwr = ep.qc.power_calc(meth_diff=0.10, coverage=20, n_per_group=5)
print(f"Power = {pwr:.2f}")
```

The model accounts for both binomial sampling noise
(`beta * (1 - beta) / coverage`) and between-replicate biological variance.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `meth_diff` | float | required | Expected methylation difference to detect |
| `coverage` | float | required | Mean per-CpG read depth |
| `n_per_group` | int | None | Samples per group (returns power) |
| `power` | float | None | Desired power (returns required n) |
| `alpha` | float | 0.05 | Significance level |
| `baseline_beta` | float | 0.5 | Baseline methylation level |
| `replicate_sd` | float | 0.05 | Between-replicate biological SD |
| `two_sided` | bool | True | Two-sided test |

Provide either `n_per_group` (to compute power) or `power` (to compute the
required sample size), but not both.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | MethylData | required | Analysis object |
| `chh_context_store` | str | None | Path to CHH-context Parquet store |
| `run_sex_check` | bool | False | Enable chrX sex inference |
| `run_contamination` | bool | False | Enable bimodality contamination score |
| `run_sample_correlation` | bool | False | Enable all-vs-all correlation |
| `correlation_method` | str | `"spearman"` | Correlation method (`"spearman"` or `"pearson"`) |
| `expected_sex_col` | str | None | Column in `md.obs` with reported sex for mismatch detection |
