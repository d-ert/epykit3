# DMC Calling

`ep.tl.dmc(md)` runs per-CpG differential methylation calling and stores the
result in `md.varm["dmc_<test>"]`. It supports 8 statistical backends, automatic
test selection, covariate adjustment, empirical FDR, and distributed execution.

## Basic Usage

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)

# Default: auto-selects lr at n>=2, fisher at n=1
ep.tl.dmc(md)

# Access results
dmc_results = md.dmc  # shorthand for md.varm["dmc_lr"]
print(f"Tested {len(dmc_results)} CpG sites")
```

## Test Backends

epykit ships 8 statistical backends for DMC calling. The `test` parameter
selects which one to use.

| Test | Description | When to use |
|------|-------------|-------------|
| `lr` (default) | Quasi-binomial LR, McCullagh-Nelder dispersion | n >= 2, general purpose |
| `score` | Pearson score test | Slightly more powerful, mildly anti-conservative |
| `glm` | Full IRLS binomial GLM | Covariate adjustment (auto-selected with `formula=`) |
| `logit_t` | Welch t on logit(beta) | Transformation-based alternative |
| `welch_t` | Welch t on raw beta | Simple mean comparison |
| `bb_lr` | True quasi-binomial LRT | Alternative to `lr` |
| `cmh` | Cochran-Mantel-Haenszel | Stratified designs |
| `fisher` | Pooled Fisher exact | n=1 fallback (auto-selected) |

### Auto-Selection Logic

When `test="auto"` (the default):

- **n >= 2 per group** -- selects `lr` (quasi-binomial likelihood-ratio).
- **n = 1 per group** -- selects `fisher` (pooled Fisher exact). Requires
  `allow_n1=True` or raises an error.
- **`formula=` specified** -- selects `glm` regardless of sample size.

```python
# Explicit test selection
ep.tl.dmc(md, test="score")

# n=1 fallback (requires opt-in)
ep.tl.dmc(md, test="auto", allow_n1=True)
```

!!! warning "Fisher exact at n=1"
    The `fisher` test pools reads across replicates and ignores
    between-replicate variance. P-values are anti-conservative and should
    not be reported as evidence of differential methylation. Use only as a
    last-resort exploratory tool.

## Dispersion Strategies

The `dispersion` parameter controls how the overdispersion parameter (phi) is
estimated for the `lr` and `score` tests.

| Strategy | Description |
|----------|-------------|
| `"eb"` (default) | Empirical-Bayes shrinkage toward a chromosome-wide inverse-Gamma prior. Better at low n, identical to `"site"` at high n. |
| `"site"` | Per-site phi from the 4-df Pearson residual sum |
| `"chrom"` | Single phi shared across all sites on a chromosome |
| `"shrink"` | James-Stein shrinkage of site-level phi toward chromosome mean |

```python
# Use per-site dispersion (no shrinkage)
ep.tl.dmc(md, dispersion="site")
```

## Empirical FDR

Permutation-based FDR adds `empirical_pvalue` and `empirical_qvalue` columns
by shuffling treatment/control labels and re-running the test.

```python
ep.tl.dmc(md, empirical_fdr=True, n_perm=100)
```

!!! note
    Empirical FDR is not supported with the covariate / multi-group path
    (`formula=`). Label shuffling invalidates the stratified design.

## Output Columns

The result DataFrame (`md.varm["dmc_<test>"]`) contains:

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | str | Chromosome |
| `pos` | int | Genomic position (0-based) |
| `strand` | str | Strand (`+` or `-`) |
| `n_case` | int | Treatment samples with coverage at this site |
| `n_control` | int | Control samples with coverage at this site |
| `mean_beta_case` | float | Mean beta in treatment group |
| `mean_beta_control` | float | Mean beta in control group |
| `meth_diff` | float | Methylation difference (case - control) |
| `meth_diff_ci_lo` | float | Lower 95% CI bound on meth_diff |
| `meth_diff_ci_hi` | float | Upper 95% CI bound on meth_diff |
| `pvalue` | float | Raw p-value |
| `qvalue` | float | BH-adjusted q-value |
| `log2_odds_ratio` | float | Log2 odds ratio |

Results are stored at `md.varm["dmc_<test>"]`, where `<test>` is the canonical
test name (e.g., `dmc_lr`, `dmc_score`, `dmc_glm`).

## Smoothed-Input Mode

When `use_smoothed=True`, the DMC test runs on pseudo-counts derived from
prior smoothing (`ep.pp.smooth(md)`):

```python
ep.pp.smooth(md, method="bsmooth")
ep.tl.dmc(md, use_smoothed=True)
# Results stored at md.varm["dmc_lr_smoothed"]
```

!!! warning "Not equivalent to DSS smoothing"
    The pseudo-count approach replaces the count signal entirely with the
    locally-averaged version. This is more aggressive than DSS's
    `DMLfit.multiFactor(smoothing=TRUE)`, which only smooths the dispersion
    step. For DSS-style behaviour, prefer `use_smoothed=False` (the default)
    with `test="lr"`.

## Resumable Computation

For long-running analyses, `resumable=True` enables checkpoint/resume. If the
same inputs and parameters have been run before, the cached result is loaded
from disk.

```python
ep.tl.dmc(md, resumable=True)
```

The fingerprint includes: methylstore path, sample lists, test, chromosomes,
dispersion, reference, and FDR parameters. Changing any of these invalidates
the cache.

## Key Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | MethylData | required | Analysis object |
| `test` | str | `"auto"` | Statistical test (see table above) |
| `dispersion` | str | `"eb"` | Dispersion estimation strategy |
| `chromosomes` | list | None | Restrict to specific chromosomes |
| `min_samples_treatment` | int | 0 | Minimum treatment samples with coverage per site |
| `min_samples_control` | int | 0 | Minimum control samples with coverage per site |
| `allow_n1` | bool | False | Allow n=1 per group (Fisher fallback) |
| `empirical_fdr` | bool | False | Run permutation-based empirical FDR |
| `n_perm` | int | 100 | Number of permutations for empirical FDR |
| `backend` | str | `"sequential"` | Execution backend (`"sequential"`, `"dask"`, `"ray"`) |
| `n_workers` | int | None | Worker pool size (None = backend default) |
| `glm_backend` | str | `"cpu"` | GLM execution target (`"cpu"`, `"gpu"`) |
| `resumable` | bool | False | Enable checkpoint/resume |
| `use_smoothed` | bool | False | Use smoothed pseudo-counts |
| `fdr_method` | str | `"fdr_bh"` | FDR correction method |
| `power_stack` | bool | False | Enable all lr+ enhancements at once |

## Distributed Execution

The `backend` parameter controls parallelism:

- `"sequential"` (default) -- One chromosome at a time on the main process.
- `"dask"` -- One task per chromosome via Dask. Requires
  `pip install 'epykit[distributed]'`.
- `"ray"` -- One task per chromosome via Ray. Requires
  `pip install 'epykit[ray]'`.

```python
ep.tl.dmc(md, backend="dask", n_workers=8)
```

The `glm_backend` parameter (for `test="glm"` only) routes the IRLS hot path:

- `"cpu"` (default) -- numpy.
- `"gpu"` -- CuPy. Requires `pip install 'epykit[gpu]'`.

## Covariate Adjustment and Multi-group Contrasts

See the dedicated [Covariate & Multi-group](covariates.md) page for
`formula=` and `contrast=` usage.

## lr+ Power Stack

See the dedicated [lr+ Power Stack](lr-plus.md) page for the four opt-in
enhancements (`fdr_method`, `neighbour_combine`, `sep_fallback`, `dispersion="eb"`).
