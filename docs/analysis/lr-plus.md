# lr+ Power Stack

The lr+ power stack is a set of four opt-in enhancements to the default `lr`
test that close the sensitivity gap to methylKit, RADMeth, and DSS without
changing the underlying statistical model. Each enhancement can be enabled
individually or all at once via `power_stack=True`.

## The Four Enhancements

### 1. Two-Stage BH / Storey q-values

```python
ep.tl.dmc(md, fdr_method="fdr_tsbh")
```

Replaces the standard Benjamini-Hochberg correction with the two-stage BH
procedure, which estimates the proportion of true nulls (pi0) and uses it to
recalibrate the q-value threshold. This is more powerful than BH when the
number of true positives is large (common in WGBS with hundreds of thousands
of differentially methylated CpGs).

### 2. Neighbour-Aware P-value Combining

```python
ep.tl.dmc(md, neighbour_combine=True, neighbour_bp=200)
```

Applies a signed Stouffer's Z combiner over nearby CpGs within `neighbour_bp`
base pairs. Adjacent CpGs in WGBS are spatially correlated -- if a CpG is
truly differentially methylated, its neighbours usually are too. Combining
their evidence increases power.

The combiner gates on sign agreement: at least 60% of CpGs in the
neighbourhood must have the same direction of effect. This prevents combining
contradictory signals.

**Column changes when enabled:**

| Column | Meaning |
|--------|---------|
| `pvalue` | Combined p-value (Stouffer Z) |
| `qvalue` | BH-corrected combined p-value |
| `pvalue_raw` | Original per-CpG p-value |
| `qvalue_raw` | BH-corrected original p-value |

!!! warning "Downstream code expecting raw p-values"
    When `neighbour_combine=True`, the `pvalue` and `qvalue` columns contain
    the **combined** values, not the per-CpG originals. If your downstream
    code (e.g., custom filtering scripts) expects raw per-CpG p-values, use
    `pvalue_raw` and `qvalue_raw` instead.

### 3. Separation-Aware Fisher Fallback

```python
ep.tl.dmc(md, sep_fallback=True, sep_threshold=0.9)
```

At CpGs where one group has beta near 0 or 1 (quasi-complete separation), the
quasi-binomial LR test can produce unstable p-values because the dispersion
estimate degenerates. The `sep_fallback` option detects these cases
(max group beta > `sep_threshold` or min group beta < 1 - `sep_threshold`) and
falls back to a Fisher exact test on pooled counts, which handles boundary
tables gracefully.

### 4. Empirical-Bayes Dispersion

```python
ep.tl.dmc(md, dispersion="eb")
```

Shrinks per-site dispersion estimates toward a chromosome-wide inverse-Gamma
prior. At low sample sizes (n = 2--4), per-site dispersion estimates are noisy
and can be wildly over- or under-estimated. Empirical-Bayes shrinkage
stabilises these estimates, producing better-calibrated p-values. At high n
(>= 10), per-site estimates are already precise and the shrinkage has negligible
effect.

This is the default dispersion strategy since version 0.7.1.

## Drop-In Recipe

Enable all four enhancements at once:

```python
ep.tl.dmc(
    md,
    fdr_method="fdr_tsbh",
    neighbour_combine=True,
    neighbour_bp=200,
    sep_fallback=True,
    sep_threshold=0.9,
    dispersion="eb",
)
```

Or use the convenience shorthand:

```python
ep.tl.dmc(md, power_stack=True)
```

`power_stack=True` sets `neighbour_combine=True` and `sep_fallback=True`. The
default `dispersion="eb"` and default `fdr_method="fdr_bh"` remain in effect.
To also switch to two-stage BH, combine:

```python
ep.tl.dmc(md, power_stack=True, fdr_method="fdr_tsbh")
```

!!! tip "Auto-engagement at low n"
    `power_stack="auto"` automatically enables `neighbour_combine` and
    `sep_fallback` when the minimum group size is <= 2. Pass
    `power_stack=False` to prevent this.

## Recipe Matrix

Choose a recipe based on your analysis goals:

| Recipe | fdr_method | neighbour_combine | sep_fallback | dispersion | Use case |
|--------|------------|-------------------|--------------|------------|----------|
| Quick run | `"fdr_bh"` | False | False | `"eb"` | Fast exploratory analysis |
| Best power | `"fdr_tsbh"` | True | True | `"eb"` | Maximise sensitivity |
| Strict FDR | `"fdr_bh"` | False | False | `"eb"` | Conservative, publication-ready |
| Low-coverage | `"fdr_tsbh"` | True | True | `"eb"` | n <= 3 or low read depth |
| Context-dependent dispersion | `"fdr_bh"` | False | False | `"chrom"` | When per-site dispersion is suspect |
| Reproduce baseline | `"fdr_bh"` | False | False | `"site"` | Match pre-0.7.1 behaviour |

## Parameters Reference

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fdr_method` | str | `"fdr_bh"` | FDR correction: `"fdr_bh"` (Benjamini-Hochberg) or `"fdr_tsbh"` (two-stage BH) |
| `neighbour_combine` | bool | False | Enable signed Stouffer Z combining over nearby CpGs |
| `neighbour_bp` | int | 500 | Maximum distance (bp) for neighbour combining |
| `sep_fallback` | bool | False | Fisher fallback for near-perfect-separation tables |
| `sep_threshold` | float | 0.9 | Beta threshold triggering the separation fallback |
| `dispersion` | str | `"eb"` | Dispersion strategy: `"eb"`, `"site"`, `"chrom"`, `"shrink"` |
| `power_stack` | bool/str | False | Convenience: `True` enables neighbour_combine + sep_fallback; `"auto"` enables them only at low n |
