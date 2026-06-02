# lr+ Power Stack

The lr+ power stack is a set of four opt-in enhancements to the default `lr`
test that close the sensitivity gap to methylKit, RADMeth, and DSS without
changing the underlying statistical model. Each enhancement can be enabled
individually, or all four at once via `power_stack="lr+"`.

`lr+` is opt-in. Out of the box -- `tl.dmc(md, test="lr")` with no other
flags -- you get bare `lr` with `dispersion="eb"` and `fdr_method="fdr_bh"`,
no neighbour combining, no separation fallback. The `power_stack` shorthand
is what flips the four knobs together.

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
ep.tl.dmc(md, neighbour_combine=True, neighbour_bp=500)
```

Applies a signed Stouffer's Z combiner over nearby CpGs within `neighbour_bp`
base pairs. Adjacent CpGs in WGBS are spatially correlated -- if a CpG is
truly differentially methylated, its neighbours usually are too. Combining
their evidence increases power.

The combiner gates on sign agreement: at least 60% of CpGs in the
neighbourhood must have the same direction of effect. This prevents combining
contradictory signals.

**Columns added when enabled (the raw `pvalue` / `qvalue` columns stay
unchanged):**

| Column | Meaning |
|--------|---------|
| `pvalue` | Raw per-CpG p-value (unchanged from bare `lr`) |
| `qvalue` | BH-corrected raw p-value (unchanged from bare `lr`) |
| `pvalue_combined` | Stouffer-combined p-value over the neighbourhood |
| `qvalue_combined` | FDR-adjusted combined p-value |
| `pvalue_combined_n_neighbours` | How many neighbours fed into each combined call |
| `qvalue_combined_reject` | Boolean: did the combined q-value pass FDR? |

!!! warning "Downstream code must opt in to combined columns"
    `pvalue` and `qvalue` always carry the **per-CpG raw** values, even with
    `neighbour_combine=True`. To filter on the combined evidence, read
    `pvalue_combined` / `qvalue_combined` explicitly. This was a deliberate
    1.0 contract change so older downstream scripts keep working unchanged
    when `neighbour_combine=True` is enabled.

### 3. Separation-Aware Fisher Fallback

```python
ep.tl.dmc(md, sep_fallback=True, sep_threshold=0.9)
```

At CpGs where one group has beta near 0 or 1 (quasi-complete separation),
the quasi-binomial LR test can produce unstable p-values because the
dispersion estimate degenerates. The `sep_fallback` option detects these
cases (max group beta > `sep_threshold` or min group beta < 1 -
`sep_threshold`) and falls back to a Fisher exact test on pooled counts,
which handles boundary tables gracefully.

### 4. Empirical-Bayes Dispersion

```python
ep.tl.dmc(md, dispersion="eb")
```

Shrinks per-site dispersion estimates toward a chromosome-wide inverse-Gamma
prior. At low sample sizes (n = 2--4), per-site dispersion estimates are
noisy and can be wildly over- or under-estimated. Empirical-Bayes shrinkage
stabilises these estimates, producing better-calibrated p-values. At high n
(>= 10), per-site estimates are already precise and the shrinkage has
negligible effect.

This is the default dispersion strategy since version 0.7.1.

## Drop-In Recipe

Enable all four enhancements at once:

```python
ep.tl.dmc(
    md,
    fdr_method="fdr_tsbh",
    neighbour_combine=True,
    neighbour_bp=500,
    sep_fallback=True,
    sep_threshold=0.9,
    dispersion="eb",
)
```

Or use the convenience shorthand:

```python
ep.tl.dmc(md, power_stack="lr+")
```

## `power_stack` modes

| Value | Behaviour |
|-------|-----------|
| `"off"` (default) / `False` | Leave every knob at whatever you passed (or the bare-engine defaults). Bare `lr`. |
| `"lr+"` (alias `"auto"` and `True`) | Engage all four knobs (`neighbour_combine=True`, `fdr_method="fdr_tsbh"`, `sep_fallback=True`, `dispersion="eb"`) at **any** sample size. |
| `"conservative"` | Engage the four knobs only when `min(n_treatment, n_control) <= 2` (the pre-1.0 behaviour, useful if you want a low-n boost without changing high-n calls). |

Unknown strings raise `ValueError`. Passing `power_stack="lr+"` together
with an explicit knob (e.g. `power_stack="lr+", fdr_method="fdr_bh"`) lets
the shorthand turn on the other three knobs but leaves your explicit choice
in place.

!!! note "What changed in 1.0"
    Before 1.0, `power_stack=True` only flipped `neighbour_combine` and
    `sep_fallback`, and `power_stack="auto"` only engaged anything at
    `min_n <= 2`. The 1.0 contract is: `"lr+"`/`"auto"`/`True` engages all
    four knobs at any n. The old behaviour is reachable via
    `power_stack="conservative"`. The **default did not change**: bare
    `tl.dmc(test="lr")` still produces bare-engine output.

## Recipe Matrix

Choose a recipe based on your analysis goals:

| Recipe | fdr_method | neighbour_combine | sep_fallback | dispersion | Use case |
|--------|------------|-------------------|--------------|------------|----------|
| Quick run | `"fdr_bh"` | False | False | `"eb"` | Fast exploratory analysis (bare `lr`) |
| Best power | `"fdr_tsbh"` | True | True | `"eb"` | Maximise sensitivity (equivalent to `power_stack="lr+"`) |
| Strict FDR | `"fdr_bh"` | False | False | `"eb"` | Conservative, publication-ready |
| Low-coverage | `"fdr_tsbh"` | True | True | `"eb"` | n <= 3 or low read depth |
| Context-dependent dispersion | `"fdr_bh"` | False | False | `"chrom"` | When per-site dispersion is suspect |
| Reproduce baseline | `"fdr_bh"` | False | False | `"site"` | Match pre-0.7.1 behaviour |

## Parameters Reference

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fdr_method` | str | `"fdr_bh"` | FDR correction: `"fdr_bh"` (Benjamini-Hochberg) or `"fdr_tsbh"` (two-stage BH) |
| `neighbour_combine` | bool | False | Enable signed Stouffer Z combining over nearby CpGs (adds `pvalue_combined`/`qvalue_combined` columns) |
| `neighbour_bp` | int | 500 | Maximum distance (bp) for neighbour combining |
| `sep_fallback` | bool | False | Fisher fallback for near-perfect-separation tables |
| `sep_threshold` | float | 0.9 | Beta threshold triggering the separation fallback |
| `dispersion` | str | `"eb"` | Dispersion strategy: `"eb"`, `"site"`, `"chrom"`, `"shrink"` |
| `power_stack` | str / bool | `"off"` | One of `"off"`, `"lr+"`, `"auto"` (alias for `"lr+"`), `"conservative"`, or `True`/`False` |
