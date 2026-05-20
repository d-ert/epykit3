# Smoothing

`ep.pp.smooth()` applies spatial smoothing to methylation values along the genome. This
can reveal regional trends in methylation that are obscured by site-level noise, and is
primarily useful for visualization and exploratory analysis.

## Function Signature

```python
ep.pp.smooth(md, method="gaussian", bandwidth=1000)
```

## Parameters

### Common parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | `MethylData` | required | The MethylData object to smooth |
| `method` | `str` | `"gaussian"` | Smoothing method: `"gaussian"` or `"bsmooth"` |

### Gaussian-specific parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `bandwidth` | `int` | `1000` | Kernel bandwidth in base pairs. Larger values produce smoother curves. |

### BSmooth-specific parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ns` | `int` | `70` | Minimum number of CpGs in each smoothing window |
| `h` | `int` | `1000` | Minimum smoothing window half-width in base pairs |

## Methods

### Gaussian kernel smoothing (default)

```python
ep.pp.smooth(md, method="gaussian", bandwidth=1000)
```

Applies a Gaussian kernel to methylation fractions along each chromosome. Each site's
smoothed value is a weighted average of nearby sites, where the weights decay with genomic
distance according to the Gaussian kernel. This method is fast and works well for
visualizing broad methylation patterns.

**Key characteristics:**

- Fixed bandwidth in base pairs.
- Computationally efficient -- processes each chromosome in a single pass.
- Best suited for visualization and exploratory analysis.

### BSmooth (B-spline smoothing)

```python
ep.pp.smooth(md, method="bsmooth", ns=70, h=1000)
```

Implements the BSmooth algorithm (Hansen et al., 2012), which fits local-likelihood
B-splines to methylation data. The smoothing window adapts to ensure that at least `ns`
CpGs are included and that the window spans at least `h` base pairs in each direction.

**Key characteristics:**

- Adaptive window size based on local CpG density.
- More statistically principled than kernel smoothing.
- Slower than Gaussian smoothing, especially on whole-genome data.

## Usage

### Gaussian smoothing for visualization

```python
import epykit as ep

md = ep.read_bismark("samplesheet.csv", treatment_group="tumor", control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)

# Smooth for plotting
ep.pp.smooth(md, method="gaussian", bandwidth=1000)
```

### BSmooth for detailed regional analysis

```python
ep.pp.smooth(md, method="bsmooth", ns=70, h=1000)
```

### Adjusting the bandwidth

```python
# Narrow bandwidth: more detail, more noise
ep.pp.smooth(md, method="gaussian", bandwidth=500)

# Wide bandwidth: smoother curves, less local detail
ep.pp.smooth(md, method="gaussian", bandwidth=3000)
```

## Important Note

This smoothing step is for **visual and exploratory purposes only**. The DMC and DMR
engines (`ep.tl.dmc()`, `ep.tl.dmr()`) apply their own statistical smoothing internally
as part of their testing procedures. You do not need to run `ep.pp.smooth()` before
differential analysis, and doing so does not affect DMC/DMR results.

## Choosing a Method

| Scenario | Recommended Method |
|----------|-------------------|
| Quick visualization of methylation landscapes | `"gaussian"` |
| Publication-quality smoothed profiles | `"bsmooth"` |
| Large datasets where speed matters | `"gaussian"` |
| Faithful reproduction of Hansen 2012 methodology | `"bsmooth"` |

## Call Order

Smoothing is optional and should be applied after the core preprocessing steps:

```python
ep.pp.filter_coverage(md)       # Required
ep.pp.normalize_coverage(md)    # Required
ep.pp.unite(md)                 # Required
ep.pp.smooth(md)                # Optional
```

## Next Steps

See [Aggregate Regions](aggregate-regions.md) for collapsing per-CpG data into
region-level summaries, or proceed directly to differential analysis with `ep.tl.dmc()`.
