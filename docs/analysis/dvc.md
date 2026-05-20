# Differential Variability

Differential variability analysis identifies CpG sites and regions where the
**between-replicate variance** differs between groups, even when mean
methylation levels do not change. This is a hallmark of cancer and aging
methylomes, where epigenetic drift increases variability at specific loci
without shifting the population mean.

## DVC: Differentially Variable CpGs

`ep.tl.dvc(md)` runs an iEVORA-style variance-equality test at every CpG site.

### Basic Usage

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)

ep.tl.dvc(md)

dvc_results = md.varm["dvc"]
n_dvc = dvc_results.filter(dvc_results["is_dvc"]).height
print(f"Found {n_dvc} differentially variable CpGs")
```

### How It Works

At each CpG site, `dvc` computes:

1. **Variance test** -- A Bartlett test for equality of variances between the
   treatment and control groups. Sites with significantly different variance
   receive a small `p_variance`.
2. **Mean filter** -- Sites are flagged as DVC only when `p_mean >
   mean_filter_alpha` (default 0.05), ensuring the variance difference is not
   simply a byproduct of a mean shift.

A CpG is classified as DVC when both conditions hold: significant variance
difference **and** non-significant mean difference.

### Output Columns

Results are stored at `md.varm["dvc"]`:

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | str | Chromosome |
| `pos` | int | Genomic position |
| `strand` | str | Strand |
| `n_treatment` | int | Treatment samples with coverage |
| `n_control` | int | Control samples with coverage |
| `var_treatment` | float | Between-replicate variance in treatment |
| `var_control` | float | Between-replicate variance in control |
| `var_log_ratio` | float | Log ratio of treatment to control variance |
| `p_variance` | float | Bartlett test p-value for variance equality |
| `q_variance` | float | BH-corrected q-value for variance test |
| `p_mean` | float | Mean-difference test p-value |
| `q_mean` | float | BH-corrected q-value for mean test |
| `is_dvc` | bool | True if the site passes both filters |

### DVC Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | MethylData | required | Analysis object |
| `test` | str | `"bartlett"` | Variance-equality test (only `"bartlett"` supported) |
| `chromosomes` | list | None | Restrict to specific chromosomes |
| `alpha` | float | 0.05 | q-value threshold for variance test |
| `mean_filter_alpha` | float | 0.05 | p-value threshold for mean filter (sites must **exceed** this) |
| `backend` | str | `"sequential"` | Execution backend |
| `n_workers` | int | None | Worker pool size |

Summary metadata is stored at `md.uns["dvc"]`:

```python
print(md.uns["dvc"])
# {'test': 'bartlett', 'alpha': 0.05, 'mean_filter_alpha': 0.05,
#  'n_sites': 22000000, 'n_dvc': 15423, 'unite': True}
```

## DVR: Differentially Variable Regions {: #dvr }

`ep.tl.dvr(md)` aggregates DVC sites into differentially variable regions
using density-based tile aggregation with per-tile binomial enrichment against
the genome-wide DVC rate.

### Basic Usage

```python
# Requires ep.tl.dvc(md) to have been run first
ep.tl.dvr(md)

dvr_results = md.uns["dvr"]
n_dvr = dvr_results.filter(dvr_results["is_dvr"]).height
print(f"Found {n_dvr} differentially variable regions")
```

### DVR Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | MethylData | required | Analysis object (must have DVC results) |
| `tile_size_bp` | int | 1000 | Tile width in base pairs |
| `min_cpgs_per_tile` | int | 5 | Minimum CpGs per tile |
| `alpha` | float | 0.05 | BH q-value threshold for the `is_dvr` flag |

Results are stored at `md.uns["dvr"]`.

## Use Case: Cancer Methylomes

In cancer, epigenetic instability often manifests as increased methylation
variability at CpG islands and shores, even at loci where the mean methylation
level is unchanged. Standard DMC analysis (mean-based) misses these loci
entirely.

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)

# Run both DMC and DVC
ep.tl.dmc(md)
ep.tl.dvc(md)
ep.tl.dvr(md)

# Compare: sites found by DMC vs DVC
dmc = md.dmc.filter(md.dmc["qvalue"] < 0.05)
dvc = md.varm["dvc"].filter(md.varm["dvc"]["is_dvc"])
print(f"DMC hits: {len(dmc)}")
print(f"DVC hits: {len(dvc)}")

# Annotate both
ep.tl.annotate(md, gtf="gencode.v44.gtf.gz")
```

!!! tip "Combine DMC and DVC for a complete picture"
    DVC captures loci with altered variability (early epigenetic drift),
    while DMC captures loci with altered mean methylation (established
    changes). Together they provide a more complete characterisation of
    the methylation landscape.
