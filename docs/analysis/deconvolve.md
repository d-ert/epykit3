# Cell-type Deconvolution

!!! warning "Experimental"
    This feature is implemented but has not been validated on real biological data.
    The API may change in future releases.

`ep.tl.deconvolve(md)` estimates cell-type proportions from bulk methylation
data using reference-based non-negative least squares (NNLS) deconvolution.
Results are stored in `md.uns["deconvolution"]` (long-form) and as
`frac_<cell_type>` columns on `md.obs` (wide-form).

## Basic Usage

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.unite(md)

ep.tl.deconvolve(
    md,
    reference="blood_reference_matrix.csv",
    manifest="HM450.hg38.manifest.tsv.gz",
)

# Long-form proportions
print(md.uns["deconvolution"])

# Wide-form: one column per cell type on md.obs
print(md.obs.select("sample_id", "frac_CD4T", "frac_CD8T", "frac_Mono"))
```

## How It Works

1. The reference matrix (cell types x CpGs) is loaded. Each column is a
   cell-type-specific methylation profile at a set of informative CpG loci.
2. The manifest maps CpG probe IDs to genomic positions.
3. For each sample, beta values at the reference CpG positions are extracted
   from the methylstore.
4. NNLS solves for the non-negative mixing coefficients that best
   reconstruct the observed bulk profile as a linear combination of the
   reference profiles.
5. Coefficients are normalized to sum to 1, giving estimated cell-type
   fractions.

## Reference Matrix Format

The reference matrix CSV has cell types as columns and CpG probe IDs as rows:

```
probe_id,CD4T,CD8T,Mono,Bcell,NK,Neutro
cg00000029,0.85,0.82,0.12,0.78,0.80,0.10
cg00000165,0.05,0.07,0.91,0.08,0.06,0.88
...
```

Standard references include:

| Reference | Cell types | Source |
|-----------|------------|--------|
| Reinius blood | 6 (CD4T, CD8T, Mono, Bcell, NK, Gran) | Reinius et al. 2012 |
| IDOL extended | 12 immune subtypes | Salas et al. 2022 |
| Houseman brain | 4 (NeuN+, NeuN-, Olig2+, Other) | Houseman et al. 2012 |

The function does not ship reference matrices. Obtain them from the original
publications or from the `FlowSorted.*` Bioconductor packages.

## Output

### Long-form: `md.uns["deconvolution"]`

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | str | Sample identifier |
| `cell_type` | str | Cell type name from the reference |
| `fraction` | float | Estimated proportion (0 to 1) |

### Wide-form: `md.obs`

Each cell type gets a column named `frac_<cell_type>`:

```python
md.obs.select("sample_id", "frac_CD4T", "frac_CD8T", "frac_Mono",
              "frac_Bcell", "frac_NK", "frac_Neutro")
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | MethylData | required | Analysis object |
| `reference` | str / Path / DataFrame | required | Reference matrix (cell types x CpGs) |
| `manifest` | str / Path / DataFrame | required | Illumina CpG manifest mapping probe IDs to positions |
| `method` | str | `"nnls"` | Deconvolution method. Currently only `"nnls"` is supported. |
| `normalize` | bool | `True` | Normalize fractions to sum to 1 |

## Use Case: Adjusting for Cell Composition

Cell-type heterogeneity is a major confounder in bulk WGBS. After
deconvolution, the estimated fractions can be included as covariates in
differential methylation testing:

```python
ep.tl.deconvolve(md, reference="blood_ref.csv",
                 manifest="HM450.hg38.manifest.tsv.gz")

# Use cell fractions as covariates in DMC testing
ep.tl.dmc(md, formula="~ group + frac_CD4T + frac_Mono + frac_Neutro")
```

## Notes

- Deconvolution accuracy depends on the overlap between the reference CpG
  set and the sites covered in your data. Check
  `md.uns["deconvolution_meta"]["cpgs_matched"]` for the number of matched
  loci.
- NNLS assumes that the bulk profile is a linear mixture of the reference
  profiles. This assumption breaks down when the true cell types are absent
  from the reference or when the methylation profiles are non-linear.
- For WGBS data (as opposed to array data), coverage variability across
  sites may reduce deconvolution precision. Consider restricting to
  well-covered sites with `ep.pp.filter_coverage(md, min_coverage=10)`
  before running deconvolution.
