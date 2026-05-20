# Analysis Overview

The `ep.tl.*` namespace contains all analysis tools in epykit. Functions operate
on a `MethylData` object and store results in `md.obs`, `md.varm`, or `md.uns`
depending on the output type.

## Analysis Functions

| Function | Description |
|----------|-------------|
| [`tl.qc()`](qc.md) | QC metrics (global methylation, coverage, clinical checks) |
| [`tl.dmc()`](dmc.md) | Per-CpG differential methylation calling |
| [`tl.dmr()`](dmr.md) | Differentially methylated region calling |
| [`tl.dvc()`](dvc.md) | Differentially variable CpG calling |
| [`tl.dvr()`](dvc.md#dvr) | Differentially variable region calling |
| [`tl.annotate()`](annotate.md) | Gene-feature and CpG-island annotation |
| `tl.age_clock()` | Epigenetic age clocks (experimental) |
| `tl.deconvolve()` | Cell-type deconvolution (experimental) |
| `tl.asm()` | Allele-specific methylation (experimental) |
| `tl.entropy()` | Methylation entropy (experimental) |
| `tl.pmd()` | Partially methylated domains (experimental) |
| `tl.hmr()` | Hypo-/low-methylated regions (experimental) |

## Typical Pipeline

The standard analysis flow follows a linear progression:

```
QC  -->  DMC  -->  DMR  -->  annotate
```

1. **QC** -- Assess data quality (coverage, conversion rate, sample identity).
   Remove or flag problematic samples before testing.
2. **DMC** -- Identify individual CpG sites with significant methylation
   differences between groups.
3. **DMR** -- Aggregate nearby significant CpGs into differentially methylated
   regions.
4. **Annotate** -- Map DMCs and DMRs to gene features (promoters, exons, etc.)
   and CpG-island context (island, shore, shelf, open-sea).

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)

ep.tl.qc(md)
ep.tl.dmc(md)
ep.tl.dmr(md)
ep.tl.annotate(md, gtf="gencode.v44.gtf.gz")
```

!!! tip "Differential variability"
    If you are studying cancer or aging methylomes where variance shifts matter
    more than mean shifts, add `ep.tl.dvc(md)` and `ep.tl.dvr(md)` after the
    DMC step. See the [Differential Variability](dvc.md) page.

## Where Results Are Stored

| Result | Location | Type |
|--------|----------|------|
| Per-sample QC metrics | `md.obs` | Polars DataFrame |
| DMC tables | `md.varm["dmc_<test>"]` | Polars DataFrame |
| DMR table | `md.uns["dmr"]` | Polars DataFrame |
| DVC table | `md.varm["dvc"]` | Polars DataFrame |
| DVR table | `md.uns["dvr"]` | Polars DataFrame |
| Annotation metadata | `md.uns["annotation"]` | dict |
| Sex check results | `md.uns["qc_sex_check"]` | Polars DataFrame |
| Sample correlations | `md.uns["qc_sample_correlation"]` | Polars DataFrame |
