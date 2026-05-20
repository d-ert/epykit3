# HMR / LMR Calling

!!! warning "Experimental"
    This feature is implemented but has not been validated on real biological data.
    The API may change in future releases.

`ep.tl.hmr(md)` identifies hypo-methylated regions (HMRs) and low-methylated
regions (LMRs) from per-CpG beta values using a MethylSeekR-style 2-state
hidden Markov model. Results are stored in `md.uns["hmr"]` and
`md.uns["lmr"]`.

## Basic Usage

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.unite(md)

ep.tl.hmr(md)

hmrs = md.uns["hmr"]
lmrs = md.uns["lmr"]
print(f"Found {len(hmrs)} HMRs and {len(lmrs)} LMRs")
```

## How It Works

1. A 2-state HMM is fitted to raw per-CpG beta values for each sample.
   The two states are:
   - **Hypo** (low methylation): beta typically < 0.3
   - **Hyper** (high methylation): beta typically > 0.7
2. The Viterbi path assigns each CpG to one of the two states.
3. Consecutive CpGs in the hypo state are merged into candidate regions.
4. Candidate regions are then split into **HMRs** and **LMRs** based on
   CpG density:
   - **HMRs**: hypo-state runs with high CpG density (above
     `lmr_max_density`). These correspond to CpG islands and active
     promoters.
   - **LMRs**: hypo-state runs with low CpG density (below
     `lmr_max_density`). These correspond to distal regulatory elements
     such as enhancers.

## HMRs vs. LMRs

| Feature | HMR | LMR |
|---------|-----|-----|
| CpG density | High | Low |
| Typical location | CpG islands, promoters | Intergenic, intronic enhancers |
| Biological role | Active promoter marks | Active distal regulatory elements |
| Size | 200 bp -- 5 kb | 200 bp -- 2 kb |
| Methylation level | Very low (< 0.2) | Low (0.2 -- 0.4) |

## Output Columns

Both `md.uns["hmr"]` and `md.uns["lmr"]` share the same schema:

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | str | Chromosome |
| `start` | int | Region start position (0-based) |
| `end` | int | Region end position |
| `sample_id` | str | Sample identifier |
| `n_cpgs` | int | Number of CpGs in the region |
| `mean_beta` | float | Mean beta value within the region |
| `cpg_density` | float | CpGs per kilobase within the region |
| `length_bp` | int | Length of the region in base pairs |

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | MethylData | required | Analysis object |
| `lmr_max_density` | float | `30.0` | CpG density threshold (CpGs/kb) to separate HMRs from LMRs. Regions above this are HMRs; below are LMRs. |
| `min_cpgs` | int | `3` | Minimum CpGs in a region to report it |
| `chromosomes` | list | None | Restrict to specific chromosomes |

## Relationship to Other Analyses

HMR/LMR calling is a **within-sample** segmentation, complementary to
between-group methods:

- **DMC/DMR** detects differences *between* groups. HMR/LMR identifies
  regulatory features *within* a single sample.
- **PMD** detects megabase-scale partial methylation. HMR/LMR detects
  kilobase-scale hypomethylation at regulatory elements.

A common workflow combines all three:

```python
ep.tl.dmc(md)
ep.tl.dmr(md)
ep.tl.hmr(md)
ep.tl.pmd(md)

# Annotate DMRs by overlap with HMR/LMR/PMD features
ep.tl.annotate(md, gtf="gencode.v44.gtf.gz")
```

## Identifying Cell-type-specific Enhancers

LMRs are strong markers of active enhancers. Comparing LMR sets between
cell types or conditions reveals condition-specific regulatory elements:

```python
import polars as pl

lmrs = md.uns["lmr"]

# Count LMRs per sample
per_sample = lmrs.group_by("sample_id").agg(
    pl.col("length_bp").count().alias("n_lmrs"),
    pl.col("length_bp").sum().alias("total_lmr_bp"),
)
print(per_sample.sort("n_lmrs", descending=True))
```

## Exporting Results

HMR and LMR regions can be exported to BED format for use in genome
browsers or downstream intersection analyses:

```python
ep.export.to_bed(md.uns["hmr"], "hmrs.bed")
ep.export.to_bed(md.uns["lmr"], "lmrs.bed")
```

## Notes

- HMR/LMR calling requires per-CpG beta values with reasonable coverage.
  WGBS at 5x or higher is recommended. RRBS data will produce valid
  results but only within the captured regions.
- The `lmr_max_density` threshold is the primary tuning knob. The default
  of 30 CpGs/kb is derived from MethylSeekR defaults. Adjust upward if
  your genome has unusually high CpG density (e.g., CpG-rich species).
- Regions are detected independently per sample. To find shared HMRs
  across samples, intersect the per-sample BED outputs.
