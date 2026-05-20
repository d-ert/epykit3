# Partially Methylated Domains

!!! warning "Experimental"
    This feature is implemented but has not been validated on real biological data.
    The API may change in future releases.

`ep.tl.pmd(md)` detects partially methylated domains (PMDs) -- megabase-scale
genomic regions with reduced methylation levels -- using a per-sample 2-state
hidden Markov model (HMM) on smoothed beta values. Results are stored in
`md.uns["pmd"]`.

## Basic Usage

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.unite(md)

ep.tl.pmd(md)

pmds = md.uns["pmd"]
print(f"Found {len(pmds)} PMD regions across all samples")
```

## How It Works

1. Per-sample beta values are smoothed using a rolling mean with a
   bandwidth of `bandwidth_bp` base pairs.
2. A 2-state HMM is fitted to the smoothed values. The two hidden states
   correspond to:
   - **Non-PMD** (high methylation): beta typically > 0.7
   - **PMD** (partial methylation): beta typically 0.5 -- 0.7
3. The Viterbi algorithm assigns each CpG to one of the two states.
4. Consecutive CpGs in the PMD state are merged into contiguous regions.
5. Regions shorter than `min_pmd_bp` are filtered out.

## PMDs vs. DMRs

PMDs and DMRs are fundamentally different features:

| Property | PMDs | DMRs |
|----------|------|------|
| Scale | Megabases (100 kb -- 10 Mb) | Kilobases (100 bp -- 10 kb) |
| Detection | Within-sample (HMM on one sample) | Between-group (statistical test) |
| Methylation level | Partial (0.5 -- 0.7 beta) | Variable |
| Biological origin | Replication-associated loss, heterochromatin | Regulatory, tissue-specific |
| Typical context | Cancer, aging, late-replicating regions | Active regulation, development |

PMDs are a within-sample structural feature. They represent large-scale
erosion of methylation and are not detected by DMC/DMR testing.

## Output Columns

The result DataFrame at `md.uns["pmd"]` contains:

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | str | Chromosome |
| `start` | int | PMD start position (0-based) |
| `end` | int | PMD end position |
| `sample_id` | str | Sample identifier |
| `n_cpgs` | int | Number of CpGs in the PMD |
| `mean_beta` | float | Mean beta value within the PMD |
| `length_bp` | int | Length of the PMD in base pairs |

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | MethylData | required | Analysis object |
| `bandwidth_bp` | int | `10000` | Smoothing bandwidth in base pairs |
| `beta_threshold` | float | `0.7` | Initial beta threshold to seed PMD vs. non-PMD states |
| `min_pmd_bp` | int | `100000` | Minimum PMD length in base pairs (100 kb) |
| `chromosomes` | list | None | Restrict to specific chromosomes |

## Use Cases

### Cancer methylomes

PMDs are a hallmark of cancer epigenomes. They often overlap
late-replicating, gene-poor regions and expand during tumour progression.
Comparing the fraction of the genome in PMDs between tumour and normal
samples quantifies global methylation erosion:

```python
import polars as pl

pmds = md.uns["pmd"]

# Total PMD coverage per sample
pmd_coverage = pmds.group_by("sample_id").agg(
    pl.col("length_bp").sum().alias("total_pmd_bp"),
    pl.col("length_bp").count().alias("n_pmds"),
)
print(pmd_coverage.sort("total_pmd_bp", descending=True))
```

### Aging

PMDs expand with age in somatic tissues. Tracking PMD boundaries across
age-stratified samples can reveal loci undergoing progressive methylation
loss.

### Cell-type identity

Different cell types have distinct PMD landscapes. PMD calls can be used
as coarse cell-type markers in heterogeneous tissue samples, complementing
CpG-level deconvolution approaches.

## Visualizing PMDs

PMD regions can be exported to BED format for visualization in genome
browsers, or plotted alongside per-CpG beta values:

```python
# Export PMDs to BED
ep.export.to_bed(md.uns["pmd"], "pmds.bed")

# Plot beta values with PMD shading for one chromosome
ep.pl.beta_track(md, chrom="chr1", sample="tumor_01", pmds=True,
                 save="chr1_pmd.png")
```

## Notes

- PMD detection requires reasonable genome-wide coverage. Sparse data
  (e.g., RRBS or targeted panels) will not produce meaningful PMD calls.
  WGBS at 5x or higher is recommended.
- The smoothing bandwidth (`bandwidth_bp`) controls sensitivity. Smaller
  values detect shorter PMDs but increase noise. The default of 10 kb
  is appropriate for typical WGBS data.
- PMDs are detected independently per sample. To identify shared PMD
  boundaries, intersect the per-sample BED outputs using standard
  genomic-interval tools.
