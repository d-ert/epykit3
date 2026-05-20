# Methylation Entropy

!!! warning "Experimental"
    This feature is implemented but has not been validated on real biological data.
    The API may change in future releases.

`ep.tl.entropy(md)` computes read-level Shannon entropy over CpG-window
methylation patterns. The output quantifies epigenetic disorder at each
genomic locus, with values ranging from fully ordered (0) to fully
disordered (1). Results are stored in `md.uns["entropy"]`.

## Prerequisites

Entropy calculation requires read-level access to BAM files. Install the
BAM extras:

```bash
pip install 'epykit[bam]'
```

!!! note "Platform restriction"
    The `[bam]` extra depends on `pysam`, which is available on Linux and
    macOS only. Windows is not supported.

## Basic Usage

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.unite(md)

ep.tl.entropy(md, bam_dir="aligned/")

entropy_results = md.uns["entropy"]
print(entropy_results.head())
```

## How It Works

1. A sliding window of `window_cpgs` consecutive CpG sites is defined
   across the genome.
2. For each window, reads spanning all CpGs in the window are extracted
   from the BAM file.
3. Each read produces a binary methylation pattern (e.g., `1010` for a
   4-CpG window where sites 1 and 3 are methylated).
4. Shannon entropy is calculated over the distribution of observed
   patterns within the window.
5. The raw entropy is divided by the theoretical maximum
   (log2 of the number of possible patterns) to yield a normalized
   entropy in [0, 1].

### Interpreting Normalized Entropy

| Value | Meaning | Biological context |
|-------|---------|--------------------|
| 0.0 | All reads share the same pattern | Homogeneous methylation (fully methylated or fully unmethylated) |
| 0.0 -- 0.3 | Low disorder | Ordered loci, e.g. stably silenced promoters |
| 0.3 -- 0.7 | Moderate disorder | Partially methylated, transitional states |
| 0.7 -- 1.0 | High disorder | Stochastic methylation, epigenetic drift |
| 1.0 | All possible patterns equally frequent | Maximum disorder |

## Output Columns

The result DataFrame at `md.uns["entropy"]` contains:

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | str | Chromosome |
| `start` | int | Window start position (0-based) |
| `end` | int | Window end position |
| `sample_id` | str | Sample identifier |
| `n_reads` | int | Reads spanning the full window |
| `n_patterns` | int | Distinct methylation patterns observed |
| `entropy` | float | Raw Shannon entropy (bits) |
| `normalised_entropy` | float | Entropy normalized to [0, 1] |

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | MethylData | required | Analysis object |
| `bam_dir` | str / Path | required | Directory containing BAM files (one per sample, indexed) |
| `window_cpgs` | int | `4` | Number of consecutive CpGs per window |
| `min_reads` | int | `10` | Minimum reads spanning the full window to report a result |
| `chromosomes` | list | None | Restrict to specific chromosomes |

## Biological Applications

### Aging drift

Methylation entropy increases with age at many loci, reflecting a gradual
loss of epigenetic maintenance fidelity. Comparing entropy distributions
between age groups can reveal loci undergoing age-related drift.

```python
import polars as pl

entropy = md.uns["entropy"]

# Mean per-sample entropy
per_sample = entropy.group_by("sample_id").agg(
    pl.col("normalised_entropy").mean().alias("mean_entropy")
)
print(per_sample.sort("mean_entropy"))
```

### Intra-tumour heterogeneity

Tumour samples often show elevated entropy at regulatory regions,
reflecting clonal diversity in methylation states. High-entropy loci in
tumours but not matched normals may indicate epigenetically heterogeneous
regions under selection.

### Stochastic methylation

Regions with intermediate mean beta (0.3 -- 0.7) can arise from two
distinct mechanisms: (a) a bimodal mix of fully methylated and fully
unmethylated alleles (low entropy), or (b) genuinely disordered
methylation across reads (high entropy). Entropy distinguishes these
cases where mean beta alone cannot.

## Notes

- Window size (`window_cpgs`) controls the resolution-sensitivity tradeoff.
  Larger windows capture longer-range coordination but require more reads
  to span the full window, reducing coverage.
- Only reads that span all CpGs in the window are used. Paired-end reads
  with large insert sizes provide better coverage of wider windows.
- BAM files must be coordinate-sorted, indexed, and contain bisulfite
  conversion tags (Bismark XM/XG format).
