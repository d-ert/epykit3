# Allele-specific Methylation

!!! warning "Experimental"
    This feature is implemented but has not been validated on real biological data.
    The API may change in future releases.

`ep.tl.asm(md)` detects allele-specific methylation (ASM) at individual CpG
sites by partitioning reads by haplotype and testing for methylation
differences between alleles. Results are stored in `md.uns["asm"]`.

## Prerequisites

ASM analysis requires read-level access to BAM files and phased genotype
information. Install the BAM extras:

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

ep.tl.asm(
    md,
    bam_dir="aligned/",
    vcf="phased_genotypes.vcf.gz",
)

asm_results = md.uns["asm"]
print(f"Tested {len(asm_results)} CpG sites for ASM")
```

## How It Works

1. For each sample, aligned reads from the BAM file are loaded.
2. Reads overlapping heterozygous SNPs (from the VCF) are assigned to
   haplotype 1 or haplotype 2 based on the phased genotype.
3. At each CpG site near a phased het SNP, methylated and unmethylated
   read counts are tallied separately per haplotype.
4. A Fisher exact test compares the methylation proportions between the
   two haplotypes.
5. Results are corrected for multiple testing (BH procedure).

## Output Columns

The result DataFrame at `md.uns["asm"]` contains:

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | str | Chromosome |
| `pos` | int | CpG position (0-based) |
| `sample_id` | str | Sample identifier |
| `n_hap1` | int | Total reads on haplotype 1 |
| `n_hap2` | int | Total reads on haplotype 2 |
| `beta_hap1` | float | Methylation fraction on haplotype 1 |
| `beta_hap2` | float | Methylation fraction on haplotype 2 |
| `meth_diff` | float | Methylation difference (hap1 - hap2) |
| `pvalue` | float | Fisher exact test p-value |
| `qvalue` | float | BH-adjusted q-value |

This schema is compatible with the DMC result family, so ASM results work
directly with `ep.pl.volcano()` and `ep.pl.manhattan()`:

```python
ep.pl.volcano(md, source="asm", save="asm_volcano.png")
ep.pl.manhattan(md, source="asm", save="asm_manhattan.png")
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | MethylData | required | Analysis object |
| `bam_dir` | str / Path | required | Directory containing BAM files (one per sample, indexed) |
| `vcf` | str / Path | required | VCF file with phased genotypes (must be bgzipped and tabix-indexed) |
| `min_reads_per_haplotype` | int | `3` | Minimum reads on each haplotype to test a site |
| `max_distance_to_snp` | int | `1000` | Maximum distance (bp) from a CpG to a phased het SNP |
| `chromosomes` | list | None | Restrict to specific chromosomes |

## Input File Requirements

### BAM files

- One BAM per sample, named `{sample_id}.bam` inside `bam_dir`.
- Must be coordinate-sorted and indexed (`.bai`).
- Must contain bisulfite-converted reads with XM/XG tags (Bismark format)
  or MD tags.

### VCF file

- Must be bgzipped (`.vcf.gz`) and tabix-indexed (`.vcf.gz.tbi`).
- Must contain phased genotypes (`0|1` or `1|0` in the GT field).
- Samples in the VCF must match the `sample_id` values in the MethylData
  object.

## Interpreting ASM Results

Sites with significant ASM (low q-value and large absolute `meth_diff`)
indicate parent-of-origin effects, cis-regulatory variation, or imprinting.
Common follow-ups include:

- Overlapping ASM sites with known imprinted regions.
- Correlating ASM with nearby eQTL or meQTL effects.
- Checking whether ASM sites cluster in regulatory elements using
  `ep.tl.annotate()`.

```python
import polars as pl

# Filter to significant ASM sites
sig_asm = md.uns["asm"].filter(
    (pl.col("qvalue") < 0.05) & (pl.col("meth_diff").abs() > 0.2)
)
print(f"Found {len(sig_asm)} significant ASM sites")
```
