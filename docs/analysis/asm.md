# Allele-specific Methylation

!!! warning "Experimental"
    This feature is implemented but has not been validated on real biological data.
    The API may change in future releases.

`ep.tl.asm(md, bam=..., vcf=...)` detects allele-specific methylation (ASM)
at individual CpG sites. For each sample it phases reads by their base at
heterozygous SNVs, tallies methylated and unmethylated calls per allele at
every CpG the phased reads cover, and tests each CpG with Fisher's exact
test. Results are stored in `md.varm["asm"]`; run metadata in
`md.uns["asm"]`.

## Prerequisites

ASM analysis requires read-level access to BAM files and a per-individual
VCF. Install the BAM extras:

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

ep.tl.asm(
    md,
    bam={"tumor_1": "aligned/tumor_1.bam", "normal_1": "aligned/normal_1.bam"},
    vcf="genotypes.vcf.gz",
)

asm_results = md.varm["asm"]
print(f"Tested {asm_results.height} (sample, CpG) pairs for ASM")
```

The lower-level `epykit.asm.call_asm(bam, vcf)` returns the same frame
without a `MethylData` object and also exposes `min_baseq` and `min_mapq`.

## How It Works

1. Per-read methylation calls are extracted from each BAM (Bismark `XM`
   tags, or MethylDackel `MM`/`ML` tags with `caller="methyldackel"`).
2. For each heterozygous biallelic SNV in the VCF, reads overlapping the
   SNV are assigned to haplotype 1 (REF base) or haplotype 2 (ALT base),
   but only when the read's bisulfite conversion strand cannot convert
   either allele (see [Bisulfite-safe phasing anchors](#bisulfite-safe-phasing-anchors)).
   A read whose assignment disagrees between SNVs is dropped, and a read
   must cover `min_phased_snvs` agreeing anchors to be phased.
3. At each CpG the phased reads cover, methylated and unmethylated calls
   are tallied per haplotype.
4. CpGs with at least `min_reads_per_haplotype` calls on both haplotypes
   are tested with Fisher's exact test.
5. p-values are BH-corrected across all samples and CpGs.

## Bisulfite-safe phasing anchors

Bisulfite conversion changes the base a read shows at a SNV. On reads
aligned to the converted top strand (Bismark tag `XG:Z:CT`) an
unmethylated C reads as T. On reads from the converted bottom strand
(`XG:Z:GA`) an unmethylated G reads as A. At a C/T anchor, a CT-strand
read that shows T is either a true T allele or an unmethylated C allele,
so assigning it by the raw base would sort reads by methylation state and
fabricate ASM.

epykit therefore reads the `XG` tag of every read and phases it only when
its conversion strand cannot convert either allele of the unordered
REF/ALT class:

| SNV class | `XG:Z:CT` reads | `XG:Z:GA` reads | Missing or unrecognised `XG` |
|-----------|-----------------|-----------------|------------------------------|
| A/T | phased | phased | phased |
| A/G | phased | dropped | dropped |
| G/T | phased | dropped | dropped |
| C/T | dropped | phased | dropped |
| A/C | dropped | phased | dropped |
| C/G | dropped | dropped | dropped |

C/G anchors, and records whose REF or ALT is not A, C, G or T, are skipped
before any read is fetched. REF and ALT are compared in upper case.

**Untagged reads.** Only Bismark writes `XG`. MethylDackel or bwa-meth
BAMs, and Bismark BAMs with the tag stripped, phase A/T anchors only.
There is no option to restore phasing by the raw base.

**What the rule costs.** The rule is deliberately conservative. It drops
every read that a methylation-aware genotyper such as Bis-SNP or BISCUIT
could still have assigned. Compared with releases before 1.2, fewer
anchors and reads phase: fabricated sites disappear, and some genuine
sites those releases reported are no longer testable. One INFO log line
per sample (`phasing summary`) reports how many anchors phased at least
one read, how many anchors were rejected by class before any read was
fetched, and how many read-anchor observations the strand rule rejected.

**Limitations.** Strand-aware exclusion is a substitute for
methylation-aware genotyping, not an implementation of it. The caller does
not fold converted bases back to their allele, does not read the `XR`
tag, and has not been validated on real biological data.

## Output Columns

`md.varm["asm"]` has one row per (sample, CpG) that reached
`min_reads_per_haplotype` calls on both haplotypes:

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | str | Sample identifier (the key in `bam`) |
| `chrom` | str | Chromosome |
| `pos` | int | CpG position (0-based) |
| `strand` | str | Strand of the reads that carry the calls (`+` or `-`) |
| `h1_meth`, `h1_unmeth` | int | Methylated and unmethylated calls on haplotype 1 (REF allele) |
| `h2_meth`, `h2_unmeth` | int | Methylated and unmethylated calls on haplotype 2 (ALT allele) |
| `n_h1`, `n_h2` | int | Total calls per haplotype |
| `meth_diff` | float | Methylation fraction difference (haplotype 1 minus haplotype 2) |
| `pvalue` | float | Fisher exact test p-value |
| `qvalue` | float | BH-adjusted q-value |
| `reject` | bool | `qvalue < 0.05` |

`md.uns["asm"]` records `n_sites`, `n_samples`, `min_reads_per_haplotype`,
`min_phased_snvs` and the `vcf` path of the run.

The `chrom`, `pos`, `meth_diff`, `pvalue` and `qvalue` columns share their
names with the DMC tables, so the same filter expressions apply.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | MethylData | required | Analysis object; every key of `bam` must appear in `md.obs["sample_id"]` |
| `bam` | mapping | required | `{sample_id: bam_path}`; each BAM coordinate-sorted and indexed |
| `vcf` | str / Path | required | Per-individual VCF (bgzipped and tabix-indexed preferred) |
| `min_reads_per_haplotype` | int | `10` | Minimum calls on each haplotype to test a CpG |
| `min_phased_snvs` | int | `1` | Minimum agreeing heterozygous SNVs a read must cover to be phased |
| `chromosomes` | list | `None` | Restrict to these chromosomes |
| `caller` | str | `"bismark"` | `"bismark"` (XM tags) or `"methyldackel"` (MM/ML tags) |

## Input File Requirements

### BAM files

- One BAM per sample, coordinate-sorted and indexed (`.bai`).
- Bismark BAMs carry the `XM` methylation string and the `XG`
  genome-conversion tag on every read. MethylDackel BAMs carry `MM`/`ML`
  tags and no `XG`, so only A/T anchors phase.
- Secondary and supplementary alignments and reads below the
  mapping-quality threshold are ignored.

### VCF file

- One individual per VCF; the first sample column is used.
- Heterozygous biallelic SNVs (`GT` of `0/1`, `0|1` or `1|0`) are the
  anchors. Phased genotypes are not required: haplotype 1 is the REF
  allele and haplotype 2 the ALT allele at every anchor, and reads that
  disagree between anchors are dropped.
- Indels and multi-allelic records are skipped.
- bgzipped (`.vcf.gz`) with a tabix index (`.vcf.gz.tbi`) is preferred.

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
sig_asm = md.varm["asm"].filter(
    (pl.col("qvalue") < 0.05) & (pl.col("meth_diff").abs() > 0.2)
)
print(f"Found {sig_asm.height} significant ASM sites")
```
