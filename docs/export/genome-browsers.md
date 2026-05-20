# Genome Browser Exports

Export methylation data and differential methylation results as track files
for genome browsers such as IGV, UCSC Genome Browser, and JBrowse.

## ep.to_bedgraph

Write a per-sample BedGraph file containing methylation values for a single sample.

```python
ep.to_bedgraph(md, sample, output, value="beta")
```

**Parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `md` | `MethylData` | MethylData object with loaded store |
| `sample` | `str` | Sample ID to export |
| `output` | `str` or `Path` | Output file path (`.bedgraph`) |
| `value` | `str` | Value column to write. `"beta"` (default) for methylation fraction, or `"coverage"` for read depth |

**Example**

```python
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")

# Export beta values for a single sample
ep.to_bedgraph(md, sample="tumor_1", output="tumor_1_beta.bedgraph")

# Export coverage track
ep.to_bedgraph(md, sample="tumor_1", output="tumor_1_cov.bedgraph", value="coverage")
```

The output is a standard four-column BedGraph (chrom, start, end, value) that can be
loaded directly into IGV or uploaded to the UCSC Genome Browser.

---

## ep.to_bigwig

Write a BigWig file for a single sample. BigWig is a compressed, indexed binary format
that loads faster than BedGraph in genome browsers and is required for large-scale
track hubs.

```python
ep.to_bigwig(md, sample, output)
```

**Parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `md` | `MethylData` | MethylData object with loaded store |
| `sample` | `str` | Sample ID to export |
| `output` | `str` or `Path` | Output file path (`.bw`) |

**Installation**

BigWig export depends on `pyBigWig`, which is included in the `export` extras group:

```bash
pip install 'epykit[export]'
```

Note: `pyBigWig` does not provide a Windows wheel. This function is available on
Linux and macOS only.

**Example**

```python
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")

ep.to_bigwig(md, sample="tumor_1", output="tumor_1.bw")
```

---

## ep.dmcs_to_bed

Export differentially methylated CpG (DMC) results as a BED file. Requires that
DMC analysis has already been run on the MethylData object.

```python
ep.dmcs_to_bed(md, output, pval_threshold=0.05, delta_threshold=0.1)
```

**Parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `md` | `MethylData` | MethylData object with DMC results attached |
| `output` | `str` or `Path` | Output file path (`.bed`) |
| `pval_threshold` | `float` | Adjusted p-value cutoff (default: `0.05`) |
| `delta_threshold` | `float` | Minimum absolute delta-beta cutoff (default: `0.1`) |

**Example**

```python
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")
ep.dmc(md, formula="~ treatment", contrast=("treatment", "tumor", "normal"))

# Export all significant DMCs
ep.dmcs_to_bed(md, output="dmcs.bed")

# Stricter filtering
ep.dmcs_to_bed(md, output="dmcs_strict.bed", pval_threshold=0.01, delta_threshold=0.2)
```

The BED file includes the CpG position, delta-beta as the score column, and
strand information, suitable for intersection with other genomic annotations.

---

## ep.dmrs_to_bed

Export differentially methylated region (DMR) results as a BED file. Requires that
DMR calling has already been run.

```python
ep.dmrs_to_bed(md, output)
```

**Parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `md` | `MethylData` | MethylData object with DMR results attached |
| `output` | `str` or `Path` | Output file path (`.bed`) |

**Example**

```python
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")
ep.dmc(md, formula="~ treatment", contrast=("treatment", "tumor", "normal"))
ep.dmr(md, method="chain_merge")

# Export DMRs
ep.dmrs_to_bed(md, output="dmrs.bed")
```

Each row represents one DMR with columns for the region coordinates, the number of
CpGs in the region, the mean delta-beta, and the region-level p-value.

---

## Typical Workflow

A common pattern is to generate all browser tracks after a differential analysis:

```python
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")
ep.pp.filter_coverage(md, min_cov=10)
ep.dmc(md, formula="~ treatment", contrast=("treatment", "tumor", "normal"))
ep.dmr(md, method="chain_merge")

# Per-sample tracks
for sample in md.samples:
    ep.to_bedgraph(md, sample=sample, output=f"{sample}.bedgraph")

# Differential results
ep.dmcs_to_bed(md, output="dmcs.bed")
ep.dmrs_to_bed(md, output="dmrs.bed")
```
