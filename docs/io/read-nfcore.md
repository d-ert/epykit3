# nf-core/methylseq

epykit provides two functions for working with nf-core/methylseq pipeline outputs:

- `ep.read_nfcore_methylseq()` -- Read methylation calls directly from a pipeline run
  directory.
- `ep.read_nfcore_methylseq_qc()` -- Extract QC metrics from the pipeline's Bismark
  and Qualimap reports.

## read_nfcore_methylseq

### Overview

`ep.read_nfcore_methylseq()` scans an nf-core/methylseq output directory, discovers
coverage files for each sample, and writes them into a partitioned Parquet methylstore.
This eliminates the need to manually build a samplesheet -- sample names and file paths
are inferred from the pipeline directory structure.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `run_dir` | `str` | required | Path to the nf-core/methylseq output directory |
| `treatment_group` | `str` | required | Group label for treatment samples |
| `control_group` | `str` | required | Group label for control samples |
| `treatment_samples` | `list[str]` | required | Sample names to assign to the treatment group |
| `control_samples` | `list[str]` | required | Sample names to assign to the control group |
| `assembly` | `str` | required | Genome assembly name (e.g. `"hg38"`) |
| `store_dir` | `str` | `"methylstore"` | Directory for the partitioned Parquet store |
| `context` | `str` | `"CpG"` | Cytosine context to retain |

### Usage

```python
import epykit as ep

md = ep.read_nfcore_methylseq(
    run_dir="/data/nfcore_methylseq/results",
    treatment_group="tumor",
    control_group="normal",
    treatment_samples=["tumor_1", "tumor_2", "tumor_3"],
    control_samples=["normal_1", "normal_2", "normal_3"],
    assembly="hg38",
    store_dir="results/methylstore",
)
```

The function automatically locates the Bismark coverage files within the standard
nf-core/methylseq directory layout:

```
results/
├── bismark/
│   ├── alignments/
│   ├── deduplicated/
│   └── methylation_extraction/
│       ├── tumor_1.bismark.cov.gz
│       ├── tumor_2.bismark.cov.gz
│       └── ...
├── multiqc/
└── pipeline_info/
```

### Output

A `MethylData` object backed by a partitioned Parquet store, identical to the output of
`ep.read_bismark()`. See [read_bismark](read-bismark.md) for details on the store layout.

---

## read_nfcore_methylseq_qc

### Overview

`ep.read_nfcore_methylseq_qc()` parses QC outputs from an nf-core/methylseq run and
returns a summary DataFrame. It extracts metrics from Bismark alignment reports, Bismark
M-bias reports, and Qualimap results.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `run_dir` | `str` | required | Path to the nf-core/methylseq output directory |

### Usage

```python
import epykit as ep

qc = ep.read_nfcore_methylseq_qc(
    run_dir="/data/nfcore_methylseq/results",
)

# qc is a DataFrame with one row per sample
print(qc.columns.tolist())
# ['sample_name', 'total_reads', 'aligned_reads', 'alignment_rate',
#  'deduplicated_reads', 'duplication_rate', 'cpg_coverage_mean',
#  'cpg_coverage_median', 'bisulfite_conversion_rate', ...]
```

### Returned Metrics

The QC DataFrame includes (where available):

| Metric | Source | Description |
|--------|--------|-------------|
| `total_reads` | Bismark | Total sequenced reads |
| `aligned_reads` | Bismark | Uniquely aligned reads |
| `alignment_rate` | Bismark | Fraction of reads aligned |
| `deduplicated_reads` | Bismark | Reads after deduplication |
| `duplication_rate` | Bismark | Fraction of duplicates removed |
| `cpg_coverage_mean` | Qualimap | Mean coverage at CpG sites |
| `cpg_coverage_median` | Qualimap | Median coverage at CpG sites |
| `bisulfite_conversion_rate` | Bismark | Estimated conversion efficiency |

### Combining Data and QC

A typical workflow reads both data and QC in sequence:

```python
import epykit as ep

# Read methylation data
md = ep.read_nfcore_methylseq(
    run_dir="/data/nfcore_methylseq/results",
    treatment_group="tumor",
    control_group="normal",
    treatment_samples=["tumor_1", "tumor_2"],
    control_samples=["normal_1", "normal_2"],
    assembly="hg38",
)

# Read QC metrics
qc = ep.read_nfcore_methylseq_qc(
    run_dir="/data/nfcore_methylseq/results",
)

# Inspect QC before proceeding
print(qc[["sample_name", "alignment_rate", "bisulfite_conversion_rate"]])

# Continue with preprocessing
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)
```

## Next Steps

After reading, proceed to [Preprocessing](../preprocessing/index.md) for coverage
filtering, normalization, and site unification.
