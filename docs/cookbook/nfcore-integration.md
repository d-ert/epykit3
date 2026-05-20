# nf-core/methylseq Integration

This guide covers using epykit to analyse output from the
[nf-core/methylseq](https://nf-co.re/methylseq) pipeline. The workflow has four
steps: parse QC metrics, ingest Bismark output, run the standard epykit
analysis, and optionally export results back to MultiQC format.

## 1. Parse QC Metrics

`ep.read_nfcore_methylseq_qc()` extracts per-sample QC metrics from the
pipeline's Bismark alignment reports, M-bias reports, and Qualimap output.

```python
import epykit as ep

qc = ep.read_nfcore_methylseq_qc(
    run_dir="/data/nfcore_methylseq/results",
)

print(qc[["sample_name", "alignment_rate", "bisulfite_conversion_rate",
           "cpg_coverage_mean"]])
```

Review these metrics before proceeding. Samples with low alignment rates,
poor bisulfite conversion, or unusually low coverage may need to be excluded.

The `samplesheet` parameter accepts the nf-core samplesheet CSV. When provided,
it is used to map pipeline sample names to your group labels:

```python
qc = ep.read_nfcore_methylseq_qc(
    samplesheet="nfcore_samplesheet.csv",
    run_dir="/data/nfcore_methylseq/results",
)
```

## 2. Ingest Bismark Output

`ep.read_nfcore_methylseq()` scans the pipeline output directory, locates the
Bismark coverage files, and writes them into a partitioned Parquet methylstore.

```python
md = ep.read_nfcore_methylseq(
    run_dir="/data/nfcore_methylseq/results",
    treatment_group="tumor",
    control_group="normal",
    treatment_samples=["tumor_1", "tumor_2", "tumor_3"],
    control_samples=["normal_1", "normal_2", "normal_3"],
    assembly="hg38",
    store_dir="methylstore",
)
```

The function expects the standard nf-core/methylseq output layout:

```
results/
├── bismark/
│   └── methylation_extraction/
│       ├── tumor_1.bismark.cov.gz
│       ├── tumor_2.bismark.cov.gz
│       └── ...
├── multiqc/
└── pipeline_info/
```

The returned `MethylData` object is identical to what `ep.read_bismark()`
produces and can be used with the full epykit pipeline.

## 3. Standard Analysis Pipeline

From this point, the workflow is the same as any other epykit analysis:

```python
# Preprocessing
ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)

# QC (epykit's own QC, complementing nf-core's upstream QC)
ep.tl.qc(md)

# DMC calling
ep.tl.dmc(md, power_stack=True, fdr_method="fdr_tsbh")

# DMR calling
ep.tl.dmr(md, method="chain_merge", preset="default")

# Annotation
ep.tl.annotate(md, gtf="gencode.v44.annotation.gtf.gz",
               cpg_islands="cpg_islands_hg38.bed")

# Visualisation
ep.pl.volcano(md)
ep.pl.manhattan(md)
ep.pl.pca(md)
```

See the [End-to-end Pipeline Walkthrough](end-to-end.md) for a detailed
explanation of each step.

## 4. Export to MultiQC

`ep.report_multiqc()` writes epykit's QC metrics and analysis summaries in
MultiQC-compatible format. This lets you combine epykit results with the
upstream nf-core QC in a single MultiQC report.

```python
ep.report_multiqc(md, output_dir="results/multiqc_epykit/")
```

This creates files that MultiQC can discover automatically:

```
results/multiqc_epykit/
├── epykit_general_stats_mqc.tsv
├── epykit_coverage_mqc.tsv
├── epykit_dmc_summary_mqc.tsv
└── epykit_dmr_summary_mqc.tsv
```

To generate a combined report, run MultiQC over both the nf-core output and the
epykit output:

```bash
multiqc /data/nfcore_methylseq/results results/multiqc_epykit/ \
    -o results/combined_multiqc/
```

## Complete Example

Putting it all together:

```python
import epykit as ep

# 1. Check upstream QC
qc = ep.read_nfcore_methylseq_qc(run_dir="/data/nfcore_methylseq/results")
print(qc[["sample_name", "alignment_rate", "bisulfite_conversion_rate"]])

# 2. Ingest
md = ep.read_nfcore_methylseq(
    run_dir="/data/nfcore_methylseq/results",
    treatment_group="tumor",
    control_group="normal",
    treatment_samples=["tumor_1", "tumor_2", "tumor_3"],
    control_samples=["normal_1", "normal_2", "normal_3"],
    assembly="hg38",
)

# 3. Analyse
ep.pp.filter_coverage(md)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)
ep.tl.qc(md)
ep.tl.dmc(md)
ep.tl.dmr(md)
ep.tl.annotate(md, gtf="gencode.v44.annotation.gtf.gz")

# 4. Export
md.report("results/report.html")
ep.report_multiqc(md, output_dir="results/multiqc_epykit/")
md.save("results/analysis")
```
