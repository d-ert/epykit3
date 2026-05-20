# MultiQC

Generate MultiQC-compatible custom content files from epykit QC data, and parse
QC output from nf-core/methylseq pipeline runs.

---

## ep.report_multiqc

Emit `*_mqc.json` custom-content files that MultiQC picks up automatically when
placed in its search directory.

```python
ep.report_multiqc(md, dir)
```

**Parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `md` | `MethylData` | MethylData object with QC data available |
| `dir` | `str` or `Path` | Output directory for the JSON files |

**Output**

The function writes one or more JSON files into the specified directory:

```
multiqc_input/
  epykit_coverage_mqc.json
  epykit_methylation_mqc.json
  epykit_dmc_summary_mqc.json    # if DMC results are present
  epykit_dmr_summary_mqc.json    # if DMR results are present
```

Each file follows the MultiQC custom content specification, including plot
configuration, section headers, and data tables.

**Example**

```python
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")
ep.pp.filter_coverage(md, min_cov=10)
ep.dmc(md, formula="~ treatment", contrast=("treatment", "tumor", "normal"))

ep.report_multiqc(md, dir="multiqc_input/")
```

---

## ep.read_nfcore_methylseq_qc

Parse QC metrics from an nf-core/methylseq pipeline run. Supports output from both
the Bismark and bwa-meth/Qualimap alignment paths.

```python
ep.read_nfcore_methylseq_qc(samplesheet, run_dir)
```

**Parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `samplesheet` | `str` or `Path` | Path to the nf-core samplesheet CSV |
| `run_dir` | `str` or `Path` | Path to the nf-core/methylseq results directory |

**Returns**

A dictionary of QC DataFrames keyed by metric type (e.g., `"alignment"`,
`"deduplication"`, `"mbias"`), which can be attached to a MethylData object or
inspected directly.

**Example**

```python
import epykit as ep

qc = ep.read_nfcore_methylseq_qc(
    samplesheet="samplesheet.csv",
    run_dir="/data/methylseq_results/"
)

# Inspect alignment rates
print(qc["alignment"])

# Attach to MethylData and include in reports
md = ep.read_methyl("samplesheet.csv", store="methylstore/")
md.qc = qc
ep.report_multiqc(md, dir="multiqc_input/")
```

---

## Integration with MultiQC

To include epykit sections in a MultiQC report, drop the generated JSON files into
the directory that MultiQC searches:

```bash
# Generate the custom content files
python -c "
import epykit as ep
md = ep.read_methyl('samplesheet.csv', store='methylstore/')
ep.report_multiqc(md, dir='multiqc_input/')
"

# Run MultiQC with the epykit files in the search path
multiqc multiqc_input/ -o multiqc_report/
```

If you already have a MultiQC run directory with other tool outputs (FastQC, Bismark,
etc.), you can write the epykit JSON files directly into that directory so they appear
alongside the other sections:

```python
ep.report_multiqc(md, dir="/data/methylseq_results/")
```

The custom content sections will appear in the MultiQC report under an "epykit"
heading, with coverage distributions, methylation summaries, and differential
methylation result tables.
