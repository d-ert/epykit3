# HTML Report

Generate a self-contained interactive HTML report from a MethylData object. The
report bundles all plots and tables into a single file with no external dependencies,
suitable for sharing with collaborators or archiving alongside publications.

## Installation

```bash
pip install 'epykit[report]'
```

This installs `Jinja2` and `Plotly` as optional dependencies.

---

## Usage

The report can be generated through the MethylData method or the top-level function:

```python
# Method on MethylData
md.report("output.html")

# Equivalent top-level function
ep.generate_report(md, "output.html")
```

**Parameters**

| Parameter | Type | Description |
|-----------|------|-------------|
| `md` | `MethylData` | MethylData object |
| `output` | `str` or `Path` | Output HTML file path |

---

## Report Sections

The report auto-populates sections based on what data is available in the MethylData
object. Sections that lack data are omitted silently.

| Section | Requires | Contents |
|---------|----------|----------|
| **Sample Overview** | Always present | Sample table, group counts |
| **Coverage QC** | Coverage data | Per-sample coverage distributions, CpG detection rates |
| **Global Methylation** | Beta values | Beta-value density plots, per-sample summaries |
| **PCA** | `pp.unite()` | PCA scatter colored by sample groups |
| **DMC Results** | `ep.dmc()` | Volcano plot, Manhattan plot, top DMC table |
| **DMR Results** | `ep.dmr()` | DMR summary table, region-level statistics |
| **Annotation** | `ep.annotate()` | Feature distribution (promoter, exon, intron, etc.), CpG island context |

---

## Example

A typical workflow that produces a full report:

```python
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")

# Preprocessing
ep.pp.filter_coverage(md, min_cov=10)
ep.pp.unite(md)

# Analysis
ep.dmc(md, formula="~ treatment", contrast=("treatment", "tumor", "normal"))
ep.dmr(md, method="chain_merge")
ep.annotate(md, gtf="gencode.v44.gtf.gz", cpg_islands="cpg_islands.bed")

# Generate report -- all sections will be populated
md.report("methylation_report.html")
```

---

## Interactive Features

All plots in the report are rendered with Plotly and support:

- **Hover** -- inspect individual data points (CpG coordinates, p-values, sample IDs)
- **Zoom** -- click and drag to zoom into regions of interest
- **Pan** -- shift-drag to pan across the plot
- **Download** -- use the Plotly toolbar to save any plot as PNG or SVG

The report is a single `.html` file with all JavaScript and CSS inlined. It can be
opened in any modern browser without a web server.

---

## CLI

The report can also be generated from the command line:

```bash
epykit report --methylstore methylstore/ --samplesheet samplesheet.csv -o report.html
```

See the [CLI Reference](../cli/index.md) for the full set of options.
