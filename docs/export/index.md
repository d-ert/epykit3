# Export Overview

epykit provides a range of export functions that convert your methylation data into
formats compatible with genome browsers, downstream analysis tools, and reporting
pipelines. All export functions are accessible through the top-level `ep` namespace.

## Export Functions

| Function | Target | Description |
|----------|--------|-------------|
| `ep.to_bedgraph()` | BedGraph | Genome browser track |
| `ep.to_bigwig()` | BigWig | Compressed genome browser track |
| `ep.dmcs_to_bed()` | BED | DMC results as BED |
| `ep.dmrs_to_bed()` | BED | DMR results as BED |
| `ep.to_anndata()` | AnnData | Sample x site matrix for scverse |
| `ep.to_mudata()` | MuData | Multi-omics bundle |
| `ep.to_methylkit_tabix()` | Tabix TSV | methylKit-compatible |
| `ep.report_multiqc()` | JSON | MultiQC custom content |
| `md.report()` | HTML | Interactive HTML report |

## Quick Start

```python
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")

# Genome browser tracks
ep.to_bedgraph(md, sample="tumor_1", output="tumor_1.bedgraph")
ep.to_bigwig(md, sample="tumor_1", output="tumor_1.bw")

# Differential methylation results
ep.dmcs_to_bed(md, output="dmcs.bed")
ep.dmrs_to_bed(md, output="dmrs.bed")

# Interoperability
adata = ep.to_anndata(md)
mdata = ep.to_mudata(md)

# methylKit cross-validation
ep.to_methylkit_tabix(md, dir="methylkit_output/")

# Reporting
ep.report_multiqc(md, dir="multiqc_input/")
md.report("report.html")
```

## Optional Dependencies

Some export targets require optional dependency groups:

```bash
pip install 'epykit[export]'    # pyBigWig (BigWig export)
pip install 'epykit[anndata]'   # anndata, mudata (scverse formats)
pip install 'epykit[report]'    # Jinja2, Plotly (HTML report)
```

## Pages in This Section

- [Genome Browser Exports](genome-browsers.md) -- BedGraph, BigWig, BED
- [AnnData / MuData](anndata.md) -- scverse integration
- [methylKit Tabix](methylkit.md) -- methylKit-compatible output
- [MultiQC](multiqc.md) -- MultiQC custom content
- [HTML Report](html-report.md) -- Interactive HTML report
