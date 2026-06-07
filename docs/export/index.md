# Export Overview

epykit provides a range of export functions that convert your methylation data into
formats compatible with genome browsers, downstream analysis tools, and reporting
pipelines. All export functions are accessible through the top-level `ep` namespace.

## Export Functions

| Function | Target | Description |
|----------|--------|-------------|
| `md.export_tables()` | TSV / CSV | **All result tables** (DMC / DMR / DVC / QC) in one call |
| `ep.export.dmc_to_tsv()` | TSV / CSV | DMC table (significant or full) |
| `ep.export.dmr_to_tsv()` | TSV / CSV | DMR table |
| `ep.export.dvc_to_tsv()` | TSV / CSV | DVC table (significant or full) |
| `ep.export.qc_to_tsv()` | TSV / CSV | Per-sample QC summary |
| `ep.to_bedgraph()` | BedGraph | Genome browser track |
| `ep.to_bigwig()` | BigWig | Compressed genome browser track |
| `ep.dmcs_to_bed()` | BED | DMC results as BED |
| `ep.dmrs_to_bed()` | BED | DMR results as BED |
| `ep.to_anndata()` | AnnData | Sample x site matrix for scverse |
| `ep.to_mudata()` | MuData | Multi-omics bundle |
| `ep.to_methylkit_tabix()` | Tabix TSV | methylKit-compatible |
| `ep.report_multiqc()` | JSON | MultiQC custom content |
| `md.report()` | HTML | Interactive HTML report |

The main analyses (`tl.dmc` / `tl.dmr` / `tl.annotate`) write a human-readable
TSV to `<analysis_root>/results/` **by default** (`tsv=False` to opt out); the
others take an opt-in `tsv=` path — see [Tabular Exports](tables.md).

## Quick Start

```python
import epykit as ep

md = ep.read_methyl("samplesheet.csv", store="methylstore/")

# Result tables (DMC / DMR / DVC / QC) as TSV -- all at once...
md.export_tables("results/tables")
# ...or one at a time with full control
ep.export.dmc_to_tsv(md, "results/dmc.tsv", alpha=0.05)

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

- [Tabular Exports](tables.md) -- DMC / DMR / DVC / QC result tables as TSV / CSV
- [Genome Browser Exports](genome-browsers.md) -- BedGraph, BigWig, BED
- [AnnData / MuData](anndata.md) -- scverse integration
- [methylKit Tabix](methylkit.md) -- methylKit-compatible output
- [MultiQC](multiqc.md) -- MultiQC custom content
- [HTML Report](html-report.md) -- Interactive HTML report
