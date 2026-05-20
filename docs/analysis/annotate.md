# Annotation

`ep.tl.annotate(md)` maps DMC and DMR results to gene features and CpG-island
context. It annotates all `dmc_*` tables in `md.varm` and the DMR table in
`md.uns["dmr"]` in a single call.

## Basic Usage

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)
ep.tl.dmc(md)
ep.tl.dmr(md)

# Annotate with GENCODE GTF and CpG islands
ep.tl.annotate(
    md,
    gtf="gencode.v44.annotation.gtf.gz",
    cpg_islands="cpg_islands.bed",
)

# Annotated DMC results
ann_dmc = md.varm["dmc_lr_annotated"]
print(ann_dmc.select(["chrom", "pos", "meth_diff", "qvalue",
                       "gene_feature", "nearest_tss_gene"]))
```

## Gene Feature Annotation

Provide a GTF file (GENCODE or Ensembl) or a UCSC `refGene.txt` file to assign
each CpG/DMR to one or more gene features.

### From GTF (GENCODE / Ensembl)

```python
ep.tl.annotate(md, gtf="gencode.v44.annotation.gtf.gz")
```

### From UCSC refGene.txt

UCSC refGene.txt is HOMER's default gene catalog. It provides high gene recall
for methylation work because it is curated and protein-coding-biased.

```python
ep.tl.annotate(md, refgene="refGene.txt")
```

!!! note "Provide one source, not both"
    Pass either `gtf=` or `refgene=`, not both. Providing both raises a
    `ValueError`.

### Gene Features

Each CpG is assigned to one of these categories based on its overlap with
gene model intervals:

| Feature | Definition |
|---------|------------|
| `promoter` | Within `promoter_upstream_bp` upstream to `promoter_downstream_bp` downstream of TSS (default: -2000 to +200) |
| `5UTR` | 5' untranslated region |
| `exon` | Exonic region |
| `intron` | Intronic region |
| `3UTR` | 3' untranslated region |
| `intergenic` | Not overlapping any gene model |

### Gene Type Filter

Restrict the gene catalog to specific biotypes before annotation:

```python
# Keep only protein-coding genes (drop lincRNAs, pseudogenes, etc.)
ep.tl.annotate(md, gtf="gencode.v44.gtf.gz", gene_type_filter="protein_coding")
```

This filters the gene catalog before building overlap intervals and the
nearest-TSS index, so all downstream columns reflect only protein-coding genes.

## CpG-Island Context

Provide a CpG-island BED file to classify each CpG by its proximity to
CpG islands:

```python
ep.tl.annotate(md, gtf="gencode.v44.gtf.gz", cpg_islands="cpg_islands.bed")
```

| Context | Definition |
|---------|------------|
| `island` | Inside a CpG island |
| `shore` | Within 2 kb flanking an island |
| `shelf` | Within 2 kb beyond the shore (2--4 kb from island) |
| `open_sea` | More than 4 kb from any CpG island |

## Multi-Annotation Mode

By default (`multi_annotation=True`), epykit adds annotatr-style columns that
capture multiple overlapping annotations per CpG:

| Column | Description |
|--------|-------------|
| `nearest_tss_gene` | Gene name of the nearest TSS |
| `nearest_tss_distance` | Distance in bp to the nearest TSS (negative = upstream) |
| `all_overlapping_genes` | Comma-separated list of all genes overlapping this position |
| `all_overlapping_features` | Comma-separated list of all features overlapping this position |

This captures the one-to-many relationship between genomic positions and gene
annotations. A CpG in a promoter region may simultaneously overlap an intron
of a neighbouring gene.

```python
# Disable multi-annotation for simpler output
ep.tl.annotate(md, gtf="gencode.v44.gtf.gz", multi_annotation=False)
```

## Significant-Only Annotation

By default (`significant_only=True`), only CpGs with `qvalue < alpha` (default
0.05) are annotated. This avoids out-of-memory issues on whole-genome datasets
where annotating all 22+ million CpGs would be prohibitive.

```python
# Annotate all CpGs regardless of significance (not recommended for WG)
ep.tl.annotate(md, gtf="gencode.v44.gtf.gz", significant_only=False)

# Custom significance threshold
ep.tl.annotate(md, gtf="gencode.v44.gtf.gz", alpha=0.01)
```

!!! warning "Memory on whole-genome data"
    Setting `significant_only=False` on a whole-genome dataset loads and
    annotates every tested CpG. This can require substantial memory. Use
    only on targeted panels or subsets.

## Where Results Are Stored

- **DMC annotations** are stored as `md.varm["dmc_<test>_annotated"]`. The
  original DMC table at `md.varm["dmc_<test>"]` is preserved unchanged.
- **DMR annotations** are applied in-place to `md.uns["dmr"]`.
- **Annotation metadata** is stored at `md.uns["annotation"]`.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `md` | MethylData | required | Analysis object |
| `gtf` | str | None | Path to GENCODE / Ensembl GTF file |
| `refgene` | str | None | Path to UCSC refGene.txt file |
| `cpg_islands` | str | None | Path to CpG-island BED file |
| `significant_only` | bool | True | Annotate only significant DMCs |
| `alpha` | float | 0.05 | Significance threshold for `significant_only` |
| `promoter_upstream_bp` | int | 2000 | Promoter region upstream of TSS |
| `promoter_downstream_bp` | int | 200 | Promoter region downstream of TSS |
| `multi_annotation` | bool | True | Add annotatr-style one-to-many columns |
| `gene_type_filter` | str/list | None | Restrict to specific gene biotypes |
| `clear_gtf_cache` | bool | True | Clear GTF cache after annotation |

## Full Example

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)
ep.tl.dmc(md)
ep.tl.dmr(md)

ep.tl.annotate(
    md,
    gtf="gencode.v44.annotation.gtf.gz",
    cpg_islands="cpg_islands.bed",
    gene_type_filter="protein_coding",
    significant_only=True,
    alpha=0.05,
)

# Inspect annotated DMCs
ann = md.varm["dmc_lr_annotated"]
print(ann.group_by("gene_feature").len().sort("len", descending=True))

# Inspect annotated DMRs
dmrs = md.uns["dmr"]
print(dmrs.select(["chrom", "start", "end", "mean_meth_diff",
                    "nearest_tss_gene", "nearest_tss_distance"]))
```
