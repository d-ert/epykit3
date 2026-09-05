# External Integrations

**Analysis Date:** 2026-06-06

## APIs & External Services

**Genomic Data Sources:**
- UCSC Genome Browser (refGene.txt) - Gene feature annotations
- GTF gene annotation files - Customizable from Ensembl, Gencode, NCBI RefSeq
- UCSC CpG Islands BED - CpG island/shore/shelf/open-sea context assignment

## Data Storage

**Parquet Store (Primary):**
- Format: Partitioned Parquet (per-sample, per-chromosome)
- Backend: `pyarrow` ≥ 11.0.0 via `polars`
- Location: `<store_dir>/.cache/<step>/sample=<sample_id>/chrom=<chrom>/part-0.parquet`
- Access pattern: Lazy scans via `pl.scan_parquet()` — never materialize full genome
- S3 support: Native via pyarrow (can read/write to S3-compatible endpoints)

**Input Formats:**

### Bismark (.cov)
- Format: Tab-delimited 6-column: `chrom`, `start`, `end`, `methyl_percent`, `N_meth`, `N_unmeth`
- Coordinate system: 0-based (BED format) — Bismark's `bismark2bedGraph` output
- Client: `src/epykit/io.py:_read_methylation_samplesheet()` + `src/epykit/convert.py`
- Processed via: `ep.read_bismark(samplesheet, treatment_group, control_group, store_dir=...)`
- Notes: 1-based files from `coverage2cytosine` must be pre-shifted by -1

### MethylDackel (.cov)
- Format: Identical 6-column layout to Bismark but with `track type="bedGraph"` header line
- Parsing: `convert.py` skips 1 header row (vs Bismark's 0)
- Client: `src/epykit/io.py:_read_methylation_samplesheet()` + `src/epykit/convert.py`
- Processed via: `ep.read_methyldackel(samplesheet, treatment_group, control_group, store_dir=...)`

### Combined-Strand BED (12-column)
- Format: Strand-collapsed CpG dyad summary: `chrom start end fwd_M fwd_T fwd_% rev_M rev_T rev_% M T methyl_percent`
- Client: `src/epykit/convert.py`
- Processed via: `ep.read_combined_strand_bed(samplesheet, treatment_group, control_group, store_dir=...)`

### nf-core/methylseq
- Format: Bismark-compatible .cov files + optional QC metadata
- Client: `src/epykit/io.py` + `src/epykit/nfcore_qc.py`
- Processed via: `ep.read_nfcore_methylseq(...)` (ingests Bismark files); `ep.read_nfcore_methylseq_qc(...)` (parses FastQC/Samtools output)

### BAM (Read-level)
- Format: Standard SAM/BAM with methylation tags
  - Bismark dialect: Per-base `XM` tag (Z/z for CpG, X/x for CHG, H/h for CHH, . for no call)
  - MethylDackel dialect: SAM standard `MM:Z` and `ML` tags (modified-base positions + likelihoods)
- Client: `src/epykit/bam_io.py:_require_pysam()` (lazy import)
- Optional dependency: `pysam` ≥ 0.22 (Linux/macOS only; `sys_platform != 'win32'`)
- Used by: `src/epykit/asm.py` (allele-specific methylation), `src/epykit/entropy.py` (methylation entropy)
- Output schema: `read_id`, `chrom`, `pos`, `methylation_status` (1/0/-1), `base_qual`, `mapq`, `mate_pair_id`, `strand`, `allele_base`

## File Storage / Exports

**Local Filesystem:**
- Parquet store writes: `<store_dir>/.cache/<step>/...` (partitioned directories)
- DMCStore handle: `src/epykit/_dmc_store.py` (per-chromosome parquet directories + `.epykit_dmc_manifest.json`)
- SmoothStore handle: `src/epykit/_smoothed_store.py` (similar structure for smoothed methylation)

**BedGraph:**
- Format: IGV/UCSC-friendly 4-column: `chrom`, `start`, `end`, `value`
- Generator: `src/epykit/export.py:to_bedgraph()` (no extra deps)
- Use case: Single-sample beta or coverage export

**BigWig:**
- Format: Binary BigWig (IGV/UCSC-friendly)
- Generator: `src/epykit/export.py:to_bigwig()` (lazy imports `pyBigWig`)
- Optional dependency: `pyBigWig` ≥ 0.3.22 (`sys_platform != 'win32'`)
- Friendly error on missing: "Install with: pip install 'epykit[export]'"

**BED Format:**
- DMC results: `src/epykit/export.py:dmcs_to_bed()` (chrom, pos, pos+1, plus DMC metadata)
- DMR results: `src/epykit/export.py:dmrs_to_bed()` (chrom, start, end, plus DMR metadata)
- No extra dependencies

**TSV/CSV:**
- DMC table: `src/epykit/export.py:dmc_to_tsv()` (supports lr+ combined p-value columns: `pvalue_combined`, `qvalue_combined`)
- DMR table: `src/epykit/export.py:dmr_to_tsv()` (full table sorted by chrom/start)
- DVC table: `src/epykit/export.py:dvc_to_tsv()` (differential variability CpG, filtered on `q_variance`)
- QC summary: `src/epykit/export.py:qc_to_tsv()` (per-sample QC summary from `md.obs`)
- Delimiter: Inferred from extension (`.csv` → comma; `.tsv` or other → tab)

**MultiQC Custom Content:**
- Format: `*_mqc.json` files (MultiQC custom-content schema)
- Generator: `src/epykit/multiqc_export.py:report_multiqc()` (no extra deps)
- Output schema: `{id, section_name, plot_type, data}` per metric
- Metrics included: bisulfite conversion rate, mean coverage, sex check, contamination, correlation heatmap
- Use case: Integrate epykit QC into nf-core pipeline MultiQC reports

## Ecosystem Interop (Heavy Deps)

**AnnData / mudata:**
- Format: HDF5-backed AnnData object with `.obs` (samples), `.var` (CpGs by chrom/pos), layers (beta, coverage, N_meth, N_unmeth)
- Generator: `src/epykit/anndata_io.py:to_anndata()` (lazy imports `anndata`, `numpy`)
- Optional dependency: `anndata` ≥ 0.10, `mudata` ≥ 0.2
- Memory strategy: Streamed per-sample, per-chromosome to avoid materializing full genome matrix
- Requirement: `pp.unite()` must run first (shared site set across samples)
- Use case: Export to scanpy, muon, or multi-omics pipelines

**methylKit Format (tabix):**
- Format: Tab-delimited methylKit methylRawDB schema: `chrBase`, `chr`, `base`, `strand`, `coverage`, `freqC`, `freqT`
- Generator: `src/epykit/methylkit_io.py:to_methylkit_tabix()` (gzip + tabix indexing via lazy `pysam` import)
- Optional dependency: `pysam` ≥ 0.22 (Linux/macOS only; Windows skips `.tbi` index silently)
- Use case: Downstream methylKit workflow (cross-platform qPCR validation, etc.)

## Annotation Sources

**GTF (Gene Feature Annotation):**
- File format: GFF3/GTF columnar text
- Parser: `src/epykit/annotate.py:annotate_features()` (custom parser, no external dep)
- Features extracted: `Chromosome`, `Start`, `End`, `Strand`, `Feature`, `gene_id`, `gene_name`
- Cache: Per-process LRU cache (default 2 slots, tunable via `EPYKIT_GTF_CACHE_SIZE` env var or `set_gtf_cache_size()`)
- Priority: Intronic-host-first (single-best gene per site); optional multi-annotation (nearest-TSS HOMER-style + all-overlapping annotatr-style)
- HOMER preset: `ep.HOMER_FEATURES` tuple (promoter, 5UTR, exon, intron, 3UTR, TTS, noncoding) for standard methylation pie charts

**UCSC CpG Islands BED:**
- File format: 3-column BED (chrom, start, end) or optional name/score columns
- Parser: `src/epykit/annotate.py:annotate_cpg_islands()` (custom parser, no external dep)
- Context assigned: island, shore (±2 kb), shelf (±2 kb further), open-sea
- Use case: Enrichment analysis, island vs. open-sea DMC counts

**Reference FASTA:**
- File format: FASTA sequence (gzipped or plain)
- Parser: `pyfaidx` ≥ 0.7
- Use case: Strand inference at CpG positions (C → forward strand, non-C → reverse strand)
- Optional: If `reference_fasta=None`, strand defaults to `*` (unknown)

**UCSC refGene.txt:**
- Legacy gene annotation format (tab-delimited, obsolete; GTF is preferred)
- Client: Not actively used in current codebase (GTF is the canonical input)

## Clustering & Imputation

**UMAP:**
- Optional dependency: `umap-learn` ≥ 0.5 (in `viz` extra)
- Use case: Sample clustering visualization (`src/epykit/pl/clustering.py`)

**KNN Imputation:**
- Client: `src/epykit/impute.py:impute_knn_beta()`, `impute_knn_anndata()`
- Backend: `scikit-learn` ≥ 1.6.1 (KNeighborsRegressor for missing beta values)
- Use case: Recover missing CpG values from sample neighbors

## Distributed Compute (Optional)

**Dask:**
- Optional dependency: `dask[distributed]` ≥ 2024.1 (in `distributed` extra)
- Use case: Chromosome-parallel DMC/DMR on multi-core or cluster
- Temp dir: Inherits `TMPDIR`/`TEMP` via `ep.set_tmp_dir()` so workers use configured staging area

**Ray:**
- Optional dependency: `ray` ≥ 2.9 (in `ray` extra)
- Alternative to Dask for distributed compute

## GPU Acceleration (Heavy, Optional)

**CuPy Backend:**
- Optional dependency: `cupy-cuda12x` ≥ 13.0 (in `gpu` extra, not in `all`)
- Module: `src/epykit/_glm_gpu.py:irls_binomial_batch_gpu()`
- Use case: GPU-accelerated binomial IRLS for GLM (large sample counts, 1000s of CpGs)
- Note: Heavy wheel (~500 MB); CUDA 12 required

**JAX Backend:**
- Optional dependency: `jax[cuda12]` ≥ 0.4.30 (in `gpu_jax` extra, not in `all`)
- Alternative GPU backend to CuPy (mutually exclusive)
- Note: Experimental; CuPy is primary GPU path

## HTML Report Generation

**Report Engine:**
- Framework: Jinja2 ≥ 3.1, Plotly ≥ 5.20 (in `report` extra)
- Generator: `src/epykit/report.py:generate_report(md, output_path)`
- Template files: `src/epykit/templates/` (`.j2` and `.css`)
- Output: Self-contained single HTML file with embedded Plotly charts
- Sections: Sample metadata, preprocessing trail, QC metrics, DMC/DMR volcano plots, annotation pie charts, optional TSS metaplot, PCA, provenance
- Conditional rendering: Unrun sections render "not run yet" stub instead of failing

**Visualization (Static):**
- Matplotlib + Seaborn: `src/epykit/pl/` modules
- Shared theme: `src/epykit/pl/_style.py`

**Visualization (Interactive):**
- Plotly twins: Charts in `report.py` (DMC/DMR scatter, annotation pie, sample correlation heatmap)

## Environmental Configuration

**Required env vars:**
- None (all critical paths have defaults or user-provided arguments)

**Optional env vars:**
- `EPYKIT_GTF_CACHE_SIZE` - LRU cache size for GTF files (default: 2)
- `EPYKIT_NO_AUTO_CSV` - Set to "1", "true", "True" to suppress automatic TSV export on CLI (via `src/epykit/cli.py`)
- `TMPDIR` / `TEMP` / `TMP` - Redirected by `ep.set_tmp_dir()` (critical for Windows)
- `PYTHONIOENCODING=utf-8` - Set in CI to ensure UTF-8 output on Windows

**Secrets:**
- No API keys, credentials, or authentication required
- All inputs are file paths, sample sheets, or local annotations

---

*Integration audit: 2026-06-06*
