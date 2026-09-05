# Codebase Structure

**Analysis Date:** 2026-06-06

## Directory Layout

```
epykit3/
├── src/epykit/              # Main library package (installable)
│   ├── __init__.py          # Public API re-exports
│   ├── methyldata.py        # MethylData central dataclass
│   ├── io.py                # Ingestion: Bismark, MethylDackel, BED, nf-core
│   ├── convert.py           # Bismark .cov → Parquet conversion
│   ├── pp.py                # Preprocessing API (filter, normalize, unite, smooth)
│   ├── tl.py                # High-level orchestrators (qc, dmc, dmr, annotate)
│   ├── filter.py            # Coverage filtering & normalization engines
│   ├── dmc.py               # DMC engines (lr, welch_t, fisher, glm)
│   ├── _dmc_store.py        # Per-chromosome DMC result store handle
│   ├── _glm.py              # Wilkinson formula → design matrix, IRLS binomial GLM
│   ├── _glm_gpu.py          # Optional CuPy/JAX GPU backend (gated by extras)
│   ├── dmr.py               # DMR engines (chain_merge, tile, sliding_window)
│   ├── _hmm.py              # HMM segmentation for DMR
│   ├── dmr_hmm.py           # DMR via HMM wrapper
│   ├── dmr_segment.py       # DSS-style segment calling wrapper
│   ├── dvc.py               # Differentially variable CpG (iEVORA-style)
│   ├── annotate.py          # Gene feature & CpG island annotation
│   ├── qc.py                # QC metrics (conversion, coverage, methylation, sex, contamination, power)
│   ├── pl/                  # Plotting (matplotlib, optional Plotly)
│   │   ├── __init__.py      # Public plotter re-exports
│   │   ├── _style.py        # Shared matplotlib theme
│   │   ├── _plotly.py       # Plotly helper utilities
│   │   ├── _utils.py        # Plotting utilities (colors, formats, etc.)
│   │   ├── _compute.py      # Compute-intensive helper (e.g. karyogram binning)
│   │   ├── qc.py            # QC plots (coverage histogram, methylation heatmap, M-bias)
│   │   ├── differential.py  # DMC plots (volcano, MA, Manhattan)
│   │   ├── genomic.py       # Genomic context plots (CpG island pie, karyogram)
│   │   ├── clustering.py    # PCA plot
│   │   ├── correlation.py   # Sample correlation heatmap
│   │   ├── metaplot.py      # TSS metaplot, gene body metaplot
│   │   ├── embedding.py     # UMAP embedding
│   │   ├── dmr_boxplot.py   # Per-DMR methylation boxplots
│   │   ├── dmr_summary.py   # DMR violin, heatmap
│   │   ├── overlap.py       # DMR overlap plots
│   │   ├── annotation.py    # Annotation-stratified plots (counts, numericals, etc.)
│   │   ├── dashboard.py     # Multi-panel QC dashboard
│   │   └── composer.py      # Grid layout for multi-plot figures
│   ├── report.py            # Interactive HTML report generator (Jinja2 + Plotly)
│   ├── templates/           # Jinja2 HTML templates for report
│   │   ├── base.html        # Main report template
│   │   └── ...              # Section templates
│   ├── export.py            # BED, BedGraph, BigWig export
│   ├── anndata_io.py        # AnnData/MuData export (requires pp.unite first)
│   ├── methylkit_io.py      # methylKit tabix export
│   ├── multiqc_export.py    # MultiQC-compatible export
│   ├── bam_io.py            # BAM file read (Linux/macOS only, optional)
│   ├── nfcore_qc.py         # nf-core methylseq QC module reader
│   ├── query.py             # Query utilities for MethylData
│   ├── impute.py            # KNN imputation (optional, biolearn)
│   ├── clocks.py            # Epigenetic age clocks
│   ├── asm.py               # Allele-specific methylation
│   ├── hmr.py               # Hypomethylated regions
│   ├── pmd.py               # Partially methylated domains
│   ├── entropy.py           # Entropy-based metrics
│   ├── dvc.py               # (duplicate entry for DVC)
│   ├── _cache.py            # Caching utilities (manifest, JSON I/O, store counting)
│   ├── _compute.py          # Low-level compute (Welford accumulators, scoring)
│   ├── _config.py           # Configuration (TMPDIR redirection)
│   ├── _smoothed_store.py   # Persistent handle for smoothed results
│   ├── _style.py            # Matplotlib style (shared with pl/_style.py)
│   └── cli.py               # argparse CLI entry point
│
├── tests/                   # Test suite
│   ├── conftest.py          # Pytest fixtures (test_data_dir, temp_cov, sample fixtures, etc.)
│   ├── fixtures/            # Test data directory (small .cov files, reference GTF, samplesheets)
│   ├── test_accuracy.py     # Accuracy validation against known results
│   ├── test_annotate_multi.py # Multi-annotation test
│   ├── test_api.py          # High-level API tests
│   ├── test_asm.py          # Allele-specific methylation tests
│   ├── test_bam_io.py       # BAM input tests
│   ├── test_bsmooth.py      # B-spline smoothing tests
│   ├── test_calibration.py  # Calibration & empirical FDR tests
│   ├── test_cli.py          # CLI integration tests
│   ├── test_combined_strand_bed.py # Combined-strand BED format tests
│   ├── test_compute_backends.py # Compute engine tests (Numba, etc.)
│   ├── test_dmc_empirical_fdr.py # DMC empirical FDR tests
│   ├── test_dmc_fisher.py   # Fisher exact test
│   ├── test_dmc_lr.py       # Likelihood-ratio test
│   ├── test_dmc_multigroup.py # Multi-group DMC & contrasts
│   ├── test_dmc_multitest.py # Multiple testing correction
│   ├── test_dmc_smooth_dispersion.py # Count smoothing & dispersion
│   ├── test_dmc_streaming_store.py # DMCStore streaming tests
│   ├── test_dmr_chain_merge.py # Chain-merge DMR tests
│   ├── test_dmr_empirical_fdr.py # DMR empirical FDR
│   ├── test_dmr_hmm.py      # HMM segmentation tests
│   ├── test_dmr_presets_and_diagnose.py # DMR preset configs
│   ├── test_dmr_segment.py  # DSS-style segment tests
│   ├── test_dmr_tile_merge.py # Tile-based DMR tests
│   ├── test_dvc.py          # Differentially variable CpG tests
│   └── ... (40+ more tests)
│
├── benchmark/               # Benchmark suite (reproducibility paper)
│   ├── data/                # Frozen parquet results committed to git
│   │   ├── MANIFEST.txt     # Data inventory & checksums
│   │   ├── seeds.json       # Piao simulator RNG seeds for reproducibility
│   │   ├── study1/          # Piao 2021 simulator results
│   │   ├── study2/          # Real methylKit comparison data
│   │   ├── study3/          # GSE263850 real data results
│   │   └── ...              # Other benchmark datasets
│   ├── paper_data/          # TSV mirrors for Excel/human inspection
│   │   ├── 01_headline_piao/
│   │   ├── 02_methods_comparison/
│   │   ├── 03_study1_sims/
│   │   ├── 04_study2_methylkit/
│   │   ├── 05_study3_gse263850/
│   │   └── 06_methodology/
│   ├── scripts/             # Benchmark orchestration
│   │   ├── regen_all.py     # Master regeneration script
│   │   ├── claims.yaml      # Claim verification specs
│   │   ├── run_epykit.py    # epykit benchmark harness
│   │   ├── run_methylkit.py # methylKit comparison
│   │   ├── run_dss.py       # DSS comparison
│   │   └── ...              # Other tool runners
│   ├── paper/               # Paper & report
│   │   ├── paper.md         # Manuscript
│   │   ├── report/
│   │   │   └── REPORT.md    # Canonical TPR/FPR/F1 results table
│   │   └── figures/         # Generated paper figures
│   └── README.md            # Benchmark setup & bootstrap instructions
│
├── docs/                    # mkdocs documentation
│   ├── index.md             # Landing page
│   ├── architecture.md      # System design documentation
│   ├── guide_*              # User guides (getting started, analysis workflow)
│   └── ...                  # Other docs
│
├── CLAUDE.md                # Project-specific Claude Code guidance
├── CHANGELOG.md             # Version history
├── README.md                # Top-level readme
├── HANDOFF.md               # Deployment & operations guide
├── pyproject.toml           # uv/pip package config (Python 3.9+, dependencies, extras)
├── mkdocs.yml               # Documentation build config
└── .planning/
    └── codebase/            # GSD codebase mapping documents (this file, ARCHITECTURE.md)
```

## Directory Purposes

**src/epykit/:**
- Purpose: Main library package installable via `pip install epykit` or `uv sync`
- Contains: All modules (API, engines, plotting, export, CLI)
- Key files: methyldata.py (central dataclass), tl.py (orchestrators), dmc.py (DMC engines), dmr.py (DMR engines), cli.py (CLI entry)

**tests/:**
- Purpose: pytest test suite
- Contains: ~50 unit, integration, and accuracy tests
- Key files: conftest.py (fixtures), test_dmc_*.py (DMC engine validation), test_dmr_*.py (DMR engine validation), test_api.py (end-to-end)

**benchmark/:**
- Purpose: Reproducibility paper benchmark suite
- Contains: Frozen parquet sources (data/), TSV human-readable mirrors (paper_data/), orchestration scripts (scripts/), paper and REPORT.md (paper/)
- Key files: data/MANIFEST.txt (inventory), paper/REPORT.md (TPR/FPR results), scripts/regen_all.py (master regenerator)

**docs/:**
- Purpose: mkdocs documentation site
- Contains: Architecture docs, user guides, API reference
- Deployed to: GitHub Pages

## Key File Locations

**Entry Points:**
- `src/epykit/__init__.py`: Public API re-exports (read_bismark, MethylData, ep.pp.*, ep.tl.*, ep.pl.*, etc.)
- `src/epykit/cli.py`: argparse CLI entry point (`epykit convert | filter | dmc | dmr | ...`)

**Configuration:**
- `pyproject.toml`: Package metadata, dependencies, extras (dev, all, report, gpu, etc.), test config (pytest, mypy, ruff)
- `src/epykit/_config.py`: TMPDIR redirection (load-bearing for Windows whole-genome temp staging)

**Core Logic:**
- `src/epykit/methyldata.py`: MethylData dataclass (obs, store, varm, uns, state properties)
- `src/epykit/tl.py`: High-level orchestrators (qc, dmc, dmr, annotate) that wire engines and mutate MethylData
- `src/epykit/dmc.py`: Per-CpG test engines (lr, welch_t, fisher, glm, dispersion shrinkage, separation fallback, neighbour combining)
- `src/epykit/dmr.py`: DMR callers (chain_merge, tile_based, sliding_window, HMM via _hmm.py)
- `src/epykit/_dmc_store.py`: Per-chromosome persistent parquet store handle for streaming memory efficiency

**Testing:**
- `tests/conftest.py`: pytest fixtures (test_data_dir, temp_cov, md_fixture, sample_samplesheet, etc.)
- `tests/test_dmc_lr.py`: Likelihood-ratio engine validation
- `tests/test_dmc_multigroup.py`: Multi-group DMC & contrasts
- `tests/test_dmr_chain_merge.py`: chain_merge DMR validation
- `tests/test_api.py`: End-to-end workflow tests

## Naming Conventions

**Files:**
- Module files: `snake_case.py` (e.g., `methyldata.py`, `_dmc_store.py`)
- Test files: `test_<feature>.py` (e.g., `test_dmc_lr.py`)
- Subpackages: `lowercase_dir/` with `__init__.py` (e.g., `pl/`, `templates/`)

**Directories:**
- Parquet caches: `.cache/<step>/` (e.g., `.cache/raw/`, `.cache/filtered/`, `.cache/dmc/lr/`)
- Private modules: Leading underscore (e.g., `_dmc_store.py`, `_glm.py`, `_compute.py`)
- Test fixtures: `fixtures/` subdirectory with descriptive names (e.g., `sample_bismark_cov.txt`)

**Functions:**
- Public API: `snake_case` (e.g., `filter_coverage()`, `read_bismark()`, `annotate_features()`)
- Private helpers: Leading underscore (e.g., `_auto_test()`, `_chrom_filename()`, `_dmc_input_signature()`)
- Tests: `test_<description>` (e.g., `test_dmc_lr_with_two_groups()`)

**Classes:**
- PascalCase (e.g., `MethylData`, `DMCStore`)
- Dataclasses frozen when immutable (e.g., `DMCStore`)

**Types:**
- Polars dtypes: `pl.Utf8`, `pl.Int32`, `pl.Float64` (used in schema dicts)

## Where to Add New Code

**New DMC Engine:**
- Primary code: `src/epykit/dmc.py` (add function like `_process_chrom_<test_name>()`)
- Register in dispatcher: Add elif branch in `dmc.py:process_chromosomes_dmc()` kwargs dispatch (line ~2000)
- CLI flag: Register in `cli.py:_cmd_dmc()` via argparse choices for --test
- Tests: Create `tests/test_dmc_<test_name>.py` with fixtures from `conftest.py`

**New DMR Engine:**
- Primary code: `src/epykit/dmr.py` or separate module if complex (e.g., `dmr_hmm.py` for HMM)
- Register in dispatcher: Add to `DMR_PRESETS` dict in `dmr.py:~line 50` and `tl.dmr()` method dispatch (~line 500)
- CLI flag: Register in `cli.py` argparse --method choices
- Tests: Create `tests/test_dmr_<method>.py`

**New Preprocessing Step:**
- API wrapper: Add function to `src/epykit/pp.py` (follow pattern of `filter_coverage()`)
- Implementation: If complex, add to `src/epykit/filter.py` or create new module
- State tracking: Call `_append_store_history(md, step, path, n_sites)` and update `md.uns` (never add new boolean flags)
- Tests: Add to `tests/test_api.py` or new `tests/test_pp_<step>.py`

**New Plotting Function:**
- Public plotter: Create in `src/epykit/pl/<category>.py` (e.g., `pl/differential.py`)
- Re-export: Add to `src/epykit/pl/__init__.py` __all__
- Tests: Add to `tests/test_plotting.py` or new `tests/test_pl_<category>.py`
- Style: Inherit theme from `src/epykit/_style.py` via import

**New Export Format:**
- Module: Create `src/epykit/<format>_io.py` (e.g., `methylkit_io.py`)
- Lazy import: Import heavy deps (scanpy, pyarrow) inside functions, not at module level
- Method on MethylData: Add to `src/epykit/methyldata.py:to_<format>()` that calls the module function
- Tests: Create `tests/test_<format>_io.py`

**New QC Metric:**
- Function: Add to `src/epykit/qc.py`
- Populate in obs: Call from `tl.qc()` orchestrator, add columns to `md.obs`
- Plotting: Add plotter in `src/epykit/pl/qc.py`, re-export in `pl/__init__.py`
- Tests: Add to `tests/test_qc.py` or new test file

**New CLI Subcommand:**
- Handler: Create `_cmd_<name>(args: argparse.Namespace)` in `src/epykit/cli.py`
- Argparse setup: Add `subparsers.add_parser('<name>', help=...)` and register arguments in `_add_<name>_parser()`
- Tests: Add to `tests/test_cli.py`

**New Annotation Context:**
- Function: Add to `src/epykit/annotate.py` (e.g., `annotate_<context>()`)
- Populate columns: Add `<context>` and `<context>_distance` (if applicable) to DMC/DMR tables
- Caching: If it requires file I/O (GTF, BED), add LRU caching in module-level OrderedDict
- Tests: Add to `tests/test_annotate_*.py`

## Special Directories

**`.cache/` (under methylstore root):**
- Purpose: Persistent caching of intermediate Parquet stores
- Generated: Yes (by filter_coverage, normalize_coverage, DMC, etc.)
- Committed: No (.gitignore whitelists only benchmark/data frozen sources)
- Layout: `.cache/raw/`, `.cache/filtered/`, `.cache/normalized/`, `.cache/dmc/<test>/`, `.cache/dmr/<method>/`

**`templates/` (in src/epykit/):**
- Purpose: Jinja2 HTML templates for report generation
- Generated: No (checked in)
- Committed: Yes
- Layout: Base template (base.html) + section templates for each report section

**`fixtures/` (in tests/):**
- Purpose: Test data (small .cov files, reference GTF, samplesheets)
- Generated: No (committed to git)
- Committed: Yes
- Layout: `sample_<format>.txt`, `test_<genome>.gtf`, `samplesheet_<scenario>.csv`

**`benchmark/data/` (frozen):**
- Purpose: Canonical parquet sources for paper reproducibility
- Generated: No (hand-curated, frozen)
- Committed: Yes (.gitignore whitelist)
- Layout: Per-study subdirs (study1, study2, study3); manifests + .parquet files; read-only

**`benchmark/paper_data/` (TSV mirrors):**
- Purpose: Human-readable TSV mirrors of benchmark/data/ for Excel inspection
- Generated: Yes (by scripts/regen_all.py)
- Committed: Yes (for paper citations)
- Layout: By paper section (01_headline_piao, 02_methods_comparison, etc.)

## Module Dependencies

**Minimal (always present):**
- polars (lazy DataFrame I/O)
- numpy, scipy (stats)
- pandas (annotation overlap via bioframe)

**Core analysis (gated by feature, not extras):**
- statsmodels (GLM, FDR methods)
- patsy (Wilkinson formula parsing)
- bioframe (interval overlap)
- numba (JIT compilation for hot paths)

**Optional extras (lazy-imported):**
- `[report]`: jinja2, plotly (HTML report generation)
- `[plot]`: matplotlib (static plotting, implicitly required for API but optional for CLI)
- `[gpu]`: cupy, jax (GPU acceleration in _glm_gpu.py)
- `[all]`: All of above + pyBigWig, pysam (for bam_io, bigwig export — Linux/macOS only)

**Development extras (pyproject.toml):**
- `[dev]`: pytest, mypy, ruff, pre-commit, sphinx (testing & linting)

---

*Structure analysis: 2026-06-06*
