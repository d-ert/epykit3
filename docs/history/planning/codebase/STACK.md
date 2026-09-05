# Technology Stack

**Analysis Date:** 2026-06-06

## Languages

**Primary:**
- Python 3.9+ (3.9, 3.10, 3.11, 3.12) - Core pipeline implementation
  - Tested on: `{ubuntu-latest, windows-latest} × {py3.9, py3.12}` (CI matrix in `.github/workflows/test.yml`)
  - Windows compatibility is load-bearing — several optional extras are Linux/macOS only

## Runtime

**Environment:**
- CPython 3.9-3.12
- **Note:** Uses `from __future__ import annotations` throughout (`src/epykit/`) to enable string-based type hints compatible with runtime 3.9 while type-checking under mypy's enforced 3.10 target

**Package Manager:**
- `uv` (primary, matches CI invocation in `.github/workflows/test.yml`)
- `pip` (supported, works as fallback)
- Lockfile: Not present (library project; `uv.lock` is gitignored)

## Frameworks

**Core Data Processing:**
- `polars` ≥ 0.20.0 - Lazy I/O over partitioned Parquet methylstore (core)
- `pyarrow` ≥ 11.0.0 - Parquet reader/writer backend

**Statistics & Numerics:**
- `numpy` ≥ 1.23 - Array operations, linear algebra
- `scipy` ≥ 1.11 - Statistical distributions (beta, chi-square), special functions
- `numba` ≥ 0.59 - JIT compilation for tight inner loops (methylation likelihood, HMM)
- `statsmodels` ≥ 0.14 - GLM inference (binomial IRLS), Benjamini-Hochberg FDR, Stouffer combining
- `patsy` ≥ 0.5 - Wilkinson formula parsing → design matrices (for GLM covariates)
- `scikit-learn` ≥ 1.6.1 - PCA (for sample QC and correlation viz), imputation

**Genomic Intervals:**
- `bioframe` ≥ 0.7 - Interval overlap operations (GTF gene-feature joins, CpG-island context)

**Utilities:**
- `pyfaidx` ≥ 0.7 - FASTA sequence lookup (strand inference from reference base at CpG position)
- `psutil` ≥ 7.2.2 - Memory monitoring, process utilities

**Visualization (base):**
- `matplotlib` ≥ 3.9.4 - Static plots (DMC/DMR density, coverage, sample correlation heatmaps)
- `seaborn` ≥ 0.13.2 - Statistical visualization, themes

**Testing:**
- `pytest` ≥ 8.4.2 - Test runner
- `pytest-cov` - Coverage reporting

**Build/Dev:**
- `ruff` ≥ 0.15.14 - Linter (pyflakes only, `F` rule set; see `pyproject.toml` [tool.ruff.lint])
- `mypy` ≥ 1.19.1 - Type checker (configured for Python 3.10)
- `setuptools` ≥ 61.0 - Build backend
- `wheel` - Wheel packaging

## Key Dependencies

**Critical:**
- `polars` - Streaming, lazy Parquet access is the core I/O contract; every pipeline step reads/writes Parquet stores
- `numpy`/`scipy` - Statistical engines (DMC likelihood-ratio, Welch t, Fisher exact, GLM binomial IRLS)
- `numba` - JIT for performance-critical paths (HMM Viterbi, large-scale likelihood calculations)
- `statsmodels` - FDR correction (Benjamini-Hochberg, two-stage variants), GLM wald/F contrasts
- `bioframe` - Genomic interval operations (GTF gene overlap, CpG island assignment)

**Infrastructure:**
- `pyarrow` - Parquet columnar I/O backend
- `patsy` - GLM design matrix construction from formula strings
- `pyfaidx` - Reference genome lookup for strand inference

## Optional Dependencies (Gated by Platform)

**report** extra:
- `jinja2` ≥ 3.1 - Templating for HTML report generation
- `plotly` ≥ 5.20 - Interactive Plotly charts in HTML reports

**anndata** extra:
- `anndata` ≥ 0.10 - AnnData format export for ecosystem interop
- `mudata` ≥ 0.2 - Companion for multi-omics interop (used with AnnData)

**viz** extra:
- `umap-learn` ≥ 0.5 - UMAP reduction (for sample clustering viz)
- `scipy` ≥ 1.11 - Already in base, but pinned for UMAP

**export** extra (Linux/macOS only):
- `pyBigWig` ≥ 0.3.22 - BigWig format export; **Windows has no PyPI wheel** (conditional in pyproject.toml: `sys_platform != 'win32'`)

**methylkit** extra (Linux/macOS only):
- `pysam` ≥ 0.22 - Tabix bgzip/indexing for methylKit-format export; **Windows has no PyPI wheel** (conditional: `sys_platform != 'win32'`)

**bam** extra (Linux/macOS only):
- `pysam` ≥ 0.22 - BAM parsing for read-level methylation analyses (ASM, entropy); **Windows incompatible** (conditional)

**distributed** extra:
- `dask[distributed]` ≥ 2024.1 - Distributed scheduler for large-scale chromosome-parallel DMC/DMR

**ray** extra:
- `ray` ≥ 2.9 - Alternative distributed compute backend

**gpu** extra (heavy, not in 'all'):
- `cupy-cuda12x` ≥ 13.0 - GPU-accelerated IRLS via CuPy; CUDA 12 wheel

**gpu_jax** extra (heavy, not in 'all'):
- `jax[cuda12]` ≥ 0.4.30 - GPU-accelerated IRLS via JAX; CUDA 12; mutually exclusive with gpu

**qc** extra:
- `diptest` - Hartigan's dip test for unimodal fallback in sex_check; improves robustness on low-coverage samples

**all** extra:
- Bundles: `report`, `anndata`, `viz`, `export` (non-Windows), `distributed`, `qc`
- **Excludes:** `methylkit`, `bam` (Windows), `gpu`, `gpu_jax` (heavy CUDA)

## Configuration

**Environment:**
- Python version: Specified via `requires-python = ">=3.9"` in `pyproject.toml`
- Mypy target: Python 3.10 (configured in `[tool.mypy] python_version = "3.10"`)
- Ruff target: Python 3.9 (configured in `[tool.ruff] target-version = "py39"`)

**Build:**
- Build backend: `setuptools.build_meta` (PEP 517/518 compliant)
- Package data: `src/` layout with `pyproject.toml` `[tool.setuptools.packages.find] where = ["src"]`
- Templates: Jinja2 templates packaged as `src/epykit/templates/*.j2` and `*.css` (declared in `[tool.setuptools.package-data]`)
- Type hints: `src/epykit/py.typed` marker file (PEP 561 compliance)

**Linting:**
- **Ruff:** Baseline `F` (pyflakes) only; `line-length = 100`; per-file ignores for test files and re-export modules
- **Mypy:** `follow_imports = "silent"`, `ignore_missing_imports = true`, `no_strict_optional = true` (permissive mode for third-party stubs)

**Testing:**
- Framework: pytest ≥ 8.4.2
- Markers: `slow` (opt-in, for tests > ~5s); registered in `[tool.pytest.ini_options] markers`
- Mode: `--strict-markers` enabled (unregistered marker = test failure)
- Filter warnings: Custom deprecation handling via `filterwarnings` — allows `DeprecationWarning:epykit` and `UserWarning:epykit` through for assertion, silences third-party noise

## Platform Requirements

**Development:**
- OS: Windows, macOS, Linux (tested on `windows-latest`, `ubuntu-latest`)
- Python: 3.9, 3.12 (CI matrix)
- uv: Automatic installation via GitHub Actions `astral-sh/setup-uv@v4`

**Production:**
- Parquet storage: Local filesystem or S3-compatible (via pyarrow's native S3 support)
- Temporary directory: Configurable via `ep.set_tmp_dir(path)` (redirects `TMPDIR`/`TEMP`/`TMS` env vars); critical for Windows where default `C:\TEMP` is often too small for whole-genome staging
- RAM: Streaming design keeps peak memory O(largest chromosome) for DMC/DMR; no whole-genome frame is ever materialized

## Version Pinning Strategy

- `numpy`, `scipy`, `polars`, `statsmodels`, `bioframe`: Locked to ≥ versions that support the known APIs; no upper bounds (forward-compatible)
- `numba` ≥ 0.59: Critical for HMM performance; pinned due to historical instability in earlier versions
- `ruff`, `mypy`: Dev tools; pinned to minimum working versions; upgraded opportunistically

---

*Stack analysis: 2026-06-06*
