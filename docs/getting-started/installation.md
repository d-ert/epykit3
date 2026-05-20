# Installation

## Requirements

- Python >= 3.9
- A C compiler is **not** required; all heavy lifting uses pre-built wheels (polars, numba, pyarrow).

## Core dependencies

The base install pulls in:

| Package | Purpose |
|---------|---------|
| polars | DataFrame engine (lazy scans, predicate pushdown) |
| pyarrow | Parquet I/O and row-group statistics |
| numpy | Numeric arrays |
| scipy | Statistical distributions |
| numba | JIT-compiled inner loops |
| bioframe | Genomic interval operations |
| pyfaidx | FASTA random access |
| statsmodels | GLM fitting, multiple-testing correction |
| patsy | Design-matrix formula language |
| psutil | Memory-aware scheduling |
| scikit-learn | KNN imputation, PCA |
| matplotlib | Static plotting backend |
| seaborn | Statistical plot themes |

## Install from source

Editable install (recommended for development):

```bash
pip install -e .
```

Or with [uv](https://github.com/astral-sh/uv):

```bash
uv pip install -e .
```

## Optional extras

Install extras with bracket syntax, e.g. `pip install -e ".[report,anndata]"`.

| Extra | Packages | Purpose |
|-------|----------|---------|
| `dev` | pytest, pytest-cov | Testing |
| `report` | jinja2, plotly | Interactive HTML reports |
| `export` | pyBigWig | BigWig export |
| `anndata` | anndata, mudata | AnnData / MuData interop |
| `viz` | umap-learn, scipy | UMAP embedding |
| `methylkit` | pysam | methylKit tabix export (Linux/macOS only) |
| `bam` | pysam | BAM ingestion for ASM / entropy (Linux/macOS only) |
| `distributed` | dask | Distributed compute via Dask |
| `ray` | ray | Distributed compute via Ray |
| `gpu` | cupy-cuda12x | GPU IRLS via CuPy (CUDA 12) |
| `gpu_jax` | jax[cuda12] | GPU IRLS via JAX (CUDA 12) |
| `all` | report + export + anndata + viz + distributed | Full feature set (no GPU extras) |

To install everything except GPU extras:

```bash
pip install -e ".[all]"
```

To include development tools:

```bash
pip install -e ".[dev]"
```

!!! note "pysam on Windows"
    The `methylkit` and `bam` extras depend on pysam, which does not ship
    Windows wheels. These extras are available on Linux and macOS only.

## CLI

Installing epykit registers a console script:

```bash
epykit --help
```

The CLI mirrors the Python API and is suitable for shell pipelines and
HPC job scripts. See `epykit --help` for subcommands and options.

## Verify installation

```python
import epykit as ep
print(ep.__version__)
```
