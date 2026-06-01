# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`epykit` is a Python-native WGBS methylation analysis pipeline (Bismark / MethylDackel `.cov` → DMC/DMR → annotation → HTML report). Core stack: `polars` lazy I/O over a partitioned Parquet methylstore, `numpy`/`scipy`/`numba` for stats, `statsmodels` + `patsy` for GLMs, `bioframe` for genomic intervals. Pre-1.0 (currently 0.7.2), MIT.

## Common commands

Use `uv` (matches CI in `.github/workflows/test.yml`). Plain `pip` works too.

```bash
uv sync --extra dev --extra all          # install dev + all optional extras
uv run pytest -m "not slow" --strict-markers -ra   # the CI invocation
uv run pytest tests/test_dmc_multigroup.py::test_name   # single test
uv run pytest -m slow                    # the slow tier (>~5s tests, opt-in)
uv run ruff check src/                   # F-only baseline (see pyproject.toml)
uv run mypy src/epykit                   # configured in pyproject.toml
uv run epykit <subcommand> --help        # CLI entry point (epykit.cli:main)
```

CI matrix is `{ubuntu-latest, windows-latest} × {py3.9, py3.12}` — Windows compatibility is load-bearing (some extras like `pyBigWig`, `pysam`-based `methylkit`/`bam` are Linux/macOS only and gated in `pyproject.toml`).

Tests rely on a `slow` marker registered in `pyproject.toml`; `--strict-markers` is enabled so any unregistered marker fails. `epykit`'s own DeprecationWarnings / UserWarnings are surfaced via `filterwarnings` so tests can assert on them — do not silence them globally.

## High-level architecture

### The methylstore is the source of truth

`read_bismark(samplesheet, store_dir=...)` converts each input to **per-chromosome, per-sample Parquet** under `<store>/.cache/raw/sample=<id>/chrom=<chr>/part-0.parquet`. Every subsequent pipeline step (`filter`, `normalize`, `unite`, `dmc`, smoothing, ...) **writes a new cached store** under `<store>/.cache/<step>/...` and **repoints `md.store`** at it. Whole-genome data (~22 M CpGs) is never loaded into RAM as a single frame — Polars lazy scans (`pl.scan_parquet`) over the partition tree are the canonical access pattern.

`MethylData` (`src/epykit/methyldata.py`) is the central dataclass:
- `obs`: per-sample metadata DataFrame (sample_id, group, plus arbitrary columns usable as GLM covariates).
- `store`: current Parquet path (mutates as pp.* steps run).
- `varm`: per-CpG result frames keyed by analysis (`"dmc_lr"`, `"dmc_glm"`, ...).
- `uns`: misc results — `uns["dmr"]`, `uns["qc"]`, `uns["_store_history"]` (the auditable list of which steps ran).
- Preprocessing state (`_filtered`, `_united`, `_smoothed`, `.state`) is **derived** from `uns["_store_history"]` rather than stored as independent booleans, so the flags can never drift from reality. When adding pipeline steps, append to `_store_history` rather than introducing a new flag.

### Scanpy-style namespaces

API mirrors scanpy: `ep.pp.*` (preprocessing — mutates `md.store`), `ep.tl.*` (tools — populates `md.obs` / `md.varm` / `md.uns`), `ep.pl.*` (plotting). The CLI (`epykit convert | filter | dmc | dmr | annotate | qc-report | smooth | report | aggregate-regions | export`) mirrors the same operations. Both API and CLI converge on the same engine functions in `dmc.py` / `dmr.py` / `annotate.py` / `qc.py` — orchestration lives in `tl.py`.

### Streaming DMC pipeline

`process_chromosomes_dmc(..., return_store=True)` returns a **`DMCStore`** handle (`src/epykit/_dmc_store.py`) — a per-chromosome parquet directory under `<methylstore>/.cache/dmc/<test>/` with a `.epykit_dmc_manifest.json`. `apply_multiple_testing_correction` and `call_dmr_sliding_window` stream from this store so peak memory is O(largest chromosome), not O(genome). When adding DMC-stage features, preserve the streaming contract: don't materialize the full per-CpG result frame in memory.

### Statistical engines (`dmc.py`)

Four per-CpG engines survive post-0.7.5: `"lr"` (quasi-binomial likelihood-ratio, default at n ≥ 2), `"welch_t"` (Welch t on raw β), `"fisher"` (pooled Fisher exact, n = 1 fallback), `"glm"` (IRLS binomial with covariates via Wilkinson formula in `_glm.py`). `"auto"` resolves to `"fisher"` at n < 2 and `"lr"` at n ≥ 2. Engines removed in 0.7.5 (raise `ValueError` with a migration hint): `"logit_t"` → use `"welch_t"`; `"bb_lr"` → use `"lr"`; `"score"` → use `"lr"`; `"cmh"` → use `formula='~ group + batch'`. Every engine outputs the same canonical schema (`chrom`, `pos`, `n_case`, `n_control`, `mean_beta_*`, `meth_diff`, `meth_diff_ci_{lo,hi}`, `pvalue`, `qvalue`, `log2_odds_ratio_pooled`) plus engine-specific extras (`coef_treatment` for GLM, `f_stat`/`df1`/`df2` for multi-group F-tests).

**`lr+` power stack** is four opt-in enhancements on top of `lr`, controlled
via the `power_stack` kwarg in `tl.dmc`. The four knobs (`neighbour_combine`,
`fdr_method="fdr_tsbh"`, `sep_fallback`, `dispersion="eb"`) are validated
against methylKit / RADMeth / DSS in `benchmark/REPORT.md`. They are
**opt-in**, not defaults — bare `lr` is what users get out of the box
(`power_stack="off"` is the 1.0 default).

`power_stack="lr+"` / `True` / `"auto"` engages all four at any n.
`power_stack="conservative"` engages only at n ≤ 2 (legacy behavior).
`power_stack="off"` / `False` leaves knobs at user-passed values.

When `neighbour_combine=True`, **`pvalue`/`qvalue` remain the raw per-CpG values; the combined values are added as `pvalue_combined`/`qvalue_combined`** (plus `pvalue_combined_n_neighbours` and `qvalue_combined_reject` as audit columns). Downstream code that wants the combined p-values must read the `_combined` columns explicitly. CLI flags for the `lr+` knobs are deferred to 1.1.

### DMR engines (`dmr.py`)

Four callers — tile-based (default, read-pooled), sliding-window with signed Stouffer combining, HMM segmentation, and DSS-compatible `chain_merge` with presets (`strict`/`default`/`permissive`). All support optional permutation empirical FDR via `tl.dmr(..., empirical_fdr=True, n_perm=N)` which shuffles labels, re-runs the engine, and adds `empirical_pvalue` / `empirical_qvalue`.

### Logging convention (load-bearing)

Library code (everything under `epykit.*` except `epykit.cli`) emits progress through the stdlib `logging` module via `logger = logging.getLogger(__name__)` and **never calls `print()`**. The CLI entry point (`epykit/cli.py`) reserves `print` for final user-facing result lines on stdout; structured progress logs go through logging and are controlled via `-v`/`-q`. This split is what lets host applications and notebooks consume epykit without stdout pollution — preserve it when adding new modules.

### `set_tmp_dir`

`ep.set_tmp_dir(path)` (`_config.py`) redirects `tempfile.tempdir` AND mirrors the value into `TMPDIR`/`TEMP`/`TMP` env vars so Dask/Ray workers inherit it. Used because the default Windows `%TEMP%` on `C:\` is often too small for whole-genome staging. Code that creates `TemporaryDirectory()` automatically honors this — don't hardcode tempdir paths.

### Benchmark directory

`benchmark/` reproduces a head-to-head against eight published DMC/DMR tools on Piao et al. 2021 simulated data; `benchmark/REPORT.md` is the canonical TPR/FPR/F1 record and is what the `lr+` knobs are validated against. Don't change `lr+` knob defaults without re-running the relevant `benchmark/scripts/ab_*.py` ablations. Raw simulated data and run caches are not bundled — see `benchmark/README.md` for the bootstrap.

## Module map (when to look where)

- `methyldata.py` — `MethylData` dataclass, save/load, `.dmc` / `.treatment_ids` / `.control_ids` properties, `region_beta()`.
- `io.py` / `convert.py` — Bismark / MethylDackel / combined-strand BED / nf-core methylseq ingestion → partitioned Parquet.
- `pp.py` — preprocessing wrappers; each function appends to `uns["_store_history"]` and repoints `md.store`.
- `dmc.py` + `_dmc_store.py` — per-CpG engines + streaming store handle.
- `_glm.py` — Wilkinson formula → design matrix, batched IRLS binomial GLM, Wald/F contrasts. `_glm_gpu.py` is a CuPy/JAX backend gated behind extras.
- `dmr.py` + `_hmm.py` + `dmr_hmm.py` — tile / sliding-window / HMM / chain-merge DMR callers + permutation FDR.
- `dvc.py` — iEVORA-style differentially variable CpG calling.
- `annotate.py` — GTF + UCSC `refGene.txt` gene features and CpG-island/shore/shelf/open-sea context.
- `qc.py` — bisulfite conversion rate, coverage uniformity, sex check, contamination estimate, sample correlation, power calc.
- `tl.py` — high-level orchestrators that wire `pp` → engines → `varm`/`uns`.
- `pl/` — matplotlib plotters with shared theme in `_style.py`; Plotly twins live in `report.py`.
- `report.py` + `templates/` — self-contained Jinja2 + Plotly HTML report.
- `export.py` / `anndata_io.py` / `mudata_io.py` / `methylkit_io.py` / `multiqc_export.py` — interop sinks.
- `cli.py` — argparse entry point exposed as the `epykit` console script.
