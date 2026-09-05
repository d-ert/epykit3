<!-- refreshed: 2026-06-06 -->
# Architecture

**Analysis Date:** 2026-06-06

## System Overview

```text
┌─────────────────────────────────────────────────────────────────────────┐
│                        User API / CLI Layer                             │
│  ep.pp.* (preprocessing) | ep.tl.* (tools) | ep.pl.* (plotting)         │
│  epykit convert | filter | dmc | dmr | annotate | qc-report | report    │
└────────────────────────┬────────────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                      MethylData Central Dataclass                        │
│  `src/epykit/methyldata.py`                                             │
│  obs (sample metadata) | store (Parquet path) | varm (results by test)  │
│  uns (pipeline state: _store_history, dmc, dmr, etc)                    │
└──────┬──────────────────────────┬──────────────────────────┬────────────┘
       │                          │                          │
       ▼                          ▼                          ▼
┌──────────────────┐    ┌──────────────────┐    ┌───────────────────────┐
│  Preprocessing   │    │   Analysis       │    │  Output & Reporting   │
│  `pp` namespace  │    │   `tl` namespace │    │  `pl`, `export`       │
│  `pp.py`         │    │  `tl.py`         │    │  `report.py`          │
│                  │    │  Orchestrators   │    │  `templates/`         │
│ filter_coverage  │    │ that wire pp →   │    │                       │
│ normalize        │    │ engines → results│    │ HTML reports, TSV,    │
│ set_unite_type   │    │                  │    │ BED, AnnData, etc     │
│ smooth           │    │ dmc()            │    │                       │
│                  │    │ dmr()            │    │                       │
│                  │    │ annotate()       │    │                       │
│                  │    │ qc()             │    │                       │
└────────┬─────────┘    └─────────┬────────┘    └───────────┬───────────┘
         │                        │                         │
         ▼                        ▼                         ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                      Statistical Engines Layer                          │
│                                                                         │
│  DMC: dmc.py (per-CpG tests)   │ DMR: dmr.py (regional tests)         │
│  • lr (default, n≥2)            │ • chain_merge (default, DSS-like)    │
│  • welch_t (Welch t on β)       │ • tile_based (read-pooled)          │
│  • fisher (n=1 fallback)        │ • sliding_window (p-value combine)  │
│  • glm (covariate-adjusted)     │ • segment (HMM)                     │
│                                 │                                      │
│  Streaming via _dmc_store.py:   │ Each DMR engine reads DMC results   │
│  • DMCStore (per-chrom parquet) │ via DMCStore.iter_chroms() or scan_ │
│  • O(largest chrom) memory      │ to keep peak memory bounded         │
│  • iter_chroms() / scan_chrom() │                                      │
│  • update_chrom() for BH        │                                      │
└────────────────────┬────────────────────────┬──────────────────────────┘
                     │                        │
                     ▼                        ▼
┌─────────────────────────────────────────────────────────────────────────┐
│              Partitioned Parquet Methylstore (Source of Truth)          │
│                     `src/epykit/io.py` (ingestion)                      │
│                                                                         │
│  <store>/.cache/raw/sample=<id>/chrom=<chr>/part-0.parquet             │
│  <store>/.cache/filtered/... (after filter_coverage)                   │
│  <store>/.cache/normalized/... (after normalize_coverage)              │
│  <store>/.cache/dmc/<test>/chrom=*.parquet (DMC streaming store)        │
│                                                                         │
│  Key design: Polars lazy scans (pl.scan_parquet) over per-sample,       │
│  per-chromosome partitions. Never materialize whole genome (~22M CpGs)  │
│  in RAM. Each preprocessing step writes a new cached store and repoints │
│  md.store at it.                                                        │
└─────────────────────────────────────────────────────────────────────────┘
```

## Component Responsibilities

| Component | Responsibility | File |
|-----------|----------------|------|
| **MethylData** | Central dataclass holding obs, store, varm (results), uns (state) | `methyldata.py` |
| **Preprocessing** | filter_coverage, normalize_coverage, set_unite_type, smooth | `pp.py` |
| **I/O & Convert** | Bismark/MethylDackel/BED ingestion → partitioned Parquet | `io.py`, `convert.py` |
| **DMC Engines** | Per-CpG likelihood-ratio, Welch t, Fisher, GLM tests | `dmc.py` |
| **DMC Store** | Persistent per-chromosome parquet handle for streaming | `_dmc_store.py` |
| **GLM** | Wilkinson formula → design matrix, batched IRLS, contrasts | `_glm.py`, `_glm_gpu.py` |
| **DMR Engines** | chain_merge, tile_based, sliding_window, HMM segmentation | `dmr.py`, `_hmm.py`, `dmr_hmm.py` |
| **Annotation** | Gene features (GTF), CpG islands (BED), overlap logic | `annotate.py` |
| **QC** | Bisulfite conversion, coverage, methylation levels, sex check | `qc.py` |
| **Plotting** | Matplotlib plotter wrappers (volcano, MA, QC, etc.) | `pl/*.py` |
| **Reporting** | Interactive HTML via Jinja2 + Plotly | `report.py`, `templates/` |
| **Export** | BED, BedGraph, BigWig, AnnData, MuData, methylKit, MultiQC | `export.py`, `anndata_io.py`, `methylkit_io.py`, etc. |
| **CLI** | argparse entry point; mirrors API via shared engine functions | `cli.py` |
| **Orchestration** | High-level wiring of pp → tl.dmc/dmr/annotate → varm/uns | `tl.py` |

## Pattern Overview

**Overall:** Scanpy-inspired modular WGBS pipeline built around an immutable partitioned-Parquet methylstore and a central MethylData dataclass. Preprocessing mutates md.store (repointing at new cached stages), while analysis populates md.varm (keyed by test name: dmc_lr, dmc_glm, etc.) and md.uns (pipeline state, results metadata). State is derived from unserialized _store_history rather than stored as independent flags, preventing drift.

**Key Characteristics:**
- **Streaming-first**: Polars lazy I/O over per-sample, per-chromosome parquets; no whole-genome in RAM
- **Immutable stores**: Each preprocessing step writes a new store under `.cache/<step>/`; md.store is repointed, not mutated in-place
- **DMC→DMR pipeline**: DMCStore (per-chrom persistent parquet) is created once by process_chromosomes_dmc(), then streamed by BH correction and DMR callers; peak memory is O(largest chromosome), not O(genome)
- **Convergent API/CLI**: Both ep.py API (ep.tl.dmc, ep.tl.dmr) and CLI (epykit dmc, epykit dmr) converge on the same engine functions; no duplication
- **Derived preprocessing state**: _filtered, _united, _smoothed are @property aliases that read from uns["_store_history"]; adding a new pp.* step means appending to _store_history

## Layers

**User-Facing API/CLI Layer:**
- Purpose: Provide scanpy-style namespaces (ep.pp.*, ep.tl.*, ep.pl.*) and argparse CLI
- Location: `tl.py`, `pp.py`, `pl/`, `cli.py`
- Contains: High-level orchestrators that wire preprocessing → engines → result storage, plotting shortcuts
- Depends on: MethylData, all engines, all I/O modules
- Used by: End users via Python notebooks or command-line

**MethylData Central Layer:**
- Purpose: Immutable dataclass holding reference to the methylstore + metadata (obs), results (varm, uns)
- Location: `methyldata.py`
- Contains: Dataclass definition, save/load, property aliases for state, convenience methods (.dmc, .dmrs_to_bed, .report, etc.)
- Depends on: polars, pathlib (for store access)
- Used by: Every API function; primary exchange format between layers

**Preprocessing Layer:**
- Purpose: Mutate md.store (repoint to filtered, normalized, united, smoothed cached stores)
- Location: `pp.py` wraps `filter.py`, `_cache.py`
- Contains: filter_coverage, normalize_coverage, set_unite_type, smooth; each appends to uns["_store_history"]
- Depends on: MethylData, Parquet I/O, count_store_rows caching
- Used by: Orchestrators in tl.py; users via ep.pp.*

**Statistical Engine Layer:**
- Purpose: Compute per-CpG (DMC) and regional (DMR) statistics, write results to MethylData.varm and .uns
- Location: `dmc.py`, `dmr.py`, `_hmm.py`, `dmr_hmm.py`, `_glm.py`
- Contains: Four DMC tests (lr, welch_t, fisher, glm), four DMR callers (chain_merge, tile, sliding_window, HMM), dispersion shrinkage, separation fallback, neighbour combining
- Depends on: numpy, scipy, statsmodels, polars lazy scans over methylstore
- Used by: Orchestrators in tl.py; users via ep.tl.dmc, ep.tl.dmr, or CLI

**DMC Streaming Store Layer:**
- Purpose: Persistent per-chromosome parquet directory handle; enables BH correction and DMR streaming without loading full table
- Location: `_dmc_store.py`
- Contains: DMCStore class with iter_chroms(), scan_chrom(), update_chrom(), to_dataframe(); manifest tracking
- Depends on: polars, pathlib
- Used by: process_chromosomes_dmc (writer), apply_multiple_testing_correction (reader/updater), DMR engines (readers)

**Annotation Layer:**
- Purpose: Overlap DMC/DMR results with genes (GTF) and CpG islands (BED), populate feature_type, gene_name, cpg_context columns
- Location: `annotate.py`
- Contains: annotate_features (HOMER-priority overlap), annotate_cpg_islands (island/shore/shelf/open-sea), caching logic for GTFs
- Depends on: bioframe (pure-Python interval overlap), pandas (for GTF parsing)
- Used by: ep.tl.annotate, plotting code that needs feature context

**QC Layer:**
- Purpose: Compute per-sample metrics (bisulfite conversion, coverage uniformity, methylation level, sex/contamination checks)
- Location: `qc.py`
- Contains: bisulfite_conversion_rate, coverage_uniformity, global_methylation_report, sex_check, contamination_estimate, power_calc
- Depends on: numpy, scipy (for stats), polars (for Parquet scans)
- Used by: ep.tl.qc, ep.pl.qc_dashboard, report.py

**Plotting Layer:**
- Purpose: Matplotlib plotters for QC, differential analysis, genomic context, clustering, metaplots
- Location: `pl/` subdirectory with _style.py (shared theme), _plotly.py (Plotly twins), _utils.py (helpers)
- Contains: Per-plot modules (qc.py, differential.py, genomic.py, etc.); Plotly counterparts in report.py
- Depends on: matplotlib, numpy, polars, optional plotly (lazy import in report.py)
- Used by: ep.pl.* shortcuts, report.py for embedded Plotly figures, user notebooks

**Reporting Layer:**
- Purpose: Generate self-contained interactive HTML report with Jinja2 templates + embedded Plotly figures
- Location: `report.py`, `templates/` directory
- Contains: generate_report orchestrator, HTML table renderers, data extraction for Jinja context
- Depends on: jinja2, plotly (lazy imports), matplotlib (for figure serialization)
- Used by: ep.report (MethylData method), CLI epykit report, users via md.report()

**Export Layer:**
- Purpose: Serialize MethylData results to external formats (BED, BedGraph, BigWig, AnnData, MuData, methylKit, MultiQC)
- Location: `export.py`, `anndata_io.py`, `mudata_io.py`, `methylkit_io.py`, `multiqc_export.py`
- Contains: Format-specific writers; each lazy-imports its optional dependency (e.g. scanpy for AnnData)
- Depends on: polars, pathlib; optional deps per module
- Used by: User convenience methods (md.to_anndata, md.to_bedgraph, etc.), CLI epykit export

## Data Flow

### Primary Request Path: Ingestion → DMC → DMR → Annotation → Report

1. **Ingestion** (`ep.read_bismark()`, `ep.convert_sample()`) — `io.py:read_bismark()`, `convert.py:convert_sample()`
   - Read Bismark .cov / MethylDackel .bedGraph / combined-strand BED
   - Write partitioned Parquet under `<store>/.cache/raw/sample=<id>/chrom=<chr>/part-0.parquet`
   - Return MethylData with md.store pointing at raw cache

2. **Preprocessing** (`ep.pp.filter_coverage()`, `ep.pp.normalize_coverage()`) — `pp.py:filter_coverage()`, `pp.py:normalize_coverage()`
   - Read from md.store (lazy Polars scans)
   - Write filtered/normalized Parquet to new `.cache/<step>/` directory
   - Repoint md.store, append to uns["_store_history"]

3. **DMC** (`ep.tl.dmc()`) — `tl.py:dmc()` → `dmc.py:process_chromosomes_dmc()`
   - Per-chromosome streaming over md.store (Polars lazy scans)
   - Accumulate per-replicate state (Welford online algorithm) for each CpG
   - Write per-chrom DMC results to `.cache/dmc/<test>/chrom=*.parquet` (DMCStore)
   - Apply multiple testing correction: `dmc.py:apply_multiple_testing_correction()` streams BH from DMCStore, writes qvalues back per-chrom
   - Populate md.varm[f"dmc_{test}"] with full table (eager load from DMCStore for in-memory analysis)

4. **DMR** (`ep.tl.dmr()`) — `tl.py:dmr()` → one of `dmc.py:call_dmr_chain_merge()`, etc.
   - Read DMC results from DMCStore (streaming via iter_chroms) or md.varm
   - Combine per-CpG p-values or pool reads per tile
   - Write DMR table to md.uns["dmr"]

5. **Annotation** (`ep.tl.annotate()`) — `tl.py:annotate()` → `annotate.py:annotate_features()` + `annotate_cpg_islands()`
   - Read GTF once, cache in LRU (bounded to 2 by default)
   - Per-chromosome overlap with DMC/DMR results (bioframe)
   - Populate feature_type, gene_name, cpg_context in-place on md.varm[dmc_*] and md.uns["dmr"]

6. **Report** (`md.report()`, `ep.report.generate_report()`) — `report.py:generate_report()`
   - Extract KPIs from md.obs, md.varm, md.uns
   - Render Jinja2 template with embedded Plotly figures
   - Write self-contained HTML to disk

### DMC Streaming Path (Memory Efficiency)

1. `process_chromosomes_dmc(..., return_store=True)` returns **DMCStore** handle
2. DMCStore.path points to `<store>/.cache/dmc/<test>/` with per-chrom files and `.epykit_dmc_manifest.json`
3. `apply_multiple_testing_correction(store)` streams via `store.iter_chroms()`, updates qvalues per-chrom with `store.update_chrom(chrom, df)`
4. DMR sliding-window streams via `store.iter_chroms()`, combines p-values per window on the fly
5. Only when materializing results for export/report do we call `store.to_dataframe()` (whole-genome eager load)

**Peak memory: O(largest chromosome) ≈ 50–200 MB per chrom, not O(genome) ≈ 3–5 GB**

### State Management

**Preprocessing state (derived, not stored):**
```python
md._filtered      # bool: True iff any h.get("step") == "filtered" in uns["_store_history"]
md._united        # bool: True iff "unite" in md.uns
md._smoothed      # bool: True iff "smooth_path" in md.uns
md.state          # list[str]: ordered steps read from _store_history + united/smoothed
```

**Analysis results (stored in varm/uns):**
```python
md.varm["dmc_lr"]              # DMC table from test="lr"
md.varm["dmc_glm"]             # DMC table from test="glm"
md.varm["dmc_lr_annotated"]    # DMC + annotation columns
md.uns["dmr"]                  # DMR table
md.uns["dmc"]                  # dict with {"last_key": "dmc_lr", "store_path": "...", "input_sig": "..."}
md.uns["_store_history"]       # list of {"step": "filtered", "path": "...", "n_sites": N} dicts
```

## Key Abstractions

**MethylData Dataclass:**
- Purpose: Immutable reference to a methylstore + sample metadata + results
- Examples: Every user-facing function accepts and mutates a MethylData instance
- Pattern: Attributes are store (str, Parquet path), obs (pl.DataFrame samples), varm (dict[str, pl.DataFrame] results), uns (dict pipeline state); most attributes are immutable, but dict contents (varm, uns) are mutated in-place

**DMCStore Handle:**
- Purpose: Opaque handle to per-chromosome parquet directory + manifest; enables streaming without loading full table
- Examples: Returned by process_chromosomes_dmc; consumed by apply_multiple_testing_correction, DMR callers
- Pattern: Frozen dataclass with path + test; iter_chroms()/scan_chrom() for streaming, update_chrom() for atomic per-chrom updates

**Statistical Test Engines:**
- Purpose: Pluggable per-CpG or regional test implementations (lr, welch_t, fisher, glm, chain_merge, tile, etc.)
- Examples: tl.dmc(..., test="lr") dispatches to dmc._process_chrom_lr; tl.dmr(..., method="chain_merge") dispatches to dmr.call_dmr_chain_merge
- Pattern: Auto-dispatcher (_auto_test) resolves test="auto" to "fisher" at n<2, "lr" at n≥2; each engine outputs canonical schema + engine-specific extras

**Caching & Manifests:**
- Purpose: Cache DMC results by input signature (SHA-256 of store path, samples, test, unite mode, etc.)
- Examples: DMC manifest includes input_sig, total_sites, per-chrom entries; rerun with same sig reads from cache
- Pattern: Each major computation step (DMC, DMR, annotation) has a cache_key function that hashes inputs; if cache hit, skip recomputation

## Entry Points

**Python API:**
- Location: `__init__.py` (public re-exports)
- Triggers: User imports ep.tl.dmc, ep.pp.filter_coverage, etc.
- Responsibilities: Wire inputs → engines → MethylData mutations

**CLI:**
- Location: `cli.py:main()`
- Triggers: User runs `epykit <subcommand> ...`
- Responsibilities: Parse argparse args, call API functions, handle result formatting/export

**Ingestion Entry Points:**
- `ep.read_bismark()` — Read Bismark .cov, build partitioned Parquet
- `ep.read_methyldackel()` — Read MethylDackel .bedGraph
- `ep.read_combined_strand_bed()` — Read combined-strand BED
- `ep.read_nfcore_methylseq()` — Read nf-core methylseq output directory
- `ep.load()` — Load saved MethylData from disk (methyldata.json + obs.parquet + varm/*.parquet + uns/*.parquet)

## Architectural Constraints

- **Threading:** Single-threaded within epykit; Polars/numpy use their internal parallelism. DMC per-chromosome loop is embarrassingly parallel (implemented for speedup in compute backends), but orchestration (pp.*, tl.dmc) is sequential
- **Global state:** GTF/built-features caches in annotate.py are module-level dicts (LRU); _FISHER_WARNED flag gates one-shot warning; TMPDIR env var mirrors set_tmp_dir (load-bearing for Windows)
- **Circular imports:** None; clean DAG: methyldata → (pp, tl, pl, export) → (dmc, dmr, annotate, qc) → (io, _cache, _compute)
- **Store mutation:** md.store is reassigned (not mutated); varm/uns dicts are mutated in-place; obs is typically mutated (columns added during QC)
- **Lazy vs. eager:** Polars lazy scans for I/O (pl.scan_parquet), eager load only when needed (md.varm tables, DMCStore.to_dataframe for export)

## Anti-Patterns

### Materializing Whole Genome in RAM

**What happens:** Old code concat'd all per-chrom parquets into one DataFrame, held it + the assembled output simultaneously
**Why it's wrong:** 22M CpGs × 8 columns × 8 bytes ≈ 1.4 GB; with concat overhead, ~3–5 GB peak. Breaks at genome scale; blocks users with limited RAM.
**Do this instead:** Use DMCStore.iter_chroms() to stream per-chrom, process chromosome-by-chromosome. See `dmc.apply_multiple_testing_correction()` in `dmc.py:2200–2300` for the BH correction example.

### Storing Preprocessing State as Independent Flags

**What happens:** Old code had _filtered, _united, _smoothed as separate attributes that could drift from reality if someone manually popped a key from uns
**Why it's wrong:** State can drift; inspection becomes unreliable; adding a new preprocessing step requires adding a new boolean
**Do this instead:** Derive state from uns["_store_history"]. See `methyldata.py:65–94` for the @property pattern. When adding pp.*, append to _store_history, never add a new attribute.

### Print Statements in Library Code

**What happens:** Library modules call print(); output gets mixed with logging and can't be controlled by the calling application
**Why it's wrong:** Notebooks and tools that import epykit get polluted stdout; can't suppress noisy output without monkeypatching
**Do this instead:** Use `logger = logging.getLogger(__name__)` and call logger.info/debug/warning. Library never calls print. CLI (cli.py) reserves print for final result lines. See `methyldata.py:14` and `cli.py:133`.

### Hardcoding Temp Directories

**What happens:** Code calls tempfile.mktemp() or writes to %TEMP% directly
**Why it's wrong:** On Windows, %TEMP% on C:\ is often too small for whole-genome staging; DSS-style smoothing allocates large temp files
**Do this instead:** Import and use _config.set_tmp_dir(path) to redirect both tempfile.tempdir and TMPDIR/TEMP/TMP env vars. See `_config.py` and `pp.py:smooth` example.

## Error Handling

**Strategy:** Explicit early checks (n≥2 per group), graceful fallbacks (allow_n1=True for Fisher), warnings for anti-conservative paths

**Patterns:**
- `_auto_test_simple()` raises ValueError if n<2 and allow_n1=False
- `_warn_fisher_once()` emits one-shot UserWarning for fisher test
- `_check_n1_and_union_footgun()` warns on risky unite+min_samples=0 combos
- GLM rank-deficiency → fallback to reduced contrasts with warning
- Separation-aware fallback (sep_fallback=True) lowers DMC p-values when effect is near-complete separation

## Cross-Cutting Concerns

**Logging:** Library modules use `logger = logging.getLogger(__name__)`; CLI controls verbosity via -v/-q flags. Library never calls print.

**Validation:** 
- Sample group membership (n_treatment, n_control) checked by _auto_test_simple
- min_samples_treatment/min_samples_control validated by tl.dmc before engine dispatch
- GTF path existence checked in annotate.py

**Authentication:** None; all I/O is file-system local (or S3 if using s3fs plugin with Polars)

**Caching:** 
- Parquet caches under .cache/<step>/ (filter, normalize, dmc, dmr, etc.)
- GTF parsed once, LRU-cached in annotate._GTF_CACHE (default 2 slots)
- Built feature index cached in annotate._BUILT_FEATURES_CACHE
- DMC results cache keyed by input signature (SHA-256); cache key includes test, unite, min_samples, dispersion, etc.

---

*Architecture analysis: 2026-06-06*
