# Codebase Concerns

**Analysis Date:** 2026-06-06

## Tech Debt

### Deprecated `log2_odds_ratio` Column (Transitional NaN-Fill)

- **Issue:** Per-CpG DMC output contains both `log2_odds_ratio` (backend-specific, transitional NaN-filled) and `log2_odds_ratio_pooled` (for `lr`/`fisher` backends). The old column is deprecated and scheduled for removal in 1.1.
- **Files:** 
  - `src/epykit/dmc.py:64` (schema definition)
  - `src/epykit/dmc.py:1796–1820` (emission logic)
  - `src/epykit/tl.py:474–475` (deprecation warning)
  - `src/epykit/dmr.py:1285–1296` (normalisation in DMR pipeline)
- **Impact:** Downstream code must read correct backend-specific columns (`log2_odds_ratio_pooled` or `coef_treatment_log2` for GLM). Existing code reading the old column gets NaN. No functional breakage in 1.0 but signals API churn ahead.
- **Fix approach:** At 1.1 release, remove column from schema and cleanup code.

### Deprecated `MethylData._analysis_root` Property

- **Issue:** Private property with a public deprecation warning directing users to use the public `analysis_root` name instead.
- **Files:** `src/epykit/methyldata.py:47, 56`
- **Impact:** Internal code still using `_analysis_root` will trigger DeprecationWarning (which is surfaced in test `filterwarnings` configuration). Users are warned not to use private API.
- **Fix approach:** Migrate internal code to public API and remove the private property in 1.1.

### Deprecated `pp.unite()` Function

- **Issue:** `pp.unite()` deprecated in favour of `pp.set_unite_type()` for explicit type control. Warning added but function still works.
- **Files:** `src/epykit/pp.py:173`
- **Impact:** Users calling `pp.unite()` see a DeprecationWarning; code still runs. CLI has no `--unite` flag; the parameter is hidden from the public surface.
- **Fix approach:** Remove `pp.unite()` in 1.1 after a stable 1.0 release.

### Deprecated `dmr_hmm` Module (Renamed to `dmr_segment`)

- **Issue:** The `epykit.dmr_hmm` module is a shim that re-exports the renamed `call_dmr_rule_segment` function from `epykit.dmr_segment` with a DeprecationWarning. The name was misleading (the engine is not Baum-Welch fitted).
- **Files:** `src/epykit/dmr_hmm.py` (entire shim module)
- **Impact:** Old code importing `from epykit.dmr_hmm import call_dmr_hmm` still works but warns. No functional loss in 1.0.
- **Fix approach:** Remove the shim module entirely in 1.1; direct users to `dmr_segment.call_dmr_rule_segment`.

## Known Issues & Research Limitations

### `lr+` Power Stack — Unvalidated on Real WGBS, Inflates Call Counts

**Critical caveat:** The `power_stack` parameter bundles four research components (`neighbour_combine`, `fdr_method="fdr_tsbh"`, `sep_fallback`, `dispersion="eb"`) into a single opt-in switch. It is **NOT validated as universally superior** to bare `lr` on realistic WGBS data.

- **Files:**
  - `src/epykit/dmc.py:2322–2470` (neighbour Stouffer combining)
  - `src/epykit/dmc.py:2567–2649` (TSBH via statsmodels)
  - `src/epykit/dmc.py:947–983` (separation fallback)
  - `src/epykit/dmc.py:775–820` (EB dispersion shrinkage)
  - `src/epykit/tl.py:498–532` (dispatcher and `power_stack` parameter)

- **Problem:** On real WGBS (GSE263850 at q=0.05), `power_stack="lr+"` inflates the DMC call count ~13× relative to bare `lr` at the same threshold. Tuning basis (Piao simulator, φ ≈ 0.4) is underdispersed compared to real WGBS (φ ≈ 1.5–5), causing FPR drift under realistic data.

- **Current mitigation:**
  - `power_stack="off"` is the default; bare `lr` is the recommended and documented engine.
  - Benchmark paper (1.0) claims around `lr`, not `lr+`.
  - EB prior (`dispersion="eb"`) is marked unvalidated against external dispersion estimators in benchmark docs.
  - Combined p-values are stored in separate `pvalue_combined`/`qvalue_combined` columns so users can opt-in explicitly.

- **Future work:**
  - CLI flags for individual `lr+` knobs deferred to 1.1 (currently only `power_stack` parameter is public).
  - Research validation needed: test against real-data cohorts at multiple dispersion regimes before promoting from "research knob" to recommended default.

- **For implementers:** When adding DMC features, do NOT change `lr+` defaults without re-running ablations and benchmarking against realistic dispersion regimes.

### Benchmark Reproducibility & Outstanding Review Items (M1–M7)

- **Files:** 
  - `benchmark/phi_sweep_export_2026-06-05/SESSION_HANDOFF.md` (latest status)
  - `benchmark/README.md` (full reproduction protocol)
  - `benchmark/scripts/renv.lock` (R environment pinned to R 4.5.3)

- **Current status (as of 2026-06-05):** The φ-sweep (dispersion robustness study) is computed but NOT YET committed. Review items addressed:
  - M4 (reproducibility): `renv.lock` + Dockerfile R4.5.3 committed + pushed.
  - M5.2 (FPR wording): "matched FPR" → "comparably negligible FPR" committed + pushed.
  - M7 (AUROC duality): 0.9999 vs 0.928 reconciliation committed + pushed.
  - M6 (EB disclaimer): Unvalidated prior caveat committed + pushed.
  - M1+M2 (φ-sweep + methylKit fairness fix): **COMPUTED but NOT committed** — found critical issue where first φ-sweep used methylKit default mode (unfair); re-run with MN mode shows epykit maintains FDR control (0.02–0.03) across φ=1–5 while methylKit MN degrades to ~0.30 at φ≈5.

- **Outstanding work (M1, M5.1, M5.3, M5.4):** Writing manuscript paragraphs, reconciling speed claim numbers, copying laptop-local files, fixing 2 stale timing claim assertions in `claims.yaml` that point to non-existent `timings_post_phase3.parquet`.

- **Impact:** The benchmark is not bit-reproducible until the φ-sweep is committed. Speed claims in the paper (§4.1) are based on 0.86s / 6.80s measurements that can't be verified from current committed data.

## Windows Compatibility

### Platform-Specific Optional Dependencies

- **Issue:** Some extras are gated to non-Windows platforms:
  - `pyBigWig` (Linux/macOS only in `pyproject.toml`)
  - `pysam` (Linux/macOS only; powers `methylkit` and `bam` extras)
  - On Windows, these extras fail to install and their modules raise helpful ImportErrors.

- **Files:**
  - `pyproject.toml:71, 83–85, 87–89` (platform gates)
  - `src/epykit/bam_io.py:24, 62–63` (pysam unavailable on Windows)
  - `src/epykit/export.py:16, 146–147, 164` (pyBigWig unavailable on Windows)
  - `src/epykit/methylkit_io.py:13, 100` (pysam gates)

- **Current mitigation:** CI tests Windows + Linux across Python 3.9–3.12. Code paths raise clear ImportError with install hints rather than cryptic failures. `@pytest.importorskip` gates tests appropriately.

- **Impact:** Users on Windows cannot use BAM ingestion or BigWig export directly. `methylkit_io` offers a pysam-free fallback for tabix-indexed TSV reading but lacks the `pysam`-based feature set.

- **For implementers:** If adding features that touch BAM/BigWig APIs, test on Windows or gate behind `pytest.importorskip`.

### Temp Directory Sizing (Critical on Windows)

- **Issue:** Windows default `%TEMP%` on `C:\` is often small/full, causing out-of-disk-space errors during DMC/DMR staging.

- **Files:** `src/epykit/_config.py` (the `set_tmp_dir` solution)

- **Current mitigation:** `ep.set_tmp_dir(path)` allows users to redirect transient files (per-chrom DMC parquets, DMR tile aggregations) to a larger disk. The function is documented and exposed as a public API.

- **Impact:** Low if users discover the config knob. High if not — whole-genome runs will OOM or fail with "no space on device."

## Streaming Contract & Memory Model

### Preservation of Streaming DMC Pipeline

- **Critical constraint:** The DMC engine (`process_chromosomes_dmc(..., return_store=True)`) returns a `DMCStore` handle rather than materializing the full per-CpG result frame in memory. Peak memory is O(largest chromosome), not O(genome). All downstream code (`apply_multiple_testing_correction`, `call_dmr_sliding_window`) streams from this store.

- **Files:** 
  - `src/epykit/_dmc_store.py` (DMCStore interface)
  - `src/epykit/dmc.py:2322–2470, 2567–2649` (combiners that stream)
  - `src/epykit/dmc.py` module docstring (architectural note)

- **Why it matters:** Whole-genome WGBS (~22M CpGs) cannot fit as a single Polars/Pandas frame on typical hardware. Breaking the streaming contract will cause catastrophic OOM on real datasets.

- **For implementers:** When adding DMC-stage features (new combiners, FDR methods, filters), ensure you read from the DMCStore iterator, not from a materialized frame. Do not call `.collect()` on lazy scans over the full partition tree.

### Whole-Genome Data Never Loaded as Single Frame

- **Architectural rule:** Every preprocessing and analysis step uses `pl.scan_parquet()` lazy iteration or chunked operations. No step reads all 22M CpGs into a single in-memory structure.

- **For implementers:** If you need to process millions of CpGs, use Polars lazy scans with pushed-down filters, or iterate per-chromosome. Never do `df = pl.read_parquet(huge_dir)` at the top level.

## Fragile Areas & Complex Modules

### Large, Interdependent Files

**Complexity hotspots:**

1. **`dmc.py` (3,048 lines, 40 functions):** The statistical engine centre with four per-CpG backends, dispersion estimation, separation fallback, EB shrinkage, Stouffer combiners, and TSBH FDR. High test coverage but tight interdependencies. Changes to dispersion or p-value logic ripple everywhere.
   - Files: `src/epykit/dmc.py` + `src/epykit/_dmc_store.py`
   - Test coverage: Good but missing coverage on `lr+` combiner weights and GLM contrast paths.

2. **`dmr.py` (1,925 lines):** Four DMR callers (chain_merge, tile, sliding-window, HMM) with optional permutation FDR, each with parameter presets. Interdependencies with the DMC store. Benchmark re-runs needed if you change scoring logic.
   - Files: `src/epykit/dmr.py` + `src/epykit/_hmm.py`
   - Test coverage: Good but permutation FDR tested sparsely (computationally expensive).

3. **`tl.py` (1,984 lines, 40+ functions):** High-level orchestrators for the analysis flow. `dmc()` is particularly complex due to the `power_stack` dispatch and covariate design routing.
   - Files: `src/epykit/tl.py`
   - Test coverage: Good but `power_stack="lr+"` branches tested minimally against real data.

4. **`_glm.py` (1,066 lines):** Batched IRLS binomial GLM with Wilkinson formula parsing, Wald/F contrasts, and numerical stability handling. Mirrored in `_glm_gpu.py` for CuPy/JAX backends.
   - Files: `src/epykit/_glm.py` + `src/epykit/_glm_gpu.py`
   - Test coverage: CPU backend well-tested; GPU backend requires CuPy (optional, tested sparsely on CI).

### GPU Backend Maturity

- **Issue:** GPU backends (`_glm_gpu.py` with CuPy, JAX variant not yet in main) are optional and experimental. They are tested only when CuPy is importable (which is rare on CI Windows runners).

- **Files:** 
  - `src/epykit/_glm_gpu.py` (CuPy IRLS port)
  - `tests/test_glm_gpu.py` (gates on CuPy availability)
  - `pyproject.toml:99–105` (gpu/gpu_jax extras, excluded from `[all]`)

- **Current state:** CuPy tests skip cleanly if the wheel isn't available. Numerical parity vs CPU (1e-6 tolerance) is enforced. JAX backend is planned but not integrated.

- **For implementers:** GPU code is load-bearing only for users who explicitly install `epykit[gpu]`. If you change IRLS, synchronize both CPU and GPU paths line-by-line.

## Test Coverage Gaps

### `lr+` Power Stack Under-Tested on Real Data

- **What's not tested:** The `lr+` power stack (`power_stack="lr+"`, `power_stack="conservative"`, or individual knobs like `neighbour_combine=True`, `dispersion="eb"`) is tested against simulated data only. Benchmark Study 3 uses bare `lr`.

- **Files:** 
  - `src/epykit/tl.py:293, 498–532` (power_stack definition and dispatch)
  - `tests/` (no real-data ablation of lr+ knobs)
  - `benchmark/paper/paper.md` (benchmark uses bare `lr`)

- **Risk:** Under real dispersion (φ ≈ 1.5–5), `lr+` inflates call counts 13× at q=0.05. The φ-sweep (pending commit) shows this; bare `lr` remains the validated engine. Users may discover inflated FPR after publication.

- **Priority:** High. Before promoting `lr+` beyond "research knob," gather multi-cohort validation showing FPR control or honest power trade-off.

### Permutation FDR Tested Sparsely

- **What's not tested:** Empirical FDR for both DMC and DMR (`tl.dmc(..., empirical_fdr=True, n_perm=N)` and `tl.dmr(..., empirical_fdr=True, n_perm=N)`) is implemented and validated on small simulated sets but not against real data at scale.

- **Files:** `src/epykit/dmc.py:2567–2649` (TSBH FDR), `src/epykit/dmr.py` (permutation dispatch)

- **Risk:** Permutation tests are computationally expensive (N re-runs of the full pipeline); most users skip them. If you run permutation FDR on real data and it diverges from BH, it's hard to debug.

- **Priority:** Medium. Recommend adding one real-data permutation validation test at slow-marker tier.

### GLM Contrast Tests Missing

- **What's not tested:** GLM contrasts (multi-group designs with `design` parameter and `formula`) are tested on simulated data but not validated end-to-end against published multi-group WGBS studies.

- **Files:** `src/epykit/_glm.py` (contrast logic), `src/epykit/tl.py:36–55` (design routing)

- **Risk:** A user applies GLM with three groups and a batch covariate; results diverge from published pipelines. Debugging requires digging into IRLS numerical stability + formula parsing.

- **Priority:** Medium. Add a published multi-group WGBS benchmark to Study 3.

## Performance Bottlenecks

### DMR Tile-Based Caller Under-Optimized

- **Issue:** The tile-based DMR caller (`call_dmr_tile_based`) reads per-tile summaries from the DMC store but rebuilds them in memory before scoring. This can be slow on high-coverage data.

- **Files:** `src/epykit/dmr.py` (tile-based caller)

- **Benchmark note:** Study 3 (real GSE263850) shows methylKit-tile 12,372 s vs epykit-tile 675 s = 18× faster, but epykit-chain_merge is 28× faster still. Tile-based is not recommended for large datasets.

- **Priority:** Low. The chain_merge caller is the paper recommendation. Tile-based is available for smaller datasets and the speed is acceptable.

### Plotly Report Generation

- **Issue:** The HTML report generation (`ep.pl.report()`, `report.py`) renders interactive Plotly figures for every DMC/DMR plot. For large datasets (millions of CpGs), this can be slow and produce huge HTML files.

- **Files:** `src/epykit/report.py`, `src/epykit/pl/_plotly.py`

- **Current state:** The report is generated on-demand and can be cached in `md.uns["_report_cache"]`. Benchmark Study 3 reports are ~50 MB uncompressed.

- **Priority:** Low in 1.0. Optimization (decimation, server-side rendering) deferred to 1.1.

## Scaling Limits

### Memory Scaling with Chromosome Size

- **Limit:** Peak memory during DMC is O(largest chromosome). For hg38 chromosome 1 (~250M bp), with deep coverage (30×+), the per-site accumulators (Welford state) can be large. Typical whole-genome run is 12–50 GB depending on coverage and sample count.

- **Current state:** All analysis uses lazy Polars scans and per-chromosome partitioning. No single CpG site holds redundant state.

- **For implementers:** When adding streaming statistics, keep per-site storage at O(1) per CpG (e.g., Welford accumulators, not full read lists).

### Benchmark Re-run Time

- **Current:** The φ-sweep (dispersion robustness) at 60 cells × 8 tools ≈ 3.5 hours end-to-end on a Linux workstation. Full re-run of all three benchmark studies (including per-tool parameter sweeps) is 48–72 hours.

- **For maintainers:** Benchmark changes (new tool, new dispersion regime) require 2–3 day turnaround. Plan accordingly; don't commit benchmark claims without test data.

## Missing Critical Features

### CLI Flags for `lr+` Knobs Deferred

- **Feature gap:** The individual `lr+` components (`neighbour_combine`, `sep_fallback`, `fdr_method="fdr_tsbh"`, `dispersion="eb"`) are exposed as `tl.dmc()` kwargs but NOT as CLI flags. Users of the `epykit dmc` command cannot tune these without Python scripting.

- **Blocks:** Parameterized research workflows and sensitivity studies via the CLI.

- **Priority:** Medium. Deferred to 1.1 with user feedback on demand.

### Multi-Cohort Validation for `lr+` & GLM Contrasts

- **Feature gap:** Benchmarks Study 3 covers one real-data cohort (GSE263850). Multi-cohort validation for `lr+` and GLM-contrast designs is missing.

- **Blocks:** Authoritative guidance on when to use `power_stack="lr+"` or GLM designs vs bare `lr`.

- **Priority:** High for future releases. Currently, users should treat `lr+` as exploratory.

## External Dependencies at Risk

### Polars Version Constraints

- **Issue:** `epykit` requires `polars>=0.20.0` (from pyproject.toml). Polars has moved fast historically; API changes (schema handling, lazy-evaluation semantics) can ripple through the codebase.

- **Current:** The code targets Polars 0.20.0+. No version pinning in the lock file, so users may install newer versions.

- **For maintainers:** Monitor Polars release notes for deprecations in lazy-scan, partition-by, and schema-handling APIs.

### statsmodels (GLM, FDR)

- **Issue:** `statsmodels>=0.14` is a hard dependency for GLM design parsing (patsy) and TSBH FDR (statsmodels.multitest). This is a large library with infrequent breaking changes.

- **Current:** stable.

- **For maintainers:** Keep an eye on statsmodels GitHub for deprecations in `sm.GLM` or `multitest` modules.

### Optional Wheel Availability (pysam, CuPy)

- **Issue:** Linux/macOS-only wheels (pysam) and GPU wheels (cupy-cuda12x) are maintained by third parties. If wheels disappear or build failures occur, users cannot install those extras.

- **Current:** pysam wheels are stable. CuPy requires CUDA 12.x (heavy ~2 GB). JAX alternative not yet integrated.

- **For maintainers:** Test optional-extra availability in CI; add fallback documentation if a wheel becomes unavailable.

## Architectural Debt

### Logging Convention Enforcement

- **Current rule (load-bearing):** Library code uses `logging.getLogger(__name__)` and never calls `print()`. CLI reserves `print()` for final user-facing output.

- **Files:** Every module under `src/epykit/` except CLI entry point, plus `epykit.cli` for final lines.

- **Risk:** If a new contributor adds a `print()` call to a library function, it pollutes stdout and breaks host applications (notebooks, servers). Tests don't catch this if logging is not checked.

- **For maintainers:** Code review must flag `print()` calls in library code. Consider a pre-commit hook or linter rule.

### Store History State Machine

- **Current rule:** Preprocessing state (`_filtered`, `_united`, `_smoothed`, `.state`) is derived from `uns["_store_history"]` rather than stored as independent booleans.

- **Why:** Flags can never drift from reality; `_store_history` is the auditable source of truth.

- **Risk:** If a new preprocessing function (e.g., `pp.normalize_batch()`) is added and stores a boolean flag instead of appending to `_store_history`, then later flag queries will be inconsistent.

- **For maintainers:** Code review must verify that every `pp.*` function appends to `_store_history`. Consider a test that validates state derivation.

### Global Deprecation Warning Gates

- **Current:** The `tl.dmc()` function uses a global `_FISHER_WARNED` flag to emit a one-shot UserWarning when `test="fisher"` is explicitly selected. This gate is session-level (not per-object), so it fires once per Python process.

- **Files:** `src/epykit/tl.py:61, 64–75`

- **Risk:** In test suites with multiple independent runs, the second run of `test="fisher"` doesn't warn (gate already fired). Reviewers may miss the warning on the first run.

- **For maintainers:** Document the gate's scope. Consider improving to per-call or per-analysis-instance warnings if test coverage demands it.

---

*Concerns audit: 2026-06-06*
