# Changelog

All notable changes to **epykit** are tracked here. Format roughly follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); the project uses
SemVer (`MAJOR.MINOR.PATCH`).

## [Unreleased]

### Changed

- **Per-engine chromosome runners.** `dmc._process_one_chromosome` now
  builds one frozen `EngineInput` record for the chromosome, dispatches to
  the engine's runner (`_run_fisher`, `_run_lr`, `_run_welch_t`, `_run_glm`,
  `_run_glm_contrast`, keyed by registry name in `_ENGINE_RUNNERS`) and
  hands the runner's reduced per-site `EngineResult` to
  `_finalise_chromosome`, which owns the effect estimates, intervals,
  minimum-sample mask and column assembly. Per-sample stacks and streaming
  accumulators end with the runner's scope. Engine output is unchanged: the
  engine hash gate holds, and every engine path compares bit-identical on a
  fixed fixture. The C901 complexity ceiling in `pyproject.toml` drops from
  38 to 32, the highest remaining function in the source tree.
- **DMC engine facts live in one registry.** `src/epykit/_dmc_engines.py`
  holds one frozen `EngineSpec` per engine (`lr`, `glm`, `welch_t`,
  `fisher` and the internal `glm_contrast`) with the facts the rest of the
  package reads: whether it is a public `test=` choice, whether the `lr+`
  power stack applies, and which effect-size column it emits. Both CLI
  `--test` choice lists come from it, in the same order with the same
  default. An unknown engine name now raises `ValueError` naming the four
  public engines from `tl.dmc` and from `dmc.process_chromosomes_dmc`
  before any DMC store directory is created; it used to reach a
  `NotImplementedError` per chromosome after the directory existed. The
  engines removed in 0.7.5 keep their migration hints, `"auto"` resolves as
  before, and engine output is unchanged (the engine hash gate holds).
- **`tl.dmc` orchestration split into stages.** The body of `ep.tl.dmc` now
  runs nine stages from `src/epykit/_dmc_stages.py` (`plan_run`,
  `run_contrast`, `lookup_resume`, `open_input_store`, `run_engine`,
  `post_process`, `publish`, `persist_resume`, `finish`), each handing a
  frozen plan or outcome record to the next; `publish` is the only writer of
  `md.uns["dmc"]`. The public signature, defaults, result keys, metadata
  record and engine output are unchanged (the engine hash gate holds). One
  observable difference: the `log2_odds_ratio` FutureWarning is now emitted
  on the `resumable=True` cache hit too, where it was silent before.
  Warnings raised by the DMC stages, including the n<2 Fisher fallback
  notice that previously pointed inside `tl.py`, now point at the caller of
  `tl.dmc`. The private `tl._run_dmc_contrast` helper is gone; its body is
  the `run_contrast` stage.
- **CI runs the BAM-backed tests on Ubuntu.** The Ubuntu legs of the test
  matrix and the slow job install the `bam` extra, so `test_asm.py`,
  `test_bam_io.py` and `test_entropy.py` execute instead of skipping.
  Windows legs are unchanged (`pysam` has no Windows wheel).

### Fixed

- **`read_methylation_calls(regions=...)` reported calls past the region
  end.** `bam_io.read_methylation_calls` fetched every read overlapping a
  requested `(chrom, start, end)` window but kept all of the read's calls, so
  positions beyond `end` (and duplicate calls for a read spanning two
  adjacent windows) leaked into the result. Calls are now clipped to the
  half-open `[start, end)` span. The existing region test in
  `tests/test_bam_io.py` catches this; CI never executed it before because
  `pysam` was not installed.

## [1.1.0] — 2026-09-05

Post-1.0 correctness and reproducibility fixes from a pre-submission code
review (Tier 1: paper-blockers + silent-wrong-science criticals), plus a
redesigned HTML report.

Upgrading from 1.0.0: epykit now requires Python 3.10 or newer. The `dev`
extra is gone and contributor tooling lives in `[dependency-groups]`, so
install it with `uv sync --group dev` or `pip install -e . --group dev`
rather than `.[dev]`. CI installs with `--locked` and runs `uv lock --check`,
so a dependency edit must land together with its `uv.lock` update. The
Bismark `.cov` coordinate fix bumps the raw-store manifest to version 2, and
a store that 1.0.0 built from real `.cov` files is rebuilt on the next
`read_bismark` (see the C1 entry under Fixed). The deprecated surfaces that
were scheduled for removal in 1.1 are retained in this release and are now
scheduled for removal in 1.2 (see Changed).

### Added

- **Redesigned HTML report — MultiQC-style dashboard.** `md.report()` /
  `epykit report` now render a fixed-sidebar dashboard (numbered table of
  contents with per-section status dots + scroll-spy, light/dark toggle)
  over a scannable main panel. New content beyond the previous report: a
  headline **"Results at a glance"** auto-narrative + KPI strip + analysis-
  completeness checklist; a preprocessing **step-flow** with per-step
  site-count deltas; per-sample **QC pass/warn/fail badges** against
  documented thresholds, a **global-methylation-per-sample bar**, and a
  clustered **sample-correlation heatmap**; a DMC **p-value histogram**
  (calibration check); a **DMR size distribution**; a **hyper/hypo-by-feature
  stacked bar**; a PCA **scree** plot; and an auto-generated, parameter-
  accurate **Methods & citations** section (with a copy button). Tables are
  sortable, searchable, and export to CSV client-side. Provenance is rendered
  as a clean key/value table with the raw JSON in a collapsible.
- **`self_contained` report flag (default `True`).** `generate_report(...,
  self_contained=True)` embeds the Plotly bundle inline so the single `.html`
  works fully offline (e.g. for emailing or archiving with a paper);
  `self_contained=False` loads Plotly from a CDN for a smaller file. Exposed
  on the CLI as `epykit report --self-contained / --no-self-contained`.
- Interactivity is implemented in inlined vanilla JS
  (`templates/report.js`) — no new runtime dependencies (still `jinja2` +
  `plotly`). Every section degrades gracefully when its upstream step has not
  been run.
- **Region-level (DMR) annotation charts.** `ep.tl.annotate()` now applies the
  **same** annotation columns to the DMR region table (`md.uns["dmr"]`) as to
  the per-CpG tables — gene features *and* CpG-island context (previously the
  DMR table received gene features only). `ep.pl.genomic_context_bar()` and
  `ep.pl.cpg_island_pie()` gained a `level="dmc"` (default) / `level="dmr"`
  argument to count the per-CpG table (`md.dmc`, density-weighted) or the
  per-region table (`md.uns["dmr"]`, the field-standard "fraction of DMRs per
  feature"). The HTML report's **Annotation** section now shows both the
  per-region (DMR) and per-cytosine (DMC) views of each pie side by side, each
  rendering only when its table is present. (`ep.pl.plot_annotation_counts()`
  was already `level`-aware, defaulting to `"dmr"`.)

### Fixed

- **Report "% significant" denominator (report redesign).** DMC summary
  statistics and the genome-wide figures (volcano / MA / Manhattan /
  p-value histogram) now use the full per-CpG result table rather than
  `md.dmc`, which prefers the annotated table and can be the
  significant-only subset when annotation was run on significant sites.
  Previously this could report "% significant = 100%" and drop the
  non-significant cloud from the volcano. The annotated table is still used
  for the annotation breakdowns and the top-DMC table.

- **Bismark `.cov` coordinate convention (C1).** `read_bismark` /
  `convert_sample` now ingest standard 1-based Bismark coverage (`start ==
  end`) correctly. A new `coordinate_base="auto"` (default) detects 1-based
  vs 0-based input from `start`/`end` and shifts 1-based positions by `-1`
  so the store is always 0-based. Previously every CpG from a real
  `*.bismark.cov.gz` was stored 1 bp downstream, mis-annotating sites and
  breaking cross-format `unite`. MethylDackel / combined-strand BED ingestion
  and the benchmark numbers are unaffected (both are genuinely 0-based). Pass
  `coordinate_base="one_based"` / `"zero_based"` to override. **Migration:**
  the raw-store manifest version bumped to 2, so stores built by an older
  epykit from real Bismark `.cov` are rebuilt on the next `read_bismark`.

- **Annotation chromosome-name mismatch (C2).** `annotate_features` /
  `annotate_cpg_islands` now raise on zero chromosome-name overlap (e.g.
  UCSC `chr1` sites vs Ensembl `1` features) and warn below 50% coverage,
  instead of silently labelling every site intergenic / open_sea.

- **Horvath epigenetic-clock transform (C3).** The negative branch of the
  anti-transform dropped the `(1 + adult_age)` factor, producing
  negative/garbage ages for samples under ~20 years; now `21*exp(x) - 1`.

- **CLI/API q-value parity (M-PKG2).** `epykit dmc` now defaults
  `dispersion="eb"` (matching `ep.tl.dmc`) instead of inheriting
  `dispersion="site"`, so the CLI and Python API produce identical q-values.
  Added `--dispersion` / `--reference` / `--fdr-method` flags.
- **DVC anti-conservative at n=2 (M-DVC1).** At n=2/group Brown-Forsythe is
  degenerate (within-group SS is a floating-point residual of 0), which
  exploded the F-statistic into spurious significance — a 2-vs-2 null cohort
  called ~75% of sites DVC. `_brown_forsythe_vectorised` now requires ≥3 finite
  obs per group (NaN otherwise) and `process_chromosomes_dvc` warns at n<3.
- **DVC mislabelled test (M-DVC2).** `tl.dvc` defaulted to `test="bartlett"`
  with a docstring claiming Bartlett runs and Brown-Forsythe is impossible —
  the reverse of reality. Default is now `"brown_forsythe"`; `"bartlett"` is a
  deprecated alias that warns. Added a `min_coverage` floor (M-DVC3).
- **Sliding-window DMR combined raw p, not q (M-DMR1).** The signed-Stouffer
  region combine used the FDR-controlled column (qvalue when present), which is
  not U(0,1) under the null; it now uses raw `pvalue` (the gate still uses
  qvalue), matching chain-merge. The `_stouffer_combine_signed` docstring's
  correlation claim was corrected (it is anti-conservative, not conservative;
  M-DMR2).
- **Annotation intron over-extension + non-determinism (M-ANN1/2).** Introns are
  now built from a running-max exon end (not `shift(1)`), so nested exons across
  transcripts no longer push introns into exonic territory; equal-priority gene
  overlaps now break ties deterministically (stable sort by gene id).
- **Epigenetic-clock / imputation hazards (M-SEC1/3).** `age_clock` warns when
  `impute_missing=False` leaves per-sample-missing CpGs (biased absolute ages);
  `impute_knn_beta(return_mask=True)` returns a `was_imputed` mask so filled
  cells can be excluded from variance analyses.
- **Power calculator (C5).** `qc.power` now uses the exact non-central-t
  two-sample power (was a z-test), adds an overdispersion (`dispersion`/φ) term
  and an `n_tests` multiple-testing adjustment, and validates against
  `statsmodels`. It no longer over-promises (told users they needed too few
  replicates).

<!-- Second pre-submission review (2026-06-07 audit) -->

- **Empirical/permutation FDR correctness (M1–M3, perm-1/2).** Stratified
  permutation now permutes *within* each stratum preserving its original
  treatment/control split (the old global split made paired/batch designs a
  degenerate, maximally-confounded null — M1). `empirical_fdr=True` on a
  non-tile DMR caller (`chain_merge`/`sliding_window`/`segment`) now raises
  `NotImplementedError` instead of silently no-opping and suppressing the
  calibration warning (M2). The DMC permutation null now reproduces the
  observed run's `sep_fallback`/`smoothing` (was anti-conservative — M3).
  Failed permutations are excluded from the empirical-p denominator (perm-1),
  and tile permutations inherit `merge_adjacent`/`backend` (perm-2).
- **`merge_strands=True` without a reference (M5).** A default
  `read_bismark("x.cov")` no longer leaves symmetric CpG dyads un-merged
  (2× sites at ½ coverage). When `merge_strands=True` and no `reference_fasta`
  is given, dyads are merged by position (C at N with its − strand partner at
  N+1); pass `reference_fasta=` for strand-aware merging. **Behavior change**
  on two-strand `.cov` input.
- **GLM over-conservative vs `lr` (M4).** The `glm` p-value path now applies
  the same `DF_PHI_FLOOR=50` to the F-reference df that `lr` uses, so
  `test="glm"` is no longer systematically less powerful than `test="lr"` on
  identical data. **Behavior change:** glm p-values decrease (toward `lr`); the
  bias was conservative so the bare-`lr` headline is unaffected.
- **`aggregate_regions` stale output + overlapping BEDs (M6, M7).** The regions
  store is cleared before re-writing (re-running with a different BED no longer
  mixes regions from two BEDs — M6), and each CpG is now assigned to *every*
  region it overlaps via a range join, fixing dropped/under-counted CpGs for
  nested or overlapping BED regions (M7).
- **methylKit export off-by-one (M9).** `to_methylkit_tabix` now writes
  1-based `base`/`chrBase` (methylKit is 1-based); previously every exported
  CpG was shifted 1 bp left versus annotations. **Behavior change** for
  existing methylKit exports.
- **`welch_t` CI used the normal quantile (M13).** The Δβ CI now uses
  `t(Satterthwaite df)`, matching its own p-value (was `z=1.96`, ~30% too
  narrow at n=3 — a CI/test disagreement).
- **Segment DMR caller unsigned Stouffer (D1).** `dmr_segment` now combines
  per-CpG p-values with the signed Stouffer Z (shared with the other callers),
  so a region's p no longer shrinks toward 0 as it grows when directions are
  mixed (was anti-conservative).
- **Neighbour combine counted untested sites (D9).** `combine_neighbour_pvalues`
  excludes NaN-p (untested) neighbours from both the Stouffer sum and the
  `_n_neighbours` audit count.
- **CLI/API parity (M10, D10–D12, contrast forwarding).** `epykit dmr
  --min-cpgs` is now honored on the default `chain_merge` caller and presets
  are no longer suppressed (M10); `epykit dmr` gained `--min-mean-qvalue`
  matching `tl.dmr`'s region q-filter across chain_merge/sliding_window/tile
  (D11); bare `epykit dmc --allow-n1` at n=1 resolves to the Fisher engine that
  actually has the n=1 fallback (D12); the `epykit dmc` formula/contrast path
  now forwards `--dispersion`/`--reference`/`--fdr-method` and `tl.dmc`'s
  contrast path honors `fdr_method` (was hardcoded `fdr_bh`).
- **Export streaming (M12).** `to_bedgraph`/`to_bigwig` stream the methylstore
  one chromosome at a time instead of materializing the whole sample as
  full-genome Python lists; peak memory is now O(largest chromosome).
- **Normalised store `N_unmeth` broke the coverage invariant.**
  `pp.normalize_coverage` scaled `N_meth` and rebuilt `coverage` from the
  scaled counts but wrote the unscaled `N_unmeth` back, so
  `coverage == N_meth + N_unmeth` did not hold on the `.cache/normalized`
  store. `N_unmeth` is now the scaled count. `N_meth` and `coverage` are
  unchanged, so every per-CpG engine (which reads only those two columns)
  returns exactly what it returned before. The one reader of the store's
  `N_unmeth` column is `pp.aggregate_regions` (and `epykit
  aggregate-regions`): on a normalised store it summed the stale column, so
  region-level `N_unmeth` and `coverage` were wrong, and so was any
  `tl.dmc` run on that regions store. Pipelines that call
  `pp.normalize_coverage` and then `pp.aggregate_regions` are affected; the
  tile DMR caller and the AnnData `N_unmeth` layer derive the count as
  `coverage - N_meth` and were correct before.
- **`tl.dmc` contrast path skipped two checks.** With `formula=` /
  `contrast=`, an unknown `power_stack` value and `materialize=False` were
  silently ignored. The contrast path now raises the same `ValueError` as
  the binary path for an unknown `power_stack`, refuses `materialize=False`
  (that path always assembles the full result onto `md.varm`), and logs one
  INFO line when a valid `power_stack` is ignored because the GLM has no
  lr+ knobs. Valid calls are unchanged.

### Changed

- **Docs and repository layout.** README, `CLAUDE.md` and the docs no longer claim
  that every DMR caller supports permutation FDR (only `method="tile"` does), the
  engine count is stated as four, the CLI `dmr` default is documented as
  `chain_merge`, and `docs/advanced/architecture.md` is the canonical engine map
  (line-number citations replaced by function names). Planning residue moved to
  `docs/history/`; the root `samplesheet.csv` moved to
  `docs/getting-started/samplesheet.example.csv`. Added `CONTRIBUTING.md` and
  GitHub issue/PR templates. Ticket-named test files (`test_p0_*`,
  `test_phase3_drops`) were folded into behaviour-named files; no test was
  removed.
- **`polars>=1.0` required (M-PKG1).** The declared floor was `>=0.20.0`, but
  the code uses `DataFrame.pivot(on=...)`, a polars-1.0 API. `numba` floor
  raised to `>=0.60` (numpy-2 support).
- **`uv.lock` is now committed (M-PKG4)** for reproducibility; documented
  `uv sync --frozen` and thread-pinning in `benchmark/README.md`.
- **Python >= 3.10 required; 3.13 added to CI.** Python 3.9 reached end of
  life in October 2025. `requires-python` moves to `>=3.10`, the CI matrix
  becomes `{ubuntu, windows} × {py3.10, py3.12, py3.13}`, and ruff targets
  `py310`. No source change is needed; the floor makes `X | None` unions
  usable outside annotations and lets future dependency floors track
  releases that already dropped 3.9 (numpy 2.1+, scipy 1.14+, statsmodels
  0.15, pandas 3, anndata 0.12+).
- **The `dev` extra is gone; contributor tooling lives in dependency groups.**
  `pytest` and `pytest-cov` moved from `[project.optional-dependencies].dev`
  into `[dependency-groups].dev` next to `mypy` and `ruff`, and a new `docs`
  group pins the mkdocs toolchain that was previously installed ad hoc via
  `uv run --with`. `uv sync` installs `dev` by default; pip users need
  `pip install -e . --group dev` (pip >= 25.1) instead of `pip install -e
  ".[dev]"`. The user-facing extras (`all`, `report`, `export`, ...) are
  unchanged. CI installs with `--locked`, checks `uv lock --check`, and pins
  the uv release; `.python-version` pins the interpreter to 3.12.
- **`epykit dmc` / `epykit dmr --method tile` default to `union` (D10/`--unite`).**
  Bare CLI DMC previously intersected sites while bare `ep.tl.dmc` unioned them
  (intersect only after an explicit `ep.pp.unite`). The CLI now defaults to
  union to match the API; `--unite` forces intersect. The benchmark scripts set
  intersect explicitly via the API, so published numbers are unaffected.
  Bare `epykit convert` now merges CpG dyads by default, matching the API
  `merge_strands=True` (D10).
- **Deprecation schedule: the 1.1 removals move to 1.2.** The transitional
  `log2_odds_ratio` column and the `epykit.dmr_hmm` import shim were
  announced as "removed in 1.1"; both are retained in 1.1.0 and their
  warnings now say 1.2. Removing a deprecated surface is a behaviour change
  and gets its own release with the owner's sign-off rather than riding on a
  version bump. The `pp.unite()` alias (scheduled for 2.0) and the `csv*`
  keyword aliases (scheduled for a future release) are also retained; their
  schedules are unchanged.

### Documentation

- Corrected overstated method-fidelity claims: `call_dmr_chain_merge` is
  DSS-callDMR-*style* (not a faithful reimplementation — it omits DSS smoothing
  and uses a non-DSS region statistic; M-DMR4); `empirical_fdr_for_dmc` is a
  Westfall-Young min-P (FWER) procedure, not a pooled-null FDR (M-STAT5);
  `contamination_estimate` is an intermediate-β fraction confounded by cell-type
  heterogeneity / ASM / imprinting / CNV, not validated contamination (M-QC2).

### Testing

- Real-engine null-calibration test (`tests/test_null_calibration.py`) asserts
  `lr` / `glm` / `welch_t` are not anti-conservative under the null; a new CI
  job runs the `slow` tier. Added coordinate-convention and CLI/API parity
  regression tests.
- Overdispersed FPR coverage extended to `welch_t` and `glm` (D18): a fast
  beta-binomial `welch_t` FPR test across φ∈{1.5,2,3,5} and a slow
  store-backed glm/welch_t null at elevated `replicate_sd`. Both come out
  not anti-conservative. Added CLI↔API parity regression tests for
  `min_cpgs`, `min_mean_qvalue`, the `--unite` default, strand-free dyad
  merging, the methylKit coordinate, the GLM `DF_PHI_FLOOR`, the `welch_t`
  CI t-quantile, overlapping-BED region assignment, and per-chromosome
  export streaming.

## [1.0.0] — 2026-06-02

First stable release. API contract is now SemVer-stable. Three targeted
breaking changes land at the major-version cutover; each ships with a
deprecation shim so 0.7.6 code continues to run with warnings.

### Breaking

- **`tl.dmc(power_stack="auto")` now engages the full lr+ stack at any
  sample size** (was: only flipped two of four knobs, and only at
  `min_n <= 2`). Specifically, when `power_stack` is `"auto"` or
  `"lr+"` (new alias), the function flips `neighbour_combine`,
  `fdr_method` (`"fdr_bh"` → `"fdr_tsbh"`), and `sep_fallback` to
  `True` regardless of `n`. The pre-1.0 conservative behavior is
  preserved under `power_stack="conservative"`. `power_stack="off"`
  (new alias) or `False` leaves knobs at user-passed values.
  `power_stack=True` aliases `"lr+"`. Unknown strings raise
  `ValueError`. **The DEFAULT value of `power_stack` remains `"off"`** —
  bare `tl.dmc(test="lr")` produces bare-engine output, no change vs.
  0.7.6. Users who explicitly set `power_stack="auto"` on data with
  `min_n > 2` will see different qvalues vs. 0.7.6.

- **`process_chromosomes_dmc`, `apply_multiple_testing_correction`,
  `empirical_fdr_for_dmc`, `fisher_exact_vectorized`, `shrink_meth_diff`
  removed from top-level `epykit.*` namespace.** Use the recommended
  `tl.dmc` wrapper, or import explicitly via `from epykit.dmc import ...`.
  A module-level `__getattr__` shim accepts the old top-level access
  pattern for 1.0 with a `DeprecationWarning`; removed in 1.2.

- **`pp.unite()` renamed to `pp.set_unite_type()`.** The old name
  suggested a verb performing a union, but the function only writes
  `md.uns["unite"]` (lazy state-marker). `pp.unite()` continues to work
  as a deprecation wrapper through 1.x; removed in 2.0.

- **`method="hmm"` removed from `tl.dmr` and the CLI `dmr` subcommand.**
  Was deprecated in 0.7.5 with `FutureWarning` ("removal in 0.8"); now
  raises `ValueError`. Use `method="segment"` (same engine, honest
  name).

- **CLI `dmr --method` default changed from `"tile"` to `"chain_merge"`**,
  matching the `tl.dmr` library default (set in 0.7.2) so the API and CLI
  converge. `chain_merge` consumes a precomputed DMC parquet via
  `--dmc-results`; scripts that relied on the implicit tile path must now
  pass `--method tile` explicitly.

### Added

- **CLI `chain_merge` DMR surface.** `epykit dmr --method chain_merge
  --dmc-results <dmc.parquet> [--preset strict|default|permissive]
  [--dis-merge-bp N] [--pct-sig F] [--minlen-bp N] [--use-q-for-sig]`
  exposes the DSS-style chain-merge caller, which was previously
  library-only — closing the `tl.dmr` / CLI parity gap.

- **`MethylData.analysis_root`** (no leading underscore) is the new
  public name for the analysis-root attribute. `MethylData._analysis_root`
  continues to work as a deprecated property alias on read AND write,
  emitting `DeprecationWarning`. Removed in 2.0.
- **`Literal[...]` type annotations** on the `value` kwarg of
  `export.to_bedgraph` and `export.to_bigwig` for IDE + mypy support.
  Runtime validation behavior unchanged.
- **Public docstring for `process_chromosomes_dmc(dispersion=...)`** now
  documents the `"eb"` option (already the default in `tl.dmc`).
- **`apply_multiple_testing_correction` docstring** now names the
  `pl.DataFrame` vs `DMCStore` code paths explicitly.

### Changed

- **README and CLAUDE.md describe `lr+` as an opt-in power stack** rather
  than as recommended defaults. Bare `lr` is the default; users opt in to
  `lr+` via `power_stack="lr+"`.
- **Project status:** classifier `Development Status :: 4 - Beta` →
  `Development Status :: 5 - Production/Stable`.

### Fixed

- **`MethylData.save()` preserves `pvalue_combined` / `qvalue_combined`
  columns** added by `neighbour_combine=True`. The DMCStore-backed
  save path hardlinks per-chromosome Parquet files, but the combined
  columns are added in-memory *after* those files were written, so the
  hardlink path was silently dropping them. `save()` now detects the
  `pvalue_combined` column on the in-memory frame and falls back to
  the single-file Parquet write so the four combined / audit columns
  survive a save/load round-trip. Regression test
  `test_save_load_preserves_neighbour_combine_columns`.

### Internal

- Public-surface audit committed at
  `docs/history/superpowers/specs/2026-06-01-public-surface-audit.md` —
  inventory of all 47 top-level exports with per-export verdicts and a
  1.1 backlog.
- Phase 4 plan checkboxes synced to reflect what actually shipped
  (Tasks 1-8 marked).

## Unreleased

## [0.7.6] — 2026-06-01

Cache-semantics correctness fix + benchmark null-calibration fixes.
Rolls up all Phase 2-3-4 work since 0.7.2 into a single release. Engine
math itself is unchanged from the `v0.7.5-phase3-engines-frozen` tag —
the only library-side change in this release is to the per-CpG DMC
cache invalidation logic. Benchmark-side changes wire the `glm` engine
correctly into the null-calibration framework and fix a sample-label
shuffling bug that had been making every null-calibration shuffle run
on the original case/control assignment.

### Fixed

- **`process_chromosomes_dmc`: weak-hit cache must not fire when
  `input_sig` differs.** Previously, when the cached manifest carried an
  `input_sig` field that did not match the current call, the function
  served the stale per-chromosome parquets and silently upgraded the
  manifest's `input_sig` in place. The weak-hit was intended to recover
  gracefully from legacy manifests written before the `input_sig` field
  existed, but the implementation also masked real input changes —
  affecting any user calling `tl.dmc` repeatedly on the same methylstore
  with different sample labels, test parameters, or covariates. The fix
  restricts the weak-hit to its documented use case (manifests with no
  `input_sig` at all) and falls through to the recompute path when
  `input_sig` is present but differs. New regression tests in
  `tests/test_dmc_streaming_store.py` cover both the recompute path and
  the legacy-manifest path. Discovered via Phase 4 null-calibration
  development; benefits all repeated-DMC use cases (cross-validation,
  sensitivity analyses, custom permutation sweeps).

- **`benchmark/scripts/_null_engines.py:_permute_md`: write
  `obs["treatment"]`, not just `obs["group"]`.** `MethylData.treatment_ids`
  reads from `obs["treatment"]` (i64), so writing only the display string
  column `obs["group"]` silently left every null-calibration shuffle
  running on the original case/control assignment. All shuffles in every
  Phase 4 cell produced identical p- and q-values; the IQR collapsed to
  the median. Fix writes both columns; the slow regression test
  `test_closure_produces_different_qvalues_per_shuffle` asserts that two
  closure invocations with different label assignments produce different
  q-value arrays.

- **`benchmark/scripts/_null_engines.py`: glm closure builds the design
  matrix.** Previously the glm registry entry called
  `ep.tl.dmc(md, test="glm")` with no `formula` or `contrast` argument,
  hitting the binary-path branch that requires explicit `design_full` /
  `design_reduced` / `coef_idx`. The fix routes glm through
  `formula="~ group", contrast="group"`, letting `build_design`
  construct the full and reduced designs from the permuted `obs.group`
  column. The glm engine now produces real per-shuffle variation in the
  null-calibration sweep.

### Added

- **`benchmark/scripts/run_phase4_null_calibration.py`** (new):
  orchestrator for the Phase 4 null-calibration sweep across every
  surviving engine × dataset × scenario cell. Builds per-dataset
  methylstores, dispatches per-cell calibration via
  `run_null_calibration.py`, aggregates per-cell parquets into a single
  `summary.parquet` with `(engine, dataset, scenario,
  observed_fdr_median, observed_fdr_q1, observed_fdr_q3,
  observed_fdr_ci_lo, observed_fdr_ci_hi)` columns. CLI accepts
  `--dataset`, `--engines`, `--skip-piao/--simulator/--gse`,
  `--keep-store`, `--skip-aggregate`, `--only-aggregate`. Companion test
  suite at `benchmark/scripts/tests/test_run_phase4_null_calibration.py`
  covers the spec registry, aggregation kernel, manifest writer, and
  per-engine parquet path layout. Final sweep produced 12 of 13 cells
  (fisher@gse263850 deferred pending parallel backend).

- **`benchmark/scripts/bug_fix_audit.py`** (new): diffs pre/post-fix
  `eval_summary.parquet` per (tool, scenario, metric); attributes each
  delta to a P0/P1 fix via `Affects: engine@scenario` trailers parsed
  from commits JSON. Unattributed changed cells cause non-zero exit.

- **`benchmark/scripts/regen_all.py`** (new): `--verify` acceptance
  gate reads `claims.yaml` and `<!-- claim: id -->` HTML comments in
  `paper.md`, asserts cited parquets match expected values to stated
  precision, exits non-zero on mismatch. Empty `claims.yaml` seed
  landed; Phase 4 populates during the locked re-run.

- **`benchmark/scripts/methylkit_stouffer_combine.R`** (new):
  adjacent-CpG Stouffer combination for methylKit output. Mirrors
  epykit's `neighbour_combine` knob for tuned-vs-tuned Phase 4
  comparisons per PROTOCOL R1. Test SKIPS without Rscript on PATH.

- **Benchmark consolidation**: the canonical benchmark suite is now
  `epykit3/benchmark/` (was previously split between `epykit3/benchmark/`
  and `benchmarkin_merges/FINAL_REPORT/`, the latter outside version
  control). One repo, one tag, reproducible via
  `uv run pytest benchmark/scripts/tests/` for the new script tests.
  Data (raw inputs, ground truth, eval_summary, per-scenario methylkit
  results) are intentionally NOT bundled — regenerated by re-running
  the scripts against the user's local `raw_sim_data/`.

- **`benchmark/scripts/simulate_piao.py`** (new): Python re-implementation
  of Piao et al. 2021's binomial DMC simulator with an intrinsic
  ``is_dmc`` flag (replaces the threshold-reconstructed truth that
  caused AUROC tautology). Outputs Piao-compatible AMP files +
  ``truth.parquet`` matching ``dmc_truth.parquet`` schema. Baseline
  beta uses ``Beta(0.75, 1.35)`` fit to Piao's marginals (1%/3% match).

- **`benchmark/scripts/wilson_bootstrap_ci.py`** (new): Wilson 95% CIs
  for proportion metrics (TPR/FPR) and percentile bootstrap CIs for
  rank/threshold metrics (AUROC, F1). Operates on the existing
  ``eval_summary.parquet`` schema without re-running engines.

- **`benchmark/scripts/evaluate.py`** (new): `--ci-only` mode appends
  Wilson 95% CIs (TPR/FPR) and bootstrap CIs (AUROC/F1, NaN without
  per-CpG cache) to `eval_summary.parquet` in place. Wires the Phase 2
  `wilson_bootstrap_ci.py` helper into the eval pipeline.

- **`benchmark/scripts/run_null_calibration.py`** (new): Label-shuffle
  empirical FDR runner. Decoupled from epykit (takes an
  ``engine_fn`` callable) so it tests with mock engines for CI and
  wraps ``ep.tl.dmc`` for real use.

- **`benchmark/scripts/_null_engines.py`** (new): real-engine
  closures (`lr`, `lr_plus`, `welch_t`, `fisher`, `glm`) for
  `run_null_calibration.py`. Each factory captures a `MethylData`
  object once and returns a closure matching the Phase 2
  `engine_fn(samples_treatment, samples_control, seed)→qvals` contract.
  `run_null_calibration.py main()` rewritten to dispatch via the
  registry with `--engine` + `--methylstore` args.

### Removed

- **DMC engine surface collapsed to `{auto, lr, welch_t, fisher, glm}`.** 
  Dropped: `logit_t` (broken near β=0/1), `bb_lr` (TPR < 8% at n ≤ 4 + 
  dispersion bug), `score` (dominated by `lr`), `cmh` (dominated by 
  `glm + batch`). All four raise `ValueError` with a one-line migration 
  hint. Removed from README, docs, CLAUDE.md, and CLI help.

- **`tl.dmc(test='logit_t')`** removed. The engine was documented by
  epykit's own source as miscalibrated near β=0/1; no paper claim
  depends on it. Calls now raise `ValueError` with a migration hint.
  Migration: `test='logit_t'` → `test='welch_t'` or `test='lr'`.

- **`tl.dmc(test='bb_lr')`** removed. TPR < 8% at n ≤ 4 in the
  published benchmark; also affected by the P1-2 dispersion-df bug
  (now closed by removal). Calls now raise `ValueError` with a
  migration hint. Migration: `test='bb_lr'` → `test='lr'`. Note:
  `irls_dispatch` stays (used by the `glm` engine).

- **`tl.dmc(test='score')`** removed. Strictly dominated by `lr` in
  finite samples; asymptotically equivalent under H₀. Migration:
  `test='score'` → `test='lr'`; output schema identical.

- **`tl.dmc(test='cmh')`** removed. Stratification semantics (one
  stratum per case-ctrl pair) were unusual; dominated by
  `tl.dmc(formula='~ group + batch')`. Migration: use the formula
  kwarg for stratified analysis. `_cmh_init/_cmh_update/_cmh_finalize`
  helpers deleted (no other callers).

### Changed (BREAKING for `log2_odds_ratio` column name)

- **`varm["dmc_lr"].log2_odds_ratio`** renamed to
  `log2_odds_ratio_pooled` (same value, clearer name). Same rename
  applies to `fisher`.
- **`varm["dmc_glm"].log2_odds_ratio`** renamed to
  `coef_treatment_log2` (it was always the logit coefficient in log₂
  units, not log₂ of an odds ratio; the old name was misleading).
- A transitional `log2_odds_ratio` column is NaN-filled in 0.7.5 with
  a `FutureWarning` on the producing call. Column removed in 0.8.

### Changed (BREAKING for `epykit.dmr_hmm` import path)

- **Renamed** `epykit.dmr_hmm` → `epykit.dmr_segment`; function
  `call_dmr_hmm` → `call_dmr_rule_segment`. The engine uses fixed
  state means / transition priors (not Baum-Welch fitted), so calling
  it an HMM was misleading. Old import path remains as a deprecated
  shim until 0.8. `tl.dmr(method='hmm')` works with `FutureWarning`.

### Fixed (P2 manifest, folded into the rename)

- **P2-4**: `call_dmr_rule_segment` now emits per-segment Stouffer-
  combined p-values (BH-corrected per chromosome) instead of NaN.
  The pre-0.7.5 implementation emitted NaN p/q-values for every
  called segment, breaking any downstream filter on qvalue.

### Fixed (P1 manifest)

- **P1-9**: `sex_check` now gates the 1D largest-gap clustering on
  Hartigan's dip test (`diptest.diptest`). On unimodal distributions
  (single-sex cohorts, dip p > 0.10), falls back to a fixed chrX-beta
  threshold (0.25) and emits `UserWarning`. Bimodal mixed-sex cohorts
  are unchanged. Clustering logic extracted into the internal helper
  `_classify_sex_from_values`. Adds `diptest` to the `qc` extra in
  `pyproject.toml`.

- **P1-7**: `tl.dvc` per-site variance test replaced from Bartlett with
  Brown-Forsythe (median-centred Levene). Bartlett assumes normality; beta
  methylation values are bounded [0,1] and U-shaped. Brown-Forsythe matches
  `scipy.stats.levene(center='median')` to 1e-6. `process_chromosomes_dvc`
  default changed to `test='brown_forsythe'`; `test='bartlett'` accepted as
  a backward-compatible alias (silently redirected). DVC not in paper; no
  headline impact.

### Fixed (P0 manifest, paper preparation)

- **P0-3 (docs)**: `tl.dmc(..., dispersion=...)` docstring previously said
  the default was `"site"`; the actual code default is `"eb"`. The
  empirical-Bayes shrinkage is the intended default. Docstring corrected;
  PROTOCOL.md and EXECUTIVE_SUMMARY downstream notes follow separately. See
  `docs/history/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md`.

- **P0-6 (docs)**: `combine_neighbour_pvalues` docstring now explicitly
  owns the Stouffer independence assumption violation under spatial
  autocorrelation. The 0.7.x FDR safety net is the sign-agreement gate,
  not the Stouffer null; Brown's-method replacement deferred to v0.8.

- **P0-5 (`dmr._merge_adjacent_tiles`)**: three Stouffer-combination bugs
  fixed: (a) input p-values are two-sided so the magnitude conversion is
  ``isf(p/2)``, not ``isf(p)``; (b) chains of length > 2 now use the
  correct ``sum_z / sqrt(n)`` denominator (was iterative pairwise
  ``/sqrt(2)`` which over-conserves on long chains); (c) running
  ``(sum_z, n)`` accumulator added. Affects tile-based DMR p-values
  in any run where `merge_adjacent=True` (the default for
  `call_dmr_tile_based`).

- **P0-2 (`dmr.empirical_fdr_for_dmr`)**: empirical p-value denominator
  changed from ``|pooled null DMRs| + 1`` to ``n_perm + 1`` (per-permutation
  tail count). The old pooled-null formula was anti-conservative -- more
  permutations grew the denominator and shrank emp_p without adding
  evidence. Affects any DMR run that used ``empirical_fdr=True``;
  empirical p-values will rise (more conservative).

- **P0-2b (`dmc.empirical_fdr_for_dmc`)**: DMC analogue of P0-2.
  Empirical p-value denominator changed from
  ``|pooled null p-values| + 1`` to ``n_perm + 1`` (per-permutation
  tail count). The old pooled-null formula was anti-conservative.
  Affects any DMC run that used ``empirical_fdr=True``; empirical
  p-values will rise (more conservative).

- **P0-4 (`dmc._score_finalize` adaptive-F branch)**: ``df_phi`` is now
  floored at 50 when the F branch fires (``phi_eff > 1.0``). Previously, in
  ``dispersion="eb"`` mode with small empirical-Bayes weight (the
  homogeneous-dispersion case), ``df_phi`` collapsed to ~4 and F(1, 4) was
  ~250x more conservative than chi^2(1) at typical test statistics --
  this drove the artifactually low FPR in eb mode in 0.7.2. The
  chi^2 branch (clamped phi) is unaffected. Expected impact on
  benchmark numbers: ``lr / eb`` FPR rises modestly, TPR essentially
  unchanged. Null-calibration smoke (n=1934, n_per_group=3, no true DMCs):
  observed FDR at q<0.05 = 0.0000 (BH conservative at low n, as expected).
  If you benchmarked under 0.7.2 with `dispersion="eb"` (the default), the previously reported FPR is artificially suppressed and should be recomputed.
  See the bug-fix audit table in Phase 2 for per-cell deltas.

### Fixed (P1 manifest)

- **P1-1**: `fisher_exact_vectorized` two-sided p now uses the mid-p
  convention (sum of hypergeometric pmf over all tables with pmf ≤
  pmf(observed)), matching `scipy.stats.fisher_exact(alternative=
  'two-sided')` to 1e-12. Previously doubled the smaller one-sided
  tail. Affects small-table cells; headline cov≥10/n≥3 unchanged.

- **P1-2** (`bb_lr` `df_resid` vs `df_phi`): **closed by removal** of
  `test='bb_lr'`. The engine incorrectly discarded `df_phi` from
  `compute_dispersion_phi` and passed `df_resid_safe` to
  `reference_pvalues` instead, causing miscalibrated F-distributed
  p-values under `dispersion != "site"`. Removed in 0.7.5 together
  with the broader TPR < 8% finding. Migration: `test='lr'`.

- **P1-3**: `lr` and `fisher` engines now emit Newcombe (1998) hybrid
  Wilson-score CIs for `meth_diff_ci_{lo,hi}`, matching
  `statsmodels.confint_proportions_2indep(method='newcombe')`.
  Previously emitted Welch-normal Wald CIs (symmetric, can violate
  [-1,1] near boundary β). `welch_t` and `glm` CIs unchanged.

- **P1-4**: `tl.dmc(formula=..., reference_level='level')` lets users
  set the categorical reference level explicitly via patsy's
  `Treatment(reference=...)`. Default behaviour unchanged (alphabetical
  reference). `_glm.build_design` now logs resolved column names and
  the chosen reference at INFO level.

- **P1-5**: IRLS non-convergence is now surfaced. Non-converged sites
  have their Wald statistics NaN-masked (previously only `separated`
  sites were NaN'd); fraction non-converged > 1% is logged at WARNING.
  Fraction ≤ 1% is logged at INFO. Separated sites that are also
  non-converged are counted only once (under the separation diagnostic).

- **P1-6**: DMR empirical FDR now (a) raises `ValueError` at
  n_treat=1, n_ctrl=1 (no valid permutations), pointing the user at
  `tl.dmc(test='fisher')`; (b) accepts `empirical_strata=<obs_column>`
  to permute within strata for paired designs.

- **P1-10**: `_storey_pi0` now clamps at `1/n` from below (Storey
  2002 standard floor). Previously could return 0 when all p-values
  fell below `lam`, causing +inf q-values via BH adjustment. Docstring
  updated to note this is the plug-in estimator, not the spline-
  smoother variant.

### Changed (breaking on the `lr+` schema)

- **P0-1**: ``tl.dmc(..., neighbour_combine=True)`` no longer overwrites
  the per-CpG ``pvalue`` column with the Stouffer-combined value. The
  raw per-CpG p-value stays in ``pvalue``; the combined value is in
  ``pvalue_combined`` and its BH q-value in ``qvalue_combined``.
  Downstream code reading ``pvalue`` (BH correction, empirical FDR,
  DMR engines) now sees the raw per-CpG value, not the combined.
  If you want the combined q-value, read ``qvalue_combined``.
  ``empirical_fdr_for_dmc`` raises ``ValueError`` if it sees a stale
  schema where ``pvalue`` was overwritten.

  Affects: anyone who consumed ``md.dmc["pvalue"]`` after
  ``neighbour_combine=True``. Migration: switch reads to
  ``md.dmc["pvalue_combined"]`` / ``md.dmc["qvalue_combined"]``.

## [0.7.2] — 2026-05-21

Eight targeted fixes identified by a full benchmark comparison against
methylKit, DSS, RADMeth, BiSeq, methylSig, and Fisher (Piao et al.
2021). Two bugs fixed, three default changes, three guardrail
improvements. Existing user code that passes explicit parameters is
unaffected; only implicit defaults change.

### Fixed

- **Fisher calibration bug (`dmc.py`).** `fisher_exact_vectorized`
  computed only the upper tail via `hypergeom.sf(meth_a - 1, ...)`.
  When group A had less methylation (hypo direction), `meth_a - 1`
  could be -1, making `sf(-1, ...)` return 1.0. Fixed to compute both
  tails and pick the correct one based on whether observed >= expected.
  Fisher TPR went from 0.000 to 0.668-0.998 across all coverage levels.

- **chain_merge DMR fragmentation (`dmr.py`).** Default `dis_merge_bp=100`
  was too narrow for real genomic CpG spacing (median 100-200 bp in
  dense regions, 500-1000 bp in intergenic/intronic). DMRs spanning
  500-2000 bp with CpGs >100 bp apart were fragmented into 4-5 pieces.
  Widened defaults: `strict` 100->250, `default` 100->500, `permissive`
  200->1000, function default 100->500, `_SIG_DEFAULTS` 100->500.
  chain_merge TPR went from 0.086 to 0.971-1.000.

### Added

- **`power_stack` parameter on `ep.tl.dmc` (default `False`).**
  Pass `power_stack=True` to enable `neighbour_combine=True` and
  `sep_fallback=True` (the lr+ stack) in one switch. Pass
  `power_stack="auto"` to auto-enable the stack at n <= 2 replicates
  per group. Default is `False` — no silent auto-engagement.

- **Adjacent tile merging in `ep.tl.dmr(method="tile")`.**
  `merge_adjacent` parameter (default `True`) on `call_dmr_tile_based`
  merges adjacent significant tiles on the same chromosome with the same
  direction via Stouffer Z-combination. P-values are re-corrected with
  BH after merging. Pass `merge_adjacent=False` for exact old behaviour.

- **bb_lr auto-shrink guardrail (`dmc.py`).** When `n_samples < 6` and
  `dispersion="site"`, bb_lr auto-promotes to `dispersion="shrink"` with
  a warning, since per-site dispersion from Pearson residuals with
  `df_resid <= 4` is extremely noisy.

- **`dispersion="eb"` support in bb_lr (`_glm.py`).** The bb_lr path
  now accepts `dispersion="eb"` (empirical-Bayes shrinkage), matching
  the lr engine. Previously bb_lr only accepted `site/chrom/shrink`
  and raised `ValueError` on `"eb"`.

### Changed

- **`dispersion` default on `ep.tl.dmc` changed from `"site"` to
  `"eb"`.** Empirical-Bayes shrinkage of per-site quasi-binomial
  dispersion toward a chromosome-wide inverse-Gamma prior. At high n
  (large per-site df), the weight on the per-site estimate dominates
  and `"eb"` reduces to `"site"`. At low n, `"eb"` shrinks toward the
  chromosome mean, stabilising noisy dispersion estimates without losing
  site-level resolution. The change improves power at n <= 3 per group
  on real WGBS data with heterogeneous overdispersion; on underdispersed
  simulation data (Piao et al. 2021) it is a no-op (per-site phi is
  clamped at 1.0 either way). Pass `dispersion="site"` to restore
  exact pre-0.7.2 behaviour.

- **`method` default on `ep.tl.dmr` changed from `"tile"` to
  `"chain_merge"`.** On the Piao et al. 2021 DMR simulation,
  chain_merge recovers 97-100% of truth DMRs at all coverage levels
  while tile recovers 54-100%. chain_merge is also more adaptive to
  irregular CpG spacing. The tile engine remains available via
  `method="tile"`. Pass `method="tile"` to restore pre-0.7.2 behaviour.

- **welch_t warning tiers (`dmc.py`).** Split `_validate_sample_size_and_warn`
  into severity tiers: CRITICAL at `min_n <= 2` (degenerate
  Welch-Satterthwaite DOF, near-zero power), softer warning at `min_n < 6`.

- **bb_lr low-n warning (`dmc.py`).** Added specific warning when
  `test="bb_lr"` and `min_n < 3`: recommends `test="lr"` instead.

- **`dis_merge_bp` default in `ep.tl.dmr`** changed from 100 to 500 to
  match the updated chain_merge defaults.

- **Docstring tuning guidance in `dmr.py`** updated: "loosen dis_merge_bp
  100 -> 200" changed to "500 -> 1000" to reflect the new defaults.

### Tests

- `test_primitives.py`: `test_fisher_reverse_separation_significant` and
  `test_fisher_symmetry` covering the dual-tail Fisher fix.
- `test_dmr_chain_merge.py`: `test_chain_merge_default_merges_300bp_gap`
  verifying sig CpGs 300 bp apart chain at the new 500 bp default.
- New `test_dmr_tile_merge.py` (6 tests): adjacent tile merging,
  direction-aware non-merging, gap handling, multi-chromosome, three-way
  merge, empty input.

---

## [0.7.1] — 2026-05-21

Targeted improvements to the `lr` DMC engine that close its
asymptotic-quasi-binomial gap to methylKit / RADMeth / DSS at low
coverage and small cohorts. All four are opt-in keyword arguments on
`ep.tl.dmc(test="lr", ...)`; defaults preserve 0.7.0 behaviour
exactly.

### Added

- **`dispersion="eb"` on `ep.tl.dmc`.** Empirical-Bayes shrinkage of
  per-site quasi-binomial dispersion toward a chromosome-wide
  inverse-Gamma prior whose pseudo-df is estimated from the
  per-site phi distribution via method-of-moments. Generalises the
  existing `dispersion="shrink"` mode (which uses a fixed
  pseudo-df = 4). No-op on data without genuine overdispersion.

- **`neighbour_combine=True`, `neighbour_bp=200` on `ep.tl.dmc`.**
  Signed-Stouffer Z combiner over neighbouring CpGs. Gated by
  `min_sign_agreement=0.6` (focal site's neighbours must agree on
  direction by ≥ 60 %) and `require_focal_signal=True` (focal raw
  p must be < `focal_p_thresh=0.5`) so spatially isolated false
  positives are not amplified. The output `pvalue` becomes the
  combined p; the raw is preserved as `pvalue_raw`. Two helper
  columns `pvalue_combined` and `pvalue_combined_n_neighbours`
  are added for audit. Exposed standalone as
  `ep.dmc.combine_neighbour_pvalues(dmc_df, ...)`.

- **`sep_fallback=True`, `sep_threshold=0.9` on `ep.tl.dmc`.**
  Separation-aware Fisher fallback inside `_score_finalize`: for
  sites where `|meth_diff| >= sep_threshold` AND the LR p-value
  failed to reject (p > 0.05), re-test with
  `scipy.stats.fisher_exact` on pooled counts and take the more
  powerful of the two. Never inflates p. Affects only sites the
  LR missed, so the overall FPR is unchanged.

- **`fdr_method="fdr_tsbh"` / `"fdr_storey"` on `ep.tl.dmc`.**
  Selects the FDR procedure run after the per-CpG test. New
  options: `"fdr_tsbh"` (Benjamini-Krieger-Yekutieli two-stage BH;
  statsmodels' adaptive variant), `"fdr_storey"` (Storey-Tibshirani
  q-values with π₀ estimated at lam = 0.5). Defaults to
  `"fdr_bh"` for back-compat. Threaded through
  `ep.apply_multiple_testing_correction` and the streaming
  `DMCStore` path; the manifest now records the method used.

- **`ep.dmc._storey_pi0`, `ep.dmc._apply_storey_qvalues`** helpers
  surfaced for re-use from custom pipelines.

- **Benchmark and write-up.** End-to-end reproduction of the Piao
  et al. 2021 (IJERPH 18:7975) simulated benchmark under
  `benchmark/`, including ground-truth
  reconstruction (35 reference DMRs / 19,999 true DMCs, both
  matching the paper's design exactly), baseline-table
  transcription, 7 figures, and a paper-style manuscript.

### Changed

- `_dmc_input_signature` now includes `sep_fallback` and
  `sep_threshold` in the SHA-256 fingerprint, so toggling them
  forces a recompute rather than serving a stale cache.

- `DMCStore.mark_bh_applied(method=...)` records the FDR method
  used so back-to-back calls with different `fdr_method` do not
  silently reuse the cached q-values.

### Tests

- New `tests/test_lr_improvements.py` (10 tests) covering Storey
  π₀ behaviour at known mixtures, fdr_method validation, the
  neighbour combiner's sign-agreement guard, the
  never-inflates-p property, and NaN preservation.

### Bug surfaced (fixed in 0.7.2)

- The pooled `fisher` backend in v0.7.0 returns `pvalue ≈ 1.0` on
  near-perfect-separation 2 × 2 tables. Fixed in 0.7.2 with a
  dual-tail hypergeometric computation.

## [0.7.0] — 2026-05-21

Adds a DSS-compatible DMR caller, DSS-style raw-count smoothing for DMC,
annotatr-style multi-annotation (nearest TSS + all overlapping genes /
features), and a UCSC `refGene.txt` parser. Existing 0.6.x defaults are
unchanged — every new capability is opt-in via an explicit kwarg.

### Added

#### DMR
- **`ep.call_dmr_chain_merge` / `ep.tl.dmr(method="chain_merge")`.** DSS
  `callDMR` semantics on top of an epykit DMC table: mark sig CpGs by
  `(pvalue < alpha) AND (|meth_diff| >= min_abs_meth_diff)`, chain
  contiguous sig CpGs within `dis_merge_bp`, then filter by `min_cpgs`,
  `pct_sig`, and `minlen_bp`. `use_q_for_sig=True` switches the gate to
  the q-value column when one is present. Cached per chrom by the same
  DMC `input_sig` fingerprint as the sliding-window caller.
- **`ep.DMR_PRESETS` parameter bundles** for chain-merge: `"strict"`
  (alpha=1e-6, min_cpgs=5, |Δβ|≥0.20 — validation-ready), `"default"`
  (alpha=1e-4, |Δβ|≥0.10, min_cpgs=3 — balanced; one order looser than
  DSS to capture real-but-moderate signal without crashing PPV), and
  `"permissive"` (alpha=1e-4, dis_merge_bp=200, |Δβ|≥0.05 — recall-
  oriented). Surfaced through `ep.tl.dmr(..., preset="...")`; any
  explicit kwarg overrides the bundle value.
- **`ep.tl.diagnose_dmr_calling(md, reference_dmrs, ...)`.** Bucket
  reference DMRs into actionable categories — `SUCCESS_OVERLAP`,
  `H1_NO_CPGS` (lost in coverage / unite), `H2_NO_SIG_CPGS` (test
  too conservative), `H3a_WEAK_ALPHA` (loosen `alpha`), `H3b_STRUCTURE`
  (loosen `dis_merge_bp` / `min_cpgs`) — so a low-recall number maps to
  a specific knob instead of guesswork.

#### DMC
- **DSS-style raw-count smoothing** in `ep.tl.dmc(..., smoothing=True,
  smoothing_span_bp=500)`. Applies DSS's uniform-box ±`span/2`
  moving-average to each sample's per-CpG `(meth, cov)` before the test
  hits them — matches `DMLfit.multiFactor(smoothing=TRUE)` semantics.
  Cached separately from the un-smoothed run (`dmc/<test>_smooth/` vs
  `dmc/<test>/`) and signature-versioned so `False → True` correctly
  invalidates a stale cache.

#### Annotation
- **annotatr-style multi-annotation** in `ep.annotate_features` (default
  on, also `ep.tl.annotate(..., multi_annotation=True)`). Adds
  `nearest_tss_gene` / `nearest_tss_distance` (HOMER-style signed
  distance, flipped on `-` strand so positive is downstream of the TSS)
  and `all_overlapping_genes` / `all_overlapping_features` (one-to-many,
  so a site that's intronic for one gene AND promoter for another is
  faithfully represented instead of collapsed by the best-pick rule).
- **UCSC `refGene.txt(.gz)` annotation source.** Pass
  `ep.tl.annotate(md, refgene=...)` (or `annotate_features(..., source=
  "refgene")`) to use HOMER's default catalog — curated, protein-coding-
  biased, and gives the highest paper-gene recall on methylation work.
  Schema-compatible with the GTF path; same downstream consumers.
- **`gene_type_filter` kwarg.** Restrict the gene catalog before
  building overlap intervals (`"protein_coding"` to drop lincRNAs /
  pseudogenes). Works on both GTF (`gene_type` / `gene_biotype`) and
  refGene (derived from accession prefix).
- **GENCODE / Ensembl GTF gene-type parity.** `_parse_gtf_streaming`
  now accepts both `gene_type` (GENCODE) and `gene_biotype` (Ensembl);
  files that omit it stay annotatable.

### Deprecated
- **`ep.tl.dmc(..., use_smoothed=True)`** (the BSmooth pseudo-count
  transform) now emits a `DeprecationWarning` pointing at the new
  `smoothing=True` path. The pseudo-count path is too aggressive
  (replaces the count signal entirely with the smoothed version, washing
  out per-CpG resolution at default BSmooth parameters) and will be
  removed in a future minor release. Existing callers keep working
  for now.

### Changed
- **Internal cleanup: ASCII docstrings.** Library docstrings, log
  messages, and inline comments have been rewritten in plain ASCII
  (`β → beta`, `μ → mu`, `Σ → Sigma`, `² → ^2`, `→ → ->`, `— → --`,
  `× → x`, `≥ → >=`, `≤ → <=`, `≈ → ~=`). No semantic changes — this
  fixes `argparse --help` crashing with `UnicodeEncodeError` on Windows
  consoles running the default `cp1252` codec and removes the need for
  the CLI's stdout/stderr `utf-8` reconfigure to be load-bearing on
  every path. The CLI still reconfigures to UTF-8 defensively.

### Tests
- New: `test_dmr_chain_merge.py`, `test_dmr_presets_and_diagnose.py`,
  `test_dmc_smooth_dispersion.py`, `test_annotate_multi.py`. Existing
  tests realigned with the ASCII docstring rewrite (no behavioural
  changes).

### Deferred to 0.8+

RefSeq / UCSC functional-element annotation tracks beyond gene models,
DMR caller for mixed designs (random-effects), Zarr storage backing,
single-cell sparse stores.

---

## [0.6.0] — HMM segmentation

One shared HMM engine, three callers on top. All three operate on the
existing methylstore (no new I/O path) and route through the 0.4
chrom-streaming dispatcher, so the distributed backend works for free.
No new dependencies — the HMM is hand-rolled in ~200 LoC of numpy.

### Added

#### Shared HMM engine
- **`epykit._hmm.segment(observations, n_states, ...)`**. Hand-rolled
  forward-backward + Viterbi for either Bernoulli or Gaussian
  emissions, with a sticky-chain transition prior (configurable
  ``self_loop`` and full ``transition_priors`` overrides).
- **`epykit._hmm.runs_of_state(viterbi, target_state, positions)`**
  extracts contiguous runs of a target state from the Viterbi path,
  optionally translated to genomic positions.

#### Callers
- **`ep.tl.pmd(md)`** — partially methylated domains. Per-sample,
  megabase-scale 2-state HMM on coverage-weighted smoothed β.
  Output in ``md.uns["pmd"]``.
- **`ep.tl.hmr(md)`** — hypo- and low-methylated regions
  (MethylSeekR-style). Per-sample 2-state HMM on raw per-CpG β.
  Tagging splits the runs into ``md.uns["hmr"]`` (dense, CpG-island-
  like) and ``md.uns["lmr"]`` (sparse, distal regulatory).
- **`ep.tl.dmr(md, method="hmm")`** — HMM-based DMR caller. Three-state
  Gaussian HMM on the per-CpG ``meth_diff`` signal from any DMC table.
  Schema-compatible with ``method="tile"``, so existing plot / export
  paths work unchanged.

### Deferred to 0.7+

Zarr storage backing, single-cell sparse stores, eQTM, motif/TFBS
enrichment, matrix-completion imputation beyond kNN, pyGenomeTracks
track plot, differential entropy CpG-window test.

---

## [0.5.0] — BAM-based read-level analyses

Adds the first analyses that need read-level methylation information.
Builds on the 0.4 engine lifts but doesn't change any existing default
behaviour. New optional `bam` extra (pysam, Linux/macOS only) gates
the new analyses; everything else stays installable on Windows.

### Added

#### BAM ingestion
- **`epykit.bam_io.read_methylation_calls(bam, ...)`**. Returns a
  long-form polars DataFrame with one row per (read, covered CpG):
  `(read_id, chrom, pos, methylation_status, base_qual, mapq,
  mate_pair_id, strand, allele_base)`. Two BAM dialects: Bismark `XM`
  tags and SAM-standard `MM`/`ML` tags (MethylDackel).
- **`bam` optional extra** in `pyproject.toml` (`pysam>=0.22`).
  Linux/macOS only — pysam has no Windows wheel.

#### Allele-specific methylation (ASM)
- **`ep.tl.asm(md, bam=..., vcf=...)`**. Per-CpG Fisher exact test of
  H1 vs H2 read methylation, with heterozygous SNVs from a VCF as
  phasing anchors. Result lands in `md.varm["asm"]` with columns
  matching the `dmc_*` family (`pvalue`, `qvalue`, `meth_diff`) so
  `pl.volcano(md, key="asm")` works without modification.
- Reuses `fisher_exact_vectorized` and `apply_multiple_testing_correction`
  from the DMC stack.

#### Methylation entropy
- **`ep.tl.entropy(md, bam=..., window_cpgs=4)`**. Per-CpG-window
  Shannon entropy over the observed read methylation patterns. Reads
  that cover all CpGs in the window contribute one binary pattern; the
  full distribution's Shannon entropy is normalised to `[0, 1]`.
  Result lands in `md.varm["entropy"]`.

### Deferred to 0.6

PMD / HMR / LMR callers and HMM-DMR all share an HMM segmentation
engine; they land together in the 0.6 release.

---

## [0.4.0] — engine lifts

Infrastructure release. Default behaviour unchanged (every 0.3 test
passes verbatim); every new capability is reached by an explicit kwarg
or optional extras install. No new analyses — those land in 0.5.

### Added

#### Compute backends
- **Distributed compute via Dask** (`tl.dmc(..., backend="dask", n_workers=4)`,
  `tl.dmr`, `tl.dvc`). Per-chromosome work submitted to a local cluster
  or an existing `dask.distributed.Client`. Requires the new
  `pip install 'epykit[distributed]'` extra. Results are bit-identical
  to the sequential path.
- **Distributed compute via Ray** (`backend="ray"`). Same surface as
  the Dask backend; uses `ray.remote` actors. Requires
  `pip install 'epykit[ray]'`.
- **GPU IRLS via CuPy** (`tl.dmc(test="glm", glm_backend="gpu")`). The
  batched binomial IRLS hot path in `_glm.py` now has a CuPy mirror
  (`_glm_gpu.py`). Requires `pip install 'epykit[gpu]'` (CUDA 12). The
  closed-form `lr` / `score` tests stay CPU-only by design.

#### Pipeline manifest + resume
- **Formal checkpoint / resume API.** Each `MethylData` analysis root
  now hosts a top-level `.epykit_manifest.json` recording completed
  pipeline stages with their input signatures and sidecar parquet
  paths. Call `ep.tl.dmc(md, ..., resumable=True)` twice with the same
  inputs and the second call loads the cached result instead of
  recomputing. `md.completed_stages` reports the recorded list;
  `md.resume_from("dmc_lr")` re-hydrates a fresh `MethylData` from the
  on-disk manifest.

#### Tabix-on-Parquet random access
- **`ep.query` module.** Three entry points —
  `query_region(store, chrom, start, end)`,
  `query_regions(store, regions_df)`, and
  `query_sites(store, sites_df)` — return long-form
  `(sample_id, chrom, pos, strand, N_meth, coverage, beta)` frames for
  arbitrary genomic loci. No new dependency: built on
  `pl.scan_parquet` predicate pushdown over the existing
  hive-partitioned store.

### Internal

- **Per-chromosome compute dispatcher** at
  `src/epykit/_compute.py:run_chrom_pipeline`. The chrom loops in
  `dmc.py`, `dmr.py` (tile), and `dvc.py` were refactored to route
  through this shared dispatcher. `backend="sequential"` (default) is
  bit-identical to the pre-0.4 in-line loop.
- **`irls_dispatch`** in `_glm.py` routes `irls_binomial_batch`
  between CPU (numpy) and GPU (CuPy via `_glm_gpu.py`).

### Deferred to 0.5+

Out of scope for 0.4: ASM, methylation entropy, PMD, HMR/LMR, HMM-DMR,
Zarr backing, single-cell sparse stores, eQTM, motif/TFBS enrichment,
matrix-completion imputation, pyGenomeTracks-style track plot.

---

## [0.3.0]

### Added

#### Visualization
- **Karyogram / chromosome painter** (`ep.pl.karyogram`). One row per
  chromosome, megabase-binned mean of any per-CpG metric (`meth_diff`,
  `-log10_qvalue`, raw β). RdBu_r by default; symmetric colour limits
  for signed metrics.
- **DMR UpSet / Venn overlap** (`ep.pl.dmr_overlap`). 2-set inputs fall
  back to a 2-circle Venn; 3-6 set inputs render an UpSet plot
  (bar chart + dot matrix + per-set totals). Matplotlib-only; no
  upsetplot or matplotlib_venn dependency.
- **Gene-body metaplot** (`ep.pl.gene_body_metaplot`). Three-zone TSS /
  body / TES plot with flanking windows; body length is normalised
  across genes so a 50 kb and a 500 kb gene contribute equally.

#### Statistical features
- **DVR — Differentially Variable Regions** (`ep.tl.dvr`,
  `ep.call_dvr_density`). Region-level aggregation of `tl.dvc` output
  via per-tile DVC-density enrichment (one-sided binomial vs the
  genome-wide rate, BH-corrected). Avoids the variance-statistic
  combining problem that defeats Fisher / Stouffer's at the
  per-CpG level.
- **Effect-size shrinkage** (`ep.shrink_meth_diff`). Empirical-Bayes
  Normal-prior James-Stein-style shrinkage of `meth_diff` toward zero.
  Adds `meth_diff_shrunk`, `meth_diff_se`, `shrinkage_factor` columns
  to a DMC table. Same spirit as `ashr` / `apeglm`; pulls
  low-coverage sites harder than well-powered ones.
- **kNN β imputation** (`ep.impute_knn_beta`, `ep.impute_knn_anndata`).
  Per-chromosome inverse-distance-weighted kNN over genomic position;
  optional `max_distance_bp` cap so cross-CGI gaps don't pull
  long-distance neighbours.

#### Clocks / deconvolution scaffolding
- **Generic linear age-clock runner** (`ep.tl.age_clock`,
  `ep.age_clock`). Takes a user-supplied `(cpg_id, coefficient)` table
  and a probe → `(chrom, pos)` manifest, computes per-sample age, and
  writes the result into `md.obs`. Supports an optional `transform`
  argument (`"horvath"` for the standard ≥20-year piecewise
  anti-transform); other published clocks (Hannum, PhenoAge,
  DunedinPACE) plug in via their own coefficient CSVs. Coefficient
  tables and manifests are **not bundled** — licensing and probe
  vendor specifics differ per clock.
- **Reference-based cell-type deconvolution** (`ep.tl.deconvolve`,
  `ep.deconvolve`). Non-negative least squares solve against a
  user-supplied reference β matrix (EpiDISH / CIBERSORT / Houseman
  style). Long-format result on `md.uns['deconvolution']`; wide
  per-cell-type columns (`frac_<celltype>`) joined onto `md.obs`.

### Tests
- New: `test_viz_new.py` (8 tests), `test_dvr.py` (6 tests),
  `test_stats_new.py` (12 tests) covering shrinkage, kNN imputation
  end-to-end on AnnData, age-clock recovery on synthetic data, and
  deconvolution NNLS round-trip.

### Notes on what's *still* left on the 0.3+ roadmap

Architectural lifts that didn't fit this round and remain unimplemented:
ASM (per-read haplotype phasing from BAM), PMD / HMR / LMR callers,
methylation entropy, single-cell methylation, distributed compute
(Dask / Ray), GPU IRLS, HMM-based DMR, Zarr backing, Tabix-on-Parquet,
formal checkpoint / resume API. These each warrant their own focused
release; the karyogram + UpSet + gene-body trio plus DVR / shrinkage /
imputation / clocks / deconvolution covers the highest-ROI subset of
the roadmap.

---

## [0.2.0]

### Added
- **MethylDackel input adapter.** `ep.read_methyldackel(samplesheet, ...)` and
  `epykit convert --format methyldackel` ingest MethylDackel `.bedGraph[.gz]`
  output through the same partitioned-Parquet pipeline as Bismark. The
  conversion cache is format-aware so a Bismark store cannot be silently
  reused for MethylDackel input.
- **Permutation empirical FDR for DMC.** `ep.tl.dmc(..., empirical_fdr=True,
  n_perm=100)` shuffles treatment / control labels, re-runs the per-CpG
  DMC engine, and emits `empirical_pvalue` / `empirical_qvalue` columns —
  parity with the existing DMR permutation FDR. Refused with `formula=` /
  `contrast=` designs (label shuffling invalidates stratified models).
- **Bismark M-bias parser and plot.** `nfcore_qc.parse_bismark_mbias(path)`
  parses the Bismark M-bias text format into a long table; `pl.mbias_plot(
  {sample_id: df_or_path})` renders percent methylation per read position
  with context / R1 / R2 lines.
- **CLI `--version` flag.** `epykit --version` prints the installed
  `__version__` from `importlib.metadata`.
- **DVC engine re-export.** `ep.process_chromosomes_dvc` is now part of the
  public surface alongside the high-level `ep.tl.dvc` orchestrator.

### Changed
- **GLM degeneracy is visible.** When the batched `(X'WX)⁻¹` solve in
  `_glm.py` falls back to a per-site solve, the helper emits a
  `logger.warning` summarising affected sites instead of failing silently.
  Per-site GLM separation is logged at `info` (or `warning` when ≥5 % of
  sites separate) rather than `debug`.
- **Bisulfite conversion rate is reported, not applied.** Doc clarification
  in `qc.bisulfite_conversion_rate` and the README: epykit follows
  `bsseq` / `methylKit` defaults and does not rescale per-CpG counts by the
  conversion rate. The rate is surfaced through QC, MultiQC export, and
  the HTML report so users can gate on it.
- **Multi-group DMC accuracy test.** `tests/test_dmc_multigroup.py` now
  asserts power (≥30 %) and FDR (≤15 %) on the 3-group joint F-test
  fixture instead of only checking column presence. Continuous-covariate
  test is now an explicit structural check (engine runs, finite
  p-values, FDR isn't catastrophic) with the fixture limitation
  documented in the test.
- **DMR sensitivity / specificity tests.** `test_accuracy.py` gains
  `test_dmr_tile_sensitivity_and_fdp` and
  `test_dmr_sliding_window_sensitivity_and_fdp`, pinning both recovery
  rate and false-discovery proportion per method. Conditional
  `pytest.skip` calls in the existing DMR direction test are now hard
  assertions to surface regressions instead of masking them.

### Fixed
- `tests/test_dmc_multigroup.py` no longer silently swallows a `ValueError`
  in the contrast-resolution test; on `ValueError` the test now asserts
  the error message is informative and on success it asserts finite
  p-values are produced.

## [0.1.0]

- Initial release. Bismark `.cov` → partitioned Parquet methylstore;
  scanpy-style `pp` / `tl` / `pl` API; 8 DMC test backends; two DMR
  engines plus permutation FDR; DVC (iEVORA-style); covariate-aware
  contrasts; AnnData / MuData / methylKit / MultiQC interop;
  self-contained HTML report.
