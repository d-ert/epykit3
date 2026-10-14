# Confirmed Epykit Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task after the user requests execution. Steps use checkbox (`- [ ]`) syntax for tracking. This task authorizes documentation only; do not execute this plan now.

**Goal:** Correct final-region seed retention and GLM batch validity, then deliver verified annotation and masked-GLS optimizations.

**Architecture:** Keep the existing statistical models and public interfaces. Deliver B1, B2, G1, G2, A1 and M1 as separate commits and benchmark variants, using independent numerical/behavioral oracles. Profile remaining GLS costs without implementing speculative rewrites.

**Tech Stack:** Python, NumPy, Numba, Polars, SciPy, statsmodels, pytest; existing dependencies only.

**Spec:** [Confirmed Epykit fixes specification](../specs/2026-10-14-confirmed-epykit-fixes-design.md).

**Status:** Draft for user review. No implementation or benchmark steps have run.

## Global Constraints

- Python compatibility remains >=3.10; retain existing dependency bounds and add no dependencies.
- Public function signatures, output column names, row ordering and dtypes remain compatible.
- Preserve the observed-cell CAR/donor covariance model and per-label support replay.
- Preserve experimental FDP/maxT semantics; make no new finite-sample FDR guarantee.
- Preserve baseline results and input stores; write every benchmark variant to a new directory.
- Run compiled Numba checks and warm timing runs; JIT-disabled results do not establish compiled correctness or speed.
- Use no simulation truth or real-paper overlap to choose algorithm parameters.
- No implementation, test creation, benchmark execution, dependency installation or source commits are authorized by this documentation task.

## Review Focus

- Trimming removes the parent's only seed while a different admissible seeded interval exists: select that interval (Task 1).
- A site has `n_eff>=p` but covered design rows are collinear: return unavailable inference and preserve neighboring valid sites (Task 3).
- An all-invalid or empty GLM batch has no active rows: terminate safely with documented shapes and NaNs (Task 4).
- TSS ties, duplicates, absent chromosomes and negative strands: preserve scalar output exactly (Task 5).
- BLAS threading or changed Gram summation perturbs a near-exact masked fit: retain degeneracy rejection and control thread counts (Task 6).

## Workspace and file map

Source paths below are relative to `development/epykit_statistical_repair` unless
marked workspace-relative. Start from clean `fe2892e`; use the worktree workflow
at execution time and record the resulting `FIX_REPO`. No worktree is created now.

| Files | Responsibility |
|---|---|
| `src/epykit/_regional_gls_native.py`, `_regional_gls_discovery.py`, `_regional_gls.py` | B1/B2, regional revision metadata and cache tests |
| `src/epykit/_glm.py`, `_dmc_config.py` | GLM validity, active iteration and resume revision |
| `src/epykit/_dmc_stages.py` | GLM result-metadata revision integration |
| `src/epykit/annotate.py` | A1 vectorized lookup |
| `src/epykit/_regional_gls_masked.py` | M1 donor products |
| `tests/test_discovery_repair.py`, `test_regional_gls_multiscale.py`, `test_regional_gls_support_repair.py`, `test_regional_gls.py` | Extend existing discovery/support/cache coverage |
| `tests/test_glm_batch_repair.py`, `test_glm_active_repair.py`, `test_glm_revision_repair.py`, `test_nearest_tss_vectorized.py`, `test_regional_gls_gram_repair.py` | New independent regression/oracle coverage |
| `benchmarks/confirmed_fixes/profile_components.py`, `benchmark_components.py`, `README.md` | Reproducible component profiling, microbenchmarks and run instructions |
| Workspace `scripts/run_statistical_repair_benchmarks.py` | Existing immutable-source pilot/real harness; use `--repo` |
| Workspace `scripts/run_statistical_repair_validation.py` | Existing repeated-null harness; use `--repo` |

Execution shell variables, defined only when implementation is authorized:

```sh
BENCH_ROOT=/Users/dert/Desktop/Coding/Epykit/wgbs_benchmark_v6
EPY_PY="$BENCH_ROOT/development/epykit_minimal_gls_trial/.venv/bin/python"
# FIX_REPO is the absolute isolated checkout path recorded by the executor.
```

For repository tests run from `FIX_REPO` with `PYTHONPATH=src`; verify
`epykit.__file__` resolves inside that checkout. The workspace harness accepts
`--repo "$FIX_REPO"` and freezes it. Commit source before a harness run, as required
by its existing source snapshot validation. Do not modify old result directories.

## Task 0: Freeze comparison inputs and establish component measurements

**Files:** Create the three `benchmarks/confirmed_fixes/` files above.
**Interfaces:** Produce a JSON measurement table containing commit/source hash,
input hash, stage, dimensions, mask density, worker/BLAS counts, cold/warm mode,
elapsed seconds and process peak RSS. Keep this schema identical across variants.

- [ ] Verify clean baseline, dependency versions and source import path; preserve `fe2892e` as the baseline checkout/ref. Record hardware and environment.
- [ ] Add an opt-in, standalone profiler using `cProfile` for Python orchestration and explicit timers around `_prepare_label_scan`, `fit_regions` and `_smooth_tricube`; time join setup with a controlled local fixture. Do not add production timing switches. Use separate subprocesses for peak-RSS measurements.
- [ ] Add fixed-seed component fixtures: GLM 100,000 sites/S=10/p=3 with zero and one invalid site; nearest TSS 1,000,000 mixed-chromosome sites and 20,000 TSS entries; masked GLS S=6/10, K=5/20/100, 0%/1%/20% dropout, with each locus still observed in both groups. Preserve fixture hashes.
- [ ] Run one warm-up plus five repeats, with workers=1 and 4 and BLAS threads fixed to 1; record compilation separately. Profile one pilot and one real chromosome using cache-disabled isolated outputs. Capture complete/masked candidate proportions and actual representative counts.
- [ ] Commit the harness and baseline measurement instructions. These measurements establish evidence; they do not assert a bottleneck or authorize deferred work.

## Task 1 (B1): Enforce the final-candidate seed invariant

**Files:** Modify `_regional_gls_native.py`, `_regional_gls.py`; extend discovery, multiscale, support and cache tests.
**Interfaces:** Extend private `_trim_and_append_subregion(..., regions, seed_signal=None, min_seed_effect=0.0)` and `_refine_and_append_region(..., boundary_min_effect=0.0, seed_signal=None)`. `discover_prepared` retains its public callable signature and passes its actual computed/supplied smoothing array to refinement. Optional defaults retain existing low-level oracle callers.

- [ ] Add `test_split_fragment_must_retain_its_own_seed`: the spec's 13-site fixture returns only `[[8,13]]` for positive and negative directions with boundary .05. Add corresponding boxcar and multiscale production-path fixtures.
- [ ] Add exhaustive-oracle tests over 100 fixed-seed small arrays: enumerate intervals satisfying span/count and a signed seed, compare optimal penalized score and specified ties. Include an optimum whose seedless alternative scores higher, a seed trimmed from an edge, exact-threshold equality, no admissible seed, and a genomic gap.
- [ ] Run the new tests on baseline; expect seedless-fragment and constrained-seed cases to fail. Existing geometry tests remain passing controls.
- [ ] Add the seed constraint to admissible prefix minima, using the most recent matching seed to limit allowed starts. Materialize the discovery signal once per scale when needed. Preserve existing B1 penalties and edge-floor rules.
- [ ] Set regional `algorithm_revision` to `spatial-gls-observed-data-multiscale-seeded-v10`. Test cache miss for baseline source digest, new-source cache reuse, complement-label symmetry, per-label eligibility replay and unchanged empty-scan accounting.
- [ ] Run compiled discovery/support/cache tests, then the complete repository suite. Commit as `fix(regional-gls): retain seeds in final candidate intervals` and save a B1-only pilot benchmark variant.

## Task 2 (B2): Decouple extension and scoring thresholds

**Files:** Modify `_regional_gls_native.py`, `_regional_gls_discovery.py`, `_regional_gls.py`; extend discovery/multiscale tests.
**Interfaces:** No new public options. Both refinement and multiscale competition use `tau=.55*min_effect`; envelope and edge-floor rules remain as specified.

- [ ] Add a parameterized threshold test for `min_effect=.1` and boundaries None/.05/.08: tau is always .055, envelope threshold is .1/.05/.08 and edge floor is .035/.05/.08. Add independently enumerated production fixtures that distinguish `.05` from `.055` scoring.
- [ ] Run on B1; expect boundary-enabled penalty cases to fail.
- [ ] Replace the two boundary-dependent penalty assignments. Set regional revision to `spatial-gls-observed-data-multiscale-seeded-decoupled-v11`; preserve source-hash cache invalidation.
- [ ] Run compiled discovery and mirror tests and the full suite. Commit as `fix(regional-gls): separate envelope extension from region scoring`.
- [ ] Benchmark B1+B2 separately against B1. If the spec's null-positive/recall guard fails, defer this commit from the combined branch and report it; do not retune the multiplier.

## Task 3 (G1): Isolate invalid and singular GLM sites

**Files:** Modify `_glm.py`, `_dmc_config.py`, `_dmc_stages.py`; create `test_glm_batch_repair.py` and `test_glm_revision_repair.py`.
**Interfaces:** Add private `_estimable_sites(X, has_cov) -> bool[n_sites]` and `_solve_stacked_systems(A, B) -> float[n_sites,p,n_rhs]`, where A has shape `(n_sites,p,p)` and B `(n_sites,p,n_rhs)`. Return NaN rows for unsolved systems. Retain `_solve_weighted_lsq(X,z,w)` and all public IRLS return signatures.

- [ ] Add tests comparing valid rows against independent scalar solves and statsmodels fits. Cover one singular site among 100 valid sites, all singular, no coverage, `n_eff>=p` with collinear covered rows, a full-rank ill-conditioned fixture and `return_cov=True` joint covariance.
- [ ] Instrument linalg calls in a regression: one bad system must not cause scalar calls for all 100 valid rows. Count p-dimensional singleton solves separately from stacked calls, and check unaffected valid rows numerically.
- [ ] Run tests on baseline; expect isolation/rank-invalid coefficient expectations to fail.
- [ ] Classify rank once per observed support pattern using the spec's SVD threshold. Solve only estimable rows. Implement 8,192-system top-level chunks with recursive splitting only after `LinAlgError`; mark failing/nonfinite singleton results invalid. Replace final batch inversion with identity-RHS solves through the same helper. Remove invalid-site ridge inference from this path.
- [ ] Add `GLM_ALGORITHM_REVISION` with the exact spec value. Include `glm_algorithm_revision` for binary `selected_test="glm"` resume signatures and GLM/contrast output metadata; test that changing the revision invalidates a cached GLM result while unrelated LR signatures remain unchanged. Do not introduce contrast resumption.
- [ ] Run focused tests plus full suite. Benchmark valid-only and one-invalid batches; inspect runtime, peak memory and NaN counts. Commit as `fix(glm): isolate rank-deficient sites in batched solves`.

## Task 4 (G2): Stop recomputing converged GLM rows

**Files:** Modify `_glm.py`; create `test_glm_active_repair.py`.
**Interfaces:** Keep `irls_binomial_batch(...)` and five/six-element results unchanged. Reuse G1 validity and stacked solves; introduce no public convergence controls.

- [ ] Write a scalar per-site reference loop with the spec's stopping rules. Test a mixed easy/slow batch, no sites, one site, all-invalid input, invalid updates, separation and max-iteration exhaustion. Compare all output arrays and NaN masks, including optional covariance and exact `n_eff`.
- [ ] Instrument `_solve_weighted_lsq` row counts to show the active batch shrinks; converged rows must receive no further iteration updates. Run on G1 and confirm this assertion fails.
- [ ] Scatter accepted updates to active indices, freeze newly converged rows, invalidate nonfinite updates and terminate on an empty active set. Calculate final diagnostics from final accepted coefficients and preserve existing separation/nonconvergence masking.
- [ ] Run G1/G2 tests and full suite. Meet `rtol=1e-6, atol=1e-8` on well-conditioned independent-oracle fixtures. Inspect any DMC selection differences near cutoffs.
- [ ] Benchmark the same mixed-convergence fixture; report row-updates avoided, median time and memory. Commit as `perf(glm): update only active IRLS sites` if acceptance gates pass.

## Task 5 (A1): Vectorize nearest-TSS annotation

**Files:** Modify `annotate.py`; create `test_nearest_tss_vectorized.py`.
**Interfaces:** Add private `_nearest_tss_arrays(site_chroms, site_mids, tss_by_chrom) -> (gene_names, signed_distances)`, returning object strings and int64 distances before existing int32 clipping/sentinel assembly.

- [ ] Freeze the old scalar lookup as a test-only oracle. Assert exact gene and distance equality for randomized unsorted chromosomes, half-integer midpoints, exact hits, midpoint ties, duplicate TSS entries, both strands, extreme distances and absent/empty TSS chromosomes.
- [ ] Add an integration fixture through `annotate_features` with `multi_annotation=True`; assert output row order/dtypes/sentinels, and that the helper is not called when false. Run new helper tests on baseline; expect missing-helper failures.
- [ ] Perform one vectorized `searchsorted(..., side="left")` per chromosome, compare left/right candidates using `<=` for left-tie preference, handle empty entries before indexing, and scatter to original indices. Preserve integer midpoint truncation and strand signs.
- [ ] Run focused tests and full suite. Compare full annotation outputs exactly. Benchmark numeric lookup separately from gene-string/output construction; apply the spec's timing/memory gates.
- [ ] Commit as `perf(annotation): vectorize chromosome nearest-TSS lookup`.

## Task 6 (M1): Remove donor-Gram copy amplification

**Files:** Modify `_regional_gls_masked.py`; create `test_regional_gls_gram_repair.py`; retain existing masked/compiled oracle suites.
**Interfaces:** Add private Numba `_donor_products(diagonal, off, donor_residual, qr) -> (gram, response)` with `(S,S)` and `(S,)` outputs. `_spectral_prepare` retains its tuple contract.

- [ ] Add a frozen scalar-product reference and dense observed-data GLS comparisons. Include S=6/10, K=5/20/100, sparse/dense dropout, the spec's covariance endpoints, an absent sample and response scales `1e-8/1/1e8` near exact fits.
- [ ] Run helper tests on baseline; expect missing-helper failures. Run existing exact-fit and dense-oracle tests as controls.
- [ ] Build the CAR-applied matrix once, then form `D.T @ QD` and `D.T @ qr` using contiguous arrays. Preserve residual sums, symmetrization and eigenvalue checks. Keep complete-data analytical code and all REML searches unchanged.
- [ ] Run with JIT enabled, workers=1/4 and BLAS threads=1; run existing mirror, parallel and exact-fit tests. Compare fixed covariance at `rtol=2e-8, atol=2e-8` and inspect all profiled-score/selection differences.
- [ ] Benchmark compiled preparation and end-to-end masked fitting. Keep this optimization only if the specified 10% gain and regression/memory gates pass. Commit accepted work as `perf(regional-gls): reuse masked donor matrix products`; otherwise retain the benchmark evidence and omit the source change.

## Task 7: Remove confirmed unused regional helpers

**Files:** Modify `_regional_gls.py`; update relevant import-only dependencies.
**Interfaces:** Remove private `_precompute_discovery_geometry`, `_fast_discover_with_geom` and `_fast_decoy_statistic` only if still unused.

- [ ] Search source, tests, benchmark scripts and package exports with `rg`; check for references to string names and saved diagnostic imports. Preserve active helpers and any independently useful oracle.
- [ ] Delete unused definitions and imports without changing callers. This task needs no mirror test of deleted implementation.
- [ ] Run full source suite and lint changed files; compare the source change to verify only dead definitions/imports were removed. Commit as `refactor(regional-gls): remove unused discovery helpers`.

## Task 8: Integrate, benchmark and record disposition

**Files:** Update `benchmarks/confirmed_fixes/README.md`; produce a versioned validation report in `docs/benchmarks/confirmed_fixes/`. No edits to existing baseline evidence.
**Interfaces:** Consume each task's frozen commit and JSON measurements. Produce a table of accepted/deferred commits, correctness results, changes in outputs, wall time and peak RSS.

- [ ] Run the full source suite and workspace simulator/harness suite with the intended checkout imported; run Ruff over changed source/tests/harness files. Required result: zero failures and zero new lint findings. Use compiled Numba for integration tests.
- [ ] Run the predeclared small pilot set (empirical null/signal seed2026, generative null/signal seed2027, empirical sparse missing signal, BLUEPRINT signal/null), all four configs, for baseline, B1 and B1+B2. Check seed violations are zero in changed variants and report recovery/null changes without claiming FDR control.
- [ ] Run the complete 18-dataset/four-config pilot suite and four full-genome real configs for baseline and the accepted combined branch. Include DMC for harness coverage, but separately exercise `irls_binomial_batch`/the public GLM path: the harness's `--include-dmc` selects depth-aware LR and does not validate the GLM changes.
- [ ] Run existing development-null datasets separately from held-out null datasets. Report completed independent datasets, missing/failed runs and binomial uncertainty. Keep holdout outcomes as validation; do not retune against them. All benchmark invocations use new output roots and source snapshots.
- [ ] For optimization-only commits compare regional coordinates/directions and selected calls against their immediate behavior-equivalent parent. Validate annotations exactly and GLM valid-site numerical tolerances; explain every changed NaN mask and threshold crossing.
- [ ] Apply spec gates to B2 and optimization commits, keeping B1/G1 correctness decisions separate from speed claims. If a gate fails, record and omit the failing optional change rather than concealing it in the combined result.
- [ ] Reprofile accepted combined code. Recommend deferred work only when a stage accounts for >=20% warmed wall time or measured peak memory. Provide the measurement and next proposed experiment; implement no follow-up architecture here.
- [ ] Publish local documentation with exact commits, parameters, input/evaluation hashes, test commands, limitations and accepted/deferred tasks. Stop for user review before external integration or publication.

### Future execution commands (not run by this task)

From the isolated source checkout:

```sh
PYTHONPATH=src "$EPY_PY" -m pytest -q tests
```

From the benchmark workspace, after each source variant has been committed:

```sh
"$EPY_PY" scripts/run_statistical_repair_benchmarks.py \
  --repo "$FIX_REPO" --kind pilot --workers 4 --include-dmc \
  --output results/confirmed_fixes_v1/combined/pilot
"$EPY_PY" scripts/run_statistical_repair_benchmarks.py \
  --repo "$FIX_REPO" --kind real --workers 4 --include-dmc \
  --output results/confirmed_fixes_v1/combined/real
"$EPY_PY" scripts/run_statistical_repair_validation.py \
  --repo "$FIX_REPO" --dataset-root data/statistical_repair_null_validation \
  --workers 4 --output results/confirmed_fixes_v1/combined/null_development
"$EPY_PY" scripts/run_statistical_repair_validation.py \
  --repo "$FIX_REPO" --dataset-root data/statistical_repair_null_holdout \
  --workers 4 --output results/confirmed_fixes_v1/combined/null_holdout
```

Run baseline/B1/B1+B2 commands with their own `--repo` and output directory.
Use `--datasets` for the named small pilot set and `--chromosomes chr1` for the
initial real profile. The repeated-null harness defaults to two configurations;
record that scope and do not describe it as validation of all four settings.

## Documentation self-review and handoff

Coverage: spec B1/B2 -> Tasks 1/2; G1/G2 and cache revision -> Tasks 3/4;
A1 -> Task 5; M1 -> Task 6; cleanup -> Task 7; measurement/acceptance/deferred
decisions -> Tasks 0/8. Each Review Focus item has an owning test task.

This plan is a reviewable draft, not implementation approval. The user explicitly
requested a stop after specification and planning. Do not create tests, edit
source, execute the listed commands, create agents or begin implementation now.
