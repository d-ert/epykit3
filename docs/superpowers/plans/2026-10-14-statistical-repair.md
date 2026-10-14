# Statistical Repair Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Repair the reviewed statistical errors, improve supported recovery, run simulation and whole-genome real benchmarks, and deliver a verified Gemini bundle.

**Architecture:** Exact observed-data masked GLS alongside the existing complete-data native path; identical eligibility and discovery rules for every label assignment; moment-based DMC EB and complete-search-family correction. All benchmark outputs carry exact source provenance and preserve prior runs.

**Tech Stack:** Python, NumPy, SciPy, Numba, Polars, PyArrow, pytest; existing local environments and count stores.

**Spec:** `docs/superpowers/specs/2026-10-14-statistical-repair-design.md`

## Global Constraints

- Python >=3.10; no new runtime dependency.
- Work in `development/epykit_statistical_repair`; preserve existing sources, input stores and results.
- Use explicit observed masks; never infer missingness from a residual equal to zero.
- Same eligibility, discovery and statistic procedure for observed and decoy labels.
- Recovery improvements must be reported with error, boundaries, runtime and source provenance.
- Keep regional pooled-permutation calibration experimental until supported by independent validation.
- User authorizes autonomous design/implementation decisions and periodic commits; no approval waits.

## Review Focus

- Missing whole sample/chromosome partitions must preserve eligible union loci.
- Zero observed residuals and exact-fit regions must have explicit, symmetric handling.
- Highly dispersed loci and infinite prior df must return calibrated finite reference p-values.
- Minimum count/span constraints must retain the globally best admissible segment.
- A resumed benchmark must use its recorded source and must not mix incomplete or old runs.

### Task 1: DMC EB and neighbour corrections

**Files:** modify `src/epykit/dmc.py`, `src/epykit/_glm.py`; create `tests/test_dispersion_repair.py`.
**Interfaces:** preserve `_smyth_eb_moderation` return tuple; `reference_pvalues` handles infinite denominator df.

- [ ] Write failing tests for constant scale across heterogeneous df, inverse-chi-square prior recovery, prior scale50 null reference, infinite-df references, invalid/small inputs and NaN singleton counts.
- [ ] Run targeted pytest and inspect expected failures.
- [ ] Implement expectation-based log-moment fitting, unclipped prior scale, post-pooling dispersion clamp, limiting references, finite-contributor counts and cache revision.
- [ ] Run targeted and full pytest; commit the independently verified milestone.

### Task 2: Constrained discovery and honest region multiplicity

**Files:** modify `_regional_gls_native.py`, `_region_search.py`, `dmr.py`, `tl.py`; create `tests/test_discovery_repair.py` and `tests/test_region_search.py`.
**Interfaces:** preserve discovery APIs and region schema; `count_intervals` counts all geometry-defined spans; `adjust_selected` uses BY.

- [ ] Write failing constrained-trim fixture and brute-force small-array oracle, independent discovery geometry checks, interval-count and BY-reference checks.
- [ ] Observe failures; implement admissible prefix minima and remove untrimmed fallback, restore complete family/BY, withdraw unsupported calibrated-default wording and arbitrary default loosening.
- [ ] Test actual batch SE against raw profiled fit, remove undocumented length modifier, and update algorithm cache revisions.
- [ ] Run targeted/full pytest and commit.

### Task 3: Exact masked GLS and exchangeable union scans

**Files:** create `src/epykit/_regional_gls_masked.py` and `tests/test_regional_gls_masked.py`; modify `_regional_gls.py`, `_regional_gls_native.py`.
**Interfaces:** add optional `obs_mask` to fixed/native fitting APIs; complete support retains numerical results; output remains five fit columns with effects, SE and REML objective.

- [ ] Write independent statsmodels/dense GLS observed-submatrix tests with irregular missingness, valid zero residuals, complementary labels, accurate residual df and complete-support agreement.
- [ ] Observe failures; implement observed covariance, nuisance-intercept elimination and identified treatment information with native numerical profiling.
- [ ] Write integration tests for relabelled support/universe, missing chromosome partitions and raw-effect reporting; observe failures then apply per-scan eligibility and label-independent universe.
- [ ] Run targeted/full pytest and performance probe; optimize only with oracle equivalence; commit.

### Task 4: Provenance-backed benchmarks and evidence-driven iteration

**Files:** create workspace `scripts/run_statistical_repair_benchmarks.py`, comparison/verification helpers and `reports/statistical_repair/`.
**Interfaces:** consume immutable committed source, existing sample sheets/count stores/truth/workbook; emit per-run manifests, candidates, selected calls, metrics, aggregate TSV/JSON/Markdown.

- [ ] Run baseline tests in isolated worktree; verify reusable count stores and source hashes.
- [ ] Add meaningful fixture tests for resumable provenance and metric aggregation before runner implementation.
- [ ] Run 18 pilot datasets for comparable coverage/boundary settings and six-sample whole-genome real benchmarks; compare with old and stage2 reports.
- [ ] Diagnose failures or loss of useful recovery, apply independently tested corrections, and rerun affected benchmarks without overwriting prior attempts.
- [ ] Run fresh-seed null/calibration validation and record limitations; commit tested source milestones and final benchmark scripts/report.

### Task 5: Independent review, integration artifact and Gemini bundle

**Files:** update source docs; create workspace `scripts/build_statistical_repair_bundle.py`, `gemini_bundle_v2/`, transferable git bundle/patch and final verification manifest.
**Interfaces:** consume exact final commits/results; produce source, tests, scripts, evidence, prompt/dossier and hashed inventory.

- [ ] Dispatch fresh whole-branch review with spec, plan, diff and verification evidence; fix material findings with regression tests.
- [ ] Run final full tests, lint on changed code, result verifiers and source consistency checks; commit.
- [ ] Create fresh bundle preserving older bundles; verify every packaged source/result hash and disclose calibration limits and benchmark reuse.
- [ ] Audit all user requirements against current files and process results; mark goal complete only when every deliverable is verified.
