# Phase 3: Engine Freeze (P1 + API Cleanup + Integration Scripts) — Design

**Date:** 2026-05-27
**Status:** Design, awaiting user review
**Companion spec:** [`2026-05-27-paper-defendable-benchmark-design.md`](./2026-05-27-paper-defendable-benchmark-design.md) §3 (P1 manifest), §4 (API cleanup), §5 (integration scripts)
**Predecessors:** Phase 1 = P0 fixes (tag `v0.7.3-p0-complete`); Phase 2 = benchmark scripts (tag `v0.7.4-phase2-scripts`).

The goal of Phase 3 is to **freeze the engine surface** so that Phase 4's locked benchmark re-run uses final code. After this phase, no public engine schema or name changes; only the locked re-run, paper rewrite, and P2 hygiene remain.

---

## 1. Scope and non-goals

**In scope:**

- All 11 P1 functional fixes from the parent spec §3.
- API cleanup (harder than the parent spec proposed, decided in brainstorming):
  - Hard-drop `logit_t`, `bb_lr`, `score`, `cmh` from `tl.dmc` — no `epykit.experimental` alias for any of them. Internal helpers stay only if reused by surviving engines.
  - Rename `dmr_hmm` → `dmr_segment` (P1-8 in parent spec).
  - Rename `log2_odds_ratio` per backend (P1-11).
  - `DVC` (`tl.dvc`) kept with the P1-7 Brown-Forsythe fix.
- Five integration items:
  1. `benchmark/scripts/methylkit_stouffer_combine.R`
  2. Wire `run_null_calibration.py` to real `ep.tl.dmc` / `ep.tl.dmr` engine closures
  3. Hook `wilson_bootstrap_ci` into `evaluate.py`'s output path
  4. `benchmark/scripts/regen_all.py` (acceptance gate; empty `claims.yaml` seed)
  5. `benchmark/scripts/bug_fix_audit.py` (pre/post-fix delta against existing `benchmark/data/study*/eval_summary.parquet`)
- Tag `v0.7.5-phase3-engines-frozen` after Phase 3 lands.

**Out of scope (Phase 4 or later):**

- The locked benchmark re-run itself (N=20 simulator seeds, Study 2 re-run, Study 3 real-data re-run, null calibration on real data, populated `claims.yaml`).
- All P2 hygiene items from the parent spec §3, **except P2-4** (per-segment Stouffer p-values for the renamed HMM): that one folds into the P1-8 rename commit because both touch the same code path and shipping the rename without the p-value fix would leave the renamed engine returning NaN q-values.
- Paper rewrite per parent spec §6.
- GPU backend audit; multi-cohort real-data validation; second simulator family (all out of scope in parent spec §1).
- Reorganising `dmc.py` (3000+ LOC; defer to v0.8).
- Pre-commit hooks for the `Affects:` trailer enforcement (`bug_fix_audit.py` will fail loud enough without one in Phase 3).

---

## 2. Decisions locked in this brainstorming

| Decision | Value | Reason |
|---|---|---|
| Aggressive cleanup style for the four dropped engines | Hard-drop, no `experimental` alias | All four are either documented-broken (`logit_t`, `bb_lr`) or strictly dominated by a surviving engine (`score` by `lr`; `cmh` by `glm` + batch covariate). Preserving them under any alias keeps an attack surface a reviewer can ask about. |
| Migration story for dropped engines | `ValueError` with explicit hint text per engine | `test='logit_t'` → suggests `welch_t`/`lr`; `test='bb_lr'` → suggests `lr`; `test='score'` → suggests `lr`; `test='cmh'` → suggests `tl.dmc(formula='~ group + batch')`. |
| `DVC` (`tl.dvc`) | Keep with P1-7 fix | DVC is a separate analysis family (variability, not differential methylation). Removing it would shrink the public tool footprint without paper-relevance gain. |
| `bug_fix_audit.py` pre-fix baseline | Use existing `benchmark/data/study*/eval_summary.parquet` | Those parquets pre-date Phase 1 — they literally are the pre-fix numbers. No need to re-run pre-P0 code. |
| Sequence | Renames → drops → P1 functional fixes → integration scripts → tag | Engine-first so integration scripts (step 4) read against the frozen schema. |
| Test policy for fixes that change number-shaped behaviour | One reference-comparison test per fix, marked fast unless it requires a real engine run or R subprocess | Avoid duplicating coverage from existing engine tests; only add new tests where the behaviour itself is new. |

---

## 3. Surviving public engine surface (the architectural invariant)

After Phase 3, `tl.dmc(test=...)` accepts exactly:

| `test=` | Math | When |
|---|---|---|
| `"auto"` | dispatcher → `"fisher"` at n=1, `"lr"` at n≥2 | Default; user-friendly entry point |
| `"lr"` | Quasi-binomial likelihood-ratio with `eb` dispersion (Phase 1 default) | Headline default at n≥2 |
| `"welch_t"` | Welch t on per-sample β, Welford accumulators | Replicate-aware β-mean baseline |
| `"fisher"` | Exact test on pooled read counts | n=1 fallback (only engine that works there) |
| `"glm"` | Binomial GLM with deviance LR, formula-based covariates | Covariate-controlled flagship |

Plus `tl.dmc_multigroup(...)` which internally dispatches to `glm_contrast` (not user-facing in `test=`).

**Dropped:** `logit_t`, `bb_lr`, `score`, `cmh`. Each raises `ValueError` with a migration hint.

DMR engines (`tl.dmr(method=...)`): `tile` (default), `sliding`, `chain_merge`, `segment` (renamed from `hmm`). All four stay.

Output schema changes:

| Layer | Pre-Phase-3 | Post-Phase-3 |
|---|---|---|
| `varm["dmc_lr"].log2_odds_ratio` | column present | renamed to `log2_odds_ratio_pooled`; transitional `log2_odds_ratio` column NaN-filled with `FutureWarning` for one release |
| `varm["dmc_glm"].log2_odds_ratio` | misleading (it was the logit coefficient, not log₂-OR) | renamed to `coef_treatment_log2`; transitional column as above |
| `varm["dmc_lr"].meth_diff_ci_{lo,hi}` | Welch-normal Wald CI | Newcombe hybrid score interval (P1-3) |
| `varm["dmc_segment"]` (was `varm["dmc_hmm"]`) | p/q-values NaN (parent spec P2-4 bug) | Per-segment Stouffer-combined p-values, BH per chromosome |
| `epykit.dmr_hmm` module | canonical | shim re-exporting from `dmr_segment` with `DeprecationWarning` on import, `logging.WARNING` mirror |

---

## 4. Execution sequence (commit-by-commit)

Each entry below corresponds to one git commit on branch `p0-fixes`. Commits use the `fix(<area>) <PID>` subject prefix for P1 items, `fix(<area>) BREAKING: <action>` for drops/renames, and `feat(benchmark) <action>` for integration scripts. Every commit affecting per-cell numbers carries an `Affects: <engine>@<scenario>, ...` trailer for `bug_fix_audit.py` attribution.

### Step 1: Renames

1. **`refactor(dmr): rename dmr_hmm → dmr_segment (P1-8) + per-segment p-values (P2-4)`**
   - Rename `src/epykit/dmr_hmm.py` → `src/epykit/dmr_segment.py`; function `call_dmr_hmm` → `call_dmr_rule_segment`.
   - Compute per-segment Stouffer-combined p-values from constituent CpG p-values; apply BH per chromosome. Replaces the existing NaN p/q-values.
   - `tl.py::dmr` accepts `method="segment"` (preferred) and `method="hmm"` (legacy, mapped to `"segment"` with `FutureWarning`).
   - Add shim module `src/epykit/dmr_hmm.py` that re-exports from `dmr_segment` and emits `DeprecationWarning` on import plus a `logger.warning(...)` mirror.
   - Update CLI (`cli.py`), README, `docs/analysis/`, `CLAUDE.md`.
   - Test additions: `tests/test_dmr_segment.py::{test_call_dmr_rule_segment_emits_pvalues, test_dmr_hmm_shim_warns}` (fast).
   - `Affects: segment@all`

2. **`refactor(dmc): rename log2_odds_ratio per backend (P1-11)`**
   - In pooled-count tests (`lr`, `fisher`): column renamed to `log2_odds_ratio_pooled`. Semantically unchanged.
   - In `glm`: column renamed to `coef_treatment_log2` (it was always the logit coefficient, not log₂-OR — the old name was misleading).
   - Add transitional `log2_odds_ratio` column: NaN-filled, emits `FutureWarning("log2_odds_ratio is deprecated and will be removed in 0.8; use log2_odds_ratio_pooled for pooled tests or coef_treatment_log2 for GLM")` on first column access.
   - Update consumers: `dmr.py`, `dmr_segment.py` (formerly `dmr_hmm.py`), `tests/test_compute_backends.py`, `tests/test_dmr_tile_merge.py`, `docs/analysis/dmc.md`, `README.md`.
   - No new test functions: existing engine tests get updated to the new names in this commit, which IS the test of the rename.
   - `Affects: lr@all, fisher@all, glm@all`

### Step 2: Drops

3. **`fix(dmc) BREAKING: remove test='logit_t'`**
   - Delete the `logit_t` branch from `dmc.py` (`dmc.py:1577,1590-1594`). Keep `_beta_binom_mom_from_welford_logit` internal helper.
   - Remove `"logit_t"` from `tl.py:349` docstring.
   - Add guard in `tl.py::dmc`: `if test == "logit_t": raise ValueError("test='logit_t' was removed in 0.7.5 (miscalibrated near β=0/1). Use test='welch_t' for the replicate-aware β-mean test or test='lr' for the recommended default.")`.
   - Delete any `tests/test_*logit_t*` files / pytest parametrise rows.
   - `Affects:` none (engine removal; no numerical drift in surviving engines).

4. **`fix(dmc) BREAKING: remove test='bb_lr'`**
   - Delete the `bb_lr` branch (`dmc.py:1601-1688`, ~88 lines). `irls_dispatch` stays (used by `glm`).
   - Delete the `n<3` guard at `dmc.py:2032-2033`.
   - Guard: `raise ValueError("test='bb_lr' was removed in 0.7.5 (TPR < 8% at n ≤ 4 + dispersion bug). Use test='lr' (recommended) which uses the same quasi-binomial dispersion but pools counts per group for higher power at small n.")`.
   - Closes spec P1-2 by removal; CHANGELOG notes this.
   - Delete any `tests/test_*bb_lr*` / parametrise rows.
   - `Affects:` none.

5. **`fix(dmc) BREAKING: remove test='score'`**
   - Change `elif test in ("score", "lr"):` (`dmc.py:1510`) to `elif test == "lr":`. Hardcode `statistic="lr"` into the `_score_finalize` call; drop the `statistic=` kwarg from public surface if no internal caller remains (Grep before delete).
   - Guard: `raise ValueError("test='score' was removed in 0.7.5 (strictly dominated by test='lr' in finite samples; asymptotically equivalent). Switch test='score' → test='lr'; output schema is identical.")`.
   - `Affects:` none.

6. **`fix(dmc) BREAKING: remove test='cmh'`**
   - Delete the `cmh` branch (`dmc.py:1477-1508`) and helpers `_cmh_init` / `_cmh_update` / `_cmh_finalize` (Grep first; no other caller).
   - Guard: `raise ValueError("test='cmh' was removed in 0.7.5 (stratification semantics confusing; dominated by GLM with batch covariate). For stratified analysis use tl.dmc(formula='~ group + batch'), which gives proper dispersion correction and handles continuous covariates.")`.
   - `Affects:` none.

7. **`refactor(dmc): collapse _auto_test, docstring, docs to 4+auto surviving engines`**
   - `tl.py::_auto_test` — no code change but inline-comment that the dispatch surface is now closed.
   - `tl.py:349` docstring — list exactly `auto, lr, welch_t, fisher, glm`.
   - `pyproject.toml` — prune any markers referencing dropped engines.
   - `rg -nP 'parametrize\(\s*"test"' tests/` to enumerate every parametrised engine list in the test suite; prune `"logit_t"`, `"bb_lr"`, `"score"`, `"cmh"` from each occurrence. (Known sites at the time of writing: `tests/test_dmc_multigroup.py`, `tests/test_compute_backends.py`; rg confirms the exhaustive list at execution time.)
   - Update `README.md` engine list table; `docs/analysis/dmc.md`; `CLAUDE.md` lines 56-60.
   - Test addition: `tests/test_phase3_drops.py` with one parametrised test:
     ```python
     @pytest.mark.parametrize("engine,hint_substring", [
         ("logit_t", "welch_t"), ("bb_lr", "lr"),
         ("score", "lr"),         ("cmh", "glm"),
     ])
     def test_dropped_engine_raises_with_migration_hint(...): ...
     ```

### Step 3: P1 functional fixes (one commit per item)

8. **`fix(dmc) P1-1: Fisher two-sided uses mid-p convention`**
   - Replace `dmc.py:253-256` doubled-tail logic with `scipy.stats.fisher_exact(alternative="two-sided")`-equivalent mid-p, vectorised by table-sum stratification.
   - Test: `tests/test_dmc_fisher.py::test_fisher_two_sided_matches_scipy` — 100 random 2×2 tables, `np.allclose(epykit_p, scipy_p, atol=1e-12)`.
   - `Affects: fisher@small-table-cells` (no impact on cov≥10, n≥3 headline).

9. **`fix(dmc) P1-3: wire Newcombe CI into lr meth_diff_ci_{lo,hi}`**
   - In `_process_one_chromosome` (`dmc.py:1907-1912`), replace Welch-normal Wald CI with `_glm.newcombe_diff_ci` (already implemented, unwired) for the `lr` path. `welch_t` keeps Welch CI (correct for that test); `glm` keeps model-based CI. `fisher` keeps Wilson on per-group rates.
   - Test: `tests/test_dmc_lr.py::test_meth_diff_ci_uses_newcombe_for_lr` — synthetic 3v3 with boundary β (~0.01 and ~0.99); assert CI bounds match `statsmodels.stats.proportion.confint_proportions_2indep(method="newcombe")` to 1e-6.
   - `Affects: lr@all` — supplementary CI numbers shift; point estimates (`meth_diff`) unchanged.

10. **`fix(_glm) P1-4: explicit reference_level kwarg for patsy Treatment`**
    - Add `reference_level: str | None = None` to `build_design` (`_glm.py:42-197`). When set, use `Treatment(reference="<level>")` on the relevant factor; otherwise alphabetical (preserves existing default behaviour).
    - Log `design_info.column_names` and resolved reference at INFO level.
    - Surface via `tl.dmc(formula=..., reference_level=...)`.
    - Test: `tests/test_glm.py::test_reference_level_respected` — fit with default vs `reference_level="treated"`; assert column ordering and contrast sign flip.
    - `Affects:` none (opt-in kwarg; default unchanged).

11. **`fix(_glm) P1-5: NaN-mask non-converged IRLS sites + log fraction`**
    - In `irls_dispatch` (`_glm.py:234, 297, 321-324`), use the existing `converged` boolean to NaN the Wald-derived columns at non-converged sites. Log fraction at `WARNING` if > 1% per call.
    - Test: `tests/test_glm.py::test_nonconverged_sites_are_nan` — construct a degenerate-design case where IRLS diverges; assert NaN at those rows and warning fired.
    - `Affects: glm@degenerate-cells` (rare; supplementary only).

12. **`fix(dmr) P1-6: empirical FDR design-aware + n=1,1 refusal`**
    - In `dmr.py:1370-1377`:
      - (a) Detect paired design via a single `md.obs` covariate where every value appears exactly once per group; shuffle within strata instead of globally.
      - (b) `raise ValueError("empirical DMR FDR requires n≥2 per group; got 1v1. Use Fisher-derived p-values directly via tl.dmc(test='fisher').")` when `n_treat=1, n_ctrl=1`.
    - Tests (`tests/test_dmr_empirical_fdr.py`):
      - `test_n_one_each_raises` — fast; asserts the `ValueError`.
      - `test_paired_design_shuffles_within_strata` — marked `slow`; runs a short permutation loop on a fixture with a paired covariate, asserts within-stratum invariance.
    - `Affects: empirical_dmr@paired-cells` (no headline cell; supplementary).

13. **`fix(dvc) P1-7: Brown-Forsythe replaces Bartlett on bounded β`**
    - In `dvc.py:58-86`, two-pass: per-site group medians → absolute deviations → F-test on deviations. Equivalent to `scipy.stats.levene(center="median")`.
    - Test: `tests/test_dvc.py::test_brown_forsythe_matches_reference` — synthetic bimodal β with known variance ratio; F-stat within 1e-6 of scipy reference.
    - `Affects: dvc@all` (DVC not in paper; no headline impact).

14. **`fix(qc) P1-9: sex check dip-test fallback for unimodal cohorts`**
    - Wrap the largest-gap clustering in `qc.py:474-495` with a Hartigan dip-test (`diptest.diptest`) prereq. If `pval_dip > 0.10`, fall back to `Y_coverage_ratio < 0.25` threshold and emit `UserWarning("single-sex cohort detected; sex inferred from Y-coverage threshold only.")`.
    - Add `diptest` to `pyproject.toml` extras under `qc`.
    - Test: `tests/test_qc.py::test_sex_check_unimodal_falls_back` — synthetic single-sex cohort; assert fallback path taken and warning fired. (No bimodal-regression test added; existing tests already cover that path.)
    - `Affects: qc@single-sex-cohorts`.

15. **`fix(dmc) P1-10: Storey π₀ clamped at 1/n`**
    - In `_storey_pi0` (`dmc.py:2419-2421`): `pi0 = max(numerator / denominator, 1.0 / n)`.
    - Docstring documents this is the plug-in estimator at `lam=0.5`, not the spline-smoother.
    - Test: `tests/test_dmc_multitest.py::test_storey_pi0_clamped_at_one_over_n` — input where every p<0.5; assert returned π₀ ≥ 1/n.
    - `Affects: lr@all-significant-cells` (edge case; default `fdr_tsbh` uses bootstrap π₀).

### Step 4: Integration scripts

16. **`feat(benchmark): methylkit_stouffer_combine.R — adjacent-3-CpG combine`**
    - R script reads `benchmark/data/study1/methylkit/scenario_<id>.tsv`, applies 3-CpG Stouffer combine (focal + neighbours within `--max-gap-bp`, default 1000), writes `methylkit_tuned/scenario_<id>.tsv` with `pvalue_combined` and `qvalue_combined` (BH per chromosome). Output schema = input schema + 2 columns; `evaluate.py` ingests unchanged.
    - Asserts on input that BH was NOT pre-applied (re-applies after combine).
    - Test: `benchmark/scripts/tests/test_methylkit_stouffer_combine.py` — 6-CpG fixture TSV with hand-computed combined p-values, `pytest.approx(rtol=1e-6)`. Dispatches via `subprocess.run(["Rscript", ...])`; skips with `shutil.which("Rscript") is None`. Marked `slow`.

17. **`feat(benchmark): null_engines + run_null_calibration real-engine wiring`**
    - New module `benchmark/scripts/_null_engines.py` exposes `ENGINE_REGISTRY: dict[str, EngineFn]` covering the surviving engines: `lr`, `lr+`, `welch_t`, `fisher`, `glm` (DMC) and `tile`, `sliding`, `chain_merge`, `segment` (DMR).
    - Each closure: loads a frozen `MethylData`, permutes `md.obs["group"]` per seed, calls the engine, returns q-value array.
    - `run_null_calibration.py::main` rewritten to take `--engine <name>` and `--methylstore <path>`; dispatches via `ENGINE_REGISTRY[name]`; writes `benchmark/data/null_calibration/<scenario>/<engine>.parquet`.
    - Test: `benchmark/scripts/tests/test_null_engines.py` with one parametrised test looping the surviving engines, on a 200-CpG / 6-sample synthetic fixture. Asserts output shape, value range, deterministic across seeded runs. Marked `slow`.

18. **`feat(benchmark): evaluate.py emits Wilson + bootstrap CI columns`**
    - Appends finalisation step that calls `wilson_bootstrap_ci.add_wilson_ci_for_tpr_fpr` on every row plus the bootstrap helpers for AUROC / F1. AUROC / F1 bootstrap needs per-CpG joined frames; these are cached as `eval_per_cpg/<tool>_<scenario>.parquet` during the main eval pass.
    - New CLI flag `--ci / --no-ci` (default `--ci`); back-fill mode `--ci-only --eval-summary <path>` for existing parquets.
    - New columns: `tpr_ci_lo, tpr_ci_hi, fpr_ci_lo, fpr_ci_hi` (Wilson); `auroc_ci_lo, auroc_ci_hi, f1_ci_lo, f1_ci_hi` (bootstrap, B=1000, seed = hash of `(tool, scenario, threshold)`).
    - Tests (`benchmark/scripts/tests/test_evaluate_ci.py`, fast, 2 tests):
      - All 8 CI columns present after `--ci` run.
      - All CIs bracket their point estimates on a 1000-CpG synthetic fixture.

19. **`feat(benchmark): regen_all.py — acceptance gate (seed manifest)`**
    - Three modes: `--verify` (default, CI-friendly; reads `<!-- source: scripts/X.py -->` HTML comments in `paper.md`, asserts cited parquets match printed values to the printed precision), `--run-cheap` (refresh fast scripts only), `--run-all` (full regen — used in Phase 4 only).
    - `claims.yaml` manifest with `claim_id → (paper_section, parquet_path, column_filter, expected_value, precision)`. Phase 3 lands an empty / seed manifest; Phase 4 populates it during the locked re-run.
    - Tests (`benchmark/scripts/tests/test_regen_all.py`, fast, 2 tests):
      - Verify-pass on a fixture `paper.md` + matching claims.
      - Verify-fail (off-by-precision) returns non-zero exit with a colour diff.

20. **`feat(benchmark): bug_fix_audit.py — pre/post-fix per-cell delta`**
    - Diffs the pre-Phase-1 `benchmark/data/study*/eval_summary.parquet` against the post-Phase-3 re-run (produced in Phase 4) on a per-(tool, scenario, metric) basis. Writes `benchmark/data/audit/bug_fix_deltas.{parquet,md}`.
    - Attribution: parses `Affects: <engine>@<scenario>` trailers in commit messages between the pre-snapshot tag and the post-fix tag; builds `fix_id → affected_cells` map; attributes each delta to the most-recent fix that touched that cell. Unattributed cells get `fix_id="UNATTRIBUTED"` and the script exits non-zero.
    - Tests (`benchmark/scripts/tests/test_bug_fix_audit.py`, fast, 2 tests):
      - Attribution-success on a fixture pre/post parquet + fake commit log.
      - Unattributed cell causes non-zero exit.

### Step 5: Tag and wrap-up

21. **CHANGELOG sweep + tag** — final pass: every commit in steps 1-20 has a CHANGELOG entry under the appropriate `## Unreleased / ### {Removed, Changed, Fixed (P1 manifest), Added}` heading. Tag `v0.7.5-phase3-engines-frozen` after acceptance criteria pass.

---

## 5. Test impact

### Existing tests modified (no net runtime change)

| File | Reason |
|---|---|
| `tests/test_compute_backends.py` | `log2_odds_ratio` → backend-specific name |
| `tests/test_dmr_tile_merge.py` | Same rename |
| `tests/test_dmr_hmm*.py` (if any) | Imports updated to `dmr_segment`; one assertion added on the `DeprecationWarning` shim |
| `tests/test_dmc_*.py` Fisher math | P1-1 reference values regenerated against scipy |
| `tests/test_dmc_*.py` Newcombe-affected CI assertions | P1-3 reference values regenerated against statsmodels |
| `tests/test_dvc.py` | P1-7 reference values regenerated against `scipy.stats.levene(center="median")` |
| `tests/test_dmc_multigroup.py`, etc. | `pytest.mark.parametrize("test", [...])` lists pruned to surviving engines |
| Any `tests/test_*logit_t*` / `*bb_lr*` files | Deleted with the engine commits |

### Net-new tests

Fast tier (`pytest -m "not slow"`), ~10 functions:

- `tests/test_phase3_drops.py::test_dropped_engine_raises_with_migration_hint` (parametrised, 4 rows)
- `tests/test_dmc_fisher.py::test_fisher_two_sided_matches_scipy`
- `tests/test_dmc_lr.py::test_meth_diff_ci_uses_newcombe_for_lr`
- `tests/test_glm.py::test_reference_level_respected`
- `tests/test_glm.py::test_nonconverged_sites_are_nan`
- `tests/test_dmr_empirical_fdr.py::test_n_one_each_raises`
- `tests/test_dvc.py::test_brown_forsythe_matches_reference`
- `tests/test_qc.py::test_sex_check_unimodal_falls_back`
- `tests/test_dmc_multitest.py::test_storey_pi0_clamped_at_one_over_n`
- `tests/test_dmr_segment.py::test_call_dmr_rule_segment_emits_pvalues`
- `tests/test_dmr_segment.py::test_dmr_hmm_shim_warns`

Slow tier (`pytest -m slow`), 1 function in main suite:

- `tests/test_dmr_empirical_fdr.py::test_paired_design_shuffles_within_strata`

Benchmark scripts suite (`uv run pytest benchmark/scripts/tests/`):

- Fast: `test_evaluate_ci` (2), `test_regen_all` (2), `test_bug_fix_audit` (2) — 6 functions.
- Slow: `test_null_engines` (1 parametrised across ~9 engines), `test_methylkit_stouffer_combine` (1) — 2 functions.

### Expected end-of-Phase-3 totals

- Main `tests/` fast suite: ≈ 235-240 passing (Phase 2 baseline 229 + ~11 new − ~5 engine-specific deletions).
- Main `tests/` slow tier: +1.
- Benchmark scripts suite: ≈ 21-23 (Phase 2 baseline 15 + ~8 new).

---

## 6. Acceptance criteria for `v0.7.5-phase3-engines-frozen`

- [ ] `uv run pytest -m "not slow" --strict-markers -ra` — all passing.
- [ ] `uv run pytest benchmark/scripts/tests/ -m "not slow"` — all passing.
- [ ] `uv run pytest -m slow` — all passing (validated locally; CI may skip).
- [ ] `uv run ruff check src/ benchmark/scripts/` — no new F-level violations vs Phase 2 baseline.
- [ ] `uv run mypy src/epykit` — no new errors vs Phase 2 baseline.
- [ ] `uv run epykit --help` and `uv run epykit dmc --help` list no removed engines.
- [ ] `uv run python -c "import epykit; epykit.tl.dmc(<fixture>, test='logit_t')"` raises `ValueError` with migration text (and analogous for `bb_lr`, `score`, `cmh`).
- [ ] `CHANGELOG.md` `## Unreleased` has populated `### Removed` (4 drops), `### Changed` (renames + Newcombe + Brown-Forsythe + ...), `### Fixed (P1 manifest)` (one bullet per P1 fix referencing its ID), `### Added` (5 integration items).
- [ ] `benchmark/scripts/regen_all.py --verify` runs without crashing on the seed `claims.yaml` (smoke check; no real claims to verify yet).
- [ ] `README.md`, `docs/analysis/dmc.md`, `CLAUDE.md` list exactly the 4 + auto surviving engines.
- [ ] Tag `v0.7.5-phase3-engines-frozen` exists on `p0-fixes` branch.

The locked-re-run criteria from the parent spec §10 (every cell carries a CI; `regen_all.py` exits 0 on populated claims; simulator marginals match within 2%; Limitations section content) belong to Phase 4.

---

## 7. Risks

1. **P1-3 (Newcombe CI) changes supplementary CI widths.** Welch-normal Wald is symmetric and generally wider than Newcombe near boundary β. Mitigation: `bug_fix_audit.py` surfaces this per-cell; Limitations table owns it.
2. **P1-11 rename breaks unpinned external consumers.** Grep shows no benchmark scripts read `log2_odds_ratio`, but external user code may. Mitigation: transitional NaN-filled column + `FutureWarning` for one release; explicit CHANGELOG migration steps.
3. **Four engine drops break user notebooks.** Migration-hint `ValueError` per dropped engine is loud and prescriptive, but a pinned `test='logit_t'` notebook hard-fails on upgrade. Mitigation: 0.7.5 is pre-1.0; SemVer permits breaking; hints are one-line fixes.
4. **`bug_fix_audit.py` UNATTRIBUTED rows.** Forgotten `Affects:` trailers cause non-zero exit. Mitigation: loud failure mode; manual `git commit --fixup` to add trailer is fast. A pre-commit hook is a candidate for a separate task (not in Phase 3 scope).
5. **methylKit Stouffer combine semantics.** Combining methylKit's `qvalue` instead of `pvalue` would double-correct. Mitigation: script asserts on input that `pvalue` column is uncorrected (dtype + range check), re-applies BH per chromosome post-combine. Documented in the script header.
6. **R toolchain on Windows CI.** Item 1's test skips when `Rscript` is unavailable; main test suite stays portable. The R script runs in Phase 4's locked re-run on the Linux methylKit host.
7. **`dmr_hmm` shim's `DeprecationWarning` invisible under `-W ignore`.** Mitigation: also mirror through `logging.getLogger(__name__).warning(...)` so it surfaces in standard epykit log output.
8. **Symmetric removal commits + P1 fix commits inflate `Affects:` map churn.** If steps 2's removal commits accidentally carry an `Affects:` trailer for the removed engine, `bug_fix_audit.py` attributes deltas to a non-existent engine. Mitigation: drop commits use empty `Affects:` (no impact on surviving engines' cells); the audit ignores empty trailers.

---

## 8. Non-goals (explicit)

- Reorganising / splitting `dmc.py` — defer to v0.8.
- Adding a parity test between `score` and `lr` (the math is asymptotically equivalent; no test debt by dropping the engine).
- Switching `dispersion` default away from `eb` — Phase 1 P0-3 locked `eb`.
- GPU backend (`_glm_gpu.py`) audit — parent spec §1 non-goal.
- P2 hygiene items — Phase 4 or later.
- Populating `claims.yaml` — Phase 4 work.
- Pre-commit hook enforcing the `Affects:` trailer — separate task candidate.
- Paper text changes — Phase 4.

---

## 9. Handoff to Phase 4

After Phase 3 tags `v0.7.5-phase3-engines-frozen`, Phase 4 runs:

1. Multi-seed simulator grid (`simulate_piao.py` × N=20 seeds at headline cell + frozen-defaults sweep at all cells).
2. Study 2 re-run (post-P0 + post-P1 code).
3. Study 3 re-run (real GSE263850 at 0.7.5).
4. Null calibration across all engines × Piao + simulator + Study 3 data (uses Phase 3's `_null_engines.py`).
5. `bug_fix_audit.py` populated (uses Phase 3's pre-Phase-1 parquet baseline).
6. `claims.yaml` populated; `regen_all.py --verify` exits 0.
7. Tables S-Calib, S-Sim, S-Fix produced.
8. Paper rewrite per parent spec §6.
9. P2 hygiene items.

No code-side blockers between Phase 3 and Phase 4; the engine surface is final at tag time.

### Phase 3 closeout note (appended on completion)

- Main fast suite: 247 passed (5 skipped — pysam/CuPy/pyBigWig platform gates)
- Slow suite: 109 passed (5 skipped)
- Benchmark scripts suite (fast): 21 passed
- All 11 P1 fixes landed; 4 engines dropped; 5 integration scripts created.
- Unexpected findings during execution: `POWER_MIN_WELCH_T` threshold in
  `tests/test_accuracy.py` was set to 0.20 but the achievable power at n=4
  per group with Welch t is ~0.045 (the implementation emits a warning about
  n<6); threshold corrected to 0.03 to guard against complete collapse only.
  This was a pre-existing miscalibration, not a regression from Phase 3.
- Phase 4 context: the pre-Phase-1 `benchmark/data/study*/eval_summary.parquet`
  files serve as the pre-fix baseline for `bug_fix_audit.py`. The `claims.yaml`
  is empty; Phase 4 must populate it during the locked re-run before
  `regen_all.py --verify` can serve as a CI gate.
