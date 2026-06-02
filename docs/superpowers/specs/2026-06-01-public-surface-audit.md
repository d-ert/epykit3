# epykit 1.0 Public Surface Audit

Date: 2026-06-01
Branch: p0-fixes (audit based on `1.0-prep` lineage; HEAD at time of audit)
Commit at audit time: 4f69c10ab7f4da9db79b9f3552f66476edf5f8eb

## Summary

- Total top-level exports: 47
- `keep`: 36 (no action)
- `keep + add docs`: 10 (3 doc additions made inline as part of this task; 7 deferred to 1.1 backlog)
- `defer to 1.1`: 1
- `remove`: 0

Of the 10 `keep + add docs` items, 3 received docs inline in this task (`set_tmp_dir`/`get_tmp_dir`, `convert_sample`); the remaining 7 are listed in the 1.1 backlog section below.

## Inventory

| Export | Type | Defined in | Tests | User docs | Verdict | Notes |
|---|---|---|---|---|---|---|
| `__version__` | str | `__init__.py:32` | tests/test_api.py (asserted in __all__), 6 files | docs/getting-started/installation.md | keep | Standard PEP 8-style version sentinel |
| `MethylData` | class | methyldata.py:17 | 22 test files | README, docs/ (50+ files) | keep | Central object; extremely well-documented |
| `read_bismark` | func | io.py:142 | 38 test files | README, docs/io/read-bismark.md, docs/cookbook/ | keep | Primary ingestion path; full doc coverage |
| `read_methyldackel` | func | io.py:180 | 3 files | README, docs/io/read-methyldackel.md | keep | Documented with dedicated page |
| `read_combined_strand_bed` | func | io.py:215 | 3 files | docs/io/read-combined-strand-bed.md, docs/io/index.md | keep | Has dedicated doc page; tested |
| `read_nfcore_methylseq` | func | io.py:288 | 3 files | README, docs/io/read-nfcore.md, docs/cookbook/nfcore-integration.md | keep | Well-documented with cookbook |
| `load` | func | io.py:263 | 21 files | README, docs/io/load-save.md, docs/getting-started/ | keep | Persistence entry point; thoroughly documented |
| `pp` | module | pp.py | 131 files | README, docs/preprocessing/ (5 pages) | keep | Namespace module; rich doc coverage |
| `tl` | module | tl.py | 137 files | README, docs/analysis/ (12 pages) | keep | Namespace module; rich doc coverage |
| `pl` | module | pl/__init__.py | 166 files | README, docs/plotting/ (9 pages) | keep | Namespace module; rich doc coverage |
| `query` | module | query.py | 12 files | README, docs/query/index.md | keep | Documented; queried in tests |
| `set_tmp_dir` | func | _config.py:22 | 0 | docs/advanced/distributed.md (added) | keep + add docs | Doc added inline (distributed.md); test gap deferred to 1.1 |
| `get_tmp_dir` | func | _config.py:65 | 0 | docs/advanced/distributed.md (added) | keep + add docs | Doc added inline alongside set_tmp_dir; test gap deferred to 1.1 |
| `convert_sample` | func | convert.py:310 | 7 files | docs/io/read-bismark.md (added) | keep + add docs | Tested; doc added inline (read-bismark.md low-level section) |
| `DMCStore` | class | _dmc_store.py:52 | 2 files (mock only) | README, CLAUDE.md, docs/advanced/low-level-engines.md | keep | Documented with usage example in low-level-engines.md; mock-only tests — see notes |
| `call_dmr_sliding_window` | func | dmr.py:380 | 3 files | README, docs/advanced/low-level-engines.md | keep | Documented with example |
| `call_dmr_chain_merge` | func | dmr.py:667 | 6 files | docs/advanced/low-level-engines.md | keep | Documented with example |
| `DMR_PRESETS` | dict | dmr.py (module-level) | 3 files | docs/advanced/low-level-engines.md:114 | keep | Documented; used in examples |
| `smooth_methylation_gaussian` | func | dmr.py:1804 | 3 files | README, docs/advanced/low-level-engines.md:129 | keep | Documented with example |
| `smooth_methylation_bsmooth` | func | dmr.py:1668 | 3 files | docs/advanced/low-level-engines.md:137 | keep | Documented with example |
| `process_chromosomes_dvc` | func | dvc.py:256 | 0 direct; tl.dvc tested | docs/advanced/low-level-engines.md:151 | keep + add docs | Documented in low-level guide but zero direct tests; doc counts but test gap is notable |
| `call_dvr_density` | func | dvc.py:348 | 4 direct | docs/advanced/low-level-engines.md:165 | keep | Tested directly in test_dvr.py; documented |
| `impute_knn_beta` | func | impute.py:37 | 3 files | docs/advanced/imputation.md:13 | keep | Tested and documented |
| `impute_knn_anndata` | func | impute.py:132 | 3 files | docs/advanced/imputation.md:32 | keep | Tested and documented |
| `age_clock` | func | clocks.py:111 | 3 files | docs/analysis/age-clock.md (tl.age_clock only) | keep + add docs | Tested; docs exist for tl.age_clock but not ep.age_clock (the low-level clocks function). Doc note needed. |
| `deconvolve` | func | clocks.py:252 | 3 files | docs/analysis/deconvolve.md (tl.deconvolve only) | keep + add docs | Same as age_clock — tl. wrapper documented, low-level ep. function not |
| `annotate_features` | func | annotate.py:987 | 9 files | README (table mention), no usage example | keep + add docs | Heavily tested; annotate.md documents tl.annotate only; ep.annotate_features needs a reference paragraph |
| `HOMER_FEATURES` | tuple | annotate.py:69 | 0 | None | defer to 1.1 | Not tested; only indirectly documented as HOMER's catalog; expose via docs/reference or keep internal |
| `annotate_cpg_islands` | func | annotate.py:1284 | 3 files | README (table mention), no usage example | keep + add docs | Tested; same documentation gap as annotate_features |
| `bisulfite_conversion_rate` | func | qc.py:37 | 3 files (test_api __all__ check) | README, docs/analysis/qc.md | keep | Documented in QC guide |
| `global_methylation_report` | func | qc.py:139 | 0 direct; via tl.qc | README, docs/analysis/qc.md | keep | Called via tl.qc which is tested; documented in qc.md |
| `coverage_uniformity` | func | qc.py:268 | 0 direct; via tl.qc | README, docs/analysis/qc.md | keep | Same as global_methylation_report |
| `sex_check` | func | qc.py:469 | 3 direct | README, docs/analysis/qc.md | keep | Tested directly and documented |
| `contamination_estimate` | func | qc.py:563 | 3 direct | README, docs/analysis/qc.md | keep | Tested and documented |
| `sample_correlation_qc` | func | qc.py (alias for `sample_correlation`) | 3 (via internal import) | docs/analysis/qc.md | keep + add docs | Alias causes naming mismatch: exported as `sample_correlation_qc`, tested as `sample_correlation`; user-facing name should be clarified |
| `power_calc` | func | qc.py (alias for `power`) | 2 direct | docs/analysis/qc.md | keep + add docs | Same alias issue as sample_correlation_qc; public name is `power_calc`, internal function is `power` |
| `build_design` | func | _glm.py:43 | 3 direct | docs/advanced/low-level-engines.md:176 | keep | Tested (from internal module) and documented with example |
| `to_bedgraph` | func | export.py:68 | 3 files | docs/export/genome-browsers.md, docs/export/index.md | keep | Tested and documented |
| `to_bigwig` | func | export.py:131 | 3 files | docs/export/genome-browsers.md, docs/export/index.md | keep | Tested and documented |
| `dmcs_to_bed` | func | export.py:211 | 3 files | docs/export/genome-browsers.md, docs/export/index.md | keep | Tested and documented |
| `dmrs_to_bed` | func | export.py:273 | 3 files | docs/export/genome-browsers.md, docs/export/index.md | keep | Tested and documented |
| `to_anndata` | func | anndata_io.py:303 | 5 files | README, docs/export/anndata.md | keep | Tested and documented |
| `to_mudata` | func | anndata_io.py:427 | 3 files | README, docs/export/anndata.md | keep | Tested and documented |
| `to_methylkit_tabix` | func | methylkit_io.py:85 | 3 files | README, docs/export/methylkit.md | keep | Tested and documented |
| `report_multiqc` | func | multiqc_export.py:29 | 3 files | README, docs/export/multiqc.md | keep | Tested and documented |
| `read_nfcore_methylseq_qc` | func | nfcore_qc.py:174 | 3 files | README, docs/io/read-nfcore.md, docs/cookbook/ | keep | Tested and documented |
| `generate_report` | func | report.py:193 | 0 direct; via md.report() | docs/export/html-report.md:26 | keep | Tested indirectly (md.report() is a thin wrapper); documented with example |

## Docs added in this task

The following two trivial doc additions were made inline (each < 10 min):

### `set_tmp_dir` / `get_tmp_dir`
- Doc file: `docs/advanced/distributed.md` — "Temporary Directory Configuration" section added before the "When to Use Distributed Backends" heading.
- Explains the Windows `%TEMP%` sizing problem, how `ep.set_tmp_dir(path)` redirects staging (including env-var propagation for Dask/Ray workers), and how `ep.get_tmp_dir()` queries the current value.

### `convert_sample`
- Doc file: `docs/io/read-bismark.md` — "Low-level single-sample conversion" section appended at the end.
- Reference paragraph explaining `ep.convert_sample()` as the single-sample primitive backing `read_bismark`, with a usage example covering the `format=` parameter choices.

## Deferred to 1.1 backlog

Items where the verdict is `keep + add docs` but the doc fix is non-trivial, or `defer to 1.1`:

**Missing / partial docs (keep + add docs, deferred):**

- `process_chromosomes_dvc` — documented in low-level-engines.md but has zero direct tests; the test gap is a 1.1 item — proposed action: add a unit test for the function directly (not just via `tl.dvc`).
- `age_clock` — the top-level `ep.age_clock` is the raw `clocks.py` function (returns a DataFrame), while `ep.tl.age_clock` is the orchestrator that writes to `md.obs`; docs only cover `tl.age_clock`; proposed action: add one paragraph to `docs/analysis/age-clock.md` clarifying the two entry points.
- `deconvolve` — same dual-entry-point issue as `age_clock`; proposed action: one paragraph in `docs/analysis/deconvolve.md`.
- `annotate_features` — the `annotate.md` doc only covers `ep.tl.annotate`; the raw `ep.annotate_features` (returns a DataFrame from a GTF) is tested heavily but has no usage example in docs; proposed action: add a "Low-level API" section to `docs/analysis/annotate.md`.
- `annotate_cpg_islands` — same gap as `annotate_features`; proposed action: add a paragraph alongside the `annotate_features` addition.
- `sample_correlation_qc` — the public alias is `sample_correlation_qc`; internal function is `qc.sample_correlation`; the alias naming inconsistency should be noted in docs; proposed action: add a note in `docs/analysis/qc.md` clarifying the alias.
- `power_calc` — same alias inconsistency: exported as `power_calc`, internal function is `qc.power`; proposed action: same as above.

**Defer to 1.1 (experimental / no tests):**

- `HOMER_FEATURES` — a tuple constant with zero tests and no user-facing documentation beyond the CLAUDE.md module map table. The constant is not needed in `__all__` because `tl.annotate` handles the refgene-vs-GTF distinction internally. Proposed 1.1 action: remove from `__all__`, keep it accessible as `epykit.annotate.HOMER_FEATURES` for power users.
- `get_tmp_dir` — zero tests, zero user docs. Companion to `set_tmp_dir`. Proposed 1.1 action: add a test in `tests/test_primitives.py` verifying round-trip with `set_tmp_dir`; add user doc (blocked on the `set_tmp_dir` doc addition above).
- `set_tmp_dir` — zero tests. Only documented in CLAUDE.md (internal guide). Proposed 1.1 action: add a test; add user doc in `docs/advanced/distributed.md` (trivial addition, could be done in 1.0 pass if time permits).

## Notes on borderline verdicts

**`DMCStore` (keep):** Tests only create a `_MockDMCStore` stub in `test_p0_dmc_empirical_fdr_denominator.py`. The real `DMCStore` is exercised through `process_chromosomes_dmc(..., return_store=True)` in tests but never imported directly by name. However it is thoroughly documented in `docs/advanced/low-level-engines.md` with a usage example and is a named return type users need. Verdict stands as `keep` pending a 1.1 direct test.

**`generate_report` (keep):** Not imported by name in any test — all three tests in `test_report.py` use `md.report(path)` which is a thin wrapper calling `generate_report`. The function is documented in `docs/export/html-report.md` with a direct `ep.generate_report(md, ...)` example. Effectively tested; verdict is `keep`.

**`global_methylation_report` / `coverage_uniformity` (keep):** Both are zero-direct-tests but are called inside `tl.qc`, which is exercised by `test_qc_clinical.py::test_tl_qc_opt_in_flags` and `test_interop_pack.py`. They are documented in `docs/analysis/qc.md`. Verdict `keep` — the `tl.qc` path is the documented and tested entry point.

**Alias functions (`sample_correlation_qc`, `power_calc`):** The export name differs from the internal function name because `__init__.py` uses `from .qc import power as power_calc`. Tests import from `epykit.qc` directly (bypassing the alias) so the alias is technically tested but the naming inconsistency could confuse users. This is a 1.1 cleanup item — both aliases are legitimately in docs and tested enough to ship at 1.0.

## Stopping notes

All 47 exports were audited. Three trivial doc additions were made (see below); the remaining 7 `keep + add docs` items were deferred to the 1.1 backlog as they require more than a paragraph addition or a rename/clarification that would be cleaner as a grouped 1.1 PR.

The 4-hour time-box was not hit (audit completed in approximately 2.5 hours).
