# Beta-binomial engine implementation plan

> **For agentic workers:** Use superpowers:executing-plans; implement inline,
> then obtain an independent statistical and code review.

**Goal:** Add a genuine replicated-count beta-binomial engine and complete
the epykit benchmark and all-recorded-tool comparison with that engine.

**Architecture:** Separate numerical likelihood/fitting from chromosome I/O
and existing result publication. Learn dispersion from counts, freeze the
inference policy after independent development simulations, and run an
immutable package snapshot with the accepted region correction unchanged.

**Tech Stack:** Python >=3.10, existing NumPy/SciPy/Numba/Polars, pytest.

**Spec:** `docs/superpowers/specs/2026-10-07-beta-binomial-engine-design.md`

## Global constraints

- Preserve all baseline results and original input/harness hashes.
- Keep accepted `_region_search.py` and `dmr.py` bytes unchanged.
- Never fit or choose inference settings using benchmark truth or paper DMRs.
- No optional permutation development or region threshold retuning.
- Reject unsupported paired/covariate, smoothing and single-replicate use.

## Review focus

- Boundary methylation: valid likelihood fits and nonzero interval uncertainty.
- Unequal coverage: retain individual counts and use biological variance.
- Missing coverage: omit it without converting it to observed zero methylation.
- Numerical or prior fitting failures: visible errors/status, no silent old engine.
- Repeated/cached runs: source revision and engine/prior policy provenance.

### Task 1: Exact likelihood and fitted means

Files: create `src/epykit/_beta_binomial.py`,
`tests/test_beta_binomial_likelihood.py`.

Interfaces: `logpmf(m,n,mu,rho)`, `fit_mean(m,n,rho)`,
`fit_fixed_dispersion(m,n,n_case,rho)` for count vectors and batch fitting.

- [x] Write SciPy likelihood/optimizer oracle tests, binomial-limit and
  boundary tests, zero-coverage/invalid-count tests; observe missing-feature failure.
- [x] Implement stable count likelihood and safeguarded mean fitting.
- [x] Verify scalar/batch equality, label swaps and coverage-imbalance oracle.

### Task 2: Dispersion and inference policy

Files: numerical module, `tests/test_beta_binomial_inference.py`,
`development/beta_binomial_20261007/{calibrate,calibrate_continuous,validate_holdout}.py`
in the benchmark project.

- [x] Test count-likelihood prior fitting and dispersion estimation against
  independent numerical objectives before implementation.
- [x] Fit and evaluate the predefined prior/reference candidates on new
  development seeds; record assumptions, power and null tails.
- [x] Freeze the selected model and implementation revision in the spec;
  run untouched evaluation seeds and retain full strata and finite-sample limits.
- [x] Implement likelihood-based boundary intervals and oracle tests.

### Task 3: Public engine integration

Files: `_dmc_engines.py`, `dmc.py`, `_dmc_config.py`, `_dmc_stages.py`, `tl.py`,
existing CLI engine choices, documentation and end-to-end tests.

- [x] Write observable API/store/CLI tests and unsupported-request tests;
  verify they fail before implementation.
- [x] Register and run the engine with fitted means/effects/intervals and
  provenance; preserve all legacy computations and region code.
- [x] Run relevant tests and the package suite; identify every baseline failure.
- [x] Obtain independent review, verify concerns, repair substantive issues.

### Task 4: Full benchmark and analysis

Files: new frozen runner and orchestration under benchmark
`development/beta_binomial_20261007`; outputs under a new results directory.

- [x] Verify old input/harness contract and snapshot candidate package+runner.
- [x] Run only epykit on full signal, full null and GSE64177; verify all finish.
- [x] Reuse frozen competitors and compare old/region-only/new epykit CpG,
  DMR, null, runtime/RSS, boundary and paper metrics.
- [x] Independently cross-check paper overlaps; produce inspectable tables,
  report and provenance with supported conclusions and remaining limitations.
- [x] Audit user requirements against actual files/status/tests before completion.

Completion evidence (8 October 2026): all three run statuses and all launcher
exit codes succeeded in `wgbs_benchmark_v4/results/beta_binomial_20261007_full`.
Its comparison directory contains the report, workbook, independent
simulation validation, GenomicRanges paper-overlap validation and artifact
validation. The original contract and frozen package hashes match. Package
verification and all 28 pre-existing suite failures are recorded separately
under `development/beta_binomial_20261007/verification.json`.

The experiment is complete; flagship readiness is not established. BB-F has
no significant CpG discoveries in the signal genome and lower region recall
than the region-only correction. Retain experimental, explicit selection.
