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

- [ ] Write SciPy likelihood/optimizer oracle tests, binomial-limit and
  boundary tests, zero-coverage/invalid-count tests; observe missing-feature failure.
- [ ] Implement stable count likelihood and safeguarded mean fitting.
- [ ] Verify scalar/batch equality, label swaps and coverage-imbalance oracle.

### Task 2: Dispersion and inference policy

Files: numerical module, `tests/test_beta_binomial_inference.py`,
`development/beta_binomial_20261007/{calibrate,validate}.py` in benchmark project.

- [ ] Test count-likelihood prior fitting and dispersion estimation against
  independent numerical objectives before implementation.
- [ ] Fit and evaluate the predefined prior/reference candidates on new
  development seeds; record assumptions, power and null tails.
- [ ] Freeze the selected model and implementation revision in the spec;
  run untouched evaluation seeds and retain full strata and finite-sample limits.
- [ ] Implement likelihood-based boundary intervals and oracle tests.

### Task 3: Public engine integration

Files: `_dmc_engines.py`, `dmc.py`, `_dmc_config.py`, `_dmc_stages.py`, `tl.py`,
existing CLI engine choices, documentation and end-to-end tests.

- [ ] Write observable API/store/CLI tests and unsupported-request tests;
  verify they fail before implementation.
- [ ] Register and run the engine with fitted means/effects/intervals and
  provenance; preserve all legacy computations and region code.
- [ ] Run relevant tests and the package suite; identify every baseline failure.
- [ ] Obtain independent review, verify concerns, repair substantive issues.

### Task 4: Full benchmark and analysis

Files: new frozen runner and orchestration under benchmark
`development/beta_binomial_20261007`; outputs under a new results directory.

- [ ] Verify old input/harness contract and snapshot candidate package+runner.
- [ ] Run only epykit on full signal, full null and GSE64177; verify all finish.
- [ ] Reuse frozen competitors and compare old/region-only/new epykit CpG,
  DMR, null, runtime/RSS, boundary and paper metrics.
- [ ] Independently cross-check paper overlaps; produce inspectable tables,
  report and provenance with supported conclusions and remaining limitations.
- [ ] Audit user requirements against actual files/status/tests before completion.
