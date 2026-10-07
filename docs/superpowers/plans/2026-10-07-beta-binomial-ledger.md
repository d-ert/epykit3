# Beta-binomial execution ledger

- Base: f73e66d in existing isolated worktree; region fix accepted and immutable.
- Preserved unfinished quasi-binomial guardrail in benchmark development archive;
  restored dmc.py to base and archived its tests. No old result changed.
- Literature/design review: supplied survey, DSS/RADmeth/methylSig sources read;
  independent reviewer supports structured count-likelihood EB and cautions
  that conditional MAP LRT + F reference is an unvalidated approximation.
- Prior design refinement: factorized mean-mixture × rho-mixture, avoiding
  unidentifiable endpoint rho allocation in an unrestricted joint mixture.
- Task 1: 12 core oracle tests observed RED (missing module), then GREEN;
  likelihood_green.log records 12 passed. Package integration not yet changed.
- Tasks 2–3 implemented: continuous log-rho MAP, conditional F reference,
  difference-profile intervals, public Python/CLI engine registry, raw-count
  guards and output-owned, hashed prior provenance. Legacy auto unchanged.
- Numerical and public tests: latest 36 passed (provenance_sparse_green.log);
  earlier 55 focused existing/new tests passed (focused_engine_tests.log).
- Independent review reproduced high-coverage boundary crash, overwritten
  priors and sparse-union sampling failure. Each observed RED and then GREEN
  after repair. `bb_df` is documented as F denominator df, including chi2 mode.
- Full suite and untouched-seed calibration running. Prior-start sensitivity
  retained as a limitation to quantify, not hidden by convergence flags.
- Dedicated benchmark runner prepared: only test and CI omission differ;
  region geometry/correction and original input/harness contract immutable.
- Task 4 incomplete. No claim of controlled FDR or flagship readiness.
- Final review also reproduced case/control label dependence when a shared
  RNG sampled case first. Independent group RNGs with the same fixed seed
  restore prior symmetry. Revision v3 records the change; the dedicated
  label-swap test observed RED before repair. Numerical p-value tolerance
  is 1e-6 near LR=0, where the F tail magnifies floating-point roundoff;
  significance decisions must also agree exactly.
- Independent reviewer confirmed all four repairs, readonly inputs and
  cache invalidation after prior tampering. No implementation blocker left.
- V2 heldout evaluation is archived/superseded; v3 uses untouched seeds
  starting 71082001. No inference/reference parameter is tuned using them.
- Final focused verification after repairs: 71 passed, 17 warnings in 128.99s
  (final_focused.log), including read-only inputs, label reversal, existing
  statistics, engine registry and metadata. CLI help lists the engine in
  both commands. The first full suite had 789 passed, 29 failed, 6 skipped:
  28 failures already existed; one obsolete rejection test was removed.
  A subsequent full suite captured an overly strict near-null tolerance
  before its correction; the final focused suite verifies that correction.
