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
- Full suite and untouched-seed calibration completed. Prior-start sensitivity
  retained as a quantified limitation, not hidden by convergence flags.
- Dedicated benchmark runner prepared: only test and CI omission differ;
  region geometry/correction and original input/harness contract immutable.
- Task 4 subsequently completed on 8 October; evidence and findings below.
  No claim of controlled FDR or flagship readiness.
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
- Full benchmark used immutable package commit b5f1135 and completed signal,
  null and all twelve GSE64177 samples. Only epykit was rerun; the real run had
  no benchmark memory cap. All benchmark/analysis/job exit codes are zero.
  An earlier interrupted attempt is retained separately and excluded from
  reported timings.
- Signal: 853 regions, 779 matches at 50% reciprocal eligible-CpG overlap,
  precision 91.3%, recall 26.9%, F1 0.416; 58.99 minutes, 6.13 GiB peak RSS.
  The region-only fix had F1 0.460 in 6.27 minutes. dmrseq had F1 0.559 and
  DSS 0.433. All five recorded competitors are included in the report.
- CpG signal: zero discoveries at BH q<=0.05, against 34,889 eligible true
  positives. AP 0.05988 exceeds original epykit's 0.03565 but remains far below
  DSS's 0.37514. Null: zero significant CpGs and zero regions. Zero discoveries
  do not establish useful power or calibrated FDR.
- Real: 48 regions, all hypo; 46 overlap the paper in agreeing direction,
  covering 47 distinct paper regions. At 50% reciprocal base-pair overlap,
  35 one-to-one matches versus 39 for the region fix and 130 for DSS.
  Runtime 41.46 minutes and peak RSS 4.83 GiB. Agreement is descriptive;
  the inherited analysis omits donor pairing and four real competitors are
  still unavailable.
- Independent verification recomputed complete CpG memberships, BH, AP and
  50%/80% region matching; maximum BH difference 4.44e-16. GenomicRanges
  confirmed paper counts for all four call sets and every cutoff. All 70
  prior hashes, source snapshots, workbook sheets and report assets verified.
  Standalone figures visually inspected. Results:
  `wgbs_benchmark_v4/results/beta_binomial_20261007_full/comparison/report.html`.
- Outcome: retain beta_binomial as experimental and opt-in. The current F
  reference, dispersion uncertainty, prior identifiability, spatial borrowing
  and runtime require further work before default promotion. The genome and
  paper results were not used to tune this frozen experiment.
