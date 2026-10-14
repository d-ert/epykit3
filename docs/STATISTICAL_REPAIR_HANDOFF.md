# Epykit statistical repair — final handoff

Date: 2026-10-14. Branch: `stats-repair-20261014`.
Base: `0bc083e`; original external checkout reviewed at `3302edf`.
Inference source last changed at `d91378f`; subsequent commits archive the harness,
fix its failure handling, clean test lint, and record this handoff.

The user stopped all three final benchmark processes, then requested this report
and committed changes instead of more benchmark execution. No benchmark restart
is authorized. Final sweep/holdout completion is superseded by that instruction;
its partial results remain preserved and explicitly incomplete.

## Relationship to the Argon conversation

The repair follows the goals of label-safe missing-data inference, accurate
boundaries, variance estimation and useful recovery. It does not reproduce every
proposed formula. Several Argon suggestions were the source of statistical errors.

| Argon topic | Delivered behavior |
|---|---|
| GLM `pearson_stat` NameError | Already corrected in the starting implementation; retained and exercised by tests. |
| Replace observed-group mean imputation | Retained label-independent centering, added explicit observation masks and an actual observed-data GLS design/covariance with correct residual degrees of freedom. A valid zero residual is never used to identify missing data. |
| Exchangeable permutations | Start from a label-independent union and replay group coverage eligibility, discovery and fitting for every label assignment; preserve chromosomes missing from the first sample. |
| Exact Kadane boundaries | Replaced unconstrained reset-and-fallback code with admissible prefix minima enforcing both CpG count and genomic span. Independent exhaustive checks cover the optimum. |
| EB log-F estimation | Replaced the unequal-df IQR standardization/median approximation with expectation-corrected log moments and safeguarded trigamma inversion; removed the prior scale20 and df120 caps; handle the infinite-prior limit. |
| Spatial variance stabilization | Removed the length-dependent multiplier that only reduced standard errors. Regional SEs now use the unmodified fitted REML variance. No new regional EB variance prior or arbitrary variance floor is claimed. |
| Anchor counting + PRDS-BH | Rejected the unsupported guarantee: count all geometry-admissible intervals and restore BY; return the conservative chain-merge seed preset. This still requires valid regional p-values and empirical model validation. |

## Additional changes

- Exact masked CAR/donor GLS uses observed sample/locus cells, locus nuisance
  intercepts, a shared treatment coefficient and an observed residual donor
  eigensystem. Independent dense GLS tests cover missingness, coefficient SE,
  degrees of freedom, profiling and complementary labels. Complete data retain
  the analytical native path.
- Pooled LR DMC moderation now accounts for depth-dependent beta-binomial
  design effects instead of pooling one constant dispersion across read depths.
  This is moment-based and its small-sample count reference remains approximate.
- Correct finite-contributor neighbour counts for NaN singletons, and report raw
  methylation differences separately from discovery support weighting.
- Add precision-weighted tricube discovery at window/4, window/2 and window
  radii, independently trim candidates, and deterministically select disjoint
  candidates. No real-paper overlap tuning was used.
- Regional default `combined_qvalue` now contains experimental complete-scan
  expected-FDP estimates, including empty scans. This changes its interpretation
  from the prior pooled-region BH calculation. Legacy calibration remains in the
  internal interface. The descriptive pooled score tail is not an exact regional
  randomization p-value.
- Add whole-scan `maxT_pvalue` and optional public
  `empirical_fdr=True, fdr_method="max_t"`. Its interpretation is joint all-null
  error under exchangeable samples, not established strong control under partial
  nulls. Balanced3+3 data have minimum exact two-sided p=2/20=0.1.
- Reject numerically exact fits symmetrically using explicit design residuals
  and a response-relative floating-point tolerance. Reject maxT selection with
  observed-pooled decoy covariance; invalidate that mode's maxT diagnostic.
- Change algorithm/cache revisions so prior inference caches cannot silently
  supply results for changed methods.
- Archive resumable benchmark runners, source/evaluation snapshots, dependency,
  count-store, truth/reference and output hashes; verify entrypoint and evaluation
  identities on resume. Correct failed null verification so it cannot enter a
  successful summary. The final active processes executed the earlier frozen
  harness; their successful outputs were independently integrity-checked. The
  failure-path correction is separately tested and included in the handoff.

## Verification

- Source suite:96 passed, with two existing `log2_odds_ratio` deprecation warnings.
- Benchmark/simulator workspace suite:86 passed in the Epykit environment.
- Ruff passes all changed source and test files after small test import/unused
  variable/pairwise cleanups.
- Original review and follow-up are included in `docs/reviews/`. The original
  three Important findings were resolved; the follow-up null-summary failure
  received an actual-main injected-verification-failure regression and repair.
- Preserved result outputs and frozen package files were checked against their
  recorded SHA-256 hashes. This is an integrity check, not proof of calibration.

## What the completed runs show

The source11163ca candidate completed72 pilot/configuration runs, four whole-genome
real configurations and80 development-null runs. Final inference source d91378f
adds numerical degeneracy/maxT safeguards and has partial reruns:39 completed
pilots, one completed real configuration and78 completed holdout/configuration
runs. There were100 held-out datasets planned across two configurations (200
caller runs);78 is not78 independent completed datasets.

Compared with the first corrected source003b8b4, candidate mean empirical recall
at coverage>=10 without boundary refinement increased14.3%→21.4%; at coverage>=5
with boundary0.05 it increased28.6%→34.5%. Corresponding generative recall increased
9.5%→15.5% and20.2%→23.0%.

On real GRCh38 human GSE263850 data at coverage>=10 without boundary refinement,
first-corrected936calls/216paper-overlap pairs/81reciprocal matches became
2341/328/102. The one completed final-source real configuration reproduces those
counts. However, the earlier imputation-era stage2 source had4940/462/162: the
repair does not beat that reference on every recovery metric. Paper concordance
is not ground truth, and reduced significance after correcting an invalid model
is not by itself evidence of a worse method.

Null evidence does not establish nominal5% FDR. For candidate empirical nulls,
GLS had discoveries in3/20datasets in each tested configuration. Depth-aware DMC
at coverage>=10 reduced empirical null-positive runs2/20→1/20 and generative
4/20→1/20; coverage>=5 generative remained3/20. Independent holdout completion
was stopped at the user's request. Repeated-seed uncertainty and geometric
recovery FDP must not be confused with a proven statistical guarantee.

## Limits and remaining work

- Full final-source performance and holdout acceptance remain unverified.
- Experimental regional FDP and maxT assumptions need further external review;
  exact masked likelihood is exact for its specified CAR/donor model, not proof
  the model captures all WGBS count covariance.
- The pre-existing adjacent-tile fixed neighbour correlation0.25 heuristic remains;
  this repair does not establish calibration for that separate tile-merging path.
- No claim that every possible issue in Epykit has been repaired.

## Use later

The branch/worktree is preserved at
`development/epykit_statistical_repair`. The external `../epykit3` checkout was
read-only in this session and has not been replaced. Use this branch directly or
transfer the supplied patch/Git bundle after reviewing it; no remote push was made.
All maintained implementation, tests, new harness entrypoints, failure fix,
verification helper, evidence summaries and this report are committed. Large raw
count stores and preserved run directories stay in the benchmark workspace, with
selected review artifacts in `gemini_bundle_v2`.
