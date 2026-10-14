# Confirmed Epykit fixes: design specification

**Status:** Draft for user review. Documentation only; implementation has not started.
**Date:** 2026-10-14.
**Target repository:** `development/epykit_statistical_repair`.
**Baseline:** `fe2892e27b83d470a19f576b32a144bcc46ebf7a`.
**Evidence:** [Verified audit](/Users/dert/Desktop/Coding/Epykit/wgbs_benchmark_v6/reports/audit_claim_verification_fe2892e.md).
**Plan:** [Implementation plan](../plans/2026-10-14-confirmed-epykit-fixes.md).

## 1. Intent and scope

The user requested a specification and plan for confirmed improvements worth
implementing, followed by a stop before implementation. The intended update fixes
candidate seed loss and GLM batch failures, makes convergence bookkeeping agree
with the iterates returned, and removes avoidable annotation and masked-GLS work.
Success means independently checked correctness, preserved model assumptions,
inspectable benchmark evidence and measured performance gains where claimed.

Three approaches were considered:

1. **Recommended: isolated corrections plus low-risk optimizations.** Deliver the
   six changes below separately, measure each, and profile the remaining costs.
   This gives attributable behavior changes and straightforward rollback.
2. Statistical redesign plus optimization: defer. Pseudocounts, t-to-z scores,
   per-cell depth covariance and regional EB need separate statistical validation.
3. Broad GLS performance rewrite: defer. Brent searches, structured masked solves
   and cross-permutation batching have not yet demonstrated accuracy and benefit.

The selected approach has independently reviewable work packages:

| ID | Change | Expected behavior |
|---|---|---|
| B1 | Enforce a seed in every final regional candidate | Remove seedless split/trim products; choose the best admissible seeded interval |
| B2 | Separate boundary extension from split/trim scoring | Use one penalty determined by `min_effect`, independently of extension threshold |
| G1 | Isolate singular GLM systems | Invalid sites cannot force every valid site into scalar linear algebra |
| G2 | Iterate only active GLM sites | Freeze converged coefficients; calculate final diagnostics from those coefficients |
| A1 | Vectorize nearest-TSS lookup by chromosome | Preserve annotation semantics and row order |
| M1 | Remove masked donor-Gram vector-copy repetition | Preserve covariance, likelihood and numerical rejection behavior |

Unused regional helper cleanup is a final maintenance task, contingent on a fresh
call-site search. Profiling produces evidence and follow-up recommendations only.

## 2. Global constraints

- Python compatibility remains >=3.10; retain existing dependency bounds and add no dependencies.
- Public function signatures, output column names, row ordering and dtypes remain compatible.
- Preserve the observed-cell CAR/donor covariance model and per-label support replay.
- Preserve experimental FDP/maxT semantics; make no new finite-sample FDR guarantee.
- Preserve baseline results and input stores; write every benchmark variant to a new directory.
- Run compiled Numba checks and warm timing runs; JIT-disabled results do not establish compiled correctness or speed.
- Use no simulation truth or real-paper overlap to choose algorithm parameters.
- No implementation, test creation, benchmark execution, dependency installation or source commits are authorized by this documentation task.

## 3. Regional candidate contract (B1/B2)

For every emitted half-open interval `[lo, hi)` with direction `d` in {-1, 1}:

- `hi-lo >= min_cpgs` and `pos[hi-1]-pos[lo]+1 >= min_length`.
- It lies within the existing gap-respecting envelope and split fragment.
- It contains an index `j` with `d * seed_signal[j] >= min_effect`.
- `seed_signal` is the actual smoothing signal used by that discovery scale,
  including the existing support-weighted boxcar signal when that mode is used.
- The seed must survive final trimming. A seed in a parent or discarded halo
  does not qualify. Raw mean effect is not an alternative seed rule.

B1 keeps current penalties and edge-floor rules, isolating the seed correction.
For each split fragment it selects the maximum penalized raw-effect sum subject
to count, span and seed constraints. For ties retain the first interval encountered
under increasing exclusive right endpoint and increasing admissible start; equal
prefix sums retain the earlier start. If no admissible seeded interval exists,
emit none. Do not merely reject an unconstrained optimum when another seeded
admissible interval exists. An admissible prefix minimum can enforce a start
bound at the most recent direction-matching seed, retaining linear scan cost.

B2 then sets `tau = 0.55 * min_effect` for valley classification, constrained
trimming and multiscale candidate competition. It adds no public tuning parameter.
Envelope thresholds remain `boundary_min_effect` when supplied, otherwise
`min_effect`; the existing edge floor remains `boundary_min_effect` when supplied,
otherwise `0.35 * min_effect`. Thus B2 changes the boundary-enabled penalty but
does not silently redesign edge trimming. Equality at the seed threshold qualifies.

Apply both changes identically to observed and decoy scans, both directions and
all discovery scales. Update regional metadata revisions and rely on the existing
source digest to invalidate inference caches. Validate caching explicitly.

The regression fixture is `pos=arange(13)*60`, raw effects `.07` repeated five
times, `.01` three times and `.20` five times; supplied smooth values are `.06`
eight times and `.15` five times. At `min_effect=.1`, boundary `.05`, minimum
five CpGs and span 50, `[0,5)` must be absent and `[8,13)` retained.

## 4. GLM validity and iteration contract (G1/G2)

G1 classifies structural estimability from the observed rows of the design.
For each distinct coverage mask, require at least p covered rows and full column
rank of `X[covered]`. Use an SVD rank threshold
`eps * max(n_covered, p) * largest_singular_value`; do not use determinant signs.
Reuse each coverage-pattern result across sites. Empty sites and rank-deficient
patterns are invalid before initialization and IRLS.

Finite, structurally valid rows use stacked linear solves. Split a failing batch
recursively to isolate numerical failures instead of looping over all sites.
Bound the top-level batches to 8,192 systems. At a failing singleton return NaNs.
Do not add a ridge or pseudoinverse to obtain ordinary unpenalized inference for
an unidentified site. Final covariance uses the same isolation strategy with an
identity right-hand side, preserving the full covariance needed by joint tests.

Invalid or unsolved sites return NaN coefficients, SEs, covariance, deviance and
Pearson statistics, while `n_eff` remains the observed coverage count. Summarize
failures once per fit rather than per failed site. Rank-deficient sites that
previously acquired finite ridge coefficients may now return NaNs intentionally.

G2 computes each IRLS update only for valid, unconverged sites. Keep existing
`max_iter`, coefficient tolerance, clipping, working weights and separation rules.
Accept the new coefficient vector before testing convergence, then freeze it.
Nonfinite updates invalidate that site; all-invalid and empty inputs terminate
without a reduction on an empty active array. Nonconverged sites remain NaN in
inference outputs. Compute diagnostics once from final accepted coefficients;
preserve the optional sixth covariance output.

Add `GLM_ALGORITHM_REVISION = "binomial-irls-estimable-active-v1"` and include it
in binary GLM resume fingerprints and GLM result metadata. The contrast path has
no resume support at baseline; retain that behavior and record the revision in
its metadata. No cache-schema redesign is required.

## 5. Annotation contract (A1)

Group midpoint indices by chromosome and search all positions for that chromosome
in one NumPy call. Preserve the scalar implementation's truncation to integer
midpoints, original site order, left-neighbour preference on distance ties,
`side="left"` behavior at exact/duplicate TSS positions, gene string selection
and strand-signed distance. Missing or empty TSS entries retain the empty-gene and
`int32.min` distance sentinel; valid distances retain the current clipping rules.
Only the existing `multi_annotation` path changes. Annotation feature overlap
assembly and coordinate conventions remain outside this work package.

## 6. Masked preparation contract (M1)

Form a contiguous CAR-applied donor-residual matrix once per sample and length
evaluation. Compute donor response and donor Gram using matrix/vector products,
without copying a donor column for each pair. Preserve the residual-quadratic
calculation that avoids cancellation, symmetrization, eigenvalue tolerances,
exact-fit rejection and observed residual df. Do not change the length or donor
search grid/tolerance, covariance structure or complete/masked dispatch rule.

This is an allocation optimization, not the proposed banded/Woodbury solver.
Keep it only if compiled benchmarks show a measurable gain without numerical or
thread-oversubscription regressions. No O(K) masked complexity claim is made.

## 7. Validation and acceptance

### Correctness

- Compare seeded trimming with an independent exhaustive small-array oracle;
  include a seeded feasible optimum that differs from the seedless optimum.
- Check complement-label symmetry, scale deduplication, partial-support replay,
  gap/span constraints and cold/warm cache behavior.
- Compare valid GLM sites with independent statsmodels binomial fits, plus a
  scalar per-site IRLS oracle using the specified stopping rules. Test single
  bad sites, all bad sites, rank deficiency despite `n_eff>=p`, separation,
  nonconvergence, empty batches and joint-contrast covariance.
- Compare TSS output exactly against a frozen scalar oracle on mixed chromosome
  order, ties, duplicates, negative-strand genes and missing/empty chromosomes.
- Compare masked objective/effect/SE with the existing dense observed-data oracle
  at length 0/10/120/10000 and donor fraction 0/.3/.98, plus near-exact residuals.

For masked fixed covariance use `rtol=2e-8, atol=2e-8`; complete-path results must
continue to meet their existing tighter tests. For well-conditioned GLM fixtures
require `rtol=1e-6, atol=1e-8` against statsmodels and the scalar stopping oracle.
Use separate expectations for separation and unavailable inference. Inspect any
profiled optimizer-basin or threshold-crossing difference rather than loosening
tolerances to conceal it.

### Benchmarks

Freeze the baseline and each variant with source/input/evaluation hashes. Use
separate B1 and B1+B2 variants before combining with output-preserving changes.
Record call counts, candidate counts, seed violations, per-scenario recovery,
FDP definitions, null-positive dataset counts, wall time and peak RSS. Enumerate
the benchmark tasks exactly in the implementation plan.

For microbenchmarks use one untimed warm-up and five measured repetitions;
report median and range with CPU, OS, versions, worker/BLAS counts and dataset
dimensions. Keep an optimization only with at least 10% median improvement on
its targeted large-input fixture and no >5% median regression or >10% peak-RSS
increase on the paired representative workload. If noise crosses these limits,
repeat once with ten measurements and retain the original results too. These
are engineering gates, not promised outcomes.

Correctness fixes may be retained without a speed gain. Any end-to-end regression
above those limits needs explicit review before integration. Changes to boundary
results must be attributed to B1/B2; optimization-only commits require identical
regional selected coordinates/directions on the comparison runs. GLM early
stopping can change last digits, so inspect threshold-crossing DMC calls.

Do not promote B2 by default if the same pilot suite gains null-positive datasets
or mean recall drops more than two percentage points in any predeclared main
empirical/generative coverage stratum. Report this as a regression and defer B2;
do not tune `.55` using the results. This pilot guard does not establish FDR.
Real-paper concordance is descriptive and never an acceptance target.

## 8. Deferred work and decision record

Profile ingestion/joins, discovery, complete and masked fitting, and Python/Numba
dispatch separately after the core fixes. Record the actual workloads: real 3+3
uses nine representative fits; fully enumerated pilot 5+5 uses 125. Recommend
further work only when a stage consumes >=20% of warmed end-to-end wall time or
is the measured peak-memory source. Profiling does not authorize those changes.

Deferred: pseudocount/FDR redesign, exact-t claims, GPD tails, regional EB,
per-cell covariance weights, Brent/coarser searches, structured masked solves,
cross-permutation batching, shared smoothing kernels, and NumPy union/scatter.
Each needs its own concrete design, equivalence/statistical validation and review.

## 9. Delivery and stop condition

This request delivers only this draft specification and its implementation plan.
Neither artifact constitutes approval to implement. Stop after document
self-review and provide both links. Implementation and all benchmark commands
remain future work until the user requests execution.
