# Fresh whole-branch statistical review

Candidate reviewed: `0bc083e..11163ca`. Verdict: **request changes**. Critical findings: 0. Important findings: 3.

Scope: supplied diff package read first, design, implementation plan, progress ledger, changed inference/discovery code and regression tests; additionally the workspace benchmark and null-validation runners because provenance is an explicit deliverable. Existing 76-test and lint evidence was accepted, without redundantly rerunning the unchanged suite. Novel targeted probes were run in the existing trial Python environment with this candidate's `src` on `PYTHONPATH`. No source files were edited.

## Important 1 — Exact fits survive the degeneracy guard and violate complementary-label equivalence

Locations: `src/epykit/_regional_gls_masked.py:149-172`, especially the strict `rss <= 0` check at line 170; complete-path analogues are `_regional_gls_native.py:300` and `:331`.

An exactly fitted treatment pattern leaves floating-point residual energy around machine precision. Because the code checks only whether the residual sum is literally nonpositive, it treats numerical residue as an estimated biological variance, produces enormous finite statistics, and can give materially different absolute statistics under complementary labels. This violates the plan's explicit exact-fit handling requirement and the identity used to avoid fitting mirrored scans.

Reproduction (real count-level public fitter, no simulated likelihood approximation):

```python
import numpy as np
from epykit._regional_gls import fit_region
p = np.arange(5) * 20.0
g = np.r_[np.zeros(3), np.ones(3)]
n = np.full((5, 6), 30.0)
m = np.tile([10, 10, 10, 20, 20, 20], (5, 1)).astype(float)
m[0, 0] = n[0, 0] = 0
print(fit_region(p, m, n, g))
print(fit_region(p, m, n, 1-g))
```

Observed original/complement results:

| Quantity | Original | Complement |
|---|---:|---:|
| Absolute statistic | 4.216583485528111e16 | 6.206648115358962e16 |
| Standard error | 1.6063363059588924e-17 | 1.0912896967929308e-17 |
| CAR length | 1.4552072822031477 | 0.0 |

The true response lies exactly in the locus-intercept plus treatment design. Fully observed versions also return artificial finite SEs (for five loci, approximately 1.80e-9), rather than an explicit degeneracy outcome. Additional probes at 3, 8, and 20 loci show the same defect; masked mirrored scores differ by large factors at 8 and 20 loci.

Required repair: detect numerical zero residuals using a justified response/design-relative tolerance, and apply one explicit symmetric policy to complete and masked fits and to all label scans. Rejecting the entire inference with a clear degeneracy error is consistent with the existing failure policy. Add exact-fit regressions for complete/masked observations and complements; retain nondegenerate low-variance fits.

## Important 2 — Pooled decoy covariance is allowed to produce/select purported maxT permutation p-values

Locations: `src/epykit/_regional_gls.py:570-571`, `:597-604`, pooled branch at `:812-838`, unconditional maximum-score attachment at `:850`, and unconditional `maxT_scope` metadata at `:905` (line numbers refer to candidate 11163ca).

`selection_method="max_t", pooled_decoy_covariance=True` is accepted. Observed candidates are fitted by per-region REML, whereas decoys are evaluated with covariance parameters estimated from observed candidates. Adding self/mirror scores does not make these different statistics exchangeable. Nevertheless the output is selected using `maxT_pvalue` and metadata describes joint whole-genome all-null error under exchangeable samples. Even the all-null interpretation needs identical label-wise statistic construction.

Novel integration probe: six temporary sample partitions, 15 loci at 20-bp intervals, depth 40, independent binomial draws with probabilities .25/.75 and RNG seed 8978; called `run_regional_gls(..., n_perm=20, pooled_decoy_covariance=True, selection_method="max_t", q_cutoff=.5, cache=False, workers=1)`. The call completed with one candidate, `selection_column="maxT_pvalue"`, and `maxT_scope="joint whole-genome all-null error under exchangeable samples; strong region control not established"`. A warning about pooled covariance is emitted, but the unsupported result is still constructed and selectable.

Required repair: reject this incompatible selector/covariance combination and omit or explicitly invalidate the maxT diagnostic whenever pooled covariance is requested. Alternatively perform the same covariance procedure for every assignment. Add regression coverage of the unsupported combination. Default per-region replay is unaffected by this finding.

## Important 3 — Resume provenance does not enforce the evaluation code or pilot truth used for metrics

Workspace locations: `scripts/run_statistical_repair_benchmarks.py:101-113`, `:138-143`, `:280-288`; `scripts/run_statistical_repair_validation.py:17-18` and its call to `snapshot_source`.

The source snapshot correctly hashes/freezes Epykit. However, on resume it compares only those package hashes and immediately returns the previous manifest. The recorded `runner_sha256` and `evaluation_sha256` are never checked. Scoring runs against live `scripts/gls_sensitivity.py`, live `wgbs_sim/evaluate.py`, and live pilot `truth.cpg.tsv.gz`/`truth.regions.tsv`; the helper and truth identities are absent from per-run identity. The real-data comparison helpers are similarly live. Also, when the null validator imports and calls `snapshot_source`, `runner_sha256=sha256(__file__)` identifies the benchmark helper module, not the validation entry point.

Novel isolated resume probe: created a temporary source snapshot from the clean candidate, replaced only its recorded runner/evaluator hashes with 64 zeroes, then called `snapshot_source` again. It returned successfully and retained both impossible hashes (`RESUME_ACCEPTED_MISMATCHED_RUNNER_AND_EVALUATOR True True`). This demonstrates the missing enforcement independently of any assumption about currently running processes.

Consequence: resuming a partial benchmark after a scoring/helper/truth change can retain old completed metrics and generate new metrics under different evaluation logic, while the common manifest still names the old evaluation hash. That fails the required reproducible, unmixed benchmark deliverable. This is not evidence that existing outputs have actually been mixed.

Required repair: record and verify the actual entry-point script, all evaluation/helper dependencies, and relevant truth/reference inputs in a frozen evaluation identity; reject changed identities on resume or rerun evaluation into a new explicitly versioned result set. Carry that identity into each run. Tests should change an evaluator/helper/truth identity while leaving the inference package unchanged and ensure old completion cannot be silently reused. Candidate tables can be reused for explicitly versioned rescoring; inference need not be rerun merely to recompute metrics.

## Checks without additional material findings

- The masked formulation eliminates the actual observed nuisance design, accounts for observed residual degrees of freedom, and retains confounded zero-eigenvalue donor contributions to coefficient variance. The dense observed-submatrix tests are independent and substantially improve confidence; no algebraic mismatch was found for nondegenerate inputs.
- Depth-relative Pearson exposure and pooled-contrast leverage match their intended beta-binomial moment quantities. The new tests exercise constant-rho depth changes, while metadata/docs retain the approximate small-sample reference limitation. Those controlled quadratic-form tests alone do not establish count-level FDR.
- Constrained trimming's increasing admissible-prefix set enforces both count and span; the exhaustive test is independent. Complete-interval counts and BY correctly include all geometry-defined endpoints.
- Union support is label independent, group eligibility is replayed per scan, and raw effects are reported. Missing first-sample chromosome partitions are retained. Observation masks are explicit rather than inferred from response zeroes.
- Multiscale candidate competition is deterministic and removes overlapping candidates; independent opposite-direction behavior tests exist. No concrete material regression was found in that procedure.
- Expected false counts use complete scan exposure, including empty scans. FDP estimates and pooled score tails are explicitly distinguished from finite-sample randomization p-values. Under the default symmetric fit procedure, maxT includes sampled self/mirror outcomes and respects the exact balanced 3+3 minimum of .1.
- Regional cache keys include changed modules, discovery method, calibration, selector, and diagnostic settings. DMC and chain-merge cache revisions changed. Python 3.10-compatible language features and existing dependencies are retained.

## Remaining acceptance gates

Tasks 4/5 are still underway. The final candidate's full 18-pilot/four-configuration and real whole-genome/four-configuration results, reserved independent null validation, final result verification, and verified Gemini bundle were not available for acceptance in this review. Their absence is an outstanding delivery gate, not an extra source-code finding. Do not infer final FDR calibration or improved recovery from this code review, earlier development null runs, or the 76 passing tests. The original plan checkboxes also remain unchecked; reconcile them with the final ledger at delivery.

After the three Important findings are resolved, run focused regressions and the required final verification, update exact source/evaluation provenance, and request review of the resulting changes before treating the branch as complete.
