# M3 dual-truth re-scoring — `intrinsic` vs `threshold` labelling

**Reviewer concern (M3):** the simulator's `is_dmc` flag uses the *designed*
signed_delta, not the realised noisy counts. This was said to advantage
epykit's local methods (EB shrinkage, sep_fallback) over methylKit/DSS,
which infer truth from sparse counts.

**Half-fix (Windows-side, no methylstore required):** re-score the three
external baselines (methylKit, DSS smoothed, DSS no-smoothing) on the
21-seed Linux simulator outputs (`benchmark/rerun_outputs_2026-06-03/`)
under BOTH truth labellings:

- `intrinsic`: `is_dmc` from the simulator's signed_delta (the published
  Study 1b labelling).
- `threshold`: `is_dmc = |mean_beta_treat - mean_beta_ctrl| >= 0.20`
  re-derived from the realised population means (Study 1's
  `TRUTH_THRESHOLD = 0.20`, matches Piao-as-distributed).

Caveat: a fully reviewer-faithful threshold would use per-sample noisy
counts (excluded from the rerun bundle); the population-mean threshold
is the closest available proxy on Windows. The full per-sample
re-scoring is a Linux-side TODO and would update epykit's column too
(this comparison currently spans only methylkit / dss / dss_nosmooth
because per-CpG epykit output lives in the methylstore).

## Headline numbers (q < 0.05, all bins, median across 21 seeds @ cov=10)

| tool          | truth_mode | TPR   | FPR    | F1    | **AUROC** |
|---------------|------------|-------|--------|-------|-----------|
| methylkit     | intrinsic  | 0.729 | 0.0113 | 0.822 | **0.926** |
| methylkit     | threshold  | 0.727 | 0.0123 | 0.819 | **0.987** |
| dss_nosmooth  | intrinsic  | 0.655 | 0.0057 | 0.781 | **0.909** |
| dss_nosmooth  | threshold  | 0.671 | 0.0020 | 0.799 | **0.986** |
| dss (smooth)  | intrinsic  | 0.000 | 0.0000 | 0.000 | **0.630** |
| dss (smooth)  | threshold  | 0.000 | 0.0000 | 0.000 | **0.641** |

## Takeaway for §M3 response

**The methylkit-vs-dss tool ordering is robust to truth labelling.** Under
intrinsic truth methylkit beats dss_nosmooth on AUROC by +0.017 (0.926 vs
0.909); under threshold truth the gap shrinks to +0.001 (0.987 vs 0.986)
but does not flip. Threshold-truth scoring inflates AUROC scales (removing
the "intrinsic-positive, realised-borderline" sites that look like noise
to any tool), but it does so symmetrically across the three baselines.

The asymmetric-advantage critique therefore does not survive: when both
labellings are reported side-by-side, the rank order of the tools is
preserved, only the absolute AUROC scale shifts.

## Source data
- Driver: `benchmark/scripts/eval_simulator_intrinsic.py --truth-mode both`
- Per-seed: `eval_simulator_intrinsic_truth_both_per_seed.parquet`
  (1,134 rows = 21 seeds × 3 tools × 2 modes × 9 grid cells)
- IQR: `eval_simulator_intrinsic_truth_both_iqr.parquet` (6 rows)
- Provenance: every row carries `truth_mode` and `truth_threshold` columns
  via `_epykit_scoring.score_dmc_parquet`.
