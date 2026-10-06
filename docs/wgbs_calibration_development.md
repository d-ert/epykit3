# WGBS calibration development — 2026-10-06

This branch starts from the exact package source used for the v4 genome benchmark. Commit `0566519` preserves the pre-existing local modifications, including optional depth-aware variance engines. Original checkout: `/scratch/epykit3`. Development checkout: `/scratch/epykit-calibration`. Branch: `codex/wgbs-calibration-20261006`.

## Diagnosis from existing whole-genome results

`benchmark/scripts/diagnose_wgbs_v4.py` reads the frozen v4 result caches, truth and original per-sample counts. Truth is used only to evaluate results, never to fit the production package.

- Epykit: 317 / 323 site-level false positives on the signal dataset are more than 1 kb from any true region. DSS: 103,628 / 112,195 are more than 1 kb away. Boundary smoothing cannot explain most of either tool's observed site-level false discoveries on this dataset.
- Epykit's null CpG raw-p rejection rate at .05 increases with observed average coverage: 4.25% at 10–20×, 6.96% at 20–40×, 9.57% at 40–80×, 16.78% above 80×. Low-coverage sites are conservative. These are associations, not proof of a unique cause.
- 224 / 239 null BH-significant sites fall in the latent tau 0.2–0.5 stratum. This is an evaluation diagnosis; latent tau is never passed to the caller.
- 339 / 1,869 null regions have mean absolute effect below 0.10. The configured 0.10 minimum gates individual seed CpGs, not the region mean; this is the existing documented behavior, not a new discovery of a missing region-level filter.

A separate chr22 null comparison reproduces the original EB p-values exactly. At raw p ≤ .001, the default EB rejection rate is 0.172% versus 0.000285% for the previously available `bb_eb` engine. The latter is strongly conservative in the tail on this development input. It must not be promoted to a default without checking signal power and held-out calibration.

## First package change

The optional count-ratio region permutation estimator previously removed self/mirror assignments and conditionally dropped failed computations. Both can distort the population of label assignments used for inference.

- Both `max_t` and `region` now require all requested scans to complete; a failed scan raises instead of returning results conditional on success.
- Every valid assignment, including self/mirror and successful empty scans, is retained for region estimation.
- Both paths reject malformed, non-finite or out-of-range decoy scores and finite observed scores outside [0,1].
- Documentation correctly distinguishes the max-T compatibility `empirical_qvalue` (already scan-adjusted p, no second BH) from approximate count-ratio q estimates.

This patch hardens the already-existing permutation path. It does not silently replace the old asymptotic default, and it cannot alone explain improvements obtained by switching inference methods. In particular, count-ratio estimates from finite permutations, including zero estimates, still need validation; they do not establish universal FDR control. Complete-null max-T control requires exchangeable labels.

## Development evaluation

`benchmark/scripts/pilot_wgbs_region_calibration.py` copies only chr22 of each saved sample's count store into private pilot directories, recomputes CpG tests, and replays DMC → BH → chain merge → region filter for 99 random label assignments. The same scans are used for count-ratio and max-T comparisons. Each observed region is scored against chr22 truth; no benchmark inputs or baseline caches are modified.

The chromosome is a development subset of the already-inspected seed, not held-out validation. Its 50% overlap metrics should not be confused with whole-genome results. The future validation plan is to repeat with new simulator seeds, compare native threshold curves, retain 80% overlap as a sensitivity check, and examine real donor-aware inference separately.

## Run

```bash
cd /scratch/epykit-calibration
export PYTHONPATH=/scratch/wgbs_benchmark_v4/vendor310:/scratch/epykit-calibration/src
export MPLCONFIGDIR=/tmp/epykit-dev-mpl
export POLARS_MAX_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 benchmark/scripts/diagnose_wgbs_v4.py
python3 benchmark/scripts/check_depth_aware_wgbs.py
python3 benchmark/scripts/pilot_wgbs_region_calibration.py --out benchmark/results/wgbs_calibration/chr22_pilot_new
```

Bulk outputs under `benchmark/results/wgbs_calibration/` are ignored by git. Small reviewed summaries are retained under `benchmark/docs/wgbs_calibration/` when the run finishes. Regression and end-to-end replay tests are described with those summaries.
