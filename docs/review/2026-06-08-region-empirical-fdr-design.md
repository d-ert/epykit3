# Region-level empirical FDR for DMRs: design note

**Status:** implemented for `method="tile"` and `method="chain_merge"` in the
1.2 development cycle, as an opt-in. `fdr_method="max_t"` stays the default
of `empirical_fdr_for_dmr`, `empirical_fdr_for_chain_merge` and `tl.dmr`, so
existing numbers do not move; `fdr_method="region"` selects the count-ratio
construction described here. `sliding_window` and `segment` still raise
`NotImplementedError`, and the CLI `--empirical-fdr` stays tile-only. The
per-CpG `empirical_fdr_for_dmc` keeps its min-P construction; region mode
exists only in the DMR API.

This note was written on 2026-06-08 on the `feat/canonical-chrom-filter`
branch, where the count-ratio estimate became the default. The port to main
(September 2026) kept the statistics, made the mode opt-in, and corrected the
chain_merge harness so that every permutation replays the observed
multiple-testing correction and streams through a private `DMCStore`. The
validation scripts and the `demo_output4` run directory named below were
never committed; the numbers are reported as recorded on the branch and were
not reproduced during the port.

**Changes landed**

- `dmr.py`: `_region_count_ratio_fdr()` (shared, caller-agnostic core),
  `_is_self_or_mirror_perm()`, `_aggregate_region_perm_results()` (the
  `max_t` / `region` dispatch both harnesses use), and the `fdr_method`
  keyword on `empirical_fdr_for_dmr`. Adds the constant `empirical_fdr_set`
  column; warns at fewer than four samples per group in region mode.
- `dmr.py`: `empirical_fdr_for_chain_merge()` and its per-permutation engine
  `_chain_merge_perm_survivors()`.
- `tl.dmr(..., fdr_method=)` passthrough for both callers; the set-level FDR
  is surfaced in `md.uns["dmr_params"]["empirical_fdr_set"]`.
- Tests: `tests/test_region_count_ratio_fdr.py` (helper units),
  `tests/test_dmr_region_fdr_mode.py` (tile dispatch),
  `tests/test_chain_merge_empirical_fdr.py` (chain_merge harness and
  replay), `tests/test_empirical_fdr_method_coverage.py` (API and CLI
  gates).

**Context.** `empirical_fdr_for_dmr` returned `empirical_qvalue == 1.0` for
every tile on GSE263850 (3 vs 3), so the 18,925 tile calls showed "0
significant". This note specifies the alternative estimate and the numbers
harness used to compare the two.

## 1. Why the shipped estimator collapses to 1.0

`empirical_fdr_for_dmr` was named and documented as an FDR but computes a
Westfall-Young max-T (min-P) statistic: for each permutation it keeps only
the genome-wide minimum null p-value, then asks how often that single best
decoy beats each observed tile.

On a genome-wide tile scan with dispersion-inflated per-tile p-values (real
WGBS phi is about 1.5 to 5), every shuffle produces a p of about 1e-100
somewhere. Evidence from the GSE263850 run: even the strongest real tile
(p = 1.3e-104) is beaten by 47 of 50 shuffles, so 18,337 of 18,925 tiles get
`empirical_pvalue = 1.0`. This is a family-wise (FWER) bar labelled and
BH-adjusted as an FDR: correct mathematics answering a different question.

The `empirical_fdr_for_dmc` twin in
[`dmc.py`](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmc.py)
does the same min-P computation and documents it as FWER. Its docstring
points users at the DMR function for a less conservative region-level
estimate; that is what region mode now provides.

## 2. The two wrong options, and the correct one

| | construction | failure mode |
|---|---|---|
| (A) max-T / min-P (`fdr_method="max_t"`) | per-permutation genome-wide min p, count permutations beating the observed value | ultra-conservative: a wall of q = 1.0 |
| (B) pooled-p + BH (the old docstring, read literally) | per-region empirical p against the pooled null, then BH | anti-conservative: the region count is inflated by dispersion and plain BH ignores it |
| **(C) count-ratio target-decoy** (`fdr_method="region"`) | **E[# null survivors] / # observed survivors** | dispersion cancels in the ratio |

(C) is the BSmooth (Hansen 2012) / SAM (Tusher 2001) construction and is the
standard for region and peak-level empirical FDR.

## 3. The estimator (C)

`call_dmr_tile_based` already BH-corrects and filters tiles on
`qvalue < alpha` and `|meth_diff| >= min_abs_meth_diff` before returning, for
the observed run and for every permutation. The chain_merge harness applies
the observed DMC correction, the same chain_merge knobs and the same
`min_mean_qvalue` filter to every permutation. Observed and decoy survivor
sets are therefore defined by an identical pipeline.

Let, over `M` usable label shuffles:

- `R` = number of observed survivors (the called regions).
- `V = (1/M) * sum over permutations of |null survivors|` = mean decoy
  survivors per shuffle.

**Set-level FDR:**

```
FDR_set = min(V / R, 1)
```

Interpretation: of the `R` regions called, about `FDR_set * R` are explained
by label noise. `V >= R` gives `FDR_set = 1`, that is, the calls are
indistinguishable from shuffled noise.

**Per-region q-values:** using the observed region p-values as the discovery
threshold `t`,

```
V(t) = (1/M) * #{ pooled null survivors with p <= t }
R(t) =         #{ observed survivors with p <= t }
fdr(t) = V(t) / R(t)
q_j    = min over t >= p_j of fdr(t)      # suffix-min: monotone; then clip to [0, 1]
```

`q` at the loosest threshold equals `FDR_set`, so the gradient and the
headline agree. Significance is `q < alpha`.

## 4. Why (C) is the right estimate

1. **Dispersion cancels.** Decoys inherit the same overdispersion, so they
   inflate `V` as the observed run inflates `R`. The ratio measures only the
   excess of real over noise, which is the quantity the dispersion drift
   breaks for the asymptotic test.
2. **A random region count is native.** Counts are the estimand; no fixed
   hypothesis set `m` is assumed, which is what makes plain BH on regions
   invalid.
3. **It estimates E[V] / R = FDR directly.** No distributional assumption.
4. **Monotone, bounded, well-defined** at the edges (section 5).

## 5. Edge cases

- **Usable permutations.** `V` divides by the number of usable
  permutations, including zero-survivor permutations: a shuffle that yields
  nothing is evidence of low noise and counts as 0. Assignments equal to the
  observed split or its mirror swap are excluded (their statistics are the
  observed ones, not null draws), and so are permutations whose engine
  failed. This differs from the `max_t` floor, which is about a min-p
  denominator, not a count mean, and which keeps self and mirror draws.
- **No usable permutation.** All three columns are NaN and a `UserWarning`
  is emitted. The estimate is undefined, not zero.
- **`R = 0`.** The empty-frame contract is unchanged: the three columns are
  added as nulls.
- **`V > R`.** `FDR_set` is capped at 1.0.
- **Non-finite observed p-values** stay NaN in both columns.
- **Monotonicity** is enforced by the suffix minimum, so `q` is
  non-decreasing in `p`.
- **Small-n caveat.** At 3 vs 3 there are only C(6, 3) = 20 (about 10
  effective) relabellings; shuffles adjacent to the truth retain real
  signal and inflate `V`, so (C) is conservative there (it over-estimates
  noise). The estimator is correct; the permutation null is weak at this
  n. A `UserWarning` fires when `min(n_treat, n_ctrl) < 4`.

## 6. API and column semantics

- `empirical_qvalue`: the per-region count-ratio q (the deliverable;
  thresholded by downstream `tl.dmr`, the report and the exports; the
  column name is unchanged).
- `empirical_pvalue`: the pooled-null tail fraction
  `#{null <= p_j} / N_null` per region. A diagnostic, not a calibrated
  individual-region p-value.
- `empirical_fdr_set`: a constant column with `FDR_set`, mirrored in
  `md.uns["dmr_params"]["empirical_fdr_set"]` (`None` when NaN).
- `fdr_method: {"max_t", "region"} = "max_t"` on `empirical_fdr_for_dmr`,
  `empirical_fdr_for_chain_merge` and `tl.dmr`. The default preserves the
  shipped FWER numbers; the tests in `tests/test_dmr_empirical_fdr.py` pin
  them.

## 7. chain_merge replay

`empirical_fdr_for_chain_merge` has to reproduce the observed analysis for
every shuffle, not only the chain_merge step. For each permutation it
streams `process_chromosomes_dmc(return_store=True)` into a private
temporary directory with the engine knobs `tl.dmc` recorded in
`md.uns["dmc"]` (`test_used`, `unite`, both minimum sample counts,
`dispersion`, `reference`, `smoothing` and its span, `sep_fallback`,
`sep_threshold`), over the chromosome universe of the observed `DMCStore`
(or the materialized DMC table), applies the observed DMC multiple-testing
method, chain-merges from the temporary store with the observed knobs and
applies the observed `min_mean_qvalue` filter. The temporary store is
removed afterwards, so the observed store and its manifest are never
touched, also with `perm_n_jobs > 1`.

`tl.dmr` rejects, before the first permutation: GLM, formula and contrast
DMCs; `use_smoothed=True` (its pseudo-count store is temporary); a DMC record
without the engine knobs; an explicit `chromosomes=` that differs from the
observed universe (rerun `ep.tl.dmc` with the restriction instead); and an
`empirical_strata` column that is missing or does not cover every sample.
Raw `pvalue` / `qvalue` drive the chain_merge gate in every permutation;
neighbour-combined columns are never substituted.

## 8. Validation harness (branch record, not reproduced)

`demo_output4/capture_null_pools.py` re-ran the GSE263850 empirical step
(n_perm = 50, seed = 42) and captured the full per-permutation null survivor
pools the library discards. `demo_output4/analyze_region_fdr.py` then
contrasted (A), (B) and (C) on those pools. Neither script nor the run
directory is part of the repository.

### Results (GSE263850, 3 vs 3, n_perm = 50, 47 usable after excluding the true split and its mirror)

Observed survivors R = 18,925. Null survivors per shuffle: min 331, median
678, **mean 2,509**, max 10,267 (nine near-true 3 vs 3 shuffles leak signal
and produce the 10,267).

- **Set-level FDR = mean(V) / R = 2,509 / 18,925 = 0.133** (about 13% of
  the calls are noise). Excluding the nine signal-leaking shuffles gives
  **0.036** (a lower bound). The truth is between about 4% and 13%.

| calls of 18,925 passing | q < 0.05 | q < 0.10 | q < 0.25 | q < 0.50 |
|---|---|---|---|---|
| (A) max-T (`max_t`) | **0** | 0 | 0 | 0 |
| (B) pooled-p + BH | 5 | 5 | 49 | 1,198 |
| (C) count-ratio (`region`) | **383** | 3,141 | 18,925 | 18,925 |

Strongest tile (p = 1.3e-104): (C) q = 0.000; (A) q = 1.000.

**Reading:** (A) is a wall of 1.0. (B) is broken the other way: its raw
per-region empirical p is anti-conservative, and BH over the inflated region
count then over-corrects to almost nothing (5). Only (C) gives a coherent
gradient: 383 high-confidence DMRs at q < 0.05, a set-level FDR of about
13%, and the whole set bounded below 0.25. The nine leaky shuffles are the
3 vs 3 small-n caveat made quantitative.

This comparison is the evidence for offering region mode. It is not a
calibration study: the engine hash gate under `benchmark/` covers selected
per-CpG `lr` output and says nothing about either permutation construction.
Users should confirm the FDR level on their own data, for example with a
null run on fully shuffled labels.
