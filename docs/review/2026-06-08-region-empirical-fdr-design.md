# Region-level empirical FDR for tile DMRs — design note

**Status:** IMPLEMENTED for `method="tile"` and `method="chain_merge"`
(2026-06-08). `sliding_window` / `segment` still raise `NotImplementedError`; the
CLI is tile-only. chain_merge reuses the shared count-ratio core via
`empirical_fdr_for_chain_merge` (recomputes the per-CpG DMC per shuffle —
expensive; simple two-group tests only). The region/max_t dispatch is now the
shared `_aggregate_region_perm_results`, used by both callers.

**Changes landed**
- `dmr.py`: `_region_count_ratio_fdr()` (shared, caller-agnostic core),
  `_is_self_or_mirror_perm()`, and `empirical_fdr_for_dmr(..., fdr_method=)`
  dispatch — default `"region"` (count-ratio), `"max_t"` preserves the old FWER
  min-P. Adds the `empirical_fdr_set` column; warns at <4/group.
- `tl.dmr(..., fdr_method="region")` passthrough; set-level FDR surfaced in
  `md.uns["dmr_params"]["empirical_fdr_set"]`.
- `dmc.py` cross-reference docstring corrected.
- Tests: `test_region_count_ratio_fdr.py` (helper units),
  `test_dmr_region_fdr_mode.py` (dispatch/integration); the P0 max-T tests
  re-pointed at `fdr_method="max_t"`.

**Original status:** design for review (no library code changed yet)
**Context:** `empirical_fdr_for_dmr` returns `empirical_qvalue ≡ 1.0` on demo4
(GSE263850, 3 vs 3), so the 18,925 tile calls show "0 significant." This note
specifies the statistically-correct replacement and the numbers harness used to
validate it.

## 1. Why the shipped estimator collapses to 1.0

`empirical_fdr_for_dmr` is named/documented as an **FDR** but computes a
**Westfall–Young max-T (min-P)** statistic ([dmr.py:1688](../../src/epykit/dmr.py)):
for each permutation it keeps only the **genome-wide minimum** null p-value, then
asks how often that single best decoy beats each observed tile.

On a genome-wide tile scan with dispersion-inflated per-tile p-values (real WGBS
φ ≈ 1.5–5), every shuffle produces a p ≈ 1e-100 tile *somewhere*. Demo4 evidence:
even the strongest real tile (p = 1.3e-104) is beaten by **47/50** shuffles, so
18,337/18,925 tiles get emp_p = 1.0. This is a family-wise (FWER) bar mislabeled
and BH-mangled as FDR; it is correct math answering the wrong question.

The `empirical_fdr_for_dmc` twin ([dmc.py:3019](../../src/epykit/dmc.py)) does the
*same* min-P but is **honestly documented** as FWER — and it explicitly tells
users to "prefer `empirical_fdr_for_dmr` … for a less conservative region-level
empirical FDR." That promise is currently false.

## 2. The two wrong options, and the correct one

| | construction | failure mode |
|---|---|---|
| (A) max-T / min-P (shipped) | per-perm genome-wide min p, count perms beating obs | ultra-conservative → wall of q=1.0 |
| (B) pooled-p + BH (docstring-literal) | per-region empirical p vs pooled null, then BH | anti-conservative — region count is inflated by dispersion; plain BH ignores it (this is what the P0 test rightly distrusts) |
| **(C) count-ratio target-decoy** | **E[# null survivors] / # observed survivors** | **correct** — dispersion cancels in the ratio |

(C) is the BSmooth (Hansen 2012) / SAM (Tusher 2001) construction and is the
standard for region/peak-level empirical FDR.

## 3. The estimator (C)

`call_dmr_tile_based` already BH-corrects and filters tiles on `qvalue < alpha`
and `|meth_diff| ≥ min_abs_meth_diff` before returning — for **both** the observed
run and every permutation. So observed and decoy survivor sets are defined by an
identical pipeline (apples-to-apples).

Let, over `M = n_perm` label shuffles:

- `R` = number of observed survivors (the called tiles).
- `V = (1/M) · Σ_perm |null survivors in perm|` = mean decoy survivors per shuffle.

**Set-level FDR (headline):**

```
FDR_set = min(V / R, 1)
```

Interpretation: "of the R tiles you called, ≈ FDR_set are explained by label
noise." `V ≥ R` ⇒ FDR = 1 ⇒ calls are indistinguishable from shuffled noise.

**Per-tile q-values (ranking gradient):** using raw tile p-values as the
discovery threshold `t`,

```
V(t) = (1/M) · #{ pooled null survivors with p ≤ t }
R(t) =          #{ observed survivors with p ≤ t }
fdr(t) = V(t) / R(t)
q_j    = min_{ t ≥ p_j } fdr(t)          # suffix-min → monotone, then clip to [0,1]
```

`q` at the loosest threshold equals `FDR_set`, so the gradient and the headline
are consistent. Significance = `q < alpha`.

## 4. Why (C) is correct and foolproof

1. **Dispersion cancels.** Decoys inherit the same overdispersion, so they inflate
   `V` exactly as the observed run inflates `R`. The ratio measures only the
   *excess* of real over noise — the precise quantity the φ-drift breaks for the
   asymptotic test.
2. **Random region count is native.** Counts are the estimand; no fixed-`m`
   hypothesis set is assumed (which is what makes plain BH on regions invalid).
3. **Directly estimates E[V]/R = FDR.** No distributional assumption.
4. **Monotone, bounded, well-defined** at the edges (see §5).

## 5. Edge cases

- **All perms in the mean.** `V` divides by `M` (all perms), *including*
  zero-survivor perms — a shuffle that yields nothing is genuine evidence of low
  noise and must count as 0, not be excluded. (This differs from the max-T
  `n_perm_used` floor, which is about a min-p denominator, not a count mean.)
- **`R = 0`** → return the empty-frame contract unchanged.
- **`V > R`** → `FDR_set` capped at 1.0 (honest "all noise").
- **Monotonicity** enforced by suffix-min so `q` is non-decreasing in `p`.
- **Small-n caveat (load-bearing for this dataset).** At 3 vs 3 there are only
  C(6,3)=20 (≈10 effective) relabelings; shuffles adjacent to the truth retain
  real signal and **inflate `V`**, so (C) is *conservative* here (over-estimates
  noise). The estimator is correct; the permutation null is intrinsically weak at
  this n. Emit a `UserWarning` when `min(n_treat, n_ctrl) < 4`. For robust
  inference at small n prefer the model-based `chain_merge` (DSS) path.

## 6. API / column semantics

- `empirical_qvalue` ← per-tile count-ratio q (the deliverable; thresholded by
  downstream `tl.dmr` / report / exports — column name unchanged).
- `empirical_pvalue` ← per-tile `V(t_j)/M`-normalised pooled tail fraction
  `(#null ≤ p_j)/N_null` (a diagnostic; the literal "fraction of null DMRs ≤
  observed" the docstring always described).
- `observed_dmr.attrs["empirical_fdr_set"]` ← the single headline `FDR_set`.
- New kwarg `fdr_method: {"region", "max_t"} = "region"`. Default = (C). `"max_t"`
  preserves the shipped FWER behaviour (honestly relabelled), so nothing the audit
  added is lost — and the existing P0 test is re-pointed at `fdr_method="max_t"`.

## 7. Validation harness (no library edits)

`demo_output4/capture_null_pools.py` re-runs demo4's exact empirical step
(n_perm=50, seed=42) tee-ing the full per-perm null survivor pools the library
discards. `demo_output4/analyze_region_fdr.py` then contrasts (A)/(B)/(C) on those
real pools.

### Results (GSE263850, 3 vs 3, n_perm=50 → 47 valid after excluding true-split + mirror)

Observed survivors R = 18,925. Null survivors/shuffle: min 331, median 678,
**mean 2,509**, max 10,267 (nine near-true 3v3 shuffles leak signal → 10,267).

- **Set-level FDR = mean(V)/R = 2,509/18,925 = 0.133** (≈13% of calls are noise).
  Excluding the 9 signal-leaking shuffles → **0.036** (lower bound). Truth ≈ 4–13%.

| # of 18,925 calls passing | q<0.05 | q<0.10 | q<0.25 | q<0.50 |
|---|---|---|---|---|
| (A) max-T (shipped) | **0** | 0 | 0 | 0 |
| (B) pooled-p+BH (docstring-literal) | 5 | 5 | 49 | 1,198 |
| (C) count-ratio (proposed) | **383** | 3,141 | 18,925 | 18,925 |

Strongest tile (p=1.3e-104): (C) q = 0.000 vs (A) max-T q = 1.000.

**Reading:** (A) is a useless wall of 1.0. (B) — the *literal docstring* — is also
broken but the other way: its raw per-region emp_p is anti-conservative, then BH
over the inflated region count m over-corrects to near-nothing (5). Only **(C)**
gives a coherent gradient: 383 high-confidence DMRs at q<0.05, an honest ~13%
set-level FDR, and the whole set bounded below 0.25. This is the deliverable.
The 9 leaky shuffles are the 3v3 small-n caveat made quantitative.
