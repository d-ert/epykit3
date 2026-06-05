# Engine architecture map

This page is the canonical reference for epykit's DMC engine dispatcher
and the `lr+` "power stack" tunable. It exists so readers can verify
what each engine name does, where it lives in the source tree, and how
the orchestration layer in `tl.py` wires them together.

## DMC engines

`ep.tl.dmc(md, test=...)` selects one of four per-CpG statistical
engines, all implemented in [`src/epykit/dmc.py`](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmc.py):

| `test=` | Engine | Recommended use |
|---|---|---|
| `"lr"` (default at n ≥ 2) | Quasi-binomial likelihood-ratio with McCullagh–Nelder dispersion | **The recommended default.** What the benchmark paper characterises. |
| `"welch_t"` | Welch t-test on raw β | Quick first pass; ignores count uncertainty |
| `"fisher"` (default at n < 2) | Pooled Fisher exact | Single-replicate fallback |
| `"glm"` | IRLS binomial GLM with Wilkinson-formula covariates | Covariate-adjusted designs, batch correction |
| `"auto"` | Resolves to `"fisher"` at n < 2 and `"lr"` at n ≥ 2 | Convenience |

All four engines emit the same canonical schema (`chrom`, `pos`,
`n_case`, `n_control`, `mean_beta_*`, `meth_diff`, `meth_diff_ci_*`,
`pvalue`, `qvalue`, `log2_odds_ratio_pooled`) plus engine-specific
extras (`coef_treatment` for GLM, `f_stat`/`df1`/`df2` for multi-group
F-tests). The output goes through
`apply_multiple_testing_correction` ([dmc.py:2567](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmc.py))
to produce the q-values, streaming per-chromosome via the
`DMCStore` handle from [`_dmc_store.py`](https://github.com/d-ert/epykit3/blob/main/src/epykit/_dmc_store.py).

## The `power_stack` dispatcher

`power_stack` is a `tl.dmc` kwarg that bundles four research-grade
extensions into a single switch. The dispatcher lives in
[`src/epykit/tl.py:498–532`](https://github.com/d-ert/epykit3/blob/main/src/epykit/tl.py)
and behaves as follows:

| `power_stack=` | Behaviour |
|---|---|
| `"off"` / `False` (default) | Leaves all four knobs at user-passed values. **Bare `lr` runs at defaults.** |
| `"lr+"` / `True` / `"auto"` | Engages all four extensions at any n |
| `"conservative"` | Engages only when n ≤ 2 (pre-1.0 behavior) |

When engaged, the dispatcher flips four downstream knobs:

| Knob set by `lr+` | Implementation site | What it does |
|---|---|---|
| `neighbour_combine=True` (default window `neighbour_bp=500`) | [`dmc.py:2322–2470`](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmc.py) | Signed Stouffer p-value combination across CpGs within `neighbour_bp` on the same chromosome. Writes new `pvalue_combined` / `qvalue_combined` columns; the raw `pvalue`/`qvalue` are preserved. |
| `fdr_method="fdr_tsbh"` | [`dmc.py:2567–2649`](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmc.py) via `statsmodels.stats.multitest.multipletests` | Two-stage Benjamini–Hochberg with a Benjamini–Krieger–Yekutieli π₀ adaptive estimator. Replaces the standard BH q-values at the DMC level. DMR-level FDR continues to use plain BH ([`dmr.py:618`](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmr.py)) — this is intentional, since region-level FDR is a separate stage. |
| `sep_fallback=True` (threshold `sep_threshold=0.9`) | [`dmc.py:947–983`](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmc.py) | When a site shows `|Δβ| ≥ 0.9` *and* the LR p-value is > 0.05 *and* the count model is quasi-complete-separated, fall back to a pooled Fisher exact test. Takes `min(p_LR, p_Fisher)` — never inflates the p-value. |
| `dispersion="eb"` (default whether or not `lr+` is engaged) | [`dmc.py:775–820`](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmc.py) | Empirical-Bayes shrinkage of the per-site Pearson-residual dispersion toward a chromosome-pooled mean via a method-of-moments inverse-Gamma fit. |

## Honest framing of `lr+`

`lr+` is positioned as a **research tunable**, not a recommended
default. The four components were tuned on the Piao 2021 simulator,
which produces under-dispersed beta-binomial counts (φ ≈ 0.4) relative
to real WGBS (φ ≈ 1.5–5).

On the GSE263850 real-data cohort, at the same q-value threshold
(q=0.05), `power_stack="lr+"` produces **≈13× more DMC calls** than
bare `lr` — consistent with FPR drift under realistic
overdispersion, not with a genuine sensitivity gain. The benchmark
paper accordingly leads with bare `lr`'s numbers; `lr+` is exposed for
users who want to experiment but is not claimed as universally
superior.

If you want to try `lr+` on your own data, do so with a null
calibration run on shuffled-label samples to confirm that the FPR is
controlled at your expected level before trusting the results.

## DMR engines

`ep.tl.dmr(md, method=...)` selects from four DMR callers in
[`src/epykit/dmr.py`](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmr.py):

| `method=` | Engine | Notes |
|---|---|---|
| `"chain_merge"` (default) | DSS-compatible directional chain-merging | Presets: `strict`, `default`, `permissive` (`DMR_PRESETS`) |
| `"tile"` | Fixed-width tiles with read-pooling | Methylkit-compatible |
| `"sliding_window"` | Sliding window with signed Stouffer combining | |
| `"hmm"` | HMM segmentation (`dmr_hmm.py`, `_hmm.py`) | |

All four support optional permutation empirical FDR via
`ep.tl.dmr(..., empirical_fdr=True, n_perm=N)`.

## Where to look in the source tree

- `src/epykit/dmc.py` — all four per-CpG engines, the dispersion
  estimators, the multiple-testing correction, and the neighbour
  combine implementation. The single biggest file in the package.
- `src/epykit/_dmc_store.py` — the per-chromosome streaming store
  handle returned by `process_chromosomes_dmc(..., return_store=True)`.
  Preserves peak memory at O(largest chromosome).
- `src/epykit/_glm.py` — Wilkinson-formula design-matrix build, batched
  IRLS binomial GLM, Wald and F contrasts.
- `src/epykit/tl.py` — orchestrators (`tl.dmc`, `tl.dmr`, `tl.qc`)
  that wire engines into `MethylData.varm` / `MethylData.uns`. The
  `power_stack` dispatcher lives here at lines 498–532.
- `src/epykit/dmr.py`, `dmr_segment.py`, `dmr_hmm.py`, `_hmm.py` —
  DMR callers.
- `src/epykit/_smoothed_store.py` — Gaussian-kernel and BSmooth
  smoothing implementations.
