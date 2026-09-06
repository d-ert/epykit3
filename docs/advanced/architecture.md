# Engine architecture map

This page is the canonical reference for epykit's DMC engine dispatcher
and the `lr+` "power stack" tunable. It exists so readers can verify
what each engine name does, where it lives in the source tree, and how
the orchestration layer in `tl.py` wires them together.

This page is the canonical architecture reference. `README.md` and
`CLAUDE.md` summarise it and link here; when they disagree, this page wins.

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
`apply_multiple_testing_correction` in [`dmc.py`](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmc.py)
to produce the q-values, streaming per-chromosome via the
`DMCStore` handle from [`_dmc_store.py`](https://github.com/d-ert/epykit3/blob/main/src/epykit/_dmc_store.py).

## The `power_stack` dispatcher

`power_stack` is a `tl.dmc` kwarg that bundles four research-grade
extensions into a single switch. The dispatcher lives in
the `power_stack` block of `tl.dmc` in [`src/epykit/tl.py`](https://github.com/d-ert/epykit3/blob/main/src/epykit/tl.py)
and behaves as follows:

| `power_stack=` | Behaviour |
|---|---|
| `"off"` / `False` (default) | Leaves all four knobs at user-passed values. **Bare `lr` runs at defaults.** |
| `"lr+"` / `True` / `"auto"` | Engages all four extensions at any n |
| `"conservative"` | Engages only when n ≤ 2 (pre-1.0 behavior) |

When engaged, the dispatcher flips four downstream knobs:

| Knob set by `lr+` | Implementation site | What it does |
|---|---|---|
| `neighbour_combine=True` (default window `neighbour_bp=500`) | `combine_neighbour_pvalues` in [`dmc.py`](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmc.py) | Signed Stouffer p-value combination across CpGs within `neighbour_bp` on the same chromosome. Writes new `pvalue_combined` / `qvalue_combined` columns; the raw `pvalue`/`qvalue` are preserved. |
| `fdr_method="fdr_tsbh"` | `apply_multiple_testing_correction` in [`dmc.py`](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmc.py) via `statsmodels.stats.multitest.multipletests` | Two-stage Benjamini–Hochberg with a Benjamini–Krieger–Yekutieli π₀ adaptive estimator. Replaces the standard BH q-values at the DMC level. DMR-level FDR continues to use plain BH (the DMR callers in [`dmr.py`](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmr.py)) — this is intentional, since region-level FDR is a separate stage. |
| `sep_fallback=True` (threshold `sep_threshold=0.9`) | the `sep_fallback` branch of `_score_finalize` in [`dmc.py`](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmc.py) | When a site shows `|Δβ| ≥ 0.9` *and* the LR p-value is > 0.05 *and* the count model is quasi-complete-separated, fall back to a pooled Fisher exact test. Takes `min(p_LR, p_Fisher)` — never inflates the p-value. |
| `dispersion="eb"` (default whether or not `lr+` is engaged) | the dispersion block of `_score_finalize` in [`dmc.py`](https://github.com/d-ert/epykit3/blob/main/src/epykit/dmc.py) | Empirical-Bayes shrinkage of the per-site Pearson-residual dispersion toward a chromosome-pooled mean via a method-of-moments inverse-Gamma fit. |

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
| `"segment"` | Rule-based 3-state segmentation over `meth_diff` (`dmr_segment.py`; `dmr_hmm.py` is a deprecated alias) | Not a fitted HMM; the HMM primitives in `_hmm.py` are separate |

Permutation empirical FDR (`ep.tl.dmr(..., empirical_fdr=True, n_perm=N)`)
is wired for `method="tile"` and `method="chain_merge"`. Both harnesses
share `_aggregate_region_perm_results` in `dmr.py`, which dispatches on
`fdr_method`: `"max_t"` (default) is the Westfall-Young min-P statistic
with a BH transform; `"region"` is the opt-in count-ratio target-decoy FDR
in `_region_count_ratio_fdr`. The tile harness re-runs
`call_dmr_tile_based` on shuffled labels. `empirical_fdr_for_chain_merge`
replays the observed DMC (engine knobs from `md.uns["dmc"]`, the observed
chromosome universe and multiple-testing method) into a private temporary
`DMCStore` per permutation, then chain-merges and filters it like the
observed run. The per-CpG `empirical_fdr_for_dmc` keeps min-P only; region
mode exists in the DMR API alone. `sliding_window` and `segment` raise
`NotImplementedError` until each gets its own label-shuffle scheme. See
[the design note](../review/2026-06-08-region-empirical-fdr-design.md).

## Canonical chromosome filter

`src/epykit/_chroms.py` is the one definition of a main-assembly
chromosome: autosomes `1` to `22`, `X`, `Y`, and the mitochondrion as `M`
or `MT`, with an optional case-insensitive `chr` prefix. It is a fixed
human-style list, not a species-aware assembly validator. Every
`canonical_only` option is opt-in (default `False`) and drops the same
contigs through `filter_canonical_logged`, which emits one INFO line per
selection naming what it dropped.

| Surface | Scope |
|---|---|
| `read_bismark`, `read_methyldackel`, `read_combined_strand_bed`, `convert_sample` | Drops non-canonical contigs before the partition write. The setting is part of the per-sample conversion manifest; a changed setting rebuilds the sample and replaces its partition directory. |
| `tl.dmc`, `process_chromosomes_dmc` | Filters the auto-detected partition list before the engine and the multiple-testing correction, on the binary and the formula / contrast path. `plan_run` resolves the list once through `resolve_dmc_chromosomes` and carries it on the `DMCPlan`, so the engine run and every `empirical_fdr` permutation test the same universe. The resolved list is part of the low-level cache signature; `canonical_only` is part of the `resumable=True` fingerprint and recorded in `md.uns["dmc"]`. An explicit `chromosomes=` list, including an empty one, is used verbatim. |
| `tl.dmr(method="tile")`, `call_dmr_tile_based` | Filters the auto-detected partition list before the tile test and the BH correction. `tl.dmr` resolves the list once and shares it with every `empirical_fdr` permutation. An explicit `chromosomes=` list, including an empty one, is used verbatim. |
| `chain_merge`, `sliding_window`, `segment` | Not supported. These callers inherit the chromosome universe of the upstream DMC run and raise `ValueError` on `canonical_only=True`; run `tl.dmc(canonical_only=True)`, filter at ingestion, or restrict `tl.dmc` with `chromosomes=`. |
| `epykit convert`, `epykit dmc`, `epykit dmr --method tile` | `--canonical-only` forwards to `convert_sample`, to `process_chromosomes_dmc` (binary) or `tl.dmc` (formula / contrast), and to the tile caller with the list resolved once for the observed run and the permutations. The other `dmr` methods exit with an error that points at `epykit dmc --canonical-only`. |
| `pl.manhattan` | Takes its axis order from `CANONICAL_CHROMS_UCSC` and hides other contigs unless `canonical_only=False` (unchanged plot behaviour). |

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
  that wire engines into `MethylData.varm` / `MethylData.uns`. The `power_stack` dispatcher lives in `tl.dmc`.
- `src/epykit/dmr.py`, `dmr_segment.py`, `_hmm.py` — DMR callers.
  `dmr_hmm.py` is the deprecated import shim for `dmr_segment.py`
  (see [Deprecations](../reference/deprecations.md)).
- `src/epykit/_smoothed_store.py` — Gaussian-kernel and BSmooth
  smoothing implementations.
- `src/epykit/_chroms.py` — the canonical chromosome predicate, the
  order-preserving filters and the UCSC order shared by ingestion, the
  tile DMR caller and the Manhattan plot.
