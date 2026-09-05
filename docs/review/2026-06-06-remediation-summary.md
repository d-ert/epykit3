# epykit remediation — autonomous session summary

**Branch:** `review-remediation` (off `main` @ `7e51dc4`). Nothing pushed; `main` untouched.
**Date:** 2026-06-06 (overnight autonomous run).
**Inputs:** [peer-review report](2026-06-06-epykit-peer-review.md) · [design spec](../history/superpowers/specs/2026-06-06-epykit-review-fixes-design.md) · [Tier 1 plan](../history/superpowers/plans/2026-06-06-epykit-review-tier1.md).

## TL;DR

Tier 1 (paper-blockers + silent-wrong-science criticals) is **complete and fully tested**. A high-value chunk of Tier 2 (correctness sweep) is **done and tested**. The remaining Tier 2 items and all of Tier 3 are **deferred with documented reasons** below. Every fix landed as an atomic commit with a test that would have caught the bug. `main` is untouched; review the branch and merge when satisfied.

**How to review:** `git log --oneline main..review-remediation`, then `git diff main..review-remediation`. Each commit is self-contained. To run: `uv run pytest -m "not slow"` (matrix tier) and `uv run pytest -m slow` (accuracy/calibration tier, now also a CI job).

## A correction worth reading first

Running the tests (rather than trusting the static review) overturned one finding and pruned two others:

- **DVC at n=2 was ANTI-conservative, not the "silent zero" the review predicted.** The review (certain) said n=2 Brown-Forsythe yields all-NaN → 0 DMRs (conservative). In fact the within-group SS is a *floating-point residual* of 0, which exploded the F-statistic so a 2-vs-2 null cohort called **~75% of sites** differentially variable. The fix (require ≥3/group) is the same direction, but the severity was the opposite — worth knowing if this caller was used.
- **M-SEC1 (clocks zero-fill) had no real code fix.** The review's proposed change (exclude missing CpGs vs zero-fill them) is mathematically identical (`0·coef` == omitting the term). The honest remediation is a warning + accurate docs, which is what landed.
- **M-STAT2 (multi-group df reporting) is more complex than a 1-line df swap.** The proposed fix (emit `df_phi`) does not make `f.sf(f_stat, df1, df2)` reproduce the stored p, because `wald_test` uses an *adaptive* F/χ² reference. Attempted, test-falsified, reverted. Deferred (see below).

## Tier 1 — complete (commit-by-commit)

| Finding | Fix | Commit |
|---|---|---|
| C3 Horvath transform | Corrected negative branch `21·exp(x)−1`; real-function test | `8671b21` |
| M-PKG1 polars floor | `polars>=1.0,<2` (pivot `on=`); numba `>=0.60`; pivot smoke test | `ebeaa13` |
| M-PKG4 lock + threads | Committed `uv.lock`; documented `uv sync --frozen` + thread pinning | `2eaf89a` |
| M-PKG2 CLI/API parity | CLI `dmc` defaults `dispersion="eb"`; new flags; subprocess parity test | `46fa31a` |
| C2 annotation chrom guard | Raise on zero chrom-name overlap, warn <50%; test | `63cefad` |
| M-PKG3 calibration + CI | Real-engine null-FPR test (lr/glm/welch_t pass); slow CI job | `c7fc734` |
| C1 Bismark coordinate | 1-based auto-detect via start/end; `coordinate_base` override; manifest v2; fixtures rewritten; end-to-end coord→annotation test | `b85323b` |

**Calibration result (validates the paper's core claim):** under the null, `lr` FPR@0.05 = 0.047 (KS not rejected), `welch_t` 0.038, `glm` 0.038 — all calibrated-to-conservative, none anti-conservative.

## Tier 2 — done (commit-by-commit)

| Finding | Fix | Commit |
|---|---|---|
| M-DVC1/2/3 | n<3 anti-conservative fix (≥3/group), `bartlett`→`brown_forsythe` rename+warn, `min_coverage` floor | `83bbc22` |
| M-DMR1/2 | Sliding-window combine uses raw p not q; correlation docstring corrected | `ee2e3d1` |
| M-ANN1/2 | Intron via running-max end; deterministic gene tie-break | `abc1e11` |
| M-SEC1/3 | Clock biased-age warning; `impute_knn_beta(return_mask=True)` | `b6c15a9` |
| C5 power | Non-central t + overdispersion + MT-aware α; statsmodels oracle | `317c036` |
| M-DMR4 / M-STAT5 / M-QC2 | Honest docstrings: DSS-*style*, min-P FWER (not FDR), intermediate-β fraction (not contamination) | `b297620` |

## Deferred — Tier 2 remainder (with reasons)

These are real findings left undone; each notes why and a pointer for the fix.

**Needs careful statistical redesign (not a safe autonomous 1-liner):**
- **M-STAT2** multi-group `df2` reporting — `wald_test` uses an adaptive F/χ² reference; to let a reviewer reconstruct p, emit the reference choice per site or force a single reference. (Attempted+reverted; test in history.)
- **M-STAT4** BH/Storey over the finite p-subset only (currently includes NaN-masked sites as p=1) — conservative-direction, but touches the core genome-wide correction; validate against the accuracy suite before changing.
- **M-DMR3** stratified permutation FDR slices treatment across whole strata — needs `(stratum, label)` threaded so labels permute *within* strata.
- **M-DMR5** region-level FDR is computed on post-filter survivors (selective inference); tile path double-BHs — correct over all candidate regions, drop the second BH.
- **M-DMR6** `dmr_segment` uses unsigned Stouffer + per-chromosome BH — sign by `meth_diff`, BH genome-wide.

**Engine internals / comparability (need a benchmark re-run to validate):**
- **M-STAT1** glm effect-size CI via the unused `delta_method_meth_diff_ci` + φ-scale `coef_se`.
- **M-STAT3** apply the `DF_PHI_FLOOR` consistently across `lr` and `glm`.

**Lower-risk but unfinished (ran out of session runway):**
- **M-DMR7** standardise DMR `start`/`end` width convention across callers; **M-DMR8** make `_MAX_DMR_BP` a parameter; **M-DMR9** NaN-safe PMD smoother.
- **M-ANN3** vectorise nearest-TSS (perf); **M-ANN4** per-transcript TSS (HOMER fidelity); **M-ANN5** treat `XM_` as protein-coding.
- **M-QC1** sex_check guards + warn when `diptest` absent; **M-QC3** pairwise-complete `sample_correlation`.
- **M-PKG5** move research-grade callers (clocks/asm/entropy/pmd/hmr) behind `epykit.experimental.*`; stop exporting `_glm.build_design`.

**Cannot be validated in this environment (Windows / no pysam / no reference data):**
- **C4 ASM bisulfite confound** — the highest-severity deferred item. The fix (restrict phasing anchors to bisulfite-safe SNV classes, or gate behind `experimental`) is well-specified in the review, but `asm.py` is pysam-gated and its tests `importorskip("pysam")` on Windows, so it cannot be exercised here. Do this on Linux with the `bam` extra installed.
- **M-SEC2** clock gold-standard-mean imputation — needs the clocks' distributed reference means.
- **M-SEC4** entropy Miller-Madow bias correction — `entropy.py` is pysam-gated (Linux-only validation).

**Tier 3 (minors/nitpicks)** — all deferred: lint ratchet (`select=["F"]`→`F,E,W,I,B`), `python -m epykit` entry, CHANGELOG ordering, save/load symmetry, samplesheet covariate dtype, `region_beta` half-open + multi-part glob, FutureWarning noise, docs-URL fix, etc. See the review's Minor/Nitpick tables.

## Test status

- Tier-1 gate: full `not slow` suite green (300 passed, 5 platform-skips); slow calibration 3/3.
- Each Tier-2 batch verified with targeted runs before commit; final full `not slow` gate run is the last validation (see the session's final report for the count).
- New tests added: coordinate convention, CLI/API parity, null calibration, DVC n2/alias, sliding-window raw-p, intron/tie-break, impute mask, power oracle.

## Recommended next steps

1. Review + merge `review-remediation` (or cherry-pick Tier 1 first — it's the paper-blocking set).
2. On a Linux box with the `bam` extra, do **C4 (ASM)** — it's the highest-severity deferred item and silently fabricates ASM on WGBS.
3. Schedule the deferred statistical redesigns (M-STAT2/4, M-DMR3/5/6) as a focused follow-up; they need the accuracy suite as a guardrail.
4. The benchmark paper itself needs no change — C1 doesn't affect its numbers, and the engine calibration claim is now backed by an executing test.
