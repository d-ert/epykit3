# Headline benchmark engine selection

## Tools in `eval_summary_post_phase3.parquet`

The Phase 4 headline benchmark (`benchmark/data/study1/eval_summary_post_phase3.parquet`) reports DMC performance for these epykit engines:

- `epykit_lr` — quasi-binomial likelihood-ratio (the recommended default at n ≥ 2)
- `epykit_lrplus` — `lr` with the power-stack on top (`fdr_method="fdr_tsbh"`, `neighbour_combine=True`, `sep_fallback=True`, `dispersion="eb"`)
- `epykit_welch_t` — Welch's t on raw β (variance-stabilising fallback)
- `epykit_fisher` — pooled Fisher exact (anti-conservative; documented as a comparator-only artefact)

Plus DMR methods (`epykit_dmr_tile`, `epykit_dmr_chain_merge`, `epykit_dmr_sliding_window`, `epykit_dmr_segment`) and the external comparators (`methylkit`, `methylkit_tuned`, `dss`, `biseq`, `methylsig`, `radmeth`, `fisher`).

## Tool *not* in the headline benchmark: `epykit_glm`

`epykit_glm` (IRLS binomial GLM with Wilkinson-formula contrast support) is **not** included in `eval_summary_post_phase3.parquet`. This is an intentional design decision, not an oversight.

### Why glm is excluded

1. **glm is designed for covariate-adjusted designs**, not for bare case-vs-control comparisons.
   Per `tl.dmc`'s own docstring and CLAUDE.md: glm is the path you take when you want `~ group + batch`, `~ group + sex + donor`, or any other formula with adjustment terms. It is reachable only via the `formula=` / `contrast=` arguments on `tl.dmc`. For a plain binary case/control split, `lr` is the closed-form, McCullagh-Nelder-dispersed quasi-binomial LR that handles the same problem with better small-n power.

2. **At n = 3v3, the Wald F reference distribution has only 4 residual degrees of freedom.** F(1, 4) at α = 0.05 has a critical value of ≈ 7.71, vs. the much-friendlier χ²(1) ≈ 3.84 that `lr` uses by default. At small n, this gap translates directly into low power: empirically, glm calls 0–1 sites at q < 0.05 even on the real treatment-vs-control split where `lr` calls ≥ 19 000 sites at the same q threshold on the same data (see `benchmark/data/null_calibration/gse263850/`).

3. **glm's per-site eligibility filter is stricter than `lr`'s.** `process_chromosomes_dmc(test="glm_contrast", ...)` requires the design matrix to be full-rank at every tested site, which excludes ~75 % of CpGs on GSE263850 (5.4 M tested vs. 22.0 M for `lr`). Even where glm fires, the comparison is over a different denominator.

4. **The Piao 3v3 simulator has no real covariates to adjust for.** Adding glm with `formula="~ group"` reduces to a slower, lower-power version of `lr` on identical inputs.

Given all four points, including glm in the headline benchmark would put a row in the table whose interpretation is "glm with bare `~ group` at n=3v3 has very low TPR" — which is statistically expected and uninformative about glm's actual purpose.

### Where glm *is* evaluated

glm is included in the **null calibration** sweep (`benchmark/data/null_calibration/summary.parquet`, `benchmark/data/null_calibration/gse263850/glm.parquet`). That table demonstrates the engine is well-calibrated under random label shuffling: median observed FDR = 0.0, with 1 site called on 1 of 10 shuffles (the shuffle that happens to recover the original 3v3 split). This is the legitimate methodological claim — "glm is correct" — without conflating it with the inappropriate-use-case claim "glm has competitive power at n=3v3 without covariates."

### Recommendations for paper Methods text

> "We evaluate epykit's `lr`, `lr+`, `welch_t`, and `fisher` DMC engines against five external comparators. `epykit_glm` is reserved for covariate-adjusted designs (Wilkinson formula syntax, e.g. `~ group + batch`) and is evaluated for calibration in the null-shuffle sweep (Table S-Calib) but not in the bare-`~ group` headline comparison, where it would offer no advantage over `lr` at small n."

### Recommended use of glm by epykit users

- **Use glm** when you have covariates that meaningfully affect methylation (batch effects, donor, sex, age, technical replicates within a contrast) and n ≥ 5 per group.
- **Use lr (or lr+)** for binary case/control without covariates at any n ≥ 2. lr's closed-form quasi-binomial LR with McCullagh-Nelder dispersion is the recommended default.
- **Use welch_t** when count-model assumptions are doubtful (very low coverage, severe overdispersion you don't want to model). Be aware it is conservative at small n.
- **Use fisher** only as a comparator artefact. It pools reads across replicates and ignores between-sample variance; documented as anti-conservative.
