# Review of the Argon conversation changes

Reviewed 2026-10-14: `argon_chat.txt`, Epykit commit `3302edf` in `../epykit3`, and the trial checkout at `development/epykit_minimal_gls_trial` commit `0bc083e`. The native GLS and DMC source files in the two current checkouts have identical SHA-256 hashes. No implementation files were changed during this review.

Assessment: the changes implement several useful ideas, especially replacing observed-group mean imputation with label-independent locus centering. However, partial-support inference remains incorrect, and several statistical claims in the supplied conversation are unsupported or contradicted by the implemented formulas. The current q-values should not be presented as validated FDR control.

## 1. High priority: missing observations still enter a complete-data GLS model

Locations: `../epykit3/src/epykit/_regional_gls.py:649`, `_regional_gls_native.py:232`, `_regional_gls_native.py:315`.

The wrapper centers on observed samples and zero-fills missing cells, but the native backend still uses the global sample contrast, its full-sample Gram, and full-sample residual degrees of freedom. It receives no observation mask. Zero-filling residuals is not equivalent to fitting only observed data.

Reproduction: positions `[0,20,40]`, unit depth weights, labels `[0,0,0,1,1,1]`, and centered observations `[-.3,-.1,0,.1,.3,0]` repeated at every locus. Samples 2 and 5 are absent. At independent fixed covariance, the six-row zero-filled fit returns effect **0.266667**, SE **0.066667**, and implied residual df **14**. Fitting the four genuinely observed rows returns effect **0.400000**, SE **0.070711**, and df **8**.

Use an explicit observation mask and a likelihood/design defined on actual observations. The conversation's proposed `yc != 0` missingness detector is also wrong: genuinely observed values can have exactly zero centered residual.

## 2. High priority: the support filter still depends on observed labels

Location: `../epykit3/src/epykit/_regional_gls.py:599`.

The chromosome universe is filtered once using minimum covered samples in each original group. All permutations subsequently reuse that universe without replaying the support filter. Thus label-independent centering does not make the whole procedure exchangeable.

For coverage mask `[1,1,0,1,1,0]`, labels `[0,0,0,1,1,1]` produce 2+2 observed samples and retain the locus. Permuted labels `[1,1,0,1,0,0]` produce 1+3, which would fail the same filter. With 5+5 samples and coverage only in one five-sample subset, an assignment can even put every observed sample in one group. `_calc_effects` then seeds no effect at that locus, despite its inclusion in the original universe.

Either use a label-independent eligible universe and apply valid group support per scan, or justify and implement a conditional randomization scheme. Do not claim exact exchangeability from the centering step alone.

## 3. High priority: anchor counting does not justify the new FDR guarantee

Locations: `../epykit3/src/epykit/_region_search.py:16`, `:36`; `dmr.py:1038`, `:1172`.

The replacement counts eligible left anchors rather than all possible intervals, then applies BH instead of BY. A single reported span per anchor does not make that span predetermined: its endpoint is selected using observed significance/effects. A fixed-span Stouffer p-value is not automatically a valid p-value for that adaptive anchor search.

Concrete check: ten positions spaced 20 bp apart, minimum five CpGs, minimum length 50 bp, gap limit 500 bp. There are **21 admissible intervals**, but the new function returns **6**. A single candidate with p=.005 gets q=.03. This arithmetic demonstrates the reduced correction; it does not establish its calibration.

The [Benjamini–Yekutieli theorem](https://www.math.tau.ac.il/~ybenja/MyPapers/benjamini_yekutieli_ANNSTAT2001.pdf) assumes the appropriate PRDS property for the tested statistics. Positive pairwise correlations alone do not establish it for this adaptive, two-sided region search, nor authorize reducing the family to anchors.

Restore a defensible search-family correction or calibrate the complete adaptive procedure. Keeping a more conservative family adjustment still requires valid regional p-values; it is not a cure for a misspecified spatial covariance.

## 4. High priority: the EB prior dispersion cap can make null tests anti-conservative

Location: `../epykit3/src/epykit/dmc.py:650`; shared by LR and GLM.

`phi0` is capped at 20, irrespective of observed dispersion, while the prior can receive 120 df. This forces substantial underestimation when true dispersion exceeds 20.

Reproduction: seed 491, 200,000 sites, df=4, residual sums `50 * chi_square(4)`. The helper estimates prior df **120**, prior scale **20**, and mean moderated dispersion **20.970**. For independent exact-null numerators `50 * chi_square(1)`, the corresponding F reference rejects **20.06%** at nominal 5%. This is a controlled statistical fixture, not a measured whole-pipeline error rate on the WGBS datasets.

Remove the arbitrary upper cap or introduce a justified, validated model that handles highly dispersed sites without forcing their scale downward.

## 5. Medium priority: heterogeneous-df EB estimation is not scale invariant

Location: `../epykit3/src/epykit/dmc.py:629`.

`z_corr / sd_sampling` divides by the sampling SD without first removing the common log-dispersion location. Multiplying every residual sum by a constant adds `log(c)/sd_sampling`, which differs across df strata and creates spurious prior heterogeneity.

Reproduction: seed 491, 200,000 residual sums, alternating df 2/8, all with one constant true dispersion. Scaling the same chi-square draws by 1, 2, 5, and 10 produces estimated prior df **120, 120, 9.626, and 4**, respectively. Biological heterogeneity is zero throughout. The same scaling with homogeneous df does not produce this effect.

More generally, a normal-IQR conversion and median substitution do not constitute an exact log-F fit. The [primary limma implementation](https://raw.githubusercontent.com/bioc/limma/master/R/fitFDist.R) uses centered log moments and expectation-based scale corrections; robust and unequal-df fits need appropriate distribution-specific treatment. Use parameter-recovery fixtures with known prior scale/df before enabling the replacement by default.

## 6. Medium priority: regional “variance moderation” only shrinks standard errors

Location: `../epykit3/src/epykit/_regional_gls_native.py:399`.

The batch adjustment multiplies SE by `sqrt((df + 6/k_scale)/(df + 6))/k_scale`. For K=5 it is exactly 1; for every K>5 it is below 1, regardless of residual variance. It supplies no positive prior variance or floor, so variance collapse remains variance collapse.

With six samples, SE multipliers at K=5/20/100 are **1.000/0.796/0.648**. A direct K=20 fit gave raw SE **0.0344965** and returned batch SE **0.0274495**. This is a length-dependent score boost, not empirical Bayes variance stabilization. Observed and default decoy fits share it, so this issue alone is not proof of permutation FDR inflation.

Implement an actual residual-variance posterior or explicitly describe and separately validate the score heuristic.

## 7. Medium priority: Kadane trimming does not solve the constrained problem

Location: `../epykit3/src/epykit/_regional_gls_native.py:44` and `:73`.

The implementation resets an unconstrained Kadane run, checks minimum length/count afterward, and falls back to the entire untrimmed interval. It can discard starts required for the optimal admissible segment.

Production discovery reproduction: positions `np.arange(7)*10`, effects `[.045,.015,.015,.015,.015,.045,.8]`, window 500, seed effect .1, minimum three CpGs, minimum span 1. The function returns `[0,7]`, penalized score **.565** at tau=.055. Valid segment `[4,7]` has higher score **.695**. With 15 bp spacing and the usual minimum span 50, `[2,7]` also beats the returned full interval: **.615** versus **.565**.

Use admissible prefix minima respecting both minimum CpG count and genomic span. Add brute-force comparisons on small arrays as an independent oracle.

## Other observations and validation limits

- `mean_meth_diff` at `_regional_gls.py:687` reports the support-shrunk discovery effect, not the ordinary mean treatment-control methylation difference. Keep the raw effect for reporting and the shrunk effect for discovery, or rename/document this change.
- Adjacent-tile combination now assumes fixed correlation .25 between immediate neighbours only (`dmr.py:1266`). This is a heuristic, not empirical calibration of arbitrary tile chains.
- An isolated untested CpG is reported as one contributing neighbour (`dmc.py:2936`) despite having NaN p-value; the correct finite-contributor count is zero.
- `reports/updated_vs_old_full_comparison.md` benchmarks trial commit **5a5e5c9**, not the current **0bc083e/3302edf** source. Those results cannot validate these latest changes.
- All **29 repository tests passed** in the repo environment. Dense GLS oracles cover complete data, and fixed-covariance checks bypass the batch SE adjustment. The discovery comparison now calls the same native implementation on both sides, so it is no longer independent.
- Ruff reported seven findings on the changed modules. Five are loop-closure warnings in a helper called immediately, so they do not by themselves demonstrate a runtime error; there is also an unused variable and import ordering issue.

Recommended priority: correct the observation model and support replay, remove the forced dispersion cap, and withdraw the anchor-BH FDR guarantee; then fix EB parameter recovery, constrained trimming, and score moderation. Validate the exact source revision with independent partial-support numerical oracles and repeated null runs before interpreting improved recovery as calibrated inference.
