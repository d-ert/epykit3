# Depth-aware LR calibration repair

The benchmark's validated development configuration is:

```python
ep.tl.dmc(md, test="lr", dispersion="bb_site", reference="F",
          smoothing=False, sep_fallback=False, fdr_method="fdr_bh")
ep.tl.dmr(md, method="chain_merge", alpha=.05,
          min_abs_meth_diff=.1, min_cpgs=5, minlen_bp=50,
          dis_merge_bp=500, pct_sig=.5, use_q_for_sig=False,
          min_mean_qvalue=None, empirical_fdr=True,
          fdr_method="max_t", n_perm=100, perm_n_jobs=1)
# Threshold empirical_qvalue <= .05; in max_t mode this is a compatibility
# alias for a scan-adjusted permutation p-value, not an additional BH q-value.
```

`bb_site` uses the squared-depth sums to account for the beta-binomial variance
of pooled proportions. Biological variation does not disappear as read depth
increases. The finite-sample F reference is retained even when a variance
estimate reaches its lower bound; the artificial degrees-of-freedom floor of
50 is removed. Effect intervals use separate depth corrections for each group.
The moment estimator and F reference are approximations, not exact finite-sample
beta-binomial likelihood inference.

The `eb` default is retained for compatibility; it is not the selected corrected
benchmark configuration and still has demonstrated calibration problems under
variable-depth simulation. `bb_eb` is an experimental moment-shrinkage ablation
and performed worse in signal recovery; it is not recommended by these results.
Smoothing, neighbour p-value combination, Fisher fallback and covariate GLM
inference have not been qualified by the new tests. Use of `bb_site` with smoothed
fractional counts is outside the validated model.

The region permutation fix includes successful empty scans in the denominator,
aborts failed scans, includes self/mirror draws in uniform Monte Carlo sampling,
and applies the plus-one correction. Max-T already adjusts for the whole scan;
it is not followed by BH again. It supports complete-null family-wise control
under exchangeability. Strong control under partial alternatives is not proved.
Permutation inference repeats the full pipeline and costs more than asymptotic
region ranking. In 5-versus-5 comparisons there are 126 distinct two-sided splits;
100 random draws can repeat assignments.

Independent compact validation: 2 of 500 complete-null datasets had CpG
rejections at BH .05 (0.4%; exact 95% interval 0.048%–1.437%). These simulations
are independent conditional on the fixed heterogeneous BLUEPRINT library. On
the already-inspected 5.9-million-CpG signal dataset, the corrected test made
zero CpG discoveries. Improved calibration therefore does not establish
competitive CpG power. See `/scratch/wgbs_benchmark_v3/results/` for source hashes,
full-scale results and comparison artifacts. No release or publication claim
of five-point noninferiority is established.
