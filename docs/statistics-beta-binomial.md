# Experimental beta-binomial count engine

Select the new engine explicitly:

```python
ep.tl.dmc(md, test="beta_binomial", dispersion="eb", reference="adaptive")
```

`test="auto"` retains the previous engine. Both command-line engine choice
lists include `beta_binomial`. The engine requires at least two covered
biological replicates in each group. It refuses smoothing, covariates,
formula/contrast inference, neighbour combination and permutation FDR.
Those features need separate statistical implementations.

The model uses each replicate's methylated and total integer counts:

```
M_gj ~ BetaBinomial(N_gj, mu_g * (1-rho)/rho, (1-mu_g) * (1-rho)/rho)
```

The group means are fitted to the actual count likelihood. Dispersion is
shared between the two groups at a site. A bounded, deterministic sample
of eligible replicate groups learns a factorized mean/dispersion mixture
by count likelihood; its log-dispersion moments define a continuous normal
shrinkage prior on log rho. Each site's log-rho MAP and conditional mean
fits use this prior. Prior estimation never reads benchmark truth or
published regions. The approximation for the binomial grid state, all
mixture weights, both optimization starts, convergence and algorithm
revision are stored in output-owned, hashed per-chromosome prior files.

The test compares separate means against one shared mean, holding the
fitted rho fixed in both fits. `adaptive` and `F` use the methylSig-inspired
F(1, n_valid_case+n_valid_control-2) small-sample approximation. It is not
an exact null distribution. `chi2` is available as an asymptotic diagnostic;
it can be liberal with few replicates and estimated dispersion. `bb_df`
records the F denominator degrees of freedom even in chi2 diagnostic mode,
whose reference degrees of freedom are always one.

Output includes fitted methylation means and their difference, `bb_rho`,
`bb_lr`, `bb_df` and `log2_odds_ratio_bb`. The odds ratio uses clipped means
only to give finite display values at exact 0/1 boundaries; fitting uses
the actual boundaries. Difference intervals invert the same conditional
likelihood-ratio test, holding fitted rho fixed. They do not integrate
dispersion or empirical-prior uncertainty. `beta_binomial_ci=False` skips
these intervals, records that omission and invalidates the corresponding
cache. Point estimates and p-values remain unchanged.

The engine is experimental. Its shared-dispersion assumption, global
mean/dispersion independence and learned prior need validation across
datasets. Convergence does not establish that a prior is well identified.
With few replicates, the F reference can have low site-level power after
genome-wide BH correction. A more explicit count likelihood alone does
not establish controlled CpG FDR, region FDR or appropriate inference for
paired donors. Region selection and spatial correlation remain separate
issues. The current comparison retains the accepted interval-family/BY
region correction unchanged.

References: [WGBS survey](https://doi.org/10.1093/bib/bbx013),
[DSS](https://pmc.ncbi.nlm.nih.gov/articles/PMC4005660/),
[RADMeth](https://pmc.ncbi.nlm.nih.gov/articles/PMC4230021/),
[methylSig](https://pmc.ncbi.nlm.nih.gov/articles/PMC4147891/).
