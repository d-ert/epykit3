# A beta-binomial statistical engine for epykit

## User objective and evidence

Implement a new flagship count engine selected from the WGBS literature,
then run epykit on the original complete signal and null simulations and
GSE64177 real data. Compare against original epykit, the accepted region-only
correction, and every recorded competitor. Preserve input and harness hashes.
The recorded real baseline has only epykit and DSS, so other real competitors
cannot be claimed complete.

The supplied survey (`/scratch/bbx013.pdf`, Shafi et al., doi:10.1093/bib/bbx013,
pp. 742–743, 750) recommends beta-binomial methods for replicated counts and
emphasizes coverage, biological variation, sample size and spatial correlation.
It does not establish a universally best method. The selected family is a
hierarchical beta-binomial count model. RADmeth supports true per-sample
likelihood regression; DSS supports dispersion borrowing; methylSig warns that
estimated dispersion makes a chi-square LRT reference too liberal at small n.
Replacing the likelihood alone therefore does not establish calibration.

Alternatives: DSS-general's transformed GLS offers speed and complex designs,
but is an approximation to the count likelihood; MACAU's mixed model models
relatedness but adds sampling and computation unsuitable for this initial
two-group genome benchmark. dmrseq is a region-inference method and remains a
competitor rather than a replacement CpG engine.

Primary references:
- https://academic.oup.com/bib/article/19/5/737/3064341
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4005660/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4230021/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4147891/
- https://bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.html
- https://github.com/haowulab/DSS/blob/master/R/BSseq_util.R

## Model and numerical core

For site i, group g, replicate j, observed methylated count M and coverage N:

    P_igj ~ Beta(mu_ig * kappa_i, (1-mu_ig) * kappa_i)
    M_igj | P_igj ~ Binomial(N_igj, P_igj)
    kappa_i = (1-rho_i)/rho_i

Thus Var(M)=N*mu*(1-mu)*(1+(N-1)*rho). Sequencing counts and biological
replicates are distinct units. Use the actual per-replicate beta-binomial
log PMF, not a scaled binomial deviance. Missing coverage contributes no
likelihood; observed M=0 is real data. Reject invalid/fractional counts.

At fixed rho, fit each group mean and the shared null mean by bounded
concave likelihood optimization, including exact all-zero/all-methylated
boundary solutions. Test the difference using the likelihood ratio with
the same estimated dispersion in the full and null mean fits. Numerical
tests compare likelihoods and optima against SciPy independently. Support
rho=0 as the binomial limit for oracles and reference diagnostics.

## Dispersion and small-sample development protocol

Dispersion information must be borrowed from counts, without simulation truth
or the published DMR list. Evaluate a likelihood-based empirical prior and
a DSS-style log-dispersion shrinkage prior on independent count simulations.
The selected empirical-prior prototype learns a factorized discrete
mean-mixture times rho-mixture using exact beta-binomial likelihoods on a
deterministic subset of replicate groups. It avoids treating noisy per-site
variance estimates as latent biological dispersion. Mean–rho independence
is a working assumption to test. An unrestricted joint mixture is rejected:
at mean=0 or1 its rho split is unidentifiable. Factorization ties these rows
to the globally identified rho prior while retaining boundary observations.
Grid discretization and fitted weights must be recorded. The selected
implementation projects the count-learned rho distribution's log-scale
moments into a continuous normal shrinkage prior on log rho. The zero-rho
grid state represents half the smallest positive rho node; this approximation
is recorded. Conditional log-rho MAP uses this continuous density, avoiding
selection of an individual mixture-grid spike. Both deterministic mixture
starts are retained for sensitivity evaluation.

The two predefined p-value references are conditional Wilks chi-square(1)
and methylSig-style F(1,n_valid_case+n_valid_control-2). The latter is a
small-sample approximation, not an exact beta-binomial null distribution.
Start with the F reference and inspect both on independent development seeds
before freezing the candidate. No invented prior degrees of freedom are added.
Use separate untouched seeds for the subsequent evaluation; never choose
priors or reference using genome-baseline truth or paper overlap.

Report raw null tail rates, probability of any BH rejection under global
null, signal precision/recall, average precision, and strata by mean, rho,
depth and replicates. Low false calls with zero power does not demonstrate
a useful flagship engine. If neither reference has acceptable calibration,
continue development and explicitly retain experimental status; do not
declare controlled FDR based on one null run.

## Package integration

Expose a first-class public `test="beta_binomial"` engine through the engine
registry, Python API and both CLI choice lists. Preserve explicit legacy
engines and original numerical behavior. Initially implement the existing
two-group flow; reject covariate, formula and contrast requests that would
silently route to the old GLM. Do not claim paired/covariate inference.
Reject count smoothing and single-replicate inference for this engine.

Emit model-fitted group means, their difference, a clearly named fitted
log2 odds ratio, rho, LR statistic, F denominator degrees of freedom and
conditional likelihood intervals for the difference. These intervals hold
the fitted rho fixed; they do not integrate dispersion uncertainty.
Intervals must retain uncertainty
at boundary means. Approximate hybrid intervals must be labeled rather than
claimed to invert the difference LRT. Store algorithm revision and prior
provenance, and invalidate old caches when their meaning changes.

Keep the accepted region family/BY correction in `_region_search.py` and
`dmr.py` fixed at commit f73e66d. No permutation path or region threshold
change belongs in this engine comparison. A CpG engine cannot by itself
establish region-level FDR or resolve within-region spatial dependence.

## Computation and verification

Python >=3.10; use existing NumPy, SciPy, Polars and Numba dependencies.
Load one chromosome at a time, fit in bounded site chunks, and record
runtime and RSS. The source snapshot used by a run must be immutable.
First verify synthetic scalar and batch oracles, missingness, label and
sample-order invariance, high coverage, boundaries and cache behavior;
then independent count calibration, end-to-end store/CLI integration,
the package suite, independent review, and the three complete genome runs.
For timed runs, `beta_binomial_ci=False` omits confidence intervals while
preserving identical fitted effects, p-values and q-values.

No run is complete until its status, expected outputs and source/input
hashes prove completion. New results must be compared with all recorded
tools using the existing scoring definitions, including CpG and region
metrics and 50%/80% overlap. Real paper agreement is descriptive, not truth.

## Authorization and execution

The user explicitly delegated method selection, implementation and benchmark
execution. Proceed with reversible work in the existing isolated worktree,
with visible decisions and inspectable artifacts, without a second consent
cycle. Preserve every previous result and archived development candidate.

## Selected candidate after independent development probes

Select true BB conditional mean/null likelihoods, a continuous normal prior
on log rho, and the methylSig-inspired F(1,n_valid-2) approximation. Learn the
latent dispersion distribution by factorized count-likelihood EM, then use its
log-scale moments as a smooth prior instead of an unstable discrete MAP spike.
Map the unresolved rho=0 grid state to half the first positive node; record
this approximation, the prior fit convergence and both-start sensitivity.
Inference bounds are rho in [1e-6,.95] and log-prior sd at least .05 for
numerical stability. These choices precede untouched hold-out evaluation.

Rejected prototypes remain in benchmark development artifacts: plugin chi2
has inflated null tails, integrated-dispersion chi2 is still too liberal in
small high-depth samples and expensive, and the experimental Bayesian
model-comparison alternative adds assumptions without solving low-depth power.
No prototype claim of exact FDR carries into this candidate.

Expose optional beta_binomial_ci=True for conditional difference-profile
intervals. The timed benchmark sets it False because primary benchmark
outputs do not use intervals; record this omission and NaN CI columns.
Hold region settings/BY correction fixed, and compare site and region power
separately. Preserve auto's legacy behavior until benchmark evidence supports
default promotion; the new engine is explicitly invoked for every new run.
