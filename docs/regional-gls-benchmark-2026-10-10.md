# Experimental regional GLS: speed and power evidence

The new engine substantially improves recovery in the tested BLUEPRINT-based
independent-residual simulations. It does **not** establish that Epykit's general
power problem is resolved. Short regions remain weak, the correlated stress test
has no recovery, and broad error-rate validation is unfinished. The engine remains
experimental; the package default remains `chain_merge`.

## Source implementation

`ep.tl.dmr(md, method="regional_gls")` reads counts directly, discovers regions
from smoothed methylation effects, and fits CpG-level observations across samples.
Individual CpG significance is no longer a prerequisite. The model retains locus
intercepts, median-depth weights, continuous autoregressive covariance and a
depth-aware biological sample-offset component. All samples use identical,
positive-coverage CpG support. The current backend supports independent binary
groups with at least two samples each. Covariate and paired permutation designs
are rejected rather than silently treated as exchangeable samples.

The numerical kernel is compiled with Numba, already an Epykit dependency.
Conditional REML estimation of donor variance is analytical: after whitening the
CAR component, the donor direction and its orthogonal complement have separate
residual variance estimates with `S-2` and `(S-1)*(K-1)` degrees of freedom.
Their ratio determines the donor fraction. Only the spatial length needs a scalar
search. Batched region fits use up to eight native threads by default;
`n_workers` controls this. Discovery prefix sums and geometry scans are compiled.

Complementary balanced permutations have identical regions and opposite signed
scores. Each pair is fitted once while both scans remain represented in the null
pool. Geometry and transformed counts are computed once per chromosome.

The earlier speed shortcut fitting observed regions individually while assigning
decoys a chromosome-wide median covariance is disabled by default. It changes the
observed and null fitting procedures. It remains an explicitly warned comparison
option in the internal runner. Nonfinite scores abort calibration; they are never
silently omitted. The source backend uses per-region fitting for both observed
and decoy regions.

Complete results are cached under `analysis_root/.cache/regional_gls`. The cache
signature covers source code, sample labels, chromosomes, analysis parameters,
and count-file paths, sizes and nanosecond modification times. Cached tables have
a content digest. Count changes invalidate reuse. Cache output is separate from
read-only count inputs when `analysis_root` is supplied.

## Measured speed

These are measurements on the local machine, not hardware-independent promises.
The numerical model and observed/decoy fitting procedure are the same working
model; analytical optimisation can improve estimates that the previous two-start
numerical optimiser missed, so threshold-near calls need not be identical.

| Measurement | Previous per-region REML | Native implementation |
|---|---:|---:|
| 100-region numerical panel | 0.334 s | 0.012 s |
| Real chr22: 242,481 common CpGs, 1,258 candidates, 18 decoy scans | 70.37 s | 1.79 s |
| Full real dataset: 15,597,046 common CpGs, complete permutation inference | — | 12.13 s |
| Repeat full inference with valid cache | — | 0.126 s |

The panel speedup is 28-fold; the chr22 speedup is approximately 39-fold.
The full result contains 818 significant regions. Initial JIT compilation can
add startup time; compilation is cached. These timings cover the DMR backend,
not ingestion, optional DMC analysis, exports or HTML report generation.

Cached `chain_merge` took 0.006 s for 262 regions from an already-computed CpG
table. That stage merges existing inference results; its timing excludes the
count-level inference and permutation calibration measured for regional GLS.
Regional GLS does not require `ep.tl.dmc` first. Running both still computes two
analyses, which is useful only when CpG results are also wanted.

## Simulation provenance and known-truth evaluation

The primary datasets use the unchanged DMRcate-inspired v4 count generator and
the supplied 206-sample BLUEPRINT chr1–3 calibration: empirical baseline means,
mean-dependent beta-binomial dispersion, CpG coordinates, depth, library factors
and dropout. New seeds were used. The independent-residual track is the original
v4 model, not a correlated biological-cohort validation.

This is not an exact reproduction of the DMRcate paper: its benchmark used real
macrophage background counts and paired samples; v4 generates background counts
and uses independent groups. BLUEPRINT parameters span heterogeneous cell types.
Independent count draws do not constitute independently calibrated biological
cohorts. The original calibration and transferred frozen results were preserved.

The n=5 signal/null datasets each contain 5,907,954 CpGs across chr1–3. The two
n=3 signal/null pairs are development subsets with 200,000 CpGs per chromosome,
600,000 total. Effects are 0.20 and nominal mean depth is 20, with empirical
variation and missing coverage. Geometry was fixed before these comparisons:
500 bp smoothing, five CpGs, 50 bp minimum span, 500 bp gaps, effect floor 0.10.
The unchanged legacy comparator uses its frozen simulation settings: raw site
p<0.05, five CpGs, effect floor 0.10, 50 bp span, 500 bp gap and 0.5 seed fraction.

All reported calls below use nominal q≤0.05. Reciprocal matching is the frozen
scorer's **CpG overlap** criterion, not 50 bp overlap. One global maximum-cardinality
matching is used; direction agreement is 100% among matched native calls.

| Dataset / seed | Truth regions | Legacy recall, 50% | Native recall, 50% | Native recall, 80% | Native calls | Wholly null calls / FDP |
|---|---:|---:|---:|---:|---:|---:|
| n=5, 9917401 | 659 | 32.93% | 59.48% | 35.96% | 443 | 21 / 4.74% |
| n=3, 9917601 | 67 | 2.99% | 32.84% | 25.37% | 24 | 1 / 4.17% |
| n=3, 9917603 | 67 | 5.97% | 41.79% | 28.36% | 29 | 1 / 3.45% |

Native F1 at 50% reciprocal CpG overlap is 0.711, 0.484 and 0.583 respectively;
legacy F1 is 0.481, 0.058 and 0.113. The paired complete-null datasets
(9917402, 9917602, 9917604) yielded zero calls with both methods. Three clean null
draws are insufficient to establish general FDR control.

Wholly null FDP counts calls containing no changed CpGs on the exact tested
support. Boundary nonmatches are reported separately. Arithmetic raw-beta
mean-null diagnostics are also retained; they are not exact nulls for the
covariance-weighted transformed GLS effect, particularly with opposite effects
that cancel. Oracle matched-FDP thresholds are descriptive and must not become
deployable cutoffs.

The separate correlated stress extension uses an unvalidated spatial fit from
mixed BLUEPRINT tissues. For n=3, 60 signal regions, seed 9917411, native recall
was zero; legacy recovered one region. Both methods had zero discoveries in the
corresponding null draw (9917412). This is an unresolved power limitation, not
evidence that within-cohort spatial calibration has succeeded.

## Where the gain occurs and where power is still lost

In the full n=5 dataset, native candidates can reciprocally match 464/659 truth
regions before significance filtering; 392 are recovered after filtering.
Observed support and geometry constrain recovery as well as statistical testing.
Support exclusions are diagnostics, not absolute overlap-recovery ceilings.

| Truth stratum, n=5 | Regions | Legacy recall | Native recall |
|---|---:|---:|---:|
| 5–9 common CpGs | 199 | 20.1% | 40.2% |
| 10–19 common CpGs | 286 | 44.8% | 76.9% |
| At least 20 common CpGs | 99 | 49.5% | 91.9% |
| Span below 100 bp | 25 | 16.0% | 8.0% |
| Span 100–249 bp | 86 | 32.6% | 24.4% |
| Span at least 500 bp | 439 | 33.3% | 71.3% |
| Lowest dispersion tertile | 220 | 55.0% | 65.0% |
| Highest dispersion tertile | 220 | 10.5% | 43.6% |

Longer patterns gain from retaining moderately informative CpGs and their sample
structure. A fixed 500 bp smoother can dilute short effects below the discovery
floor. Low baseline methylation also remains a weakness: for mean baseline below
0.20, recall decreased from 50.6% to 29.2% in this n=5 draw. These findings argue
for separately validated multiscale/adaptive discovery; they do not justify
tuning settings against the same truth and claiming a new calibrated result.

## Real AKAP11 study

Against the 813 published reference regions, the corrected intersection legacy
baseline reports 262 calls and recovers 17.59% at at least 50 bp overlap. Native
GLS reports 818 calls and recovers 25.58%. Reference precision decreases from
54.96% to 24.94%, and reference F1 from 0.2665 to 0.2526. Hypo-reference recovery
decreases from 17.71% to 8.00%; hyper recovery increases from 17.55% to 30.41%.
This is improved overall reference recovery with a material boundary/precision
tradeoff, not a biological FDR estimate. The published list is not exhaustive
truth. The data support further evaluation rather than default promotion.

## Validation and reproduction

Numerical checks compare fixed-covariance estimates with independent dense GLS,
analytical donor profiling with numerical REML, serial with parallel fits, mirror
scores, and compiled discovery with the NumPy reference across window widths and
gaps. Cache tests verify reuse, invalidation and output isolation. The Epykit
suite has 23 passing tests; the benchmark project's canonical suite has 41.
The changed engine, runner and tests pass Ruff. Known legacy odds-ratio
deprecation warnings remain. A bare recursive benchmark collection also reaches
transferred repository snapshots and fails with `ImportPathMismatchError` for
their duplicated `tests.conftest`; the canonical `wgbs_benchmark_v4/tests` suite
passes. No claim is made that those transferred suites were validated.

For an isolated real-data run:

```sh
uv run python run_power_regional.py --samplesheet YOUR_SHEET.csv --treatment-group CASE --control-group CONTROL --input-store YOUR_FILTERED_COUNT_STORE --output-dir FRESH_OUTPUT --n-perm 100
```

The existing local `test.py` already calls `regional_gls` and uses the updated
source. DMC analysis is optional for this engine. Three-versus-three designs have
20 balanced assignments; excluding self and mirror leaves 18 scans, only nine
unique two-sided partitions. Increasing the requested permutation count beyond
20 does not create more assignments. Genome pooling provides tail resolution
under candidate-score comparability assumptions, not 18 independent null
experiments or exact per-region randomization guarantees. Empty scans are retained
and recorded, but contribute no candidate scores to the pooled tail denominator.

Aggregate results and strata are included under `docs/benchmarks/regional_gls`.
Local simulation scripts and complete outputs remain in
`wgbs_benchmark_v4/development/power_prototype_20261009`.

The architecture is informed by the
[dmrseq implementation](https://github.com/kdkorthauer/dmrseq), independently
implemented and extended with the depth-aware sample-effect component.
Further promotion requires untouched mixed-signal and null evaluations with
validated within-cohort spatial dependence, heterogeneous depths/dispersion,
outliers, missing coverage and valid nuisance designs.
