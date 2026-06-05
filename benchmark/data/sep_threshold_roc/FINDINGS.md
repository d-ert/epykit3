# M11(b) sep_threshold ROC — findings (2026-06-04, corrected 2026-06-05)

**Result: `sep_threshold` is inert on the Piao-style simulator at every available
coverage, so a simulator ROC cannot justify the 0.9 default.** This is a real
property of the data, confirmed by cold-store recompute (NOT a cache artifact —
see correction below).

## What was tested

The lr+ power stack's separation fallback (`sep_fallback=True`, default
`sep_threshold=0.9`) is meant to rescue quasi-separated CpGs (a group is ~all-
or ~none-methylated) where the bare quasi-binomial LR misbehaves. M11(b) asked
for a ROC justifying 0.9.

A ROC needs ground truth (TPR **and** FPR), so this was run on the intrinsic-
truth simulator (not null-permuted real data, which has no true positives and
would give only an FPR axis).

## Result (cold store per (seed, sep_threshold), full recompute each step)

`run_sep_threshold_roc_simulator.py`, 3 headline seeds, cov=10, 3v3. Every
`sep_threshold ∈ {0.7, 0.8, 0.9, 0.95}` gives identical metrics:

| sep_threshold | TPR | FPR | F1 | AUROC | n_called |
|--------------:|----:|----:|---:|------:|---------:|
| 0.7 / 0.8 / 0.9 / 0.95 (median over 3 seeds) | 0.7494 | 0.0652 | 0.7453 | 0.9055 | 20222 |

Cold-store cross-check (fresh store per call, sep_fallback toggled fully off too):
output is byte-identical with `sep_fallback=False` and across sep_threshold at
both cov=10 and cov=5. The fallback never engages: the simulator at 3v3,
coverage ∈ {5,10,15,20,25} contains no quasi-separated sites, so the threshold
(or disabling the fallback) changes nothing. The flat ROC is a property of the
data, **not** evidence for or against 0.9.

## Correction re: the cache (2026-06-05)

An earlier version of this note claimed the first sweep's flat numbers were a
DMC-cache artifact. **That was wrong.** `_dmc_input_signature` in `dmc.py`
(lines 146-148) already hashes `sep_fallback` and `sep_threshold`, so the
DMCStore engine cache — the layer used on the `resumable=False` path this sweep
takes — correctly invalidates when the threshold changes. The cold-store rerun
reproduces the original numbers exactly (e.g. n_called=20069 for seed 2026000),
confirming the original sweep was already correct.

A genuinely separate cache gap *did* exist and was fixed independently in
`tl.py`: the `resumable=True` resume_sig omitted the lr+ knobs
(`power_stack`/`sep_fallback`/`sep_threshold`/`neighbour_combine`/`neighbour_bp`/
`fdr_method`). That matters only for `resumable=True` sweeps (which neither this
script nor the M5 calibration use), but it is the right defensive fix. Note the
DMCStore-layer signature still omits the *post-processing* knobs
(`neighbour_combine`/`neighbour_bp`/`fdr_method`/`power_stack`) — correctly, since
those are applied downstream of the cached per-site engine output, not inside it.

## Real-data confirmation (GSE263850, cov>=5, 2026-06-05)

`run_sep_prevalence_gse.py` on the real contrast (3 Het_AKAP11_KO vs 3 WT,
21,993,377 CpGs) shows the fallback is inert here too: at every
sep_threshold ∈ {0.7,0.8,0.9,0.95} there are **0 candidate sites** and the call
count is unchanged at **48,185** DMCs (q<0.05).

| sep_threshold | candidates | n_called | rescued |
|--------------:|-----------:|---------:|--------:|
| baseline (off)| –          | 48,185   | –       |
| 0.7/0.8/0.9/0.95 | 0       | 48,185   | +0      |

Reason: the fallback fires only for sites with |meth_diff| >= threshold (>=0.7)
that the LR test *missed* (p>0.05). At coverage>=5 (the pipeline's filter), the
LR test already rejects every such large-effect site, so the candidate set is
empty. The fallback is a safeguard for pathologically low coverage (<= ~2 reads),
which the standard cov>=5 filter removes.

**Paper-ready conclusion:** sep_threshold=0.9 has no effect on any reported
result; the separation fallback never activates at coverage>=5 on either the
simulator or real WGBS. Its exact value therefore needs no tuning. (To *see* it
fire, one would rerun without the coverage filter, cov>=1 — not done; the
pipeline always filters at cov>=5.)

Output: benchmark/data/sep_threshold_roc/sep_prevalence_gse263850.csv

## How to actually justify sep_threshold=0.9

A simulator ROC is the wrong instrument here. Options:
1. **Count separation events vs threshold on real data.** On GSE263850 (or a
   low-coverage cohort) report how many CpGs are flagged quasi-separated at each
   `sep_threshold`, and the effect on their p-values — the FPR/specificity
   evidence the reviewer can use (no truth needed).
2. **Generate a low-coverage separation regime** (cov ≤ 2, or higher group
   imbalance) where separation actually occurs, then run the truth-based ROC
   there. The current simulator grid (cov ≥ 5) does not reach that regime.
3. State plainly that on realistic-coverage data the fallback is a rare-event
   safeguard, not a routinely-active knob — which is why the exact threshold has
   little effect at typical coverage.
