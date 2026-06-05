# M5 null calibration report — engine `lr`

**Dataset:** GSE263850 (label-permuted, exhaustive k=10, 3v3 n=6)
**DMC config:** `ep.tl.dmc(test='lr')  [default dispersion='eb']`
**Verdict:** **calibrated (slightly conservative)**

## Calibration vs conservatism

| quantity | observed | nominal/uniform |
|----------|---------:|----------------:|
| mean null p-value | 0.5055 | 0.5 |
| frac p<0.01 | 0.0096 | 0.01 |
| frac p<0.05 | 0.0472 | 0.05 |
| frac p<0.10 | 0.0987 | 0.10 |
| KS D vs Uniform | 0.0506 | 0 |
| mean(observed − expected) | +0.0055 | 0 |

`mean(observed − expected) > 0` ⇒ p-values run slightly *larger* than uniform ⇒
mildly **conservative**: FDR control is valid, with a small (<6% relative) power
cost at the tails. Not anti-conservative.

## Observed FDR under the null (q<0.05)

median = 1.764e-05, IQR = [1.604e-05, 2.028e-05]
across the 10 exhaustive label permutations (n=6 ⇒ C(6,3)/2 = 10 unique).

This replaces the old single-point "1.53e-5" (which could not distinguish
calibrated from conservative) with the distribution **and** the directional
Q-Q/KS verdict above. Q-Q figure: `benchmark/data/null_calibration/gse263850/lr_qq.png`.
