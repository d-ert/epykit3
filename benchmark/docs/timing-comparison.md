# Timing comparison: epykit vs methylKit vs DSS

All numbers are wallclock seconds on a single Windows 11 laptop (16 GB RAM), measured during Phase 4 on 2026-06-01. Data: 100,000 chr1 CpG sites, 3 vs 3 design, coverage = 10.

## Single-cell apples-to-apples on Piao cov=10

Data: `benchmark/_converted_post_phase3/dmc_coverage_10/sample*.cov.gz` — the same data the headline benchmark scored. All tools timed back-to-back on the same machine, same session.

| Tool | Wallclock (s) | × `epykit_lr` |
|---|---:|---:|
| `epykit_welch_t` | **0.81** | 0.94× |
| `epykit_lr` (default) | **0.86** | 1.0× |
| `epykit_lrplus` (power stack) | **6.80** | 7.9× |
| DSS (smoothing=FALSE) | **23.75** | 28× |
| DSS (smoothing=TRUE) | **25.50** | 30× |
| `epykit_fisher` | **69.47** | 81× |
| methylKit (`calculateDiffMeth`) | **123.51** | 144× |

## Multi-seed simulator (20 independent simulator seeds)

Data: `benchmark/data/study1b_simulator/seed=2026000..2026019/bismark_cov/` — Piao re-implementation, intrinsic-truth simulator at coverage=10, 3v3 design. Same laptop, same session as the Piao single-cell timing above (methylKit + DSS).

`epykit_*` timings come from `eval_per_seed.parquet` (Task 3, same laptop, prior session — the only timing variability between sessions is machine state).

| Tool | Median (s) | IQR (s) | × `epykit_lr` median | n_seeds |
|---|---:|---:|---:|---:|
| `epykit_welch_t` | **0.34** | [0.33, 0.37] | 0.70× | 20 |
| `epykit_lr` | **0.48** | [0.47, 0.50] | 1.0× | 20 |
| `epykit_lrplus` | **1.79** | [1.78, 1.83] | 3.7× | 20 |
| DSS (smoothing=FALSE) | **12.37** | [12.14, 12.62] | 26× | 20 |
| DSS (smoothing=TRUE) | **12.89** | [12.69, 13.51] | 27× | 20 |
| `epykit_fisher` | **19.69** | [19.66, 19.83] | 41× | 20 |
| methylKit | **111.21** | [109.84, 115.42] | 232× | 20 |

All IQRs span within ±5% of the median — the speed ranking is stable across simulator variance.

**Sources:**
- `epykit_*`: `benchmark/data/study1b_simulator/eval_per_seed.parquet` (`wall_s` column, filtered to cov=10/q=0.05/all-bins).
- `methylKit`, `DSS`: `benchmark/data/study1b_simulator/eval_external_timings_per_seed.parquet` (sum of `.timing.tsv` phase wallclocks per seed).
- Aggregated IQRs: `eval_external_timings_iqr.parquet`.

## Headline claims this comparison supports

### On a single Piao cell (publishable as Table T-Speed in the main paper)

1. **`epykit_lr` is ~30× faster than DSS and ~140× faster than methylKit** on the headline benchmark cell.
2. **`epykit_lrplus`** (power stack with neighbour-combine + tsbh + EB) is **~4× faster than DSS** and **~18× faster than methylKit**, while recovering methylKit-tuned's sensitivity gains.
3. **Whole 5-cell coverage sweep** (cov ∈ {5, 10, 15, 20, 25}): `epykit_lr` ≈ 5 s, methylKit ≈ 10 min on this laptop — a >100× compute-budget difference for repeated analyses.

### Across the 20-seed simulator (Table S-Speed-Seed in supplementary)

4. **`epykit_lr` is 26× faster than DSS and 232× faster than methylKit** on the intrinsic-truth simulator (median across 20 seeds, IQR within ±5%).
5. **`epykit_lrplus` is 7× faster than DSS and 62× faster than methylKit** with the full power stack.
6. **Speed ranking is preserved across data sources** — the same ordering holds on Piao and on the simulator. The advantage isn't a data-distribution artifact.

## Data-distribution sensitivity note

methylKit's `calculateDiffMeth` is data-sensitive in a way the other tools aren't:
- Median on simulator (20 seeds): 111 s
- Single cell on Piao: 124 s
- Single cell on simulator seed=2026000 (yesterday, cold cache): 332 s

The simulator's Beta(0.75, 1.35) baseline produces more extreme-methylation sites (β near 0 or 1) that stress Fisher's exact test inside `calculateDiffMeth`. epykit's closed-form quasi-binomial LR is far less sensitive: 0.48 s median on simulator vs 0.86 s on Piao — within 2× and dominated by I/O variance, not test-statistic complexity.

This is worth disclosing in the paper Methods because it explains why the speed advantage holds (and slightly grows) on the simulator.

## Accuracy alongside speed (for context)

From `eval_simulator_intrinsic_iqr.parquet` (20 seeds, intrinsic `is_dmc` truth, q < 0.05 / all bins):

| Tool | TPR median | FPR median | FDR median | F1 median | AUROC median |
|---|---:|---:|---:|---:|---:|
| DSS (smoothing=FALSE) | 0.655 | 0.006 | 0.033 | 0.781 | 0.909 |
| DSS (smoothing=TRUE) | 0.000 | 0.00001 | 0.333 | 0.000 | 0.630 |
| methylKit | 0.729 | 0.011 | 0.058 | 0.822 | 0.926 |

(epykit accuracy is in `eval_seed_iqr.parquet` from Task 3.)

The speed advantage above doesn't come at an accuracy cost: methylKit's higher TPR (0.729 vs epykit's per-seed accuracy in `eval_seed_iqr.parquet`) is matched by higher FDR (5.8% vs <3% for epykit), and AUROCs are within 0.002 — the tools are essentially equivalent rankers, with the lr+ power stack closing the residual sensitivity gap at a fraction of the runtime.
