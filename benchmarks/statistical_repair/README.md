# Repair benchmark harness archive

These files archive the exact workspace entry points and their regression tests.
They run from the benchmark workspace layout (scripts/, wgbs_sim/, data/, results/),
not directly from this archive directory. The Gemini delivery includes that layout,
the transitive helper scripts and simulator. Source/evaluation signatures and count,
truth, reference and output hashes are recorded alongside every final result.

Commands in the workspace:

```sh
PYTHONDONTWRITEBYTECODE=1 development/epykit_minimal_gls_trial/.venv/bin/python scripts/run_statistical_repair_benchmarks.py --kind pilot --include-dmc --output results/statistical_repair_final_full
PYTHONDONTWRITEBYTECODE=1 development/epykit_minimal_gls_trial/.venv/bin/python scripts/run_statistical_repair_benchmarks.py --kind real --include-dmc --output results/statistical_repair_final_real
PYTHONDONTWRITEBYTECODE=1 development/epykit_minimal_gls_trial/.venv/bin/python scripts/run_statistical_repair_validation.py --dataset-root data/statistical_repair_null_holdout --output results/statistical_repair_final_null_holdout
```

Experimental regional FDP estimates are not finite-sample FDR guarantees. Whole-scan
maxT is a joint all-null diagnostic requiring exchangeable samples and the same
per-assignment statistic; strong region control under partial nulls is unproved.
Real data have 3+3 samples, so the exact balanced two-sided randomization minimum is
2/20=0.1. Reference-paper overlaps are concordance, not ground truth.
