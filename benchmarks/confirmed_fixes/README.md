# Confirmed-fix measurements

`benchmark_components.py` runs one fixed-seed fixture per process. Invoke with
`--repo`, `--component glm|tss|masked`, and a new JSON `--output` path. Use
`--size 1000000` for TSS, `--invalid` for the single-invalid GLM fixture, and
`--samples 6|10 --k 5|20|100 --dropout 0|.01|.2 --inner 100` for masked preparation.
Run workers 1 and 4 with OPENBLAS_NUM_THREADS, VECLIB_MAXIMUM_THREADS and
OMP_NUM_THREADS all set to 1. Each subprocess warms once, measures five repeats,
and records source/input hashes, environment and process peak RSS. Set
`--repeats 10` for the specification's noisy-result repeat. Warm-up includes
compilation/import costs, reported separately. These fixtures measure component
costs and do not establish statistical calibration.

`profile_components.py --workspace <benchmark-root> --repo <checkout> --output
<new-profile-directory> -- --kind pilot --datasets empirical_signal_seed2026
--configs cov10_boundNone --workers 4` wraps the existing frozen-source harness.
Use `--kind real --chromosomes chr1` for a real chromosome. Preserve python.prof
and stages.json; inclusive nested timings must not be summed. Profile cold and
warm processes separately using separate output directories. cProfile identifies
Polars join time in its Python profile; the residual outer runtime also includes
count loading, validation and calibration. No production instrumentation is added.

End-to-end variant sweeps may overlap to shorten correctness/recovery validation.
Those times are explicitly unsuitable for performance claims. Performance gates
use dedicated component subprocesses without concurrent sweeps; preserve all
measurements and record any repeat. Existing results are never overwritten.
