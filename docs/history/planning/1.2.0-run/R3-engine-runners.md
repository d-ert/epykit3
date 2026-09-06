# R3: extract the chromosome engine runners

Start from `refactor-2-engine-registry` after its implementation and gates pass.
Use `refactor-3-engine-runners` and target R2.
Read the [run rules](README.md), issue [15](https://github.com/d-ert/epykit3/issues/15), and the complete current `_process_one_chromosome`.

Own `src/epykit/dmc.py`, the complexity limit in `pyproject.toml`, and the relevant tests.

## Implement

1. Introduce a frozen `EngineInput` record. Carry the full `canonical_df`, including pos and strand, both sample lists, the selected test, and both minimum-sample thresholds. Carry every engine and design argument consumed by the original function. Derive `canonical_pos` and site count from that frame. The original brief omitted strand, selected test, and the minimum-sample thresholds.
2. Introduce `EngineResult` for reduced per-site outputs: p-values, effect coefficients, both Welford triples, extras, pooled counts and dispersion inputs for the Newcombe interval, and the existing multigroup outputs.
3. Move the five existing branches into `_run_fisher`, `_run_lr`, `_run_welch_t`, `_run_glm`, and `_run_glm_contrast`. Give them the common `EngineInput -> EngineResult` signature. Define a module-level `_ENGINE_RUNNERS` table.
4. Keep `_process_one_chromosome` as the empty guard, input construction, validated runner lookup, and `_finalise_chromosome(inp, result)`.
5. Move effect estimates, intervals, minimum-sample masking, column construction, and sorting into the finalizer without changing their order. Keep the GLM adjusted `meth_diff` and dispersion-scaled `coef_se`, not just its interval bounds. Keep empty and multigroup schemas.
6. Extract a smaller effect helper only if needed to meet the complexity target. If extracted, return every value the effect block changes. Do not lose the adjusted difference or standard error by returning only the interval.
7. Keep temporary sample arrays inside each runner. Fisher and Welch must retain their current sample-by-sample accumulation. The existing LR and GLM array requirements may remain, but no runner may retain arrays for another chromosome.
8. Set `max-complexity = 32` only after the whole source tree meets it. Measure the highest remaining functions on the integrated code, rather than assuming the earlier numbers still hold.

## Accept when

- Runner keys equal registry engine names.
- Existing numerical tests cover all five paths, minimum-sample masking, empty input, strand preservation, single contrasts, and multigroup output.
- Compare pre-refactor and post-refactor results on a fixed small fixture for each path, including NaNs, dtypes, effect sizes, intervals, and GLM extras. Reuse existing fixtures.
- The LR hash gate remains unchanged. That gate does not cover Fisher, Welch, or all GLM behavior.
- No runner returns a sample stack or materializes another chromosome.
- Report measured complexity. Target at most 6 for orchestration, 15 for the largest runner, and 12 for the finalizer. The enforced source-tree ceiling is 32.
- All code-layer gates pass.

PR title: `Extract per-engine chromosome runners and lower the complexity ceiling`.
