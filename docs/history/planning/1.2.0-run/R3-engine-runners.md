# R3: decompose _process_one_chromosome and set the ceiling

Branch `refactor-3-engine-runners` from `refactor-2-engine-registry`. PR against it. Decision: map ticket 15. Read `_process_one_chromosome` in `dmc.py` end to end before editing; it is about 600 lines, complexity 38, the current ceiling.

You own `dmc.py`, `pyproject.toml` (the mccabe setting only), and the tests you add.

## Commits, in this order

1. `refactor(dmc): EngineInput and EngineResult`
   - Two frozen dataclasses in `dmc.py`. `EngineInput`: `methylstore_path`, `chrom`, `canonical_pos`, `n_sites`, `samples_case`, `samples_control`, the engine knobs (`dispersion`, `reference`, `smoothing`, `smoothing_span_bp`, `sep_fallback`, `sep_threshold`, `glm_backend`) and the design bundle (`design_full`, `design_reduced`, `coef_idx`, `contrast_matrix`, `contrast_label`, `samples_all_ordered`, `group_labels_per_sample`). `EngineResult`: `pvals`, `log2_ors`, the two Welford triples, `extras`, an optional Newcombe block (`meth_a`, `cov_a`, `meth_b`, `cov_b`, `phi`, `df`) and an optional multigroup block (`level_mean_beta`, `f_stat`, `df1`, `df2`).
2. `refactor(dmc): one runner per engine`
   - `_run_fisher`, `_run_lr`, `_run_welch_t`, `_run_glm`, `_run_glm_contrast`, each `(inp: EngineInput) -> EngineResult`, bodies moved from the five branches without edits beyond the record plumbing. A runner returns reduced per-site arrays only, never a sample stack; the per-branch `del` statements become scope exits. `_ENGINE_RUNNERS: dict[str, Callable[[EngineInput], EngineResult]]` at the bottom of `dmc.py`.
   - `_process_one_chromosome` becomes: the empty guard, build `EngineInput`, look up the runner (a `KeyError` cannot happen: the registry validated the name), call `_finalise_chromosome`.
3. `refactor(dmc): _finalise_chromosome and _effect_ci`
   - Everything after the branches moves: the min-samples NaN pass, `_effect_ci(res)` (Newcombe when the pooled block is present, delta method when `coef_treatment` is in extras, Welch otherwise, NaN in multigroup mode), column assembly, the sort.
4. `test(dmc): runner keys match the registry`
   - One test: `set(_ENGINE_RUNNERS) == set(ENGINES)`.
5. `chore(lint): lower the complexity ceiling to 32`
   - `[tool.ruff.lint.mccabe] max-complexity = 32`. Report in the commit body the new numbers for `_process_one_chromosome`, each runner and `_finalise_chromosome`, and the three highest functions in the tree (expected: `_build_features_index` 32, `call_dmr_chain_merge` 28, `process_chromosomes_dmc` 27). Targets: `_process_one_chromosome` at most 6, `_run_lr` and `_run_glm_contrast` at most 15, `_finalise_chromosome` at most 12. If a target is missed, restructure; do not raise the ceiling.

## Contract

No results change; this is the layer where the hash gate matters most. Run the regen script after every commit. Memory stays per chromosome: no runner may hold more than one chromosome's arrays.

## Deliver

PR title: `Per-engine runners for the chromosome loop, complexity ceiling to 32`. Then `worker_done`.
