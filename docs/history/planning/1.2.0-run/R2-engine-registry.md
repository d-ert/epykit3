# R2: the DMC engine registry

Branch `refactor-2-engine-registry` from `refactor-1-stage-split`. PR against it. Decision: map ticket 14. Read `_dmc_config.py`, the top of `dmc.py` (`_GLM_BACKENDS`, `_canonicalise_test_name`, `_TEST_RECOMMENDATIONS`), `_validate_sample_size_and_warn`, `tl.py` `_auto_test_simple`, and the `--test` choices in `cli.py`.

You own `_dmc_engines.py` (new), `_dmc_config.py`, `dmc.py`, the `--test` choices in `cli.py`, and the tests you add.

## Commits, in this order

1. `feat(engines): declarative engine registry`
   - `src/epykit/_dmc_engines.py`: a frozen `EngineSpec` with exactly these fields: `name`, `public`, `supports_n1`, `power_stack_applies`, `uses_glm_backend`, `effect_column`, `needs_design`. `ENGINES: dict[str, EngineSpec]` for `lr`, `welch_t`, `fisher`, `glm`, `glm_contrast` (`public=False`). `REMOVED_ENGINES: dict[str, str]` moved from `_dmc_config.py` with the hint text unchanged. `PUBLIC_ENGINES: tuple[str, ...]` in today's CLI order. No callables, no aliases. The module imports nothing from `dmc.py`.
2. `refactor: consult the registry instead of comparing names`
   - `DMCConfig.validate`: removed names raise their hint; a name not in `ENGINES` (and not `"auto"`) raises `ValueError` naming the four public engines. `DMCConfig.apply_power_stack`: `ENGINES[selected_test].power_stack_applies` replaces `selected_test != "lr"`. `dmc.py`: `effect_column` replaces the `_GLM_BACKENDS` check for the log2 column; `uses_glm_backend` replaces the set elsewhere; delete `_GLM_BACKENDS` and `_canonicalise_test_name` with its four call sites. `cli.py`: both `--test` `choices` lists come from `PUBLIC_ENGINES`. `_auto_test_simple`, `_warn_fisher_once` and `_validate_sample_size_and_warn` are not changed.
3. `test(engines): registry contract`
   - Three tests, fast tier: public names equal the CLI choices and the documented set; the removed names raise their hints through `tl.dmc`; `tl.dmc(test="nope")` raises `ValueError` before any store directory is created (assert on the filesystem under the analysis root).
4. `docs(dmc): empirical_fdr_for_dmc docstring`
   - The hand-off from S1: correct the docstring to describe both DMR `fdr_method` values with `max_t` as the default.

## Contract

No results change. The only behaviour change is for invalid engine names (error type and timing). Regen hashes unchanged.

## Deliver

PR title: `One registry for the DMC engines`. Then `worker_done`.
