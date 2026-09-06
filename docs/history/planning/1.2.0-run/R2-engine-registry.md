# R2: centralize DMC engine facts and validation

Start from `refactor-1-stage-split` after its implementation and gates pass.
Use `refactor-2-engine-registry` and target R1.
Read the [run rules](README.md) and issue [14](https://github.com/d-ert/epykit3/issues/14).

Own the new `_dmc_engines.py`, `_dmc_config.py`, DMC stages, `dmc.py`, CLI engine choices, and their tests.
R1 has already moved two name-normalization call sites from `tl.py` into stages. Update those too.

## Implement

1. Add a dependency-light registry with frozen `EngineSpec` records for lr, glm, welch_t, fisher, and internal glm_contrast.
2. Give each record only the consumed fields: `name`, `public`, `power_stack_applies`, and `effect_column`. Drop the proposed unused `supports_n1`, `needs_design`, and `uses_glm_backend` fields. Preserve the existing sample-size and design validation.
3. Keep `REMOVED_ENGINES` as data with existing migration messages. Set `PUBLIC_ENGINES` to `("lr", "glm", "welch_t", "fisher")`, which is the current CLI order.
4. Make `DMCConfig.validate` reject removed and unknown names. Keep auto selection in its current high-level path. Do not expose glm_contrast as a CLI choice.
5. Validate resolved names at `process_chromosomes_dmc` too. The CLI, tile DMR path, and public low-level callers bypass `DMCConfig`. Raise a useful ValueError before opening or creating an output store, rather than allowing a runner-map KeyError.
6. Replace concrete registry-backed checks. Keep `_auto_test_simple`, sample-size warnings, and the numerical algorithms unchanged.
7. Delete the identity function `_canonicalise_test_name` and every reference, including imports and fields in R1's stages. Remove `canonical_test` from `DMCPlan` when it duplicates `selected_test`.
8. Use `PUBLIC_ENGINES` for both CLI choice lists without changing their order or defaults.
9. Correct `empirical_fdr_for_dmc`'s cross-reference after S1's contract. Its own algorithm remains max-T style. Only the DMR function accepts `fdr_method="region"`.

## Accept when

- Both CLI engine lists equal the registry's public names in the original order.
- Each removed name retains its migration error through `tl.dmc`.
- An unknown name fails before output creation through both `tl.dmc` and `process_chromosomes_dmc`.
- Auto selection and valid GLM contrast calls retain their previous behavior.
- Registry import does not import `dmc.py`; no import cycle exists.
- No `_canonicalise_test_name` reference remains.
- Code-layer gates pass. Only invalid-name error type and timing may change.

PR title: `Centralize DMC engine facts and validate low-level calls`.
