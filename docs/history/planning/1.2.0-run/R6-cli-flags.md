# R6: CLI flags and the DMC side of canonical_only

Branch `refactor-6-cli-flags` from `refactor-5-removals`. PR against it. Top of the refactor stack. Start only after S3 has merged and the coordinator has rebased the refactor stack onto `main`, because this layer uses `convert_sample(canonical_only=)` from S2 and the helper from L. Decisions: map tickets 27 and 28.

You own `tl.dmc`, `_run_dmc_contrast` (now `run_contrast` in `_dmc_stages.py`), `_dmc_config.py`, `dmc.py`, the `cli/` package, and their tests.

## Commits, in this order

1. `feat(dmc): canonical_only on process_chromosomes_dmc and tl.dmc, default off`
   - `process_chromosomes_dmc(canonical_only: bool = False)`: when on and `chromosomes` is auto-detected, filter with `filter_canonical_logged` and log the dropped contigs; an explicit list is honoured verbatim. `DMCConfig` gains the field; `plan_run` passes it; `run_contrast` too. `to_uns` records it. Carry the DMC tests from the branch's `tests/test_canonical_chrom_filter.py` with the default inverted (they were left out of S2 for you).
2. `feat(cli): --canonical-only on convert, dmc and dmr`
   - Default off. `convert` forwards to `convert_sample`; `dmc` and `dmr` forward to the Python calls. No `--all-contigs`.
3. `feat(cli): --smoothing and --smoothing-span-bp on dmc and dmr`
   - Mapped to the existing Python parameters, defaults unchanged. Add the parity cases to `tests/test_cli_api_parity.py`.
4. `docs: the new flags`
   - The cli page, README CLI table, CHANGELOG under Added.

## Contract

Default off everywhere; no existing number moves. Regen hashes unchanged.

## Deliver

PR title: `CLI flags for canonical_only and DSS-style smoothing, canonical_only on the DMC engine`. Then `worker_done`.
