# S3: the additive parts of bench-tune

Branch `salvage-3-bench-tune-additive` from `salvage-2-canonical-filter`. PR against it. Top of the salvage stack. Decision: map ticket 28. Read `origin/bench-tune` commit `66ee2ef` (which includes `fine-tune`'s `47d7878`) for the diff; re-type by hand, cherry-pick nothing.

## Commits

1. `feat(dmr): record the smoothing knobs in dmr_params`
   - `tl.dmr` records `smoothing` and `smoothing_span_bp` in `md.uns["dmr_params"]` exactly as the branch did (the span only when smoothing is on).
2. `docs(dmr): docstring fix from bench-tune`
   - Re-type the docstring correction from the branch if it is in a file you own; if it is in `tl.dmc` or `dmc.py`, hand it to the refactor stack in `worker_done` instead.

## Not here

The `min_cpgs` 5 to 3 and permissive `pct_sig` 0.5 to 0.4 default changes (maintainer question, out of scope), `AGENTS.md`, the calibration plan, the lock change, and the CLI flags (R6 in the refactor stack).

## Contract

No behaviour change; one metadata dict gains two keys. Regen hashes unchanged.

## Deliver

PR title: `Record the smoothing knobs in dmr_params`. Then `worker_done`.
