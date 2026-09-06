# R5: the 1.2 removals

Branch `refactor-5-removals` from `refactor-4-cli-package`. PR against it. Decision: map ticket 22. Every removal is announced by a `FutureWarning` or `DeprecationWarning` naming 1.2; grep for `1.2` in `src/` to find them all.

You own every file a removal touches; name any file outside the refactor stack's list in `worker_done` (there should be none except `pp.py` and `dmr_hmm.py`, which no other stack touches).

## Commits, one per removal, in this order

1. `feat(dmc)!: remove the transitional log2_odds_ratio column`
   - The NaN-filled column, its `_EMPTY_SCHEMA` entry, the FutureWarning in the `finish` stage, and any reader that special-cases it. `log2_odds_ratio_pooled` and `coef_treatment_log2` stay.
2. `feat(pp)!: remove pp.unite`
   - The function and its DeprecationWarning; the documented replacement stays.
3. `feat!: remove the dmr_hmm shim`
   - `src/epykit/dmr_hmm.py`; `dmr_segment` is the name.
4. `feat(tl)!: remove the csv keyword aliases`
   - `csv`, `csv_full`, `csv_alpha` on `tl.dmc` and `tl.dmr` (and `DMCConfig`); `_resolve_auto_tsv` loses its csv arguments.
5. `docs(changelog): the 1.2 removals`
   - A `### Removed` section under `[Unreleased]` with one line per item and the replacement.

For each removal, delete the tests that asserted on the warning text and update `tests/test_dmc_metadata.py` only where it checks the version string. Update `README.md`, `CLAUDE.md` and the docs pages that mention the removed names.

## Contract

Results are unchanged for every call that did not use a removed name. Regen hashes unchanged. The PR is the maintainer's review point; if they object to one item, that commit is dropped and the others stay.

## Deliver

PR title: `Remove the four items deprecated for 1.2`. Then `worker_done`.
