# R5: document compatibility and retain deprecated APIs

Start from R4 only after S3 merges and the refactor stack is updated onto main.
Use `refactor-5-compatibility` and target `refactor-4-cli-package`.
Read the [run rules](README.md) and issue [22](https://github.com/d-ert/epykit3/issues/22).

The filename is retained so existing links resolve. This brief replaces the proposed removal layer.
The previous resolution treated deprecation messages as approval for a breaking change. It was also factually wrong:
`pp.unite` names 2.0, while the CSV aliases say only "a future release".
No maintainer-approved removal is established by those messages.

Own current compatibility docs, README and CLAUDE migration guidance, and deprecation wording in `tl.py`, R1's `finish` stage, `dmc.py`, and `dmr_hmm.py`.
This layer starts after the other owners of those files finish.

## Implement

1. Keep `pp.unite`, `epykit.dmr_hmm`, the transitional `log2_odds_ratio` column, and all CSV keyword aliases.
2. Keep the existing pp.unite 2.0 promise and generic CSV deprecation wording.
3. Replace the unapproved 1.2 removal promise for the transitional column and dmr_hmm shim with "a future major release". Keep warning categories and behavior. Update tests only where they assert that release text.
4. Document the actual replacements: `pp.set_unite_type`; `epykit.dmr_segment.call_dmr_rule_segment`; `log2_odds_ratio_pooled` for pooled effects and `coef_treatment_log2` for GLM effects; and the matching TSV keyword aliases.
5. Keep file format distinct from keyword naming. A TSV keyword with a .csv path still requests comma-delimited output.
6. Check CSV aliases on `tl.qc`, `tl.dmc`, `tl.dmr`, `tl.dvc`, and `tl.annotate`. Do not partially remove the shared resolvers or break these other callers.
7. Update current docs and CHANGELOG to state that 1.2 retains compatibility. Preserve historical release notes.

## Accept when

All existing signatures, aliases, transitional columns, and warning categories remain.
The metadata test keeps its schema assertions, including `log2_odds_ratio`; it is not a version-string test.
Existing alias, export, and deprecation tests pass with only planned wording updates.
The code-layer gates and strict docs build pass.

A future removal requires a separate maintainer decision on version and migration.
It is out of scope for this run and does not block these compatibility changes.

PR title: `Retain deprecated APIs and correct the 1.2 migration guidance`.
