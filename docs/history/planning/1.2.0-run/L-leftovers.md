# L: prepare the shared inputs

Start now from current `main` on `leftovers-1.1.0`. Target the PR at main.
This layer must merge before either stack starts.
Read the [run rules](README.md). Planning sources are issues [12](https://github.com/d-ert/epykit3/issues/12) and [27](https://github.com/d-ert/epykit3/issues/27).

Own `README.md`, `CLAUDE.md`, the archive index and field guide, `src/epykit/_chroms.py`, and `tests/test_chroms.py`.

## Implement

1. Correct current claims that the lr+ CLI options are deferred to 1.1. State that the lr+ options are Python API only, with no promised release. Preserve dated historical plans and release notes.
2. Copy [the supplied field guide](inputs/field-guide.html) to `docs/history/2026-09-05-field-guide.html` on the implementation branch. Add its date and purpose to `docs/history/README.md`. The input is archived on the planning branch so a worker does not need access to the original untracked file. Its SHA-256 is `484fe917697579c0527e41e07c0649d6174e7b85dc727ab99f9e44a6e4be5dfd`. Preserve the original file in the main checkout. `mkdocs.yml` already excludes `history/`; do not add a redundant exclusion.
3. Port only the chromosome helper and its tests from commit `cd9f89b3fd7a56bd3f9c7fbd2bbdf78aac5f1676`. Add `is_canonical_chrom`, `filter_canonical`, `filter_canonical_logged`, `CANONICAL_CHROM_CORES`, and `CANONICAL_CHROMS_UCSC`. Use current typing and lint conventions.
4. Keep the existing helper semantics: case-insensitive optional chr prefix, identifiers 1 through 22, X, Y, M, and MT, and order-preserving filtering. Document that this is a fixed human-style chromosome list, not a species-aware assembly validator.
5. Correct the donor helper's default-on wording and `--all-contigs` advice. Filtering is opt-in. The INFO message reports the dropped contigs and says to omit `canonical_only=True` to retain them.

No branch deletion is part of this layer. Retain `fine-tune`, `optimize/autonomous-v1`, `bench-tune`, and `feat/canonical-chrom-filter`.

## Accept when

- Helper tests cover both naming styles, prefix case, mitochondrial names, excluded contigs, input order, and the INFO message.
- The field guide copy has the recorded hash and remains excluded from the site.
- The helper has no callers yet. Existing defaults and numerical output are unchanged.
- The code-layer gates in the run rules pass.

PR title: `Prepare the 1.2 shared chromosome helper and archive the field guide`.
