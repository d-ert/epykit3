# L: leftovers from the 1.1.0 run

Branch `leftovers-1.1.0` from `main`. Standalone PR against `main`, merged before both stacks start. Decisions: map ticket 12 (leftovers) and ticket 27 (the canonical helper, item 3 below). Read `README.md` in this directory first.

## Commits, in this order

1. `docs: the lr+ knobs are Python-API only`
   - `README.md` (the `dmc` row of the CLI table), `CLAUDE.md` (the neighbour_combine paragraph) and `docs/history/planning/codebase/CONCERNS.md` (the two "deferred to 1.1" statements) say the `lr+` knobs are Python-API only, with no version promise. Grep `rg -n "deferred to 1\.1|in 1\.1\b" README.md CLAUDE.md docs/` and fix every remaining forward-looking statement; historical changelog text stays.
2. `docs(history): add the 2026-09-05 field guide`
   - The map driver supplies the HTML file; it is not in the repository. Commit it as `docs/history/2026-09-05-field-guide.html`. Add one line to `docs/history/README.md` describing it as a dated onboarding presentation. Add `2026-09-05-field-guide.html` to `exclude_docs` in `mkdocs.yml` so the site does not ship it; `mkdocs build --strict` must stay green.
3. `feat(chroms): canonical chromosome helper`
   - Re-type `src/epykit/_chroms.py` from `origin/feat/canonical-chrom-filter` commit `cd9f89b`: `is_canonical_chrom`, `filter_canonical`, `filter_canonical_logged`, `CANONICAL_CHROM_CORES`, `CANONICAL_CHROMS_UCSC`. Accepts chr1 to chr22, X, Y, M and MT under UCSC and Ensembl naming. Fix it under the current lint set (the branch used `typing.Iterable`). Carry `tests/test_chroms.py` from the branch. Nothing calls the helper yet; both stacks will. Say so in the commit body.

## Not a commit

After the PR is open, delete the two remote branches the triage report drops, first confirming their tips: `fine-tune` at `47d7878` and `optimize/autonomous-v1` at `99b2b4f`. `git push origin --delete fine-tune optimize/autonomous-v1`. Record the two SHAs in the PR body. Leave `bench-tune` and `feat/canonical-chrom-filter` alone.

## Contract

No behaviour change. The helper is unused. Regen hashes unchanged.

## Deliver

PR title: `Leftovers from 1.1.0: lr+ wording, the field guide, and the canonical-chromosome helper`. Body per the run README, plus the deleted branch tips. Then `worker_done`.
