# REL: cut 1.2.0

Branch `release-1.2` from `main` after every other PR in this run has merged. Standalone PR against `main`, the last one. Model: the 1.1.0 release PR (#7) and its changelog section.

You own `pyproject.toml` (version line only), `CHANGELOG.md`, `README.md` (badge and BibTeX), and `docs/` statements of the current version.

## Commits, in this order

1. `build: bump version to 1.2.0`
   - `version = "1.2.0"`; README badge and BibTeX; grep `docs/` for `1.1.0` as the current version. `uv lock --check` must stay green (the lock records the project by path).
2. `docs(changelog): cut 1.2.0`
   - Insert `## [1.2.0] — the merge date` under an empty `## [Unreleased]` and move the body under it. Add an upgrade paragraph, four to six sentences: the four removals and their replacements; the ASM fix and which inputs it affects; the new opt-in features (`fdr_method="region"`, chain_merge empirical FDR, `canonical_only`, the CLI flags); the pysam tests now running in CI; the module layout changes (`_dmc_stages`, `_dmc_engines`, the `cli` package) that only matter to anyone importing private modules.

## Contract

No behaviour change. All gates green.

## Deliver

PR title: `Prepare the 1.2.0 release`. Then `worker_done`. After the squash merge the map driver tags `v1.2.0`, waits for the full-matrix main CI, and creates the GitHub release from the changelog section.
