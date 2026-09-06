# The 1.2.0 dev run

Planned on the wayfinder map, d-ert/epykit3 issue 11. Every decision below links to the ticket that holds its detail; this document only fixes the order and the shared rules. Each brief in this directory is dispatch-ready: hand it to one worker in its own worktree.

## Shape

Two stacks plus three standalone PRs. Squash merges, bottom to top.

| Order | Brief | Branch | Base | Stack |
|---|---|---|---|---|
| 1 | `L-leftovers.md` | `leftovers-1.1.0` | `main` | standalone, merges first |
| 2 | `S1-region-fdr.md` | `salvage-1-region-fdr` | `main` (after L) | salvage, bottom |
| 3 | `S2-canonical-filter.md` | `salvage-2-canonical-filter` | `salvage-1-region-fdr` | salvage |
| 4 | `S3-bench-tune-additive.md` | `salvage-3-bench-tune-additive` | `salvage-2-canonical-filter` | salvage, top |
| 5 | `R1-stage-split.md` | `refactor-1-stage-split` | `main` (after L) | refactor, bottom |
| 6 | `R2-engine-registry.md` | `refactor-2-engine-registry` | `refactor-1-stage-split` | refactor |
| 7 | `R3-engine-runners.md` | `refactor-3-engine-runners` | `refactor-2-engine-registry` | refactor |
| 8 | `R4-cli-package.md` | `refactor-4-cli-package` | `refactor-3-engine-runners` | refactor |
| 9 | `R5-removals.md` | `refactor-5-removals` | `refactor-4-cli-package` | refactor |
| 10 | `R6-cli-flags.md` | `refactor-6-cli-flags` | `refactor-5-removals` | refactor, top; written after the salvage stack merged |
| 11 | `C4-asm-anchors.md` | `fix-asm-anchor-classes` | `main` (after L) | standalone, any time |
| 12 | `REL-1.2.0.md` | `release-1.2` | `main` (after everything) | standalone, last |

Parallelism: L runs alone and merges first. Then the salvage stack (S1 to S3) and the refactor stack (R1 to R5) run in parallel, one worktree each, sequential inside a stack. C4 runs whenever a Linux or macOS worker is free. R6 starts only after S3 has merged and the refactor stack has been rebased onto `main`. REL is last.

Merge order between the stacks: the salvage stack merges first (it is short). Then rebase the refactor stack onto `main`, prove each layer's tree is unchanged (`git diff --stat <old-tip> <new-tip>` is empty apart from the rebase), and merge it. Then C4 if not already merged, then REL, then the tag.

## Ownership

| Owner | Files |
|---|---|
| salvage stack | `src/epykit/dmr.py`, the `tl.dmr` function in `src/epykit/tl.py`, `src/epykit/convert.py`, `src/epykit/io.py`, `src/epykit/pl/_compute.py`, their tests, `docs/` pages for dmr, read-bismark and architecture |
| refactor stack | `tl.dmc` and `_run_dmc_contrast` in `src/epykit/tl.py`, `src/epykit/_dmc_config.py`, `src/epykit/_dmc_stages.py`, `src/epykit/_dmc_engines.py`, `src/epykit/dmc.py`, `src/epykit/cli.py` and the `cli/` package, `.github/workflows/test.yml`, `pyproject.toml` lint settings, their tests |
| C4 | `src/epykit/asm.py`, `tests/test_asm.py` |
| L | `README.md`, `CLAUDE.md`, `docs/history/`, `mkdocs.yml`, `src/epykit/_chroms.py`, `tests/test_chroms.py` |
| REL | `pyproject.toml` version line, `CHANGELOG.md`, `README.md` badge and BibTeX |

An edit outside your files goes to the owning stack as a hand-off, named in `worker_done`.

## Rules for every brief

- No results change except where a brief says so, and then only on opt-in paths. `uv run --frozen python benchmark/scripts/regen_small.py` must print `OK: all 2 hashes match` after every commit; if it does not, stop and report failed, never re-snapshot.
- Everything lands opt-in. Defaults do not move. The four maintainer questions on the map are closed as out of scope; do not reopen them by changing a default.
- Use `uv run --frozen` and `uv sync --locked`; never commit a `uv.lock` change. If the rtk shell wrapper refuses `git status`, `git log` or `git diff`, call `/usr/bin/git`.
- Commit messages in the repo's `type(scope): summary` style, one concern per commit. No AI attribution anywhere, in commits or PRs.
- Read `CLAUDE.md` and `CONTRIBUTING.md` before the first edit.

## Gates, all must pass after each commit

```bash
uv sync --locked --group dev --extra all            # add --extra bam on Linux or macOS for C4 and R1
uv lock --check
uv run --frozen ruff check src/ tests/
uv run --frozen ruff format --check src/
uv run --frozen mypy src/epykit
uv run --frozen pytest -m "not slow" --strict-markers -q
uv run --frozen pytest -m slow --strict-markers -q
uv run --frozen python benchmark/scripts/regen_small.py     # must print: OK: all 2 hashes match
uv run --locked --only-group docs mkdocs build --strict
```

## PR conventions

Title as given in the brief. Body in plain prose with sentence-case headings, no em dashes, no bold-label bullet lists: what and why in three or four sentences; the behaviour contract (what can change, what cannot); a table of commits; a table of gate results with numbers; a section "deliberately not here". Open stacked PRs against the branch below (`gh pr create --base <branch>`); the bottom of each stack and the standalone PRs go against `main`.

## Worker protocol

Send `worker_done` with `--outcome succeeded`, the PR URL, the pushed HEAD sha and the gate numbers. If a gate cannot go green without a results change, do not push; send `worker_done --outcome failed` with details. Use `orca orchestration ask` for anything the brief does not cover; the coordinator decides.

## After the last merge

The map driver tags `v1.2.0` on the squash of REL, waits for the full-matrix main CI, and creates the GitHub release from the `## [1.2.0]` section of the changelog.
