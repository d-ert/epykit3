# The 1.2.0 implementation plans

Reviewed on 2026-09-06 against main commit `900ea6cfc22f521b1a708ae43c664c09c563f318`.
GitHub had no open issues and 18 closed planning issues at review time.
[Issue 11](https://github.com/d-ert/epykit3/issues/11) records the planning map.
A closed planning issue means its planning work ended. It does not mean its implementation shipped.

These reviewed briefs replace conflicting instructions in the original issue resolutions and prototype notes.
[The review record](REVIEW.md) explains each correction and its evidence.
The branch `planning/1.2.0-run` stores plans and their inputs. Do not merge that branch into main.

## Start in dependency order

Each brief is ready for implementation when its listed start condition holds.
The table gives branch bases for a stacked review. The maintainer controls merges.

| Layer | Brief | Branch | Base | Start condition |
|---|---|---|---|---|
| L | [Leftovers](L-leftovers.md) | `leftovers-1.1.0` | `main` | Can start now |
| S1 | [Region FDR](S1-region-fdr.md) | `salvage-1-region-fdr` | `main` after L | L merged |
| S2 | [Canonical filter](S2-canonical-filter.md) | `salvage-2-canonical-filter` | `salvage-1-region-fdr` | S1 implementation and gates complete |
| S3 | [Preset documentation](S3-bench-tune-additive.md) | `salvage-3-bench-tune-additive` | `salvage-2-canonical-filter` | S2 implementation and gates complete |
| R1 | [DMC stages](R1-stage-split.md) | `refactor-1-stage-split` | `main` after L | L merged |
| R2 | [Engine registry](R2-engine-registry.md) | `refactor-2-engine-registry` | `refactor-1-stage-split` | R1 implementation and gates complete |
| R3 | [Engine runners](R3-engine-runners.md) | `refactor-3-engine-runners` | `refactor-2-engine-registry` | R2 implementation and gates complete |
| R4 | [CLI package](R4-cli-package.md) | `refactor-4-cli-package` | `refactor-3-engine-runners` | R3 implementation and gates complete |
| R5 | [Compatibility guidance](R5-removals.md) | `refactor-5-compatibility` | `refactor-4-cli-package` | S3 merged and refactor stack updated onto main |
| R6 | [CLI options and DMC filter](R6-cli-flags.md) | `refactor-6-cli-flags` | `refactor-5-compatibility` | R5 complete on the integrated base |
| C4 | [ASM anchors](C4-asm-anchors.md) | `fix-asm-anchor-classes` | `main` after R1 | R1 merged so CI runs the BAM tests |
| REL | [Release preparation](REL-1.2.0.md) | `release-1.2` | `main` | All implementation PRs merged |

L runs first. S1 through S3 and R1 through R4 can then run in parallel, sequentially within each stack.
Merge the salvage stack first. Update the refactor stack onto that main, then implement R5 and R6.
C4 can run alongside the remaining refactor work after R1 merges. REL runs last.

After a squash merge, replay only each remaining layer's own commits onto its new base.
Record full old and new base and tip SHAs. Compare the old and new patch series with `git range-diff`.
Inspect the diff from the new base and run the affected gates.
An old-tip versus new-tip tree diff normally includes changes from the new base. It is not an identity proof.
Retarget the next PR before deleting its former base branch.

## Keep ownership explicit

| Owner | Files and responsibility |
|---|---|
| L | `README.md`, `CLAUDE.md`, the archived field guide and its index, `src/epykit/_chroms.py`, `tests/test_chroms.py` |
| S1 to S3 | `src/epykit/dmr.py`, `tl.dmr` in `src/epykit/tl.py`, `convert.py`, `io.py`, `pl/_compute.py`, their tests, DMR and ingestion docs |
| R1 to R4 | `tl.dmc`, `_run_dmc_contrast`, DMC-only helpers, `_dmc_config.py`, `_dmc_stages.py`, `_dmc_engines.py`, `dmc.py`, `cli.py` and `cli/`, their tests, `.github/workflows/test.yml`, the complexity setting in `pyproject.toml` |
| R5 | Compatibility docs and deprecation wording listed in its brief, after S3 and R4 |
| R6 | DMC config, stages, engine and cache tests, CLI package, DMC and CLI docs, README |
| C4 | `src/epykit/asm.py`, `tests/test_asm.py`, `docs/analysis/asm.md` |
| REL | Project version metadata, the matching root package entry in `uv.lock`, README version text, CHANGELOG and current-version docs |

The salvage and refactor workers share `tl.py` by function. Preserve the other worker's edits.
Shared imports can conflict. Resolve them on integration rather than copying an older file.
Leave helpers used by both DMC and DMR in `tl.py` during R1. Import them locally from stages when needed.
Do not create a module import cycle or move DMR's helpers without coordinating with its owner.

Every implementation layer may add its own entry under `[Unreleased]` in CHANGELOG.
Keep those entries separate and retain both entries during integration. C4 follows the same rule.
S1 to S3 own shared architecture and DMR documentation during parallel work.
R5 and R6 update shared README and architecture text only after S3 merges.

## Preserve the approved scope

- Preserve defaults and supported APIs in 1.2. The proposed removals require a separate maintainer decision.
- Keep `canonical_only=False`, empirical DMR `fdr_method="max_t"`, and the current DMR preset values.
- Add DSS-style smoothing flags only to DMC. There is no existing DMR smoothing API to expose.
- Preserve numerical output on unchanged paths. C4 intentionally changes affected ASM output. New FDR and filtering behavior requires explicit opt-in.
- Do not delete old source branches as part of this run. Their full commits remain useful review inputs.
- Read `CLAUDE.md`, `CONTRIBUTING.md`, and any local instructions before editing. Use the personal GitHub identity for this repository.
- Snapshot Git status before edits. Preserve unrelated files, including the untracked field guide in the main checkout.
- Use `uv sync --locked` and `uv run --frozen`. Only REL may update the root project's version in the lockfile. Dependency versions stay pinned.
- Keep one concern per commit. Use `type(scope): summary`. Do not add AI attribution.

## Verify the completed layer

Run focused checks while editing. Run the relevant full gates on each finished code layer and after integration changes its behavior.
Do not repeat the entire suite after each intermediate commit.

```bash
uv sync --locked --group dev --extra all
uv lock --check
uv run --frozen ruff check src/ tests/
uv run --frozen ruff format --check src/
uv run --frozen mypy src/epykit
uv run --frozen pytest -m "not slow" --strict-markers -ra
uv run --frozen pytest -m slow --strict-markers -ra
uv run --frozen python benchmark/scripts/regen_small.py
uv run --locked --only-group docs mkdocs build --strict
```

On Linux or macOS, R1 and C4 add `--extra bam` to the sync command.
Check that ASM, BAM ingestion, and entropy tests execute instead of skipping.
The engine gate must report `OK: all 2 hashes match`. Do not update the reference hashes in this run.
Those hashes cover selected LR outputs only. Each brief also names checks for behavior outside that slice.

For a documentation-only layer, check links, names, defaults, placeholders, `git diff --check`, and the strict docs build.
Full numerical tests are not required for text-only changes.
A failing baseline check must be recorded with its command and output. Do not report an unrun check as passed.

## Deliver a reviewable PR

Open each stacked PR against the branch below it. Standalone and bottom-layer PRs target main.
Describe the problem, the behavior after the change, and the actual validation.
Include a limitation only when it affects review. No fixed number of paragraphs or tables is required.

Return the PR URL, full pushed HEAD SHA, changed files, test results, and any dependency still waiting to merge.
Use the active session's handoff mechanism. A standalone worker does not need an Orca coordinator or a `worker_done` event.
A brief authorizes its implementation when dispatched. It does not authorize merges, branch deletion, tags, or releases.

After REL merges, the maintainer can publish 1.2.0 once CI for that exact main commit passes.
See [REL](REL-1.2.0.md) for version and lockfile checks.
