# REL: prepare the 1.2.0 release

Start from current main after L, S1, S2, any needed S3 change, R1 through R6, and C4 have merged.
Use `release-1.2` and target main.
Read the [run rules](README.md). Verify merged PRs and main CI live before preparing the release.

Own `pyproject.toml`, the root epykit package entry in `uv.lock`, CHANGELOG, README version text, and current-version docs.

## Implement

1. Change the project version to 1.2.0 and update current README badge and citation version text. Leave historical version references intact.
2. Use the uv version pinned in CI to refresh the lockfile after the version change. The lock records `epykit = 1.1.0` at the reviewed baseline, even though its source is editable. A pyproject-only bump leaves the lock stale.
3. Inspect the lock diff. Accept only the root epykit package version change and metadata strictly required by that change. Do not upgrade dependencies or rewrite unrelated lock entries.
4. Run `uv lock --check` and `uv sync --locked --group dev --extra all`. Verify the installed distribution reports 1.2.0.
5. Move the completed Unreleased entries into a new 1.2.0 section. Keep an empty Unreleased section. Use the actual preparation date; update the date before merge if the release date changes.
6. Write upgrade guidance for what actually merged. State that deprecated APIs remain in 1.2. Explain the ASM anchor policy, opt-in region FDR and canonical filtering, DMC-only smoothing flags, and Ubuntu BAM-test coverage.
7. Include the module moves only as a note for consumers of private imports. Do not claim new DMR smoothing support or any unmerged feature.

## Accept when

- Only the intended package version changes in the lockfile. All dependency versions remain pinned.
- Installed package metadata, pyproject, README, and the changelog agree.
- All code-layer gates and the strict docs build pass on the final release commit.
- The changelog contains no incomplete entries or promises about unmerged work.

PR title: `Prepare the 1.2.0 release`.

## Publish after review

The maintainer merges REL and verifies full-matrix CI for that exact main commit.
Before a separate tag or release action, check whether `v1.2.0` and the GitHub release already exist.
The tag must point to the approved release commit.
Preparing this PR does not authorize its merge, a tag, a package upload, or a GitHub release.
