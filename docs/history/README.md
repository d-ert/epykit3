# Project history

Frozen planning and review artefacts from the 2026 development campaigns.
They are kept for provenance and are excluded from the published docs site
(`exclude_docs` in `mkdocs.yml`). Nothing here is maintained; the current
architecture reference is [`../advanced/architecture.md`](../advanced/architecture.md)
and the current debt list is tracked in issues.

| Path | Contents |
|---|---|
| `2026-09-05-field-guide.html` | The epykit Field Guide, a self-contained HTML onboarding guide to the whole system, read on 2026-09-05 against commit `1f501e0`: the domain, the pipeline, the data model, the quality gates, and a sequenced cleanup plan. It was the input to the 1.2 planning run. Open it in a browser; it is a snapshot and is not maintained. |
| `planning/` | The publication-hardening campaign charter (`PROJECT.md`) and a codebase map dated 2026-06-06 (architecture, concerns, conventions, integrations, stack, structure, testing). `codebase/CONCERNS.md` is the most useful file: the maintainers' own itemised debt list with file pointers. |
| `superpowers/specs/` | Design specs, 2026-05-27 to 2026-06-06: benchmark design, engine freeze, 1.0 and paper finish, public-surface audit, CSV output, review fixes, report redesign. |
| `superpowers/plans/` | Dated implementation plans that executed those specs. |

The peer review that drove the June 2026 fixes, and its remediation summary,
live in [`../review/`](../review/) because they are still referenced from the
changelog.
