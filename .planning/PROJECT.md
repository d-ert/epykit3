# epykit — Publication-Hardening Campaign

## What This Is

A thorough, adversarial LLM-driven review-and-fix campaign to get the `epykit`
WGBS methylation analysis package — and its target manuscript
`benchmark/paper/paper.md` — over the bar for a peer-reviewed software/methods
paper. Claude (plus subagents/workflows, run directly in this repo) reviews the
package across four dimensions, surfaces weak points, and **fixes them** until
the package and paper are publication-ready. The audience is the manuscript's
reviewers and the WGBS/bioinformatics community who would adopt the tool.

## Core Value

When a hostile peer reviewer attacks the package or the paper, every weak point
has already been found and either fixed or defensibly documented — nothing real
is left for them to discover.

## Requirements

### Validated

<!-- Existing epykit capabilities — the artifact under review, inferred from the codebase map. -->

- ✓ Python-native WGBS pipeline: Bismark/MethylDackel `.cov` → DMC/DMR → annotation → HTML report — existing
- ✓ Streaming Polars/Parquet methylstore; peak memory O(largest chromosome), never whole-genome-in-RAM — existing
- ✓ Four per-CpG DMC engines (`lr`, `welch_t`, `fisher`, `glm`) + streaming `DMCStore` — existing
- ✓ Four DMR callers (chain_merge, tile, sliding-window, HMM) + permutation empirical FDR — existing
- ✓ Scanpy-style `pp`/`tl`/`pl` API mirrored by a CLI; interop sinks (AnnData, MuData, methylKit, MultiQC) — existing
- ✓ Head-to-head benchmark vs 8 published tools on Piao 2021 simulated + GSE263850 real cohort — existing
- ✓ CI on `{ubuntu, windows} × {py3.9, py3.12}`; Windows compatibility is load-bearing — existing

### Active

<!-- The campaign. Each is a hypothesis until the work ships and the "done" bar is met. -->

**Statistical correctness**
- [ ] Adversarial review of all DMC engines, DMR callers, FDR/dispersion logic; every HIGH/CRITICAL statistical-validity issue fixed or defensibly documented
- [ ] `lr+` power stack re-litigated from scratch (ship / reframe / drop) with evidence — treat the "research knob" framing as unproven
- [ ] Every statistical claim in `paper.md` traced to re-runnable evidence (no orphan or unverifiable claims)

**Benchmark rigor & reproducibility**
- [ ] Comparison fairness audited end-to-end (e.g. methylKit mode parity); any unfair comparison corrected or disclosed
- [ ] `claims.yaml` ↔ committed-results chain verified; stale/broken claim assertions fixed (e.g. references to non-existent `timings_post_phase3.parquet`)
- [ ] Benchmark reproducible from committed data (φ-sweep committed; repro protocol — README/Docker/renv — validated)

**Code quality & API**
- [ ] Correctness-bug sweep across core modules (`dmc.py`, `dmr.py`, `tl.py`, `_glm.py`, …); HIGH/CRITICAL bugs fixed
- [ ] API surface, deprecation shims, and 1.0 stability story reviewed for a clean public contract
- [ ] Lint / type / cross-platform (Windows) / CI hygiene green

**Docs & usability**
- [ ] Public API docs, tutorials, and examples audited; a new user can onboard end-to-end from documentation
- [ ] Manuscript prose hardened for defensibility and clarity

**Cross-cutting (defines "done")**
- [ ] Hostile-reviewer objection register: every plausible objection has a documented answer or mitigation
- [ ] No open HIGH/CRITICAL findings across all four dimensions + green CI
- [ ] Venue decision made (or explicitly deferred with documented criteria)

### Out of Scope

- Full benchmark re-execution from raw data — benchmark numbers are treated as mostly frozen; review targets rigor, provenance, and claim traceability, not re-running 48–72h studies
- A bug-fix changelog in the manuscript — fixes live in the repo's git history, not in `paper.md` (software-package papers present the current stable version)
- New 1.1 analysis features (CLI flags for `lr+` knobs, JAX GPU backend, report decimation) — unless one is a genuine publication blocker
- Standing up a separate/re-runnable review tool or harness — the review is run here with Claude + subagents

## Context

- **Brownfield, mature.** `epykit` is at 1.0 (stable API, MIT). The codebase has been mapped — see `.planning/codebase/` (ARCHITECTURE, STACK, CONCERNS, CONVENTIONS, INTEGRATIONS, STRUCTURE, TESTING).
- **`CONCERNS.md` is a pre-built weak-points inventory** and a primary input: it already flags the `lr+` FPR-drift risk, benchmark reproducibility gaps (M1/M5 outstanding, φ-sweep computed-but-uncommitted, stale timing claims in `claims.yaml`), test-coverage gaps (`lr+` on real data, permutation FDR, GLM contrasts), and architectural debt.
- **Known landmine — `lr+` power stack.** On GSE263850 at q=0.05, `power_stack="lr+"` inflates DMC calls ~13× vs bare `lr` (FPR drift; tuned on an underdispersed φ≈0.4 simulator vs real WGBS φ≈1.5–5). Currently positioned as a "research knob," default off, paper claims around bare `lr`. The campaign re-litigates whether that framing survives a hostile reviewer.
- **Publication history signal.** Recent git history (`gb-resubmission-scaffolding`) suggests a prior/ongoing Genome Biology track, but venue is not committed.
- **Working files in flight.** `benchmark/paper/paper.md` and `benchmark/scripts/claims.yaml` are currently modified; uncommitted φ-sweep export exists at `benchmark/phi_sweep_export_2026-06-05/`.

## Constraints

- **Reproducibility**: Raw benchmark data isn't bundled; full studies take 48–72h — Benchmark numbers are mostly frozen; verify on-paper rigor and committed-data reproducibility, not by re-execution.
- **Platform**: Windows compatibility is load-bearing — `pysam`/`pyBigWig`-based extras are Linux/macOS only and gated; fixes must not break the Windows CI leg.
- **Process**: Mechanical fixes executed inline, not over-orchestrated — reserve subagent/workflow fan-out for genuine review judgment and context-heavy analysis.
- **Manuscript hygiene**: No bug-fix history in the paper — the manuscript presents the current stable version; remediation is recorded in the repo only.
- **Streaming contract**: DMC/DMR fixes must preserve O(largest-chromosome) memory — no materializing the whole-genome frame.

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Find **and** fix (not review-only) | Goal is a publication-ready artifact, not a critique document | — Pending |
| Re-litigate everything (adversarial) | A thorough review must attack even "settled" decisions, incl. whether `lr+` ships | — Pending |
| Benchmark numbers mostly frozen | No bundled raw data; 48–72h re-runs; review rigor/provenance instead | — Pending |
| Review runs here (Claude + subagents) | No appetite for separate review tooling; capture methodology only if cheap | — Pending |
| `paper.md` is the target manuscript | Harden the existing benchmark/methods paper, not a new artifact | — Pending |
| Done = findings-clean+traceable **and** reviewer-proof | Both an objective bar (no HIGH/CRITICAL, claims traceable, green CI) and a defensibility bar | — Pending |
| Venue undecided | Keep options open (GB resubmission vs fresh methods journal vs JOSS) until review clarifies the strongest framing | — Pending |

## Evolution

This document evolves at phase transitions and milestone boundaries.

**After each phase transition** (via `/gsd-transition`):
1. Requirements invalidated? → Move to Out of Scope with reason
2. Requirements validated? → Move to Validated with phase reference
3. New requirements emerged? → Add to Active
4. Decisions to log? → Add to Key Decisions
5. "What This Is" still accurate? → Update if drifted

**After each milestone** (via `/gsd-complete-milestone`):
1. Full review of all sections
2. Core Value check — still the right priority?
3. Audit Out of Scope — reasons still valid?
4. Update Context with current state

---
*Last updated: 2026-06-06 after initialization*
