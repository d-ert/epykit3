# Default canonical-chromosome filtering for DMC/DMR

- **Date:** 2026-06-08
- **Status:** Approved (design); pending spec review
- **Scope:** Library only (benchmark/paper regeneration is a separate follow-up)

## Context

epykit's annotation step warns when a large share of site chromosomes are
absent from the gene model. Investigation traced this to **non-canonical
contigs** (`*_random`, `chrUn_*`, `GL000*`, `KI270*`) flowing through the whole
pipeline. They are not mis-annotated at random — unmatched sites get a
deterministic `intergenic`/`open_sea` fallback ([annotate.py:609](../../../src/epykit/annotate.py),
[annotate.py:616](../../../src/epykit/annotate.py)) and per-chromosome
annotation is isolated, so canonical-chromosome signal is never corrupted.
The real problem is upstream: scaffold CpGs are **tested for differential
methylation at all**, despite poor mappability on unplaced/alt contigs.

Evidence from the benchmark real cohort (GSE263850):

- `benchmark/data/study3/dmr_significant_lenient.csv`: **34 of 3,433** lenient
  DMRs (~1%) sit on scaffolds.
- `chrUn_KI270742v1` is the **single most significant DMR** in that set — a
  mappability red flag presented as the top real-cohort hit.

The simulator (Piao) benchmark is canonical-only, so headline TPR/FPR/F1 are
unaffected by filtering; only the real-data section would shift, and shifting
it would *improve* defensibility (drops a likely artifact).

## Decisions (locked)

1. **Placement:** `canonical_only=True` **default** on the DMC + DMR analysis
   paths, plus an **optional opt-in** (`canonical_only=False` default) at
   ingestion.
2. **Canonical definition includes the mitochondrial contig** (`chrM` / `MT`).
3. **This branch is library-only.** Benchmark/paper regeneration is tracked
   separately.

## Design

### 1. Single source of truth — `src/epykit/_chroms.py` (new)

Dependency-free module so analysis, ingestion, and plotting share one
definition:

```python
CANONICAL_CORES = {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}

def is_canonical_chrom(name: str) -> bool:
    """True for the main mammalian assembly chromosomes under either the
    UCSC (`chr1`, `chrM`) or Ensembl (`1`, `MT`) naming convention."""
    core = name[3:] if name[:3].lower() == "chr" else name
    return core.upper() in CANONICAL_CORES

def filter_canonical(chroms: Iterable[str]) -> list[str]:
    """Order-preserving filter to canonical chromosomes."""
    return [c for c in chroms if is_canonical_chrom(c)]
```

Predicate behaviour (kept → dropped):

| Input | Result | Why |
|-------|--------|-----|
| `chr1`, `chr22`, `chrX`, `chrY`, `chrM` | keep | UCSC canonical |
| `1`, `22`, `X`, `Y`, `MT` | keep | Ensembl canonical |
| `chr14_KI270722v1_random` | drop | core `14_KI270722v1_random` ∉ set |
| `chrUn_KI270742v1` | drop | core `Un_KI270742v1` ∉ set |
| `GL000216v2`, `KI270722.1` | drop | not a canonical core |
| `chr23`, `chr0` | drop | no such canonical chromosome |

No denylist — the allowlist rejects every scaffold pattern intrinsically.

The inline canonical list at [pl/_compute.py:614](../../../src/epykit/pl/_compute.py)
is replaced by an import from this module, so the plotter and the analysis
path can never drift apart.

### 2. DMC + DMR: `canonical_only: bool = True` (the behavioural change)

Applied **only when `chromosomes` is left `None`** (auto-detect). Each
*independent* enumeration point routes its detected list through
`filter_canonical`:

- **`process_chromosomes_dmc`** ([dmc.py:2160](../../../src/epykit/dmc.py)) —
  covers `tl.dmc` and the CLI `dmc` subcommand. Filter applied right after
  `_detect_chromosomes()` ([dmc.py:1169](../../../src/epykit/dmc.py)).
- **`call_dmr_tile_based`** ([dmr.py:1194](../../../src/epykit/dmr.py)) — has
  its own `chrom=*` glob ([dmr.py:1291](../../../src/epykit/dmr.py)); filter
  applied there.
- **`tl.dmc`** ([tl.py:399](../../../src/epykit/tl.py)) and **`tl.dmr`**
  ([tl.py:1208](../../../src/epykit/tl.py)) gain a `canonical_only=True`
  passthrough.

**Inherited (no direct change needed):** `call_dmr_chain_merge` (default DMR
caller) and `call_dmr_sliding_window` consume an already-built
`DMCStore`/table ([dmr.py:739](../../../src/epykit/dmr.py),
[dmr.py:386](../../../src/epykit/dmr.py)). If the upstream DMC ran
canonical-only, their output is canonical automatically. Consequently
`tl.dmr`'s `canonical_only` is **meaningful only for `method="tile"`** — for
`chain_merge`/`sliding_window` it is inherited from the DMC table. This is
documented and consistent with the codebase's existing per-method parameter
scoping (e.g. `empirical_fdr` is tile-only).

### 3. Ingestion opt-in: `read_bismark(..., canonical_only: bool = False)`

Default **off** — current ingestion is unchanged. When `True`, non-canonical
chromosome partitions are **not written** to the methylstore, so every
downstream step (QC, smoothing, DVC, DMC, DMR) is canonical and no compute is
spent on scaffolds. Injected at the partition-writing path in the convert core
([io.py:144](../../../src/epykit/io.py) and the shared convert helper). If
`read_methyldackel` / other ingestion entry points share that convert core,
they gain the same parameter; otherwise this branch wires `read_bismark` only
and leaves a note.

### 4. Precedence (escape hatches)

| Situation | Behaviour |
|-----------|-----------|
| `chromosomes=[...]` passed explicitly | Used **verbatim**; `canonical_only` ignored (you asked for exactly those) |
| `canonical_only=False` | All detected contigs (pre-change behaviour) |
| `canonical_only=True` (default), `chromosomes=None` | Canonical subset of detected contigs |
| CLI `epykit dmc/dmr --all-contigs` | Sets `canonical_only=False` |
| CLI `epykit convert --canonical-only` | Opt into ingestion-time filtering |

### 5. Audit trail

When the filter drops ≥1 contig, emit a single `logger.info` (never `print` —
the library logging convention):

> `canonical_only: testing N chromosomes; dropped M non-canonical contigs (chrUn_KI270742v1, …). Pass canonical_only=False (or --all-contigs) to include them.`

This also defuses the non-mammalian footgun (a yeast/plant user sees
`chrI, chrII…` being dropped and knows to flip the flag). `canonical_only` is
recorded in the DMC manifest (`.epykit_dmc_manifest.json`) for provenance.

## Canonical definition

`chr1`–`chr22`, `chrX`, `chrY`, `chrM` (UCSC) and the bare `1`–`22`, `X`, `Y`,
`MT` (Ensembl). Mitochondrial contig **kept** (decision 2), matching the
existing plotter list. Note: the GSE263850 real-data DMR set contains no
`chrM` rows, so the chrM choice does not affect the benchmark.

## Behavioural change & versioning

This **changes default output**: `tl.dmc(md)` / `tl.dmr(md)` will return fewer
calls than before (scaffold calls removed). The **API stays
backward-compatible** (only new keyword args with defaults are added), but the
*semantics* change. Mitigations: prominent CHANGELOG **"Changed"** entry, the
per-run audit log, and docs. Version bump is the maintainer's call; flagged
here, not decided.

## Risks

- **Mammalian assumption.** The canonical set is human/mouse numeric naming.
  Non-mammalian assemblies (roman-numeral or named chromosomes) would be
  over-dropped. Mitigation: audit log makes it visible + `canonical_only=False`
  restores all contigs. No organism auto-detection (out of scope, YAGNI).
- **Silent surprise for existing users.** Addressed by the audit log + CHANGELOG.

## Testing plan

CI-tier (`not slow`), `--strict-markers`, Windows-safe.

- **Unit (`_chroms`):** both naming conventions kept; every scaffold pattern
  and `chr23`/`chr0` dropped; order preserved.
- **DMC integration:** tiny store with `chr1` + one scaffold →
  - default → scaffold absent from DMC results;
  - `canonical_only=False` → scaffold present;
  - explicit `chromosomes=[scaffold]` → scaffold present (override wins).
- **DMR tile integration:** same three cases via `tl.dmr(method="tile")`.
- **DMR inheritance:** `method="chain_merge"` after a default `tl.dmc` →
  output has no scaffold DMRs (inherited).
- **Ingestion opt-in:** `read_bismark(canonical_only=True)` omits the scaffold
  partition; default keeps it.
- **CLI:** `--all-contigs` keeps scaffolds; `convert --canonical-only` drops
  them at ingest.
- **Audit log:** assert the `logger.info` fires with the dropped-contig count
  (surfaced via the existing `filterwarnings`/caplog pattern).

## Docs / CHANGELOG

- CHANGELOG **"Changed"** entry describing the new default + escape hatches.
- Update the DMC/DMR analysis docs (and a note in `docs/analysis/annotate.md`
  explaining the annotation warning is now mostly avoided upstream).

## Out of scope (this branch)

- Benchmark/paper regeneration of the GSE263850 real-data section (separate
  tracked follow-up).
- Site-weighted rewrite of the annotation chromosome-overlap warning.
- Extending `canonical_only` to DVC / entropy / ASM (trivial later via the
  shared helper).

## Acceptance criteria

1. `ep.tl.dmc(md)` and `ep.tl.dmr(md)` exclude non-canonical contigs by default.
2. `canonical_only=False`, `--all-contigs`, and explicit `chromosomes=`
   each restore/override as specified in the precedence table.
3. `read_bismark(canonical_only=True)` filters at ingest; default unchanged.
4. One shared canonical definition; plotter imports it.
5. Audit log fires when contigs are dropped.
6. All new tests pass on the CI matrix; CHANGELOG + docs updated.
