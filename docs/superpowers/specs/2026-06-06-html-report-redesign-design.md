# HTML Report Redesign — Design

**Date:** 2026-06-06
**Status:** Approved (brainstorming complete)
**Scope:** Redesign `epykit`'s single-file HTML report into a MultiQC-style dashboard with new analytical content, while preserving the dependency budget, Windows compatibility, and graceful-skip robustness of the current implementation.

---

## 1. Goal

Make `md.report()` / `epykit report` produce a report that is **more informative, better looking, more professional, more scientific, and more useful** than the current single-column card stack.

The current report (`src/epykit/report.py`, `templates/report.html.j2`, `templates/report.css`, Plotly twins in `pl/_plotly.py`) is clean but generic: one scrolling column, light-only, no navigation, raw-path-overflowing tables, plain KPIs, CDN-only Plotly, a raw-JSON provenance dump, and no scientific framing.

## 2. Decisions (locked with user)

| Decision | Choice |
|---|---|
| Scope | Full redesign — presentation overhaul **plus** new analytical content. |
| Offline | `self_contained` flag on `generate_report`, **default `True`** (embed Plotly inline); `False` keeps CDN. |
| Aesthetic | **MultiQC-style dashboard** — fixed sidebar TOC with per-section status dots, scannable main panel. |
| Palette | **Refine + stay consistent**: keep `hyper`/`hypo`/feature/CpG-context semantics from `_style.PALETTE`, modernize dashboard accents (indigo accent, donut colors). **Matplotlib plots untouched.** |
| Robustness | **Graceful skips** preserved; every section conditional on its data; partial pipelines render a "not run" state + hint and a hollow sidebar dot. Add tests for partial-pipeline rendering. |

## 3. Non-goals (YAGNI)

- No JS framework, no build step, no new runtime dependencies (still `jinja2` + `plotly` only).
- No changes to the matplotlib plotters in `pl/*.py` or to `_style.PALETTE`'s existing keys (may *add* dashboard-only color constants in the report layer, not mutate the shared palette).
- No changes to compute semantics of DMC/DMR/QC engines — the report only *reads* `md`.
- No server/interactive backend — output remains a single static `.html`.
- CLI flags for every new knob are not required for 1.0; expose `self_contained` and keep existing flags.

## 4. Architecture

The existing layering is sound and is preserved:

```
md (MethylData)
  │  read-only
  ▼
_compute.py  ── small picklable result objects (numpy/polars), heavy scans, cached in uns["_report_cache"]
  │
  ▼
_plotly.py   ── render adapters → plotly.graph_objects.Figure (interactive twins)
  │
  ▼
report.py    ── orchestrator: builds the Jinja context dict, calls compute/plotly via _safe(), writes HTML
  │
  ▼
templates/report.html.j2 + report.css  ── presentation
```

**What changes per layer:**

- **`_compute.py`** — add small compute helpers for new sections (summary stats, preprocessing deltas, QC pass/warn/fail evaluation, p-value-histogram bins, DMR-size bins, scree, global-methylation-per-sample, correlation matrix densification, methods-text facts). Reuse existing computes wherever possible (`compute_categorical_proportions`, volcano/MA/manhattan/PCA/metaplot, `qc_sample_correlation`). Same caching + single-scan discipline.
- **`_plotly.py`** — add Plotly twins for the new figures: p-value histogram, DMR-size histogram, global-methylation bar, sample-correlation heatmap, scree, hyper/hypo-by-feature stacked bar. All on white backgrounds (theme-agnostic), captioned by the template. Apply the refined dashboard accent constants while keeping `hyper`/`hypo`/feature colors.
- **`report.py`** — rewritten orchestration: new context keys for the new sections; `self_contained` parameter controlling `include_plotlyjs` (embed once vs CDN); helper to elide file paths; helper to build the auto-generated methods text and narrative summary from provenance/uns; status computation per section (ok/warn/skip). Backward-compatible signature: all new params keyword-only with defaults; existing params unchanged.
- **`templates/report.html.j2`** — rewritten to the sidebar+main dashboard layout with numbered sections, status pills, captions, and the new sections. Inlined vanilla JS block (scroll-spy, table sort/search, CSV download, theme toggle, copy-methods).
- **`templates/report.css`** — rewritten as a design-token system (CSS variables for light + dark), sidebar/main grid, components (cards, KPI, badges, chips, step-flow, tables), print stylesheet.

## 5. Section inventory

Each section is conditional on its data and shows a status (ok / warn / skip). `*` = new content.

| # | Section | Source | New? |
|---|---|---|---|
| 00 | Results at a glance — narrative + KPI strip + completeness checklist | derived from uns/obs/dmc/dmr | * |
| 01 | Samples — chips + formatted, path-elided, sortable/searchable table | `md.obs` | improved |
| 02 | Preprocessing — step-flow with site-count deltas + params | `uns["_store_history"]`, `uns["filter"]` | improved |
| 03 | Quality control — per-sample pass/warn/fail badges, coverage hist, global-meth bar*, correlation heatmap* | `md.obs`, `uns["qc_*"]` | improved + * |
| 04 | DMC — enriched KPIs, volcano, p-value histogram*, MA, Manhattan, sortable+CSV top table, engine provenance line | `md.dmc`, `uns["dmc"]` | improved + * |
| 05 | DMR — KPIs, size distribution*, sortable+CSV table | `uns["dmr"]`, `uns["dmr_params"]` | improved + * |
| 06 | Annotation — feature donut, CpG-context donut, hyper/hypo-by-feature stacked bar* | `md.dmc` (annotated) | improved + * |
| 07 | TSS metaplot | `compute_tss_metaplot` (gtf_path) | unchanged |
| 08 | Sample similarity — PCA + scree* | `compute_pca` | improved + * |
| 09 | Methods & citations — auto methods paragraph (copy), how-to-cite, parameters table | provenance/uns | * |
| 10 | Provenance — key/value table + raw JSON collapsible | provenance payload | improved |

### QC thresholds (status badges)
Threshold-based, surfaced explicitly in the UI. Defaults (tunable via report kwargs, sensible for WGBS):
- `bisulfite_conversion_rate` ≥ 0.99 → pass, ≥ 0.98 → warn, else fail (only when CHH store supplied; otherwise "not assessed").
- `frac_ge_10x` ≥ 0.70 → pass, ≥ 0.50 → warn, else fail.
- `mean_coverage` ≥ 10× → pass, ≥ 5× → warn, else fail.
- `low_coverage_flag` true → warn.
- `min_pairwise_corr` ≥ 0.85 → pass, ≥ 0.75 → warn, else fail (only when `run_sample_correlation` was used).
- `sex_mismatch` true → fail; `contamination_score` above its method threshold → warn.
A sample's row status is the worst of its assessed metrics. The section pill aggregates: any fail → warn pill with count; otherwise ok. Thresholds rendered in a small legend so the report is self-documenting.

## 6. Interactivity (inlined vanilla JS, ~60–90 lines)

- **Scroll-spy**: `IntersectionObserver` highlights the active TOC entry.
- **Table sort**: click a header → numeric-aware sort (strips `, % ×`), toggles asc/desc.
- **Table search**: per-table filter box hides non-matching rows.
- **CSV download**: client-side `Blob` from the rendered table (top-N is what downloads; full-table export remains `epykit export`).
- **Theme toggle**: light/dark via `data-theme` on `<html>`, persisted in `localStorage`. Figures live on white plot-cards so Plotly never needs re-theming.
- **Copy methods**: copies the methods paragraph to clipboard.

All JS is defensive (null-guarded) so missing sections never throw.

## 7. `self_contained` behavior

- `generate_report(..., self_contained: bool = True)`.
- `True`: the first emitted figure uses `include_plotlyjs=True` (full inline bundle), subsequent figures use `include_plotlyjs=False`; if **zero** figures render (e.g. nothing run), no Plotly is needed. File is portable/offline (~3–4 MB).
- `False`: figures use `include_plotlyjs="cdn"` (current behavior); smaller file, needs internet.
- The CSS and JS are always inlined regardless of this flag (they are tiny).
- CLI `epykit report` gains `--self-contained/--no-self-contained` (default self-contained).

## 8. Backward compatibility

- `generate_report` keeps its current positional/keyword params; new params (`self_contained`, optional QC threshold overrides) are keyword-only with defaults.
- `MethylData.report(...)` wrapper passes through.
- Output is still a single self-contained `.html` at the same path.
- Existing tests for report generation must still pass; new tests added.

## 9. Testing

- **Full pipeline**: synthetic fixture run through filter→unite→qc→dmc→dmr→annotate; assert the HTML contains each section's anchor, the KPI numbers, the methods text, and that `self_contained=True` embeds Plotly (no `cdn` script src) while `False` references the CDN.
- **Partial pipelines** (new): (a) only import+filter (no DMC/DMR/QC) → DMC/DMR/QC sections render "not run", sidebar dots hollow, no crash; (b) DMC but no DMR; (c) QC run without `run_sample_correlation` → no heatmap, no crash.
- **No-figure path**: a report with nothing plottable still writes valid HTML and embeds no Plotly when `self_contained=True`.
- **Determinism**: report renders without network access when `self_contained=True` (Plotly embedded).
- Keep tests Windows-safe and within the existing `pytest -m "not slow"` budget; heavy paths (metaplot, PCA) reuse cached computes / small fixtures.

## 10. Risks & mitigations

- **File size** with embedded Plotly (~3–4 MB): acceptable for self-contained default; `self_contained=False` available. Document in the docstring.
- **Plotly inline-once ordering**: the "first figure inlines the bundle" trick depends on render order; mitigate by having `report.py` track an `_plotly_embedded` flag and pass `include_plotlyjs` explicitly per figure rather than relying on template order.
- **Dark-mode + Plotly**: avoided entirely by keeping figures on white plot-cards.
- **Scope creep**: the new sections are bounded to data the pipeline already produces; no new analyses are computed beyond cheap derivations (bins, deltas, thresholds).

## 11. File-level work summary

- `src/epykit/report.py` — rewrite orchestration; add `self_contained`, status/methods/narrative/path-elide helpers.
- `src/epykit/pl/_compute.py` — add summary/QC-status/p-hist/dmr-size/scree/global-meth/corr-matrix/methods-facts computes (reusing existing where possible).
- `src/epykit/pl/_plotly.py` — add p-value histogram, DMR-size hist, global-meth bar, correlation heatmap, scree, stacked-by-feature twins; apply refined accents.
- `src/epykit/templates/report.html.j2` — rewrite to sidebar+main dashboard; inline JS.
- `src/epykit/templates/report.css` — rewrite as token-based design system + dark + print.
- `src/epykit/cli.py` — add `--self-contained/--no-self-contained` to `report`.
- `tests/` — extend report tests; add partial-pipeline + self-contained tests.

A throwaway visual reference of the target look lives at `demo_output/mockup_report.html` (not shipped).
