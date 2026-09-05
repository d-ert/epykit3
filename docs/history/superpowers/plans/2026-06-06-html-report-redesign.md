# HTML Report Redesign — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (inline) to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Rebuild `epykit`'s single-file HTML report as a MultiQC-style dashboard (fixed sidebar TOC + scroll-spy, dark/light, scannable sections) with new scientific content, a `self_contained` offline flag, and preserved graceful-skip robustness.

**Architecture:** Keep the existing `_compute.py → _plotly.py → report.py → Jinja template` layering. Add small compute helpers + Plotly twins for new figures, rewrite `report.py` orchestration (status/narrative/methods/path-elide + `self_contained` embed control), and rewrite the Jinja template + CSS into a token-based dashboard with ~80 lines of inlined vanilla JS. No new dependencies.

**Tech Stack:** Python, polars, numpy, jinja2, plotly; vanilla JS/CSS (no framework).

**Spec:** `docs/history/superpowers/specs/2026-06-06-html-report-redesign-design.md`

**Compatibility contract (existing `tests/test_report.py` must keep passing):** the rendered HTML must still contain the substrings `ep.tl.dmc` and `ep.tl.dmr` (skip hints), `Provenance`, a case-insensitive `volcano`, `plotly`, the report title, and each sample_id. Size > 30 kB for a full report.

---

## File Structure

- `src/epykit/pl/_compute.py` — **modify**: add `compute_sample_correlation_matrix`, `compute_pvalue_histogram`, `compute_dmr_size_distribution`, `compute_global_methylation`, `compute_scree`. New result dataclasses where useful.
- `src/epykit/pl/_plotly.py` — **modify**: add `pvalue_histogram_plotly`, `dmr_size_hist_plotly`, `global_methylation_bar_plotly`, `sample_correlation_plotly`, `scree_plotly`, `feature_direction_stacked_plotly`. Add a shared `_dashboard_layout()` helper (white bg, refined accents).
- `src/epykit/report.py` — **rewrite orchestration**: add `self_contained` param + per-figure embed control; helpers `_section_status`, `_summary_facts`, `_methods_text`, `_qc_rows`, `_elide_path`, `_preproc_flow`. Build the expanded Jinja context.
- `src/epykit/templates/report.html.j2` — **rewrite**: sidebar + main dashboard, numbered sections w/ status pills, captions, inlined JS.
- `src/epykit/templates/report.css` — **rewrite**: design-token system (light+dark), components, print stylesheet.
- `src/epykit/cli.py` — **modify** `_cmd_report` + `report` subparser: add `--self-contained/--no-self-contained`.
- `tests/test_report.py` — **extend**: self_contained on/off, new section markers, QC badge states, partial-pipeline robustness.

---

## Task 1: Compute helpers for new figures

**Files:**
- Modify: `src/epykit/pl/_compute.py`
- Test: `tests/test_report_compute.py` (create)

- [ ] **Step 1: Write failing tests** in `tests/test_report_compute.py`:

```python
"""Tests for new report compute helpers."""
from __future__ import annotations
import numpy as np
import polars as pl
import pytest


def test_pvalue_histogram(synth_md_filtered):
    import epykit as ep
    from epykit.pl._compute import compute_pvalue_histogram
    ep.tl.dmc(synth_md_filtered, test="lr")
    counts, edges = compute_pvalue_histogram(synth_md_filtered, bins=20)
    assert counts.sum() > 0
    assert len(edges) == len(counts) + 1
    assert edges[0] == 0.0 and abs(edges[-1] - 1.0) < 1e-9


def test_dmr_size_distribution(synth_md_filtered):
    import epykit as ep
    md = synth_md_filtered
    ep.tl.dmc(md, test="lr")
    ep.tl.dmr(md, method="tile", tile_size_bp=500, min_cpgs_per_tile=3, min_mean_qvalue=1.0)
    sizes = compute_dmr_size_distribution(md)
    from epykit.pl._compute import compute_dmr_size_distribution
    assert sizes.ndim == 1


def test_global_methylation(synth_md_filtered):
    import epykit as ep
    from epykit.pl._compute import compute_global_methylation
    ep.tl.qc(synth_md_filtered)
    samples, values, groups = compute_global_methylation(synth_md_filtered)
    assert len(samples) == len(values) == synth_md_filtered.n_samples
    assert all(0.0 <= v <= 1.0 for v in values if v == v)


def test_correlation_matrix(synth_md_filtered):
    import epykit as ep
    from epykit.pl._compute import compute_sample_correlation_matrix
    ep.tl.qc(synth_md_filtered, run_sample_correlation=True)
    mat, labels = compute_sample_correlation_matrix(synth_md_filtered)
    assert mat.shape == (len(labels), len(labels))
    assert np.allclose(np.diag(mat), 1.0, atol=0.05)


def test_scree(synth_md_filtered):
    from epykit.pl._compute import compute_scree
    ev = compute_scree(synth_md_filtered, n_sites=2000, max_components=4)
    assert ev.ndim == 1 and len(ev) >= 1
    assert (ev >= 0).all()
```

- [ ] **Step 2: Run, expect ImportError/fail**

Run: `uv run pytest tests/test_report_compute.py -x -q`
Expected: FAIL (functions undefined).

- [ ] **Step 3: Implement** in `src/epykit/pl/_compute.py` (append before `__all__`, and extend `__all__`):

```python
def compute_pvalue_histogram(md, *, bins: int = 30):
    """Histogram of raw per-CpG p-values (calibration check)."""
    dmc = md.dmc
    if dmc is None or "pvalue" not in dmc.columns:
        raise ValueError("Run ep.tl.dmc(md) first")
    p = dmc["pvalue"].to_numpy()
    p = p[~np.isnan(p)]
    counts, edges = np.histogram(p, bins=bins, range=(0.0, 1.0))
    return counts, edges


def compute_dmr_size_distribution(md) -> np.ndarray:
    """Per-DMR CpG counts (or genomic width if n_cpgs absent)."""
    dmr = md.uns.get("dmr")
    if dmr is None or not isinstance(dmr, pl.DataFrame) or dmr.is_empty():
        raise ValueError("Run ep.tl.dmr(md) first")
    if "n_cpgs" in dmr.columns:
        return dmr["n_cpgs"].drop_nulls().to_numpy()
    if {"start", "end"} <= set(dmr.columns):
        return (dmr["end"] - dmr["start"]).to_numpy()
    raise ValueError("DMR table lacks n_cpgs and start/end")


def compute_global_methylation(md):
    """Return (samples, global_methylation, groups) for the per-sample bar."""
    obs = md.obs
    if "global_methylation" not in obs.columns:
        raise ValueError("Run ep.tl.qc(md) first (no global_methylation in obs)")
    samples = obs.get_column("sample_id").to_list()
    values = [float(v) if v is not None else float("nan")
              for v in obs.get_column("global_methylation").to_list()]
    _, groups = _resolve_group_col(md, None)
    return samples, values, groups


def compute_sample_correlation_matrix(md):
    """Densify uns['qc_sample_correlation'] (long form) into a clustered matrix.

    Returns (matrix, labels). Hierarchical clustering order applied when
    scipy is available; otherwise sample order is preserved.
    """
    corr_df = md.uns.get("qc_sample_correlation")
    if corr_df is None or not isinstance(corr_df, pl.DataFrame) or corr_df.is_empty():
        raise ValueError("Run ep.tl.qc(md, run_sample_correlation=True) first")
    labels = list(dict.fromkeys(
        corr_df.get_column("sample_a").to_list()
        + corr_df.get_column("sample_b").to_list()
    ))
    idx = {s: i for i, s in enumerate(labels)}
    n = len(labels)
    mat = np.full((n, n), np.nan)
    for row in corr_df.iter_rows(named=True):
        i, j = idx[row["sample_a"]], idx[row["sample_b"]]
        mat[i, j] = mat[j, i] = row["correlation"]
    np.fill_diagonal(mat, 1.0)
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform
        filled = np.where(np.isnan(mat), np.nanmean(mat), mat)
        dist = 1.0 - filled
        np.fill_diagonal(dist, 0.0)
        order = leaves_list(linkage(squareform(dist, checks=False), method="average"))
        mat = mat[np.ix_(order, order)]
        labels = [labels[i] for i in order]
    except Exception:
        pass
    return mat, labels


def compute_scree(md, *, n_sites: int = 10_000, max_components: int = 6,
                  seed: int = 42) -> np.ndarray:
    """Explained-variance ratios for the first ``max_components`` PCs."""
    try:
        from sklearn.decomposition import PCA
    except ImportError as exc:
        raise ImportError("scikit-learn is required for the scree plot") from exc
    matrix, _samples, _n = compute_sample_site_matrix(md, n_sites=n_sites, seed=seed)
    k = min(max_components, matrix.shape[0], matrix.shape[1])
    pca = PCA(n_components=k)
    pca.fit(matrix)
    return np.asarray(pca.explained_variance_ratio_, dtype=float)
```

- [ ] **Step 4: Run, expect PASS**

Run: `uv run pytest tests/test_report_compute.py -x -q`
Expected: PASS (5 tests).

- [ ] **Step 5: Commit** — `feat(report): compute helpers for p-value hist, DMR size, global meth, correlation matrix, scree`.

---

## Task 2: Plotly twins for new figures

**Files:**
- Modify: `src/epykit/pl/_plotly.py`
- Test: `tests/test_report_compute.py` (extend with smoke tests)

- [ ] **Step 1: Add smoke tests** asserting each twin returns a Plotly Figure with ≥1 trace (skip if compute raises). Pattern:

```python
def test_plotly_twins_smoke(synth_md_filtered):
    import epykit as ep
    from epykit.pl import _plotly as P
    md = synth_md_filtered
    ep.tl.qc(md, run_sample_correlation=True)
    ep.tl.dmc(md, test="lr")
    ep.tl.dmr(md, method="tile", tile_size_bp=500, min_cpgs_per_tile=3, min_mean_qvalue=1.0)
    for fn in (P.pvalue_histogram_plotly, P.global_methylation_bar_plotly,
               P.sample_correlation_plotly, P.scree_plotly, P.dmr_size_hist_plotly):
        fig = fn(md)
        assert fig is not None and len(fig.data) >= 1
```

- [ ] **Step 2: Run, expect fail.** `uv run pytest tests/test_report_compute.py::test_plotly_twins_smoke -x -q`

- [ ] **Step 3: Implement** in `_plotly.py`. Import the new computes. Add a shared layout + the twins:

```python
# refined dashboard accents (figures only; matplotlib PALETTE untouched)
_ACCENT = "#2563eb"
_VIOLET = "#7c3aed"

def _dash_layout(go, **over):
    base = dict(paper_bgcolor="white", plot_bgcolor="white",
                font=dict(family="system-ui,-apple-system,Segoe UI,sans-serif",
                          size=12, color="#1a2230"),
                margin=dict(l=54, r=18, t=18, b=44), template="simple_white")
    base.update(over)
    return base

def pvalue_histogram_plotly(md):
    go = _require_plotly()
    counts, edges = compute_pvalue_histogram(md, bins=30)
    centers = (edges[:-1] + edges[1:]) / 2
    fig = go.Figure([go.Bar(x=centers, y=counts, width=(edges[1]-edges[0])*0.95,
                            marker_color=_ACCENT, opacity=0.85)])
    fig.update_layout(**_dash_layout(go, xaxis_title="p-value", yaxis_title="CpG count",
                                     height=320, bargap=0.02))
    return fig

def dmr_size_hist_plotly(md):
    go = _require_plotly()
    sizes = compute_dmr_size_distribution(md)
    fig = go.Figure([go.Histogram(x=sizes, marker_color=_VIOLET, opacity=0.85)])
    fig.update_layout(**_dash_layout(go, xaxis_title="CpGs per DMR",
                                     yaxis_title="DMR count", height=300, bargap=0.05))
    return fig

def global_methylation_bar_plotly(md):
    go = _require_plotly()
    samples, values, groups = compute_global_methylation(md)
    uniq = list(dict.fromkeys(groups))
    cmap = {uniq[0]: PALETTE["treatment"]} if uniq else {}
    if len(uniq) > 1:
        cmap = {g: (PALETTE["treatment"] if i == 0 else PALETTE["control"]
                    if i == 1 else PALETTE["neutral"]) for i, g in enumerate(uniq)}
    colors = [cmap.get(g, PALETTE["neutral"]) for g in groups]
    fig = go.Figure([go.Bar(x=samples, y=values, marker_color=colors)])
    fig.update_layout(**_dash_layout(go, yaxis_title="global methylation",
                                     height=300))
    return fig

def sample_correlation_plotly(md):
    go = _require_plotly()
    mat, labels = compute_sample_correlation_matrix(md)
    fig = go.Figure([go.Heatmap(z=mat, x=labels, y=labels,
                                colorscale=[[0, "#ffffff"], [1, _ACCENT]],
                                zmin=float(min(0.7, np.nanmin(mat))), zmax=1.0,
                                colorbar=dict(thickness=12, len=0.75))])
    fig.update_layout(**_dash_layout(go, height=360, margin=dict(l=60, r=20, t=18, b=50)))
    return fig

def scree_plotly(md):
    go = _require_plotly()
    ev = compute_scree(md, max_components=6) * 100.0
    labels = [f"PC{i+1}" for i in range(len(ev))]
    fig = go.Figure([go.Bar(x=labels, y=ev, marker_color=_ACCENT)])
    fig.update_layout(**_dash_layout(go, yaxis_title="% variance explained", height=340))
    return fig

def feature_direction_stacked_plotly(md):
    """Hyper vs hypo proportion within each gene feature (annotated DMC)."""
    go = _require_plotly()
    dmc = md.dmc
    if dmc is None or "feature_type" not in dmc.columns or "meth_diff" not in dmc.columns:
        return None
    import polars as _pl
    work = dmc.with_columns(
        _pl.when(_pl.col("meth_diff") > 0).then(_pl.lit("hyper"))
        .otherwise(_pl.lit("hypo")).alias("dmr_type")
    )
    prop = compute_categorical_proportions(work, group_col="dmr_type",
                                           annot_col="feature_type",
                                           include_all_group=False, normalize=True)
    feats = prop["feature_type"].unique().to_list()
    fig = go.Figure()
    for direction, color in (("hyper", PALETTE["hyper"]), ("hypo", PALETTE["hypo"])):
        sub = prop.filter(_pl.col("dmr_type") == direction)
        ymap = dict(zip(sub["feature_type"].to_list(), sub["proportion"].to_list()))
        fig.add_trace(go.Bar(x=feats, y=[ymap.get(f, 0.0) for f in feats],
                             name=direction, marker_color=color))
    fig.update_layout(**_dash_layout(go, barmode="stack", height=300,
                                     yaxis_title="proportion", showlegend=True,
                                     legend=dict(orientation="h")))
    return fig
```

Add the new compute imports at the top (`compute_pvalue_histogram, compute_dmr_size_distribution, compute_global_methylation, compute_sample_correlation_matrix, compute_scree, compute_categorical_proportions`) and `import numpy as np` already present. Extend `__all__`.

- [ ] **Step 4: Run, expect PASS.**

- [ ] **Step 5: Commit** — `feat(report): plotly twins for new dashboard figures`.

---

## Task 3: `report.py` orchestration rewrite

**Files:**
- Modify: `src/epykit/report.py`
- Test: covered by Task 7.

- [ ] **Step 1: Add `self_contained` + embed control.** Replace `_fig_html` with a closure/counter so the first embedded figure inlines Plotly:

```python
def _make_fig_renderer(self_contained: bool):
    state = {"embedded": False}
    def render(fig):
        if fig is None:
            return None
        if self_contained:
            inc = True if not state["embedded"] else False
            if inc:
                state["embedded"] = True
        else:
            inc = "cdn"
        try:
            return fig.to_html(include_plotlyjs=inc, full_html=False, default_height=None)
        except Exception as exc:
            logger.warning("Failed to render Plotly figure: %s", exc)
            return None
    return render
```

- [ ] **Step 2: Add helpers** to `report.py`:
  - `_elide_path(p, keep=1)` → show basename (and parent), full path retained for `title=`.
  - `_section_status(...)` → returns `"ok" | "warn" | "skip"` per section (sidebar dot + pill).
  - `_summary_facts(md, dmc_stats, dmr_stats)` → dict for the narrative + KPI strip (n_samples, group breakdown, n_raw/n_filtered + pct, hyper/hypo, % sig, median |Δβ| of sig, global meth per group).
  - `_qc_rows(md, thresholds)` → list of dicts per sample with formatted metrics + a `status` (`pass|warn|fail`) computed from thresholds (defaults from spec §5; allow override via `qc_thresholds` kwarg).
  - `_preproc_flow(md)` → list of `{step, n_sites, delta_pct}` from `_store_history`.
  - `_methods_text(md, alpha, min_abs_diff)` → parameter-accurate paragraph string built from `uns["dmc"]`, `uns["dmr_params"]`, `uns["filter"]`, version.

- [ ] **Step 3: Update `generate_report` signature** — add keyword-only params (defaults preserve behavior):

```python
def generate_report(md, output, *, title=None, gtf_path=None,
                    alpha=0.05, min_abs_diff=0.1, dmc_top_n=50, dmr_top_n=50,
                    metaplot_max_genes=5000, pca_n_sites=10_000,
                    coverage_max_points=200_000, clear_cache=False,
                    self_contained: bool = True,
                    qc_thresholds: Optional[dict] = None) -> str:
```

- [ ] **Step 4: Build the expanded context dict** — add keys consumed by the new template: `report_title`, chips, `summary` (facts), `kpis`, `completeness` (checklist items w/ done flag), `sections` statuses, `samples_rows` (pre-elided), `preproc_flow`, `qc_rows`, `qc_legend`, the new figure HTML (`pvalue_hist`, `dmr_size_hist`, `global_meth_bar`, `corr_heatmap`, `scree`, `feature_direction_stacked`), `methods_text`, `cite_text`, `params_rows`, `provenance_rows` (clean k/v) + `provenance` (raw JSON). Keep `volcano_plot`, `ma_plot`, `manhattan_plot`, `coverage_plot`, `feature_pie`, `cpg_pie`, `metaplot`, `pca_plot`. Wrap every figure call in `_safe(...)` and render via the new renderer. **Preserve compat substrings**: skip-state text must include `ep.tl.dmc` / `ep.tl.dmr`; keep a `Provenance` heading; volcano title remains "DMC volcano".

- [ ] **Step 5: Smoke-run** `uv run python -c "import epykit.report"` to catch syntax errors.

- [ ] **Step 6: Commit** — `feat(report): dashboard orchestration + self_contained embed control`.

---

## Task 4: Template rewrite (`report.html.j2`)

**Files:**
- Modify: `src/epykit/templates/report.html.j2`

- [ ] **Step 1:** Rewrite to the sidebar+main structure validated in `demo_output/mockup_report.html`, driven by the Task 3 context. Sections 00–10 per spec §5. Sidebar nav loops the section list with status dots; each section is conditional (`{% if ... %}`) and renders a "not run" state otherwise. Plotly figure HTML strings are inserted with `| safe`. Include the inlined `<script>` block (scroll-spy, table sort/search, CSV download, theme toggle, copy methods) and `<style>{{ css_inline }}</style>`. Keep skip hints `ep.tl.dmc` / `ep.tl.dmr` and a `Provenance` heading.

- [ ] **Step 2:** Verify Jinja renders without undefined-variable errors via the Task 7 tests.

- [ ] **Step 3: Commit** — `feat(report): MultiQC-style dashboard template + inlined JS`.

---

## Task 5: CSS rewrite (`report.css`)

**Files:**
- Modify: `src/epykit/templates/report.css`

- [ ] **Step 1:** Replace with the token-based system from the mockup: `:root` light vars + `html[data-theme="dark"]` overrides; sidebar/main layout; components (`.card`, `.kpi`, `.badge`, `.chip`, `.grp`, `.flow/.step`, `table.df`, `.plot-card`, `.callout`, `.params`); `@media print` (hide sidebar, expand main, avoid break-inside on cards). Keep figures on white plot-cards in both themes.

- [ ] **Step 2: Commit** — `feat(report): token-based dashboard CSS with dark mode + print`.

---

## Task 6: CLI `--self-contained` flag

**Files:**
- Modify: `src/epykit/cli.py`

- [ ] **Step 1:** In the `report` subparser (≈line 899) add:

```python
p_rep.add_argument("--self-contained", dest="self_contained",
                   action="store_true", default=True,
                   help="Embed Plotly inline so the HTML works offline (default)")
p_rep.add_argument("--no-self-contained", dest="self_contained",
                   action="store_false",
                   help="Load Plotly from CDN (smaller file, needs internet)")
```

In `_cmd_report` add `kwargs["self_contained"] = args.self_contained`.

- [ ] **Step 2:** `uv run epykit report --help` shows the flags. Commit — `feat(cli): report --self-contained/--no-self-contained`.

---

## Task 7: Report integration tests

**Files:**
- Modify: `tests/test_report.py`

- [ ] **Step 1: Keep the 3 existing tests passing** (they assert compat substrings). Add:

```python
def test_report_self_contained_embeds_plotly(synth_md_filtered, tmp_path):
    import epykit as ep
    ep.tl.dmc(synth_md_filtered, test="lr")
    out = tmp_path / "sc.html"
    synth_md_filtered.report(str(out), self_contained=True)
    html = out.read_text(encoding="utf-8")
    # Embedded bundle defines Plotly inline; no CDN <script src=...plotly...>
    assert "https://cdn.plot.ly" not in html
    assert "Plotly" in html
    assert out.stat().st_size > 1_000_000  # full bundle inlined


def test_report_cdn_mode_uses_cdn(synth_md_filtered, tmp_path):
    import epykit as ep
    ep.tl.dmc(synth_md_filtered, test="lr")
    out = tmp_path / "cdn.html"
    synth_md_filtered.report(str(out), self_contained=False)
    html = out.read_text(encoding="utf-8")
    assert "cdn.plot.ly" in html


def test_report_dashboard_sections_and_qc_badges(synth_md_filtered, tmp_path):
    import epykit as ep
    md = synth_md_filtered
    ep.tl.qc(md, run_sample_correlation=True)
    ep.tl.dmc(md, test="lr")
    out = tmp_path / "dash.html"
    md.report(str(out), title="dash")
    html = out.read_text(encoding="utf-8")
    for anchor in ('id="summary"', 'id="qc"', 'id="dmc"', 'id="methods"', 'id="prov"'):
        assert anchor in html, anchor
    assert "Results at a glance" in html
    assert "Methods" in html
    # QC badge classes present
    assert "badge" in html


def test_report_no_figures_still_self_contained(synth_md_filtered, tmp_path):
    # Nothing plottable beyond coverage; must not crash and must be valid html
    out = tmp_path / "min.html"
    synth_md_filtered.report(str(out), self_contained=True)
    assert out.exists()
    assert "<html" in out.read_text(encoding="utf-8").lower()
```

- [ ] **Step 2: Run full report suite** — `uv run pytest tests/test_report.py tests/test_report_compute.py -q`. Expected: all pass.

- [ ] **Step 3:** Update the existing `test_report_partial_pipeline` CDN comment/assertion only if it conflicts (it asserts `"plotly" in html.lower()`, which holds for embedded too). Keep as-is if green.

- [ ] **Step 4: Commit** — `test(report): self_contained, dashboard sections, QC badges, robustness`.

---

## Task 8: Full-suite verification + real render

- [ ] **Step 1:** `uv run pytest -m "not slow" --strict-markers -q` — entire suite green (no regressions in plot/CLI tests).
- [ ] **Step 2:** `uv run ruff check src/` and `uv run mypy src/epykit` — no new errors.
- [ ] **Step 3: Render the real demo** by re-running the report on the existing synthetic store and visually verify in the browser (sidebar, scroll-spy, dark toggle, tables, all sections populated). Regenerate `demo_output/demo_report.html`.
- [ ] **Step 4:** Update `CHANGELOG.md` with the report redesign entry. Commit — `docs(changelog): HTML report redesign`.
- [ ] **Step 5:** Remove the throwaway `demo_output/mockup_report.html` and the temporary `.claude/launch.json` if not wanted (leave demo_output gitignored).

---

## Self-Review

**Spec coverage:** §5 sections → Tasks 3–5 (all 11 sections in template); new figures (p-hist, dmr-size, global-meth, corr, scree, stacked) → Tasks 1–2; QC thresholds §5 → `_qc_rows` (Task 3 Step 2) + Task 7 badge test; `self_contained` §7 → Task 3 Step 1 + Task 6 + Task 7; backward-compat §8 → compat-substring contract + existing tests retained (Task 7); robustness §5 → conditional template (Task 4) + partial-pipeline tests (Task 7). Methods/citations §5 → `_methods_text` + template §09. Interactivity §6 → inlined JS (Task 4).

**Placeholder scan:** compute + plotly code is complete; template/CSS reference the validated mockup (`demo_output/mockup_report.html`) as the concrete source, copied and parameterized — not "TBD". Acceptable because the artifact exists in-repo.

**Type consistency:** `compute_global_methylation` returns `(samples, values, groups)` consumed identically by `global_methylation_bar_plotly`. `compute_sample_correlation_matrix` → `(mat, labels)` consumed by `sample_correlation_plotly`. `_make_fig_renderer` returns a callable used for every figure. Context keys produced in Task 3 match those consumed in Task 4 (single author, verified at render via Task 7).
