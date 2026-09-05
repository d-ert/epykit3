"""HTML report generator for a MethylData object.

Renders a single-file, MultiQC-style dashboard: a fixed sidebar table of
contents (with per-section status dots and scroll-spy) over a scannable
main panel covering a headline summary, sample metadata, the preprocessing
trail, QC, DMC/DMR, annotation, optional TSS metaplot, PCA, auto-generated
methods & citations, and provenance.

Each section is conditional on the relevant data being populated on ``md``;
unrun sections render a short "not run yet" notice (and a hollow sidebar
dot) instead of crashing.

Interactivity (scroll-spy, sortable/searchable tables, CSV download, theme
toggle, copy-methods) is implemented in a small block of vanilla JS inlined
into the template -- no JS framework, no new runtime dependencies.

Optional deps: ``jinja2`` and ``plotly`` (install via
``pip install 'epykit[report]'``). Both are imported lazily so the rest of
the package keeps working with the base install.
"""

from __future__ import annotations

import datetime
import html
import json
import logging
from pathlib import Path
from typing import Any, Callable

import numpy as np
import polars as pl

from .methyldata import MethylData

logger = logging.getLogger(__name__)

_TEMPLATE_DIR = Path(__file__).parent / "templates"

# Default QC pass/warn/fail thresholds (WGBS-sensible). Override per call via
# ``generate_report(..., qc_thresholds={...})``.
_DEFAULT_QC_THRESHOLDS: dict = {
    "bisulfite_conversion_rate": (0.99, 0.98),  # >= pass, >= warn, else fail
    "frac_ge_10x": (0.70, 0.50),
    "mean_coverage": (10.0, 5.0),
    "min_pairwise_corr": (0.85, 0.75),
}


# ---------------------------------------------------------------------------
# Dependency shims
# ---------------------------------------------------------------------------


def _require_jinja():
    try:
        from jinja2 import Environment, FileSystemLoader, select_autoescape

        return Environment, FileSystemLoader, select_autoescape
    except ImportError as exc:
        raise ImportError(
            "jinja2 is required for HTML report generation. "
            "Install with: pip install 'epykit[report]'"
        ) from exc


# ---------------------------------------------------------------------------
# Small formatting helpers
# ---------------------------------------------------------------------------


def _esc(v: Any) -> str:
    return html.escape("" if v is None else str(v))


def _fmt_value(v: Any) -> str:
    if v is None:
        return "--"
    if isinstance(v, float):
        if v != v:  # NaN
            return "--"
        if abs(v) >= 1000 or (abs(v) < 0.001 and v != 0):
            return f"{v:.4g}"
        return f"{v:.4f}".rstrip("0").rstrip(".") or "0"
    if isinstance(v, int):
        return f"{v:,}" if abs(v) >= 1000 else str(v)
    return str(v)


def _fmt_int(v: Any) -> str:
    try:
        return f"{int(v):,}"
    except (TypeError, ValueError):
        return "--"


def _fmt_pct(v: Any, digits: int = 1) -> str:
    """Render a 0-1 fraction as a percentage. Values > 1 are treated as
    already-percent and shown verbatim."""
    if v is None or (isinstance(v, float) and v != v):
        return "--"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "--"
    if f <= 1.0:
        f *= 100.0
    return f"{f:.{digits}f}%"


def _fmt_sci(v: Any) -> str:
    if v is None or (isinstance(v, float) and v != v):
        return "--"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "--"
    if f == 0:
        return "0"
    if f < 1e-3 or f >= 1e4:
        return f"{f:.2e}"
    return f"{f:.4g}"


def _elide_path(p: Any, *, keep: int = 1) -> tuple[str, str]:
    """Return ``(short, full)`` for a file path -- short is the basename
    (optionally prefixed with ``keep`` parent dirs), full is the original."""
    if p is None:
        return "--", ""
    full = str(p)
    parts = full.replace("\\", "/").rstrip("/").split("/")
    short = "/".join(parts[-(keep + 1) :]) if len(parts) > keep else full
    return short, full


# ---------------------------------------------------------------------------
# Generic, spec-driven HTML table (works with the inlined sort/search/CSV JS)
# ---------------------------------------------------------------------------


def _html_table(rows: list[dict], columns: list[dict], table_id: str) -> str:
    """Render ``rows`` as a ``table.df``.

    ``columns`` is a list of dicts: ``{"key", "label", "num"?, "fmt"?}``
    where ``fmt(value, row) -> html`` (already-trusted HTML; plain strings
    are escaped by the default formatter).
    """
    if not rows:
        return '<div class="skipped">no rows</div>'
    head = "".join(
        f'<th class="{"num" if c.get("num") else ""}">{_esc(c["label"])}'
        f'<span class="ar">↕</span></th>'
        for c in columns
    )
    body_rows = []
    for row in rows:
        cells = []
        for c in columns:
            fmt: Callable = c.get("fmt") or (lambda v, r: _esc(_fmt_value(v)))
            cls = "num" if c.get("num") else ""
            cells.append(f'<td class="{cls}">{fmt(row.get(c["key"]), row)}</td>')
        body_rows.append("<tr>" + "".join(cells) + "</tr>")
    return (
        f'<table class="df" id="{table_id}"><thead><tr>{head}</tr></thead>'
        f"<tbody>{''.join(body_rows)}</tbody></table>"
    )


def _group_class_map(md: MethylData) -> dict:
    """Map each group label to a CSS class: 't' (treatment), 'c' (control),
    or 'n' (neutral) for colouring chips consistently with the figures."""
    out: dict = {}
    if "group" not in md.obs.columns:
        return out
    try:
        tids = set(md.treatment_ids)
        cids = set(md.control_ids)
    except Exception:
        tids, cids = set(), set()
    for row in md.obs.iter_rows(named=True):
        g = row.get("group")
        sid = row.get("sample_id")
        if g in out:
            continue
        if sid in tids:
            out[g] = "t"
        elif sid in cids:
            out[g] = "c"
    # Anything still unclassified cycles t/c by first-seen order.
    seen = [g for g in dict.fromkeys(md.obs.get_column("group").to_list())]
    for i, g in enumerate(seen):
        out.setdefault(g, "t" if i % 2 == 0 else "c")
    return out


def _chip(label: Any, cls: str) -> str:
    return f'<span class="grp {cls}">{_esc(label)}</span>'


# ---------------------------------------------------------------------------
# Section tables
# ---------------------------------------------------------------------------


def _samples_table(md: MethylData) -> str:
    gmap = _group_class_map(md)
    has = set(md.obs.columns)
    rows = list(md.obs.iter_rows(named=True))
    cols: list[dict] = [{"key": "sample_id", "label": "sample_id"}]
    if "group" in has:
        cols.append(
            {"key": "group", "label": "group", "fmt": lambda v, r: _chip(v, gmap.get(v, "n"))}
        )
    if "global_methylation" in has:
        cols.append(
            {
                "key": "global_methylation",
                "label": "global meth",
                "num": True,
                "fmt": lambda v, r: _esc(_fmt_pct(v)),
            }
        )
    if "mean_coverage" in has:
        cols.append(
            {
                "key": "mean_coverage",
                "label": "mean cov",
                "num": True,
                "fmt": lambda v, r: _esc(f"{v:.1f}×") if v is not None else "--",
            }
        )
    if "frac_ge_10x" in has:
        cols.append(
            {
                "key": "frac_ge_10x",
                "label": "≥10×",
                "num": True,
                "fmt": lambda v, r: _esc(_fmt_pct(v)),
            }
        )
    # Source file (elided, full path on hover) if a path-like column exists.
    path_col = next((c for c in ("path", "file", "source", "filepath") if c in has), None)
    if path_col:

        def _pf(v, r):
            short, full = _elide_path(v)
            return f'<span class="path" title="{_esc(full)}">{_esc(short)}</span>'

        cols.append({"key": path_col, "label": "source file", "fmt": _pf})
    return _html_table(rows, cols, "tSamples")


def _qc_status(value: Any, thresholds: tuple, higher_is_better: bool = True) -> str | None:
    if value is None or (isinstance(value, float) and value != value):
        return None
    pass_t, warn_t = thresholds
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    if higher_is_better:
        if v >= pass_t:
            return "pass"
        if v >= warn_t:
            return "warn"
        return "fail"
    if v <= pass_t:
        return "pass"
    if v <= warn_t:
        return "warn"
    return "fail"


def _qc_rows(md: MethylData, thresholds: dict) -> tuple[list[dict], list[str]]:
    """Build per-sample QC rows with a worst-of-metrics status badge, plus a
    human-readable legend of the thresholds actually applied."""
    has = set(md.obs.columns)
    if "global_methylation" not in has and "mean_coverage" not in has:
        return [], []
    gmap = _group_class_map(md)
    out: list[dict] = []
    applied: list[str] = []

    def note(metric, label):
        if metric in thresholds and any(metric in c for c in [has]):
            pass

    for row in md.obs.iter_rows(named=True):
        statuses: list[str] = []
        rec = {"sample_id": row.get("sample_id"), "group": row.get("group")}
        # --- absolute, meaningful per-sample gates: depth, breadth, conversion ---
        if "mean_coverage" in has:
            rec["mean_coverage"] = row.get("mean_coverage")
            s = _qc_status(row.get("mean_coverage"), thresholds["mean_coverage"])
            if s:
                statuses.append(s)
        if "frac_ge_1x" in has:
            rec["frac_ge_1x"] = row.get("frac_ge_1x")
        if "frac_ge_10x" in has:
            rec["frac_ge_10x"] = row.get("frac_ge_10x")
            s = _qc_status(row.get("frac_ge_10x"), thresholds["frac_ge_10x"])
            if s:
                statuses.append(s)
        if "global_methylation" in has:
            rec["global_methylation"] = row.get("global_methylation")
        if "bisulfite_conversion_rate" in has:
            rec["bisulfite_conversion_rate"] = row.get("bisulfite_conversion_rate")
            s = _qc_status(
                row.get("bisulfite_conversion_rate"), thresholds["bisulfite_conversion_rate"]
            )
            if s:
                statuses.append(s)
        if "min_pairwise_corr" in has:
            rec["min_pairwise_corr"] = row.get("min_pairwise_corr")
        if row.get("low_coverage_flag"):
            statuses.append("warn")
        if row.get("sex_mismatch"):
            statuses.append("fail")
        rec["_statuses"] = statuses
        rec["_group_class"] = gmap.get(row.get("group"), "n")
        out.append(rec)

    # --- sample correlation: outlier-RELATIVE advisory, never an absolute gate ---
    # Absolute pairwise r depends heavily on the methylation level (a highly
    # methylated genome saturates beta near 1, so even perfect replicates score
    # a low rank-correlation) and on the study design (a real case/control
    # cohort legitimately has moderate cross-group r -- that difference *is* the
    # signal). So we flag a sample only when its min pairwise r is a clear LOW
    # OUTLIER versus the cohort (robust MAD rule) or is absolutely broken.
    corr_floor = None
    if "min_pairwise_corr" in has:
        vals = [
            r["min_pairwise_corr"]
            for r in out
            if r.get("min_pairwise_corr") is not None
            and r["min_pairwise_corr"] == r["min_pairwise_corr"]
        ]
        if len(vals) >= 4:
            med = float(np.median(vals))
            mad = float(np.median([abs(v - med) for v in vals])) * 1.4826
            corr_floor = (med - 3.0 * mad, 0.9 * med)
        for r in out:
            v = r.get("min_pairwise_corr")
            if v is None or v != v:
                continue
            is_outlier = corr_floor is not None and v < corr_floor[0] and v < corr_floor[1]
            if is_outlier or v < 0.2:
                r["_statuses"].append("warn")

    # finalize per-sample status as the worst of its flags
    order = {"pass": 0, "warn": 1, "fail": 2}
    for r in out:
        st = r.pop("_statuses")
        r["_status"] = max(st, key=lambda s: order[s]) if st else "pass"

    # Legend (only the metrics that actually gate pass/warn/fail)
    if "mean_coverage" in has:
        applied.append("mean cov ≥ 10× pass · ≥ 5× warn")
    if "frac_ge_10x" in has:
        applied.append("≥10× fraction ≥ 70% pass · ≥ 50% warn")
    if "bisulfite_conversion_rate" in has:
        applied.append("conversion ≥ 99% pass · ≥ 98% warn")
    if "min_pairwise_corr" in has:
        applied.append(
            "min pairwise r is advisory (flags only cohort outliers, not the heatmap below)"
        )
    return out, applied


def _qc_table(md: MethylData, qc_rows: list[dict]) -> str:
    if not qc_rows:
        return ""
    has = {k for r in qc_rows for k in r}
    cols: list[dict] = [{"key": "sample_id", "label": "sample"}]
    if "mean_coverage" in has:
        cols.append(
            {
                "key": "mean_coverage",
                "label": "mean cov",
                "num": True,
                "fmt": lambda v, r: _esc(f"{v:.1f}×") if v is not None else "--",
            }
        )
    if "frac_ge_1x" in has:
        cols.append(
            {
                "key": "frac_ge_1x",
                "label": "≥1×",
                "num": True,
                "fmt": lambda v, r: _esc(_fmt_pct(v)),
            }
        )
    if "frac_ge_10x" in has:
        cols.append(
            {
                "key": "frac_ge_10x",
                "label": "≥10×",
                "num": True,
                "fmt": lambda v, r: _esc(_fmt_pct(v)),
            }
        )
    if "global_methylation" in has:
        cols.append(
            {
                "key": "global_methylation",
                "label": "global meth",
                "num": True,
                "fmt": lambda v, r: _esc(_fmt_pct(v)),
            }
        )
    if "bisulfite_conversion_rate" in has:
        cols.append(
            {
                "key": "bisulfite_conversion_rate",
                "label": "conversion",
                "num": True,
                "fmt": lambda v, r: _esc(_fmt_pct(v, 2)),
            }
        )
    if "min_pairwise_corr" in has:
        cols.append(
            {
                "key": "min_pairwise_corr",
                "label": "min pairwise r",
                "num": True,
                "fmt": lambda v, r: _esc(f"{v:.2f}") if v is not None else "--",
            }
        )
    cols.append(
        {
            "key": "_status",
            "label": "status",
            "fmt": lambda v, r: f'<span class="badge {v}">{v}</span>',
        }
    )
    return _html_table(qc_rows, cols, "tQC")


def _dbeta_cell(v, r):
    if v is None or (isinstance(v, float) and v != v):
        return "--"
    cls = "cell-hyper" if v > 0 else "cell-hypo"
    sign = "+" if v > 0 else ""
    return f'<span class="{cls}">{sign}{v:.2f}</span>'


def _q_cell(alpha):
    def f(v, r):
        if v is None or (isinstance(v, float) and v != v):
            return "--"
        cls = "cell-sig" if v < alpha else ""
        return f'<span class="{cls}">{_fmt_sci(v)}</span>'

    return f


def _dir_cell(gmap_dir):
    def f(v, r):
        diff = r.get("meth_diff")
        if diff is None:
            return "--"
        return _chip("hyper" if diff > 0 else "hypo", "t" if diff > 0 else "c")

    return f


def _top_dmc_table(md: MethylData, n: int, alpha: float) -> str | None:
    dmc = md.dmc
    if dmc is None:
        return None
    p_col = "qvalue" if "qvalue" in dmc.columns else "pvalue"
    sub = dmc.filter(pl.col(p_col).is_not_nan()).sort(p_col).head(n)
    rows = list(sub.iter_rows(named=True))
    cols: list[dict] = [
        {"key": "chrom", "label": "chrom"},
        {"key": "pos", "label": "pos", "num": True, "fmt": lambda v, r: _esc(_fmt_int(v))},
        {"key": "meth_diff", "label": "Δβ", "num": True, "fmt": _dbeta_cell},
        {
            "key": p_col,
            "label": p_col[0] if p_col == "qvalue" else "p",
            "num": True,
            "fmt": _q_cell(alpha),
        },
        {"key": "_dir", "label": "dir", "fmt": _dir_cell(None)},
    ]
    for opt, lab in (
        ("feature_type", "feature"),
        ("gene_name", "gene"),
        ("cpg_context", "CpG ctx"),
    ):
        if opt in sub.columns:
            cols.append(
                {
                    "key": opt,
                    "label": lab,
                    "fmt": lambda v, r: (
                        _esc(v) if v not in (None, "") else '<span class="muted">—</span>'
                    ),
                }
            )
    return _html_table(rows, cols, "tDMC")


def _top_dmr_table(md: MethylData, n: int) -> str | None:
    dmr = md.uns.get("dmr")
    if dmr is None or not isinstance(dmr, pl.DataFrame) or dmr.is_empty():
        return None
    q_col = next(
        (
            c
            for c in ("qvalue", "combined_qvalue", "combined_pvalue", "empirical_qvalue")
            if c in dmr.columns
        ),
        None,
    )
    sub = dmr.sort(q_col).head(n) if q_col else dmr.head(n)
    rows = list(sub.iter_rows(named=True))

    def region(v, r):
        s, e = r.get("start"), r.get("end")
        c = r.get("chrom", "")
        if s is None or e is None:
            return _esc(c)
        return _esc(f"{c}:{int(s):,}-{int(e):,}")

    cols: list[dict] = [{"key": "chrom", "label": "region", "fmt": region}]
    if "n_cpgs" in sub.columns:
        cols.append({"key": "n_cpgs", "label": "CpGs", "num": True})
    diff_col = (
        "meth_diff"
        if "meth_diff" in sub.columns
        else ("mean_meth_diff" if "mean_meth_diff" in sub.columns else None)
    )
    if diff_col:
        cols.append({"key": diff_col, "label": "Δβ", "num": True, "fmt": _dbeta_cell})
    if q_col:
        cols.append({"key": q_col, "label": "q", "num": True, "fmt": _q_cell(1.0)})
    if "dmr_type" in sub.columns:
        cols.append(
            {
                "key": "dmr_type",
                "label": "type",
                "fmt": lambda v, r: _chip(v, "t" if v == "hyper" else "c"),
            }
        )
    for opt, lab in (("gene_name", "gene"), ("feature_type", "feature")):
        if opt in sub.columns:
            cols.append(
                {
                    "key": opt,
                    "label": lab,
                    "fmt": lambda v, r: (
                        _esc(v) if v not in (None, "") else '<span class="muted">—</span>'
                    ),
                }
            )
    return _html_table(rows, cols, "tDMR")


# ---------------------------------------------------------------------------
# Preprocessing / history
# ---------------------------------------------------------------------------


def _preproc_flow(md: MethylData) -> list[dict]:
    out: list[dict] = []
    prev = None
    for h in md.uns.get("_store_history", []) or []:
        n = h.get("n_sites")
        delta_pct = None
        delta_cls = "flat"
        if isinstance(n, int) and isinstance(prev, int) and prev > 0:
            delta_pct = (n - prev) / prev * 100.0
            delta_cls = "dn" if delta_pct < -0.05 else ("up" if delta_pct > 0.05 else "flat")
        out.append(
            {
                "step": h.get("step", "?"),
                "n_sites_str": _fmt_int(n) if isinstance(n, int) else "--",
                "delta": (
                    f"{delta_pct:+.1f}% sites" if delta_pct is not None else h.get("note", "")
                ),
                "delta_cls": delta_cls,
                "path": h.get("path", "?"),
            }
        )
        if isinstance(n, int):
            prev = n
    return out


# ---------------------------------------------------------------------------
# Summary facts + narrative + KPIs + completeness
# ---------------------------------------------------------------------------


def _dmc_stats(
    md: MethylData, alpha: float, min_abs_diff: float, dmc: pl.DataFrame | None = None
) -> dict:
    if dmc is None:
        dmc = md.dmc
    if dmc is None:
        return {"available": False}
    p_col = "qvalue" if "qvalue" in dmc.columns else "pvalue"
    diff = dmc["meth_diff"].to_numpy()
    pval = dmc[p_col].to_numpy()
    valid = ~np.isnan(pval) & ~np.isnan(diff)
    sig = valid & (pval < alpha) & (np.abs(diff) >= min_abs_diff)
    n_sig = int(sig.sum())
    med_abs = float(np.median(np.abs(diff[sig]))) if n_sig else float("nan")
    n_total = int(valid.sum())
    return {
        "available": True,
        "n_total": n_total,
        "n_sig": n_sig,
        "n_hyper": int((sig & (diff > 0)).sum()),
        "n_hypo": int((sig & (diff < 0)).sum()),
        "pct_sig": (n_sig / n_total) if n_total else 0.0,
        "median_abs_diff": med_abs,
        "p_col": p_col,
    }


def _dmr_stats(md: MethylData) -> dict:
    dmr = md.uns.get("dmr")
    if dmr is None or not isinstance(dmr, pl.DataFrame) or dmr.is_empty():
        return {"available": False}
    stats: dict = {"available": True, "n_total": len(dmr)}
    if "dmr_type" in dmr.columns:
        types = dmr["dmr_type"].to_list()
        stats["n_hyper"] = sum(1 for t in types if t == "hyper")
        stats["n_hypo"] = sum(1 for t in types if t == "hypo")
    if "n_cpgs" in dmr.columns:
        vals = dmr["n_cpgs"].drop_nulls().to_numpy()
        stats["median_cpgs"] = int(np.median(vals)) if len(vals) else None
    return stats


def _group_counts(md: MethylData) -> dict:
    if "group" not in md.obs.columns:
        return {}
    counts: dict = {}
    for g in md.obs.get_column("group").to_list():
        counts[g] = counts.get(g, 0) + 1
    return counts


def _global_meth_by_group(md: MethylData) -> dict:
    if "global_methylation" not in md.obs.columns or "group" not in md.obs.columns:
        return {}
    df = md.obs.select(["group", "global_methylation"]).drop_nulls()
    if df.is_empty():
        return {}
    agg = df.group_by("group").agg(pl.col("global_methylation").mean().alias("m"))
    return {r["group"]: r["m"] for r in agg.iter_rows(named=True)}


def _summary_narrative(md, facts, dmc_stats, dmr_stats) -> str:
    """Build a colour-coded HTML sentence describing the run."""
    gc = facts["group_counts"]
    grp_phrase = (
        ", ".join(f"{n} {_esc(g)}" for g, n in gc.items()) if gc else f"{md.n_samples} samples"
    )
    parts = [
        f"This report summarises a {_esc(md.context or 'methylation')} comparison of "
        f"<b>{md.n_samples} samples</b> ({grp_phrase}) across "
        f"<b>{facts['n_sites_str']} sites</b> retained after filtering."
    ]
    if dmc_stats.get("available"):
        parts.append(
            f" Differential testing ({_esc(facts['dmc_engine'])}, "
            f"{_esc(facts['fdr_method'])}-corrected) identified "
            f"<b>{dmc_stats['n_sig']:,} differentially methylated cytosines</b> "
            f"(q&lt;{facts['alpha']}, |Δβ|≥{facts['min_abs_diff']}) — "
            f"<span class='hy'>{dmc_stats['n_hyper']:,} hyper-</span> and "
            f"<span class='ho'>{dmc_stats['n_hypo']:,} hypo-methylated</span>."
        )
    if dmr_stats.get("available"):
        parts.append(f" These merged into <b>{dmr_stats['n_total']:,} DMRs</b>.")
    gm = facts.get("global_meth_by_group") or {}
    if len(gm) >= 1:
        gm_phrase = "; ".join(f"{_esc(g)} {_fmt_pct(v)}" for g, v in gm.items())
        parts.append(f" Genome-wide mean methylation: {gm_phrase}.")
    return "".join(parts)


def _build_facts(md, dmc_stats, dmr_stats, alpha, min_abs_diff) -> dict:
    n_sites = (
        md.uns.get("n_sites_filtered") or md.uns.get("n_sites_regions") or md.uns.get("n_sites_raw")
    )
    hist = md.uns.get("_store_history", []) or []
    n_raw = next((h.get("n_sites") for h in hist if isinstance(h.get("n_sites"), int)), None)
    n_final = (
        n_sites
        if isinstance(n_sites, int)
        else (
            next(
                (h.get("n_sites") for h in reversed(hist) if isinstance(h.get("n_sites"), int)),
                None,
            )
        )
    )
    dmc_uns = md.uns.get("dmc") or {}
    return {
        "n_sites_str": _fmt_int(n_final) if isinstance(n_final, int) else "?",
        "n_sites_raw_str": _fmt_int(n_raw) if isinstance(n_raw, int) else None,
        "pct_retained": (n_final / n_raw)
        if (isinstance(n_final, int) and isinstance(n_raw, int) and n_raw)
        else None,
        "group_counts": _group_counts(md),
        "global_meth_by_group": _global_meth_by_group(md),
        "dmc_engine": str(
            dmc_uns.get("test_used") or dmc_uns.get("test_requested") or "per-CpG test"
        ),
        "fdr_method": str(dmc_uns.get("fdr_method") or "BH").replace("fdr_", "").upper(),
        "alpha": alpha,
        "min_abs_diff": min_abs_diff,
    }


def _build_kpis(md, facts, dmc_stats, dmr_stats) -> list[dict]:
    kpis: list[dict] = [
        {
            "big": str(md.n_samples),
            "label": "samples",
            "cls": "",
            "ctx": ", ".join(f"{n} {_esc(g)}" for g, n in facts["group_counts"].items()) or None,
        }
    ]
    if facts["n_sites_str"] != "?":
        ctx = None
        if facts["pct_retained"] is not None and facts["n_sites_raw_str"]:
            ctx = f"of {facts['n_sites_raw_str']} raw ({facts['pct_retained']:.0%})"
        kpis.append({"big": facts["n_sites_str"], "label": "sites retained", "cls": "", "ctx": ctx})
    if dmc_stats.get("available"):
        total = max(dmc_stats["n_hyper"] + dmc_stats["n_hypo"], 1)
        kpis.append(
            {
                "big": f"{dmc_stats['n_hyper']:,}",
                "label": "hyper DMCs",
                "cls": "hyper",
                "spark": {"pct": dmc_stats["n_hyper"] / total * 100, "color": "var(--hyper)"},
            }
        )
        kpis.append(
            {
                "big": f"{dmc_stats['n_hypo']:,}",
                "label": "hypo DMCs",
                "cls": "hypo",
                "spark": {"pct": dmc_stats["n_hypo"] / total * 100, "color": "var(--hypo)"},
            }
        )
        kpis.append(
            {
                "big": _fmt_pct(dmc_stats["pct_sig"], 2),
                "label": "% sig of tested",
                "cls": "",
                "ctx": f"{dmc_stats['n_sig']:,} significant",
            }
        )
    if dmr_stats.get("available"):
        ctx = None
        if "n_hyper" in dmr_stats:
            ctx = f"{dmr_stats['n_hyper']} hyper · {dmr_stats['n_hypo']} hypo"
        kpis.append({"big": str(dmr_stats["n_total"]), "label": "DMRs", "cls": "", "ctx": ctx})
    return kpis


def _completeness(md, dmc_stats, dmr_stats) -> list[dict]:
    has = set(md.obs.columns)
    annotated = md.dmc is not None and "feature_type" in (
        md.dmc.columns if md.dmc is not None else []
    )
    return [
        {"label": "Import & methylstore built", "done": bool(md.uns.get("_store_history"))},
        {
            "label": "Filtering / normalisation",
            "done": md._filtered if hasattr(md, "_filtered") else bool(md.uns.get("filter")),
        },
        {"label": "DMC calling", "done": dmc_stats.get("available", False)},
        {"label": "DMR calling", "done": dmr_stats.get("available", False)},
        {"label": "Genomic annotation", "done": annotated},
        {"label": "QC metrics", "done": "global_methylation" in has},
        {
            "label": "Bisulfite conversion",
            "done": "bisulfite_conversion_rate" in has,
            "note": "no CHH store" if "bisulfite_conversion_rate" not in has else "",
        },
        {"label": "Sample correlation", "done": md.uns.get("qc_sample_correlation") is not None},
    ]


# ---------------------------------------------------------------------------
# Methods & citations
# ---------------------------------------------------------------------------


def _methods_text(md, facts, alpha, min_abs_diff) -> str:
    try:
        from . import __version__ as version
    except Exception:
        version = "unknown"
    dmc_uns = md.uns.get("dmc") or {}
    dmr_uns = md.uns.get("dmr_params") or {}
    filt = md.uns.get("filter") or {}
    sentences = [
        f"Methylation calls were imported and stored as a partitioned Parquet "
        f"methylstore using epykit v{version}."
    ]
    if filt:
        mc = filt.get("min_coverage") or filt.get("lo_count")
        if mc is not None:
            sentences.append(f"CpGs were retained at a minimum coverage of {mc}.")
    if dmc_uns:
        engine = dmc_uns.get("test_used") or dmc_uns.get("test_requested") or "a per-CpG test"
        fdr = str(dmc_uns.get("fdr_method") or "fdr_bh")
        sentences.append(
            f"Differential methylation was tested per CpG using the '{engine}' engine; "
            f"p-values were corrected with the {fdr} procedure. Cytosines with "
            f"q < {alpha} and |Δβ| ≥ {min_abs_diff} were called differential."
        )
    if dmr_uns:
        method = dmr_uns.get("method") or "the configured"
        sentences.append(f"Regions were called with the '{method}' DMR caller.")
    if md.dmc is not None and "feature_type" in md.dmc.columns:
        sentences.append(
            "Differential sites were annotated against gene features and CpG-island context."
        )
    sentences.append("Analyses used polars, NumPy/SciPy, statsmodels and bioframe.")
    return " ".join(sentences)


def _params_rows(md, alpha, min_abs_diff) -> list[dict]:
    dmc_uns = md.uns.get("dmc") or {}
    dmr_uns = md.uns.get("dmr_params") or {}
    rows = []
    if dmc_uns.get("test_used"):
        rows.append({"k": "DMC engine", "v": dmc_uns["test_used"]})
    if dmc_uns.get("fdr_method"):
        rows.append({"k": "Multiple testing", "v": dmc_uns["fdr_method"]})
    rows.append({"k": "Significance threshold", "v": f"q < {alpha}, |Δβ| ≥ {min_abs_diff}"})
    if dmr_uns.get("method"):
        rows.append({"k": "DMR caller", "v": dmr_uns["method"]})
    return rows


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------


def _serialisable(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, pl.DataFrame):
        return f"<DataFrame: {len(value)} rows x {value.width} cols>"
    if isinstance(value, dict):
        return {k: _serialisable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serialisable(v) for v in value]
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


# ---------------------------------------------------------------------------
# Figure rendering with self-contained embed control
# ---------------------------------------------------------------------------


def _make_fig_renderer(self_contained: bool) -> Callable:
    """Return a ``render(fig) -> html|None`` that, when ``self_contained``,
    inlines the full Plotly bundle into the *first* figure only and references
    it from the rest -- so the file works offline without duplicating ~3 MB
    per figure. When not self-contained, every figure pulls Plotly from a CDN.
    """
    state = {"embedded": False}

    def render(fig):
        if fig is None:
            return None
        if self_contained:
            inc: Any = not state["embedded"]
            if inc:
                state["embedded"] = True
        else:
            inc = "cdn"
        try:
            # responsive=True makes each chart size to its container (and resize
            # with the window) instead of a fixed ~700px width that overflows the
            # two-column grid and overlaps the neighbouring card.
            return fig.to_html(
                include_plotlyjs=inc,
                full_html=False,
                default_height=None,
                config={"responsive": True},
            )
        except Exception as exc:  # pragma: no cover - plotly version drift
            logger.warning("Failed to render Plotly figure: %s", exc)
            return None

    return render


def _safe(fn, *args, **kwargs):
    """Call a Plotly-twin function; swallow exceptions so a single failing
    section (e.g. PCA needs >= 2 samples) never aborts the whole report."""
    try:
        return fn(*args, **kwargs)
    except Exception as exc:
        logger.info("report: section skipped (%s): %s", getattr(fn, "__name__", fn), exc)
        return None


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def generate_report(
    md: MethylData,
    output: str,
    *,
    title: str | None = None,
    gtf_path: str | None = None,
    alpha: float = 0.05,
    min_abs_diff: float = 0.1,
    dmc_top_n: int = 50,
    dmr_top_n: int = 50,
    metaplot_max_genes: int | None = 5000,
    pca_n_sites: int = 10_000,
    coverage_max_points: int = 200_000,
    dmc_max_points: int = 200_000,
    clear_cache: bool = False,
    self_contained: bool = True,
    qc_thresholds: dict | None = None,
) -> str:
    """Render a single-file, MultiQC-style HTML dashboard report.

    Parameters
    ----------
    md : MethylData
        The data object. Whatever subset of pp / tl steps has been run is
        what gets rendered; unrun sections render a short notice and a hollow
        sidebar status dot.
    output : str
        Output HTML file path. Parent directory is created if needed.
    title : str, optional
        Report title. Defaults to ``"epykit report -- <assembly>"``.
    gtf_path : str, optional
        If supplied, a TSS metaplot section is rendered using this GTF.
    alpha : float
        DMC q-value threshold for KPI counts. Default 0.05.
    min_abs_diff : float
        Minimum |meth_diff| for KPI counts. Default 0.1.
    self_contained : bool
        When True (default) the full Plotly bundle is embedded inline so the
        resulting HTML works offline (larger file, ~3-4 MB). When False,
        Plotly is loaded from a CDN (smaller file, needs internet).
    qc_thresholds : dict, optional
        Override the default QC pass/warn/fail thresholds. Keys:
        ``bisulfite_conversion_rate``, ``frac_ge_10x``, ``mean_coverage``,
        ``min_pairwise_corr``; each value is a ``(pass, warn)`` tuple.
    pca_n_sites : int
        Cap for sites entering PCA. Lower (e.g. 5_000) to shave RAM/time on
        huge stores; the visual is virtually unchanged below ~20_000.
    coverage_max_points : int
        Cap for points entering the coverage histogram.
    dmc_max_points : int
        Cap for points rendered in the volcano / MA / Manhattan scatter
        figures. ALL significant CpGs are always kept; only the
        non-significant background is subsampled. Bounds memory and HTML size
        on whole-genome tables (tens of millions of CpGs) -- without this a
        genome-wide self-contained report can exhaust RAM. Default 200_000.
    clear_cache : bool
        If True, drop cached compute results on ``md.uns['_report_cache']``
        before rendering. Use after re-running upstream steps so a stale PCA /
        metaplot doesn't survive into the new report.

    Returns
    -------
    str
        Absolute path of the written HTML file.
    """
    Environment, FileSystemLoader, select_autoescape = _require_jinja()

    from .pl._plotly import (
        coverage_histogram_plotly,
        cpg_island_pie_plotly,
        dmr_size_hist_plotly,
        feature_direction_stacked_plotly,
        feature_pie_plotly,
        global_methylation_bar_plotly,
        ma_plot_plotly,
        manhattan_plotly,
        pca_plotly,
        pvalue_histogram_plotly,
        sample_correlation_plotly,
        scree_plotly,
        tss_metaplot_plotly,
        volcano_plotly,
    )

    if clear_cache:
        from .pl._compute import clear_report_cache

        clear_report_cache(md)

    thresholds = dict(_DEFAULT_QC_THRESHOLDS)
    if qc_thresholds:
        thresholds.update(qc_thresholds)

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template("report.html.j2")
    css_inline = (_TEMPLATE_DIR / "report.css").read_text(encoding="utf-8")
    js_inline = (_TEMPLATE_DIR / "report.js").read_text(encoding="utf-8")

    try:
        from . import __version__ as version
    except Exception:
        version = "unknown"

    render = _make_fig_renderer(self_contained)

    # Resolve the *full* per-CpG DMC table for genome-wide stats and figures.
    # ``md.dmc`` prefers the annotated table, which may have been annotated on
    # only the significant subset -- using it as the denominator for "% sig" or
    # as the volcano cloud would be wrong. The annotated table is still used
    # for the annotation pies and the top-DMC table.
    try:
        dmc_full = md.get_dmc(annotated=False)
    except Exception:
        dmc_full = None
    if dmc_full is None or (hasattr(dmc_full, "is_empty") and dmc_full.is_empty()):
        dmc_full = md.dmc

    dmc_stats = _dmc_stats(md, alpha=alpha, min_abs_diff=min_abs_diff, dmc=dmc_full)
    dmr_stats = _dmr_stats(md)
    facts = _build_facts(md, dmc_stats, dmr_stats, alpha, min_abs_diff)

    qc_rows, qc_legend = _qc_rows(md, thresholds)
    qc_has_fail = any(r["_status"] == "fail" for r in qc_rows)
    qc_has_warn = any(r["_status"] == "warn" for r in qc_rows)

    # Figures (each guarded; None -> section shows skip state).
    coverage_plot = render(_safe(coverage_histogram_plotly, md, max_points=coverage_max_points))
    global_meth_bar = render(_safe(global_methylation_bar_plotly, md))
    corr_heatmap = render(_safe(sample_correlation_plotly, md))
    volcano_plot = (
        render(
            _safe(
                volcano_plotly,
                md,
                alpha=alpha,
                min_abs_diff=min_abs_diff,
                dmc=dmc_full,
                max_points=dmc_max_points,
            )
        )
        if dmc_stats.get("available")
        else None
    )
    pvalue_hist = (
        render(_safe(pvalue_histogram_plotly, md, dmc=dmc_full))
        if dmc_stats.get("available")
        else None
    )
    ma_plot = (
        render(
            _safe(
                ma_plot_plotly,
                md,
                alpha=alpha,
                min_abs_diff=min_abs_diff,
                dmc=dmc_full,
                max_points=dmc_max_points,
            )
        )
        if dmc_stats.get("available")
        else None
    )
    manhattan_plot = (
        render(_safe(manhattan_plotly, md, alpha=alpha, dmc=dmc_full, max_points=dmc_max_points))
        if dmc_stats.get("available")
        else None
    )
    dmr_size_hist = render(_safe(dmr_size_hist_plotly, md)) if dmr_stats.get("available") else None
    # Annotation pies at BOTH levels. Per-CpG (DMC) answers "where do
    # differential cytosines fall?" (density-weighted); per-region (DMR) gives
    # the field-standard "what fraction of DMRs hit each feature?". Each returns
    # None when its table/column is absent, so missing levels just don't render.
    feature_pie = render(_safe(feature_pie_plotly, md, level="dmc"))
    cpg_pie = render(_safe(cpg_island_pie_plotly, md, level="dmc"))
    feature_pie_dmr = render(_safe(feature_pie_plotly, md, level="dmr"))
    cpg_pie_dmr = render(_safe(cpg_island_pie_plotly, md, level="dmr"))
    feature_stacked = render(_safe(feature_direction_stacked_plotly, md))
    metaplot = (
        render(_safe(tss_metaplot_plotly, md, gtf_path, max_genes=metaplot_max_genes))
        if gtf_path
        else None
    )
    pca_plot = render(_safe(pca_plotly, md, n_sites=pca_n_sites))
    scree_plot = render(_safe(scree_plotly, md))

    # Sidebar section statuses (ok / warn / skip).
    def st(ok, warn=False):
        return "warn" if warn else ("ok" if ok else "skip")

    sections = [
        {"id": "summary", "num": "00", "label": "Summary", "status": "ok"},
        {"id": "samples", "num": "01", "label": "Samples", "status": "ok"},
        {
            "id": "preproc",
            "num": "02",
            "label": "Preprocessing",
            "status": st(bool(md.uns.get("_store_history"))),
        },
        {
            "id": "qc",
            "num": "03",
            "label": "Quality control",
            "status": st(bool(qc_rows), warn=(qc_has_fail or qc_has_warn)),
        },
        {
            "id": "dmc",
            "num": "04",
            "label": "Diff. methylation (DMC)",
            "status": st(dmc_stats.get("available", False)),
        },
        {
            "id": "dmr",
            "num": "05",
            "label": "Diff. regions (DMR)",
            "status": st(dmr_stats.get("available", False)),
        },
        {
            "id": "annot",
            "num": "06",
            "label": "Annotation",
            "status": st(bool(feature_pie or cpg_pie or feature_pie_dmr or cpg_pie_dmr)),
        },
        {"id": "metaplot", "num": "07", "label": "TSS metaplot", "status": st(bool(metaplot))},
        {"id": "pca", "num": "08", "label": "Sample similarity", "status": st(bool(pca_plot))},
        {"id": "methods", "num": "09", "label": "Methods & citations", "status": "ok"},
        {"id": "prov", "num": "10", "label": "Provenance", "status": "ok"},
    ]

    provenance_payload = {
        "epykit_version": version,
        "assembly": md.assembly,
        "context": md.context,
        "state": list(md.state),
        "store": md.store,
        "n_samples": md.n_samples,
        "treatment_ids": md.treatment_ids,
        "control_ids": md.control_ids,
        "dmc_uns": _serialisable(md.uns.get("dmc")),
        "dmr_params": _serialisable(md.uns.get("dmr_params")),
        "annotation": _serialisable(md.uns.get("annotation")),
        "regions": _serialisable(md.uns.get("regions")),
        "filter": _serialisable(md.uns.get("filter")),
        "unite": _serialisable(md.uns.get("unite")),
        "alpha": alpha,
        "min_abs_diff": min_abs_diff,
    }
    provenance_rows = [
        {"k": "epykit version", "v": version},
        {"k": "assembly / context", "v": f"{md.assembly} / {md.context}"},
        {"k": "store", "v": md.store},
        {"k": "state", "v": " → ".join(md.state) if md.state else "raw"},
        {"k": "samples", "v": md.n_samples},
    ]

    now = datetime.datetime.now()
    ctx = {
        "title": title or f"epykit report -- {md.assembly}",
        "css_inline": css_inline,
        "js_inline": js_inline,
        "epykit_version": version,
        "generated_at": now.isoformat(timespec="seconds"),
        "generated_short": now.strftime("%Y-%m-%d %H:%M"),
        "assembly": md.assembly,
        "context": md.context,
        "n_samples": md.n_samples,
        "n_sites": facts["n_sites_str"],
        "state": list(md.state),
        "alpha": alpha,
        "min_abs_diff": min_abs_diff,
        "sections": sections,
        # Summary
        "summary_narrative": _summary_narrative(md, facts, dmc_stats, dmr_stats),
        "kpis": _build_kpis(md, facts, dmc_stats, dmr_stats),
        "completeness": _completeness(md, dmc_stats, dmr_stats),
        # Samples
        "samples_table": _samples_table(md),
        "group_counts": facts["group_counts"],
        "group_class_map": _group_class_map(md),
        # Preprocessing
        "preproc_flow": _preproc_flow(md),
        "filter_params": _pretty_dict(md.uns.get("filter")),
        "regions_params": _pretty_dict(md.uns.get("regions")),
        # QC
        "qc_table": _qc_table(md, qc_rows),
        "qc_legend": qc_legend,
        "coverage_plot": coverage_plot,
        "global_meth_bar": global_meth_bar,
        "corr_heatmap": corr_heatmap,
        # DMC
        "dmc_available": dmc_stats.get("available", False),
        "dmc_n_total": _fmt_int(dmc_stats.get("n_total")),
        "dmc_n_sig": _fmt_int(dmc_stats.get("n_sig")),
        "dmc_n_hyper": _fmt_int(dmc_stats.get("n_hyper")),
        "dmc_n_hypo": _fmt_int(dmc_stats.get("n_hypo")),
        "dmc_median_abs": (
            _fmt_value(dmc_stats["median_abs_diff"])
            if dmc_stats.get("available")
            and dmc_stats["median_abs_diff"] == dmc_stats["median_abs_diff"]
            else None
        ),
        "dmc_engine": facts["dmc_engine"],
        "dmc_fdr": facts["fdr_method"],
        "volcano_plot": volcano_plot,
        "pvalue_hist": pvalue_hist,
        "ma_plot": ma_plot,
        "manhattan_plot": manhattan_plot,
        "dmc_top_table": _top_dmc_table(md, dmc_top_n, alpha)
        if dmc_stats.get("available")
        else None,
        "dmc_top_n": dmc_top_n,
        # DMR
        "dmr_available": dmr_stats.get("available", False),
        "dmr_n_total": _fmt_int(dmr_stats.get("n_total")),
        "dmr_n_hyper": _fmt_int(dmr_stats.get("n_hyper")) if "n_hyper" in dmr_stats else None,
        "dmr_n_hypo": _fmt_int(dmr_stats.get("n_hypo")) if "n_hyper" in dmr_stats else None,
        "dmr_median_cpgs": dmr_stats.get("median_cpgs"),
        "dmr_caller": (md.uns.get("dmr_params") or {}).get("method"),
        "dmr_size_hist": dmr_size_hist,
        "dmr_top_table": _top_dmr_table(md, dmr_top_n) if dmr_stats.get("available") else None,
        "dmr_top_n": dmr_top_n,
        # Annotation
        "feature_pie": feature_pie,
        "cpg_pie": cpg_pie,
        "feature_pie_dmr": feature_pie_dmr,
        "cpg_pie_dmr": cpg_pie_dmr,
        "feature_stacked": feature_stacked,
        # Metaplot / PCA
        "metaplot": metaplot,
        "pca_plot": pca_plot,
        "scree_plot": scree_plot,
        # Methods
        "methods_text": _methods_text(md, facts, alpha, min_abs_diff),
        "params_rows": _params_rows(md, alpha, min_abs_diff),
        # Provenance
        "provenance_rows": provenance_rows,
        "provenance": json.dumps(provenance_payload, indent=2, default=str),
    }

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(template.render(**ctx), encoding="utf-8")
    logger.info("Report written: %s", out)
    return str(out.resolve())


def _pretty_dict(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, dict):
        return json.dumps(_serialisable(value), indent=2, default=str)
    return str(value)


__all__ = ["generate_report"]
