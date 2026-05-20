"""HTML report generator for a MethylData object.

Renders an interactive single-file HTML report covering sample metadata,
preprocessing trail, QC, DMC/DMR, annotation, optional TSS metaplot, PCA
and provenance. Each section is conditional on the relevant data being
populated on ``md``; unrun sections render a short "not run yet" note
instead of crashing.

Optional deps: ``jinja2`` and ``plotly`` (install via
``pip install 'epykit[report]'``). Both are imported lazily so the rest
of the package keeps working with the base install.
"""

from __future__ import annotations

import datetime
import json
import logging
from pathlib import Path
from typing import Any, Optional

import numpy as np
import polars as pl

from .methyldata import MethylData

logger = logging.getLogger(__name__)

_TEMPLATE_DIR = Path(__file__).parent / "templates"


def _require_jinja():
    try:
        from jinja2 import Environment, FileSystemLoader, select_autoescape
        return Environment, FileSystemLoader, select_autoescape
    except ImportError as exc:
        raise ImportError(
            "jinja2 is required for HTML report generation. "
            "Install with: pip install 'epykit[report]'"
        ) from exc


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


def _df_to_html_table(df: pl.DataFrame, max_rows: int = 50) -> str:
    """Render a Polars DataFrame as a compact HTML table."""
    if df is None or len(df) == 0:
        return '<div class="skipped">no rows</div>'
    head = df.head(max_rows)
    cols = head.columns
    header_html = "".join(f"<th>{c}</th>" for c in cols)
    rows_html = []
    for row in head.iter_rows(named=True):
        cells = "".join(f"<td>{_fmt_value(row.get(c))}</td>" for c in cols)
        rows_html.append(f"<tr>{cells}</tr>")
    return (
        '<table class="df"><thead><tr>'
        + header_html
        + "</tr></thead><tbody>"
        + "".join(rows_html)
        + "</tbody></table>"
    )


def _samples_table(md: MethylData) -> str:
    return _df_to_html_table(md.obs, max_rows=200)


def _qc_summary(md: MethylData) -> Optional[str]:
    """Concise per-sample QC table, if QC has been run."""
    qc_cols = [
        c for c in ("sample_id", "group", "treatment",
                    "global_methylation", "mean_coverage",
                    "frac_ge_1x", "frac_ge_10x",
                    "bisulfite_conversion_rate", "low_coverage_flag")
        if c in md.obs.columns
    ]
    if not qc_cols or "global_methylation" not in md.obs.columns:
        return None
    return _df_to_html_table(md.obs.select(qc_cols), max_rows=200)


def _history_entries(md: MethylData) -> list[dict]:
    out = []
    for h in md.uns.get("_store_history", []) or []:
        out.append({
            "step": h.get("step", "?"),
            "path": h.get("path", "?"),
            "n_sites_str": (
                f"{h.get('n_sites'):,}" if isinstance(h.get("n_sites"), int)
                else "--"
            ),
        })
    return out


def _dmc_stats(md: MethylData, alpha: float, min_abs_diff: float) -> dict:
    """Compute KPI counts on the current DMC table."""
    dmc = md.dmc
    if dmc is None:
        return {"available": False}
    p_col = "qvalue" if "qvalue" in dmc.columns else "pvalue"
    diff = dmc["meth_diff"].to_numpy()
    pval = dmc[p_col].to_numpy()
    valid = ~np.isnan(pval) & ~np.isnan(diff)
    sig = valid & (pval < alpha) & (np.abs(diff) >= min_abs_diff)
    hyper = int((sig & (diff > 0)).sum())
    hypo = int((sig & (diff < 0)).sum())
    return {
        "available": True,
        "n_total": int(valid.sum()),
        "n_sig": int(sig.sum()),
        "n_hyper": hyper,
        "n_hypo": hypo,
    }


def _top_dmcs(md: MethylData, n: int = 50) -> Optional[pl.DataFrame]:
    dmc = md.dmc
    if dmc is None:
        return None
    p_col = "qvalue" if "qvalue" in dmc.columns else "pvalue"
    cols_pref = [
        "chrom", "pos", "strand", "meth_diff", "mean_beta_case",
        "mean_beta_control", p_col, "feature_type", "gene_name", "cpg_context",
    ]
    cols = [c for c in cols_pref if c in dmc.columns]
    return (
        dmc.filter(pl.col(p_col).is_not_nan())
        .sort(p_col)
        .head(n)
        .select(cols)
    )


def _top_dmrs(md: MethylData, n: int = 50) -> tuple[Optional[pl.DataFrame], dict]:
    dmr = md.uns.get("dmr")
    if dmr is None or not isinstance(dmr, pl.DataFrame) or len(dmr) == 0:
        return None, {"available": False}
    q_col = (
        "qvalue" if "qvalue" in dmr.columns
        else "combined_qvalue" if "combined_qvalue" in dmr.columns
        else "combined_pvalue" if "combined_pvalue" in dmr.columns
        else None
    )
    stats: dict = {"available": True, "n_total": len(dmr)}
    if "dmr_type" in dmr.columns:
        types = dmr["dmr_type"].to_list()
        stats["n_hyper"] = sum(1 for t in types if t == "hyper")
        stats["n_hypo"] = sum(1 for t in types if t == "hypo")
    if q_col is None:
        return dmr.head(n), stats
    pref = [
        "chrom", "start", "end", "n_cpgs", "meth_diff", "mean_meth_diff",
        q_col, "dmr_type", "feature_type", "gene_name",
    ]
    cols = [c for c in pref if c in dmr.columns]
    return dmr.sort(q_col).head(n).select(cols), stats


def _fig_html(fig) -> Optional[str]:
    if fig is None:
        return None
    try:
        return fig.to_html(include_plotlyjs="cdn", full_html=False, default_height=None)
    except Exception as exc:  # pragma: no cover - plotly version drift
        logger.warning("Failed to render Plotly figure: %s", exc)
        return None


def _safe(fn, *args, **kwargs):
    """Call a Plotly-twin function; swallow exceptions to keep the report
    rendering even if one section fails (e.g. PCA needs >=2 samples)."""
    try:
        return fn(*args, **kwargs)
    except Exception as exc:
        logger.info("report: section skipped (%s): %s", fn.__name__, exc)
        return None


def generate_report(
    md: MethylData,
    output: str,
    *,
    title: Optional[str] = None,
    gtf_path: Optional[str] = None,
    alpha: float = 0.05,
    min_abs_diff: float = 0.1,
    dmc_top_n: int = 50,
    dmr_top_n: int = 50,
    metaplot_max_genes: Optional[int] = 5000,
    pca_n_sites: int = 10_000,
    coverage_max_points: int = 200_000,
    clear_cache: bool = False,
) -> str:
    """Render a single-file HTML report.

    Parameters
    ----------
    md : MethylData
        The data object. Whatever subset of pp / tl steps has been run is
        what gets rendered; unrun sections render a short notice.
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
    pca_n_sites : int
        Cap for sites entering PCA. Lower (e.g. 5_000) to shave RAM/time
        on huge stores; the visual is virtually unchanged below ~20_000.
    coverage_max_points : int
        Cap for points entering the coverage histogram. The store is
        subsampled deterministically beyond this limit.
    clear_cache : bool
        If True, drop any cached compute results on ``md.uns['_report_cache']``
        before rendering. Use after re-running upstream steps so a stale
        PCA / metaplot doesn't survive into the new report.

    Returns
    -------
    str
        Absolute path of the written HTML file.
    """
    Environment, FileSystemLoader, select_autoescape = _require_jinja()

    from .pl._plotly import (
        coverage_histogram_plotly,
        volcano_plotly,
        ma_plot_plotly,
        manhattan_plotly,
        feature_pie_plotly,
        cpg_island_pie_plotly,
        pca_plotly,
        tss_metaplot_plotly,
    )
    if clear_cache:
        from .pl._compute import clear_report_cache
        clear_report_cache(md)

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template("report.html.j2")

    css_inline = (_TEMPLATE_DIR / "report.css").read_text(encoding="utf-8")

    try:
        from . import __version__ as version
    except Exception:
        version = "unknown"

    n_sites = (
        md.uns.get("n_sites_filtered")
        or md.uns.get("n_sites_regions")
        or md.uns.get("n_sites_raw")
        or "?"
    )

    dmc_stats = _dmc_stats(md, alpha=alpha, min_abs_diff=min_abs_diff)
    top_dmc_df = _top_dmcs(md, n=dmc_top_n)
    top_dmr_df, dmr_stats = _top_dmrs(md, n=dmr_top_n)

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

    ctx = {
        "title": title or f"epykit report -- {md.assembly}",
        "css_inline": css_inline,
        "epykit_version": version,
        "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "assembly": md.assembly,
        "context": md.context,
        "n_samples": md.n_samples,
        "n_sites": f"{n_sites:,}" if isinstance(n_sites, int) else str(n_sites),
        "state": ", ".join(md.state) if md.state else None,
        "alpha": alpha,
        "samples_table": _samples_table(md),
        "history": _history_entries(md),
        "filter_params": _pretty_dict(md.uns.get("filter")),
        "regions_params": _pretty_dict(md.uns.get("regions")),
        "qc_summary": _qc_summary(md),
        "coverage_plot": _fig_html(_safe(coverage_histogram_plotly, md, max_points=coverage_max_points)),
        "dmc_available": dmc_stats.get("available", False),
        "dmc_n_total": _fmt_value(dmc_stats.get("n_total")),
        "dmc_n_sig": _fmt_value(dmc_stats.get("n_sig")),
        "dmc_n_hyper": _fmt_value(dmc_stats.get("n_hyper")),
        "dmc_n_hypo": _fmt_value(dmc_stats.get("n_hypo")),
        "volcano_plot": _fig_html(_safe(volcano_plotly, md, alpha=alpha, min_abs_diff=min_abs_diff))
            if dmc_stats.get("available") else None,
        "ma_plot": _fig_html(_safe(ma_plot_plotly, md, alpha=alpha, min_abs_diff=min_abs_diff))
            if dmc_stats.get("available") else None,
        "manhattan_plot": _fig_html(_safe(manhattan_plotly, md, alpha=alpha))
            if dmc_stats.get("available") else None,
        "dmc_top_table": _df_to_html_table(top_dmc_df, max_rows=dmc_top_n) if top_dmc_df is not None else None,
        "dmc_top_n": dmc_top_n,
        "dmr_available": dmr_stats.get("available", False),
        "dmr_n_total": _fmt_value(dmr_stats.get("n_total")),
        "dmr_n_hyper": _fmt_value(dmr_stats.get("n_hyper")) if "n_hyper" in dmr_stats else None,
        "dmr_n_hypo": _fmt_value(dmr_stats.get("n_hypo")) if "n_hyper" in dmr_stats else None,
        "dmr_top_table": _df_to_html_table(top_dmr_df, max_rows=dmr_top_n) if top_dmr_df is not None else None,
        "dmr_top_n": dmr_top_n,
        "feature_pie": _fig_html(_safe(feature_pie_plotly, md)),
        "cpg_pie": _fig_html(_safe(cpg_island_pie_plotly, md)),
        "metaplot": (
            _fig_html(_safe(tss_metaplot_plotly, md, gtf_path,
                            max_genes=metaplot_max_genes))
            if gtf_path else None
        ),
        "pca_plot": _fig_html(_safe(pca_plotly, md, n_sites=pca_n_sites)),
        "provenance": json.dumps(provenance_payload, indent=2, default=str),
    }

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(template.render(**ctx), encoding="utf-8")
    logger.info("Report written: %s", out)
    return str(out.resolve())


def _serialisable(value: Any) -> Any:
    """Render a value into something json-friendly."""
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


def _pretty_dict(value: Any) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, dict):
        return json.dumps(_serialisable(value), indent=2, default=str)
    return str(value)


__all__ = ["generate_report"]
