"""MultiQC custom-content emitter .

Writes ``*_mqc.json`` files into ``output_dir`` so MultiQC's
``--config``-driven custom-content scanner picks them up. Each metric
becomes its own file with the
``{id, section_name, plot_type, data}`` schema. No new dependencies --
just JSON serialisation.

Reference: https://multiqc.info/docs/custom_content/#a-custom-content-file-with-static-data
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import polars as pl

logger = logging.getLogger(__name__)


def _write_json(out: Path, payload: dict) -> Path:
    out.write_text(json.dumps(payload, indent=2, default=str))
    return out


def report_multiqc(md, output_dir: str) -> str:
    """Write per-metric ``*_mqc.json`` files for MultiQC pickup.

    Returns the output directory. Files written depend on what's
    populated on ``md.obs`` / ``md.uns`` -- missing metrics are simply
    skipped (no empty stubs).
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    obs = md.obs
    written: list[Path] = []

    samples = obs.get_column("sample_id").to_list()

    if "bisulfite_conversion_rate" in obs.columns:
        vals = obs.get_column("bisulfite_conversion_rate").to_numpy()
        payload: dict[str, Any] = {
            "id": "epykit_conversion_rate",
            "section_name": "epykit: Bisulfite conversion rate",
            "description": "1 - mean CHH methylation per sample.",
            "plot_type": "bargraph",
            "data": {
                s: {"conversion_rate": (float(v) if v is not None else None)}
                for s, v in zip(samples, vals)
            },
        }
        written.append(_write_json(
            out / "epykit_conversion_rate_mqc.json", payload,
        ))

    if "mean_coverage" in obs.columns:
        cov = obs.get_column("mean_coverage").to_numpy()
        payload = {
            "id": "epykit_coverage",
            "section_name": "epykit: Mean coverage",
            "description": "Genome-wide mean per-CpG coverage per sample.",
            "plot_type": "bargraph",
            "data": {
                s: {"mean_coverage": float(v) if v is not None else None}
                for s, v in zip(samples, cov)
            },
        }
        written.append(_write_json(
            out / "epykit_coverage_mqc.json", payload,
        ))

    if "global_methylation" in obs.columns:
        gm = obs.get_column("global_methylation").to_numpy()
        payload = {
            "id": "epykit_global_methylation",
            "section_name": "epykit: Global methylation",
            "description": "Genome-wide CpG mean beta per sample.",
            "plot_type": "bargraph",
            "data": {
                s: {"global_methylation": float(v) if v is not None else None}
                for s, v in zip(samples, gm)
            },
        }
        written.append(_write_json(
            out / "epykit_global_methylation_mqc.json", payload,
        ))

    if "qc_sample_correlation" in md.uns and isinstance(
        md.uns["qc_sample_correlation"], pl.DataFrame
    ):
        df = md.uns["qc_sample_correlation"]
        nodes = sorted(set(
            df.get_column("sample_a").to_list()
            + df.get_column("sample_b").to_list()
        ))
        data: dict[str, dict[str, float]] = {n: {} for n in nodes}
        for row in df.iter_rows(named=True):
            data[row["sample_a"]][row["sample_b"]] = float(row["correlation"])
        payload = {
            "id": "epykit_sample_correlation",
            "section_name": "epykit: Sample correlation",
            "description": "All-vs-all sample correlation matrix.",
            "plot_type": "heatmap",
            "data": data,
        }
        written.append(_write_json(
            out / "epykit_sample_correlation_mqc.json", payload,
        ))

    # DMC / DMR summary
    if "dmc" in md.uns and isinstance(md.uns["dmc"], dict):
        info = md.uns["dmc"]
        payload = {
            "id": "epykit_dmc_summary",
            "section_name": "epykit: DMC summary",
            "description": "DMC engine settings and call counts.",
            "plot_type": "table",
            "data": {
                "cohort": {
                    "test_used": info.get("test_used"),
                    "n_sites": info.get("n_sites"),
                    "unite": info.get("unite"),
                },
            },
        }
        written.append(_write_json(
            out / "epykit_dmc_summary_mqc.json", payload,
        ))

    logger.info("MultiQC custom files: wrote %d file(s) to %s", len(written), out)
    return str(out)


__all__ = ["report_multiqc"]
