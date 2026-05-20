"""nf-core/methylseq run-dir QC ingestion .

Picks up Bismark alignment reports, Qualimap, and preseq outputs from
the run directory and returns a per-sample ``pl.DataFrame`` that
``MethylData.obs`` can be left-joined against.

Parsers here are small targeted regex passes (the same approach MultiQC
uses). No MultiQC dependency.
"""

from __future__ import annotations

import csv
import logging
import re
from pathlib import Path
from typing import Optional

import polars as pl

logger = logging.getLogger(__name__)


_BISMARK_PATTERNS = {
    "bismark_aligned_reads": re.compile(
        r"Number of alignments analysed in total:\s*(\d+)", re.MULTILINE
    ),
    "bismark_unique_alignments": re.compile(
        r"Number of paired-end alignments with a unique best hit:\s*(\d+)",
        re.MULTILINE,
    ),
    "bismark_mapping_efficiency": re.compile(
        r"Mapping efficiency:\s*([0-9.]+)%", re.MULTILINE
    ),
    "bismark_pct_meth_cpg": re.compile(
        r"C methylated in CpG context:\s*([0-9.]+)%", re.MULTILINE
    ),
    "bismark_pct_meth_chg": re.compile(
        r"C methylated in CHG context:\s*([0-9.]+)%", re.MULTILINE
    ),
    "bismark_pct_meth_chh": re.compile(
        r"C methylated in CHH context:\s*([0-9.]+)%", re.MULTILINE
    ),
}

_QUALIMAP_PATTERNS = {
    "qualimap_mean_coverage": re.compile(
        r"mean coverageData =\s*([0-9.]+)X", re.MULTILINE
    ),
    "qualimap_std_coverage": re.compile(
        r"std coverageData =\s*([0-9.]+)X", re.MULTILINE
    ),
}


def _parse_bismark_report(path: Path) -> dict:
    txt = path.read_text(errors="replace")
    out = {}
    for key, pat in _BISMARK_PATTERNS.items():
        m = pat.search(txt)
        if m:
            try:
                out[key] = float(m.group(1))
            except ValueError:
                out[key] = m.group(1)
    return out


def parse_bismark_mbias(path: str) -> pl.DataFrame:
    """Parse a Bismark M-bias report into a long DataFrame.

    Bismark's ``--mbias_only`` (or the M-bias section of the standard
    methylation extraction report) emits per-position methylation tallies
    grouped into 1-6 context+read panels, e.g.::

        CpG context (R1)
        ================
        position    count methylated    count unmethylated    % methylation    coverage
        1    123    45    73.21    168
        ...

        CpG context (R2)
        ================
        ...

    Single-end runs only have ``(R1)`` blocks. CHG / CHH contexts may be
    absent for libraries that didn't emit them.

    Returns
    -------
    pl.DataFrame
        Columns ``position`` (Int64), ``context`` (Utf8: CpG / CHG / CHH),
        ``read`` (Utf8: R1 / R2), ``n_meth`` (Int64), ``n_unmeth`` (Int64),
        ``percent`` (Float64), ``coverage`` (Int64). Sorted by
        (context, read, position).
    """
    header_re = re.compile(
        r"^(CpG|CHG|CHH)\s+context\s*(?:\(R(\d)\))?\s*$", re.IGNORECASE
    )
    txt = Path(path).read_text(errors="replace")
    lines = txt.splitlines()
    rows: list[dict] = []
    i = 0
    while i < len(lines):
        m = header_re.match(lines[i].strip())
        if not m:
            i += 1
            continue
        context = m.group(1).upper().replace("CPG", "CpG")
        read = f"R{m.group(2) or '1'}"
        # Skip the underline + column header lines.
        i += 1
        while i < len(lines) and (
            lines[i].startswith("=")
            or lines[i].lower().startswith("position")
            or lines[i].strip() == ""
        ):
            i += 1
        # Consume data rows until a blank line or the next header.
        while i < len(lines):
            line = lines[i].strip()
            if not line or header_re.match(line):
                break
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 4:
                i += 1
                continue
            try:
                pos = int(parts[0])
                n_meth = int(parts[1])
                n_unmeth = int(parts[2])
                pct = float(parts[3]) if parts[3] not in ("nan", "NA") else float("nan")
                cov = int(parts[4]) if len(parts) >= 5 else n_meth + n_unmeth
            except ValueError:
                i += 1
                continue
            rows.append({
                "position": pos,
                "context": context,
                "read": read,
                "n_meth": n_meth,
                "n_unmeth": n_unmeth,
                "percent": pct,
                "coverage": cov,
            })
            i += 1
    if not rows:
        return pl.DataFrame(schema={
            "position": pl.Int64, "context": pl.Utf8, "read": pl.Utf8,
            "n_meth": pl.Int64, "n_unmeth": pl.Int64,
            "percent": pl.Float64, "coverage": pl.Int64,
        })
    return pl.DataFrame(rows).sort(["context", "read", "position"])


def _parse_qualimap(path: Path) -> dict:
    txt = path.read_text(errors="replace")
    out = {}
    for key, pat in _QUALIMAP_PATTERNS.items():
        m = pat.search(txt)
        if m:
            try:
                out[key] = float(m.group(1))
            except ValueError:
                out[key] = m.group(1)
    return out


def _resolve_sample_ids(samplesheet: str) -> list[str]:
    with open(samplesheet, newline="") as fh:
        return [row["sample_id"] for row in csv.DictReader(fh) if "sample_id" in row]


def read_nfcore_methylseq_qc(
    samplesheet: Optional[str],
    run_dir: str,
    *,
    sample_ids: Optional[list[str]] = None,
) -> pl.DataFrame:
    """Walk an nf-core/methylseq run dir and pull per-sample QC metrics.

    Searches for ``<sample>*_PE_report.txt`` / ``<sample>*_SE_report.txt``
    under ``run_dir/`` (recursively) and Qualimap ``genome_results.txt``
    in a Qualimap output dir keyed by sample.

    Parameters
    ----------
    samplesheet : str, optional
        Path to the samplesheet used in the pipeline run. Required when
        ``sample_ids`` isn't provided.
    run_dir : str
        nf-core/methylseq output directory.
    sample_ids : list[str], optional
        Explicit list of sample IDs to look for. Overrides samplesheet.

    Returns
    -------
    pl.DataFrame
        One row per sample with the union of parsed metrics. Missing
        metrics are NaN. The DataFrame can be left-joined onto
        ``md.obs`` via ``obs.join(qc_df, on="sample_id", how="left")``.
    """
    run = Path(run_dir)
    if sample_ids is None:
        if samplesheet is None:
            raise ValueError("pass either samplesheet or sample_ids")
        sample_ids = _resolve_sample_ids(samplesheet)

    rows: list[dict] = []
    for sample in sample_ids:
        record: dict = {"sample_id": sample}
        # Bismark report(s)
        for pattern in (
            f"**/{sample}*_PE_report.txt", f"**/{sample}*_SE_report.txt"
        ):
            for hit in run.glob(pattern):
                try:
                    record.update(_parse_bismark_report(hit))
                    break
                except Exception as exc:
                    logger.warning(
                        "failed to parse Bismark report %s: %s", hit, exc
                    )
            if any(k.startswith("bismark_") for k in record):
                break
        # Qualimap
        for hit in run.glob(f"**/{sample}*/genome_results.txt"):
            try:
                record.update(_parse_qualimap(hit))
                break
            except Exception as exc:
                logger.warning(
                    "failed to parse Qualimap output %s: %s", hit, exc
                )
        rows.append(record)
    return pl.DataFrame(rows)


__all__ = ["read_nfcore_methylseq_qc"]
