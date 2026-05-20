from __future__ import annotations

import csv
import logging
from pathlib import Path

import polars as pl

from .convert import ensure_converted_sample
from .methyldata import MethylData

logger = logging.getLogger(__name__)


def _count_store_rows(store_dir: str) -> int | None:
    try:
        import pyarrow.parquet as pq

        total = 0
        for path in Path(store_dir).rglob("part-*.parquet"):
            total += pq.read_metadata(str(path)).num_rows
        return total
    except Exception:
        return None


def _build_obs_from_samplesheet(
    samplesheet: str,
    treatment_group: str | None,
    control_group: str | None,
    groups: list[str] | None,
) -> tuple[list[dict], list[tuple[str, str]]]:
    """Parse the samplesheet CSV into obs rows and (path, sample_id) pairs.

    Shared by every read_* entry point so format-specific code can focus
    on the conversion step.
    """
    with open(samplesheet) as handle:
        rows = list(csv.DictReader(handle))

    required = {"sample_id", "group", "path"}
    if not rows:
        raise ValueError("samplesheet contains no rows")
    missing_cols = required - set(rows[0].keys())
    if missing_cols:
        raise ValueError(f"samplesheet missing required columns: {sorted(missing_cols)}")

    if groups is not None:
        allowed = set(groups)
    elif treatment_group is not None and control_group is not None:
        allowed = {treatment_group, control_group}
    else:
        raise ValueError(
            "Either pass (treatment_group, control_group) for binary mode "
            "or groups=[...] for multi-group mode."
        )

    obs_rows: list[dict] = []
    files: list[tuple[str, str]] = []
    for row in rows:
        group = row["group"]
        if group not in allowed:
            continue
        obs_row = {
            "sample_id": row["sample_id"],
            "group": group,
            "path": row["path"],
        }
        if treatment_group is not None:
            obs_row["treatment"] = 1 if group == treatment_group else 0
        for key, value in row.items():
            if key not in {"sample_id", "group", "path"}:
                obs_row[key] = value
        obs_rows.append(obs_row)
        files.append((row["path"], row["sample_id"]))

    if not obs_rows:
        raise ValueError(
            "No samples matched the requested groups from samplesheet. "
            f"groups={sorted(allowed)}"
        )
    return obs_rows, files


def _read_methylation_samplesheet(
    samplesheet: str,
    *,
    pipeline: str,
    source_format: str,
    treatment_group: str | None,
    control_group: str | None,
    groups: list[str] | None,
    assembly: str,
    store_dir: str,
    context: str,
    reference_fasta: str | None,
) -> MethylData:
    """Backend shared by ``read_bismark`` and ``read_methyldackel``."""
    obs_rows, files = _build_obs_from_samplesheet(
        samplesheet, treatment_group, control_group, groups,
    )

    analysis_root = Path(store_dir)
    cache_store_dir = str(analysis_root / ".cache" / "raw")
    Path(cache_store_dir).mkdir(parents=True, exist_ok=True)
    for path, sample_id in files:
        converted = ensure_converted_sample(
            path,
            sample_id,
            cache_store_dir,
            context=context,
            reference_fasta=reference_fasta,
            format=source_format,
        )
        status = "converted" if converted else "cached"
        logger.info("  %s: %s", sample_id, status)

    n_sites_raw = _count_store_rows(cache_store_dir)
    uns: dict = {
        "samplesheet": samplesheet,
        "pipeline": pipeline,
        "source_format": source_format,
        "epykit_version": "0.1.0",
    }
    if n_sites_raw is not None:
        uns["n_sites_raw"] = n_sites_raw
    uns["_store_history"] = [
        {"step": "raw", "path": cache_store_dir, "n_sites": n_sites_raw}
    ]

    md = MethylData(
        obs=pl.DataFrame(obs_rows),
        store=cache_store_dir,
        assembly=assembly,
        context=context,
        uns=uns,
    )
    md._analysis_root = str(analysis_root)
    return md


def read_bismark(
    samplesheet: str,
    treatment_group: str | None = None,
    control_group: str | None = None,
    assembly: str = "unknown",
    store_dir: str = "methyl_store",
    context: str = "CpG",
    reference_fasta: str | None = None,
    groups: list[str] | None = None,
) -> MethylData:
    """Read a samplesheet of Bismark ``.cov[.gz]`` files into a MethylData.

    Expected samplesheet columns: sample_id, group, path

    Two operating modes:

    * Binary (legacy, default): pass ``treatment_group`` and
      ``control_group``. Only those two groups are kept; ``obs.treatment``
      is set to 1 for treatment-group samples and 0 for control.
    * Multi-group : pass ``groups=[g1, g2, ...]`` to load any
      subset of groups. ``obs.treatment`` is added only if
      ``treatment_group`` is also supplied; otherwise the column is
      omitted and downstream code falls back to formula-based contrasts.
    """
    return _read_methylation_samplesheet(
        samplesheet,
        pipeline="bismark",
        source_format="bismark",
        treatment_group=treatment_group,
        control_group=control_group,
        groups=groups,
        assembly=assembly,
        store_dir=store_dir,
        context=context,
        reference_fasta=reference_fasta,
    )


def read_methyldackel(
    samplesheet: str,
    treatment_group: str | None = None,
    control_group: str | None = None,
    assembly: str = "unknown",
    store_dir: str = "methyl_store",
    context: str = "CpG",
    reference_fasta: str | None = None,
    groups: list[str] | None = None,
) -> MethylData:
    """Read a samplesheet of MethylDackel ``.bedGraph[.gz]`` files into a
    MethylData.

    MethylDackel's ``extract`` output uses the same 6-column layout as
    Bismark ``.cov`` (``chrom, start, end, percent, M, U``) with one
    ``track`` header line at the top; that header is skipped automatically.

    Parameters are identical to :func:`read_bismark`; ``path`` entries in
    the samplesheet point at MethylDackel ``.bedGraph[.gz]`` files instead
    of Bismark ``.cov[.gz]``.
    """
    return _read_methylation_samplesheet(
        samplesheet,
        pipeline="methyldackel",
        source_format="methyldackel",
        treatment_group=treatment_group,
        control_group=control_group,
        groups=groups,
        assembly=assembly,
        store_dir=store_dir,
        context=context,
        reference_fasta=reference_fasta,
    )


def read_combined_strand_bed(
    samplesheet: str,
    treatment_group: str | None = None,
    control_group: str | None = None,
    assembly: str = "unknown",
    store_dir: str = "methyl_store",
    context: str = "CpG",
    reference_fasta: str | None = None,
    groups: list[str] | None = None,
) -> MethylData:
    """Read a samplesheet of 12-column strand-collapsed methylation BEDs.

    File schema (one row per CpG, tab-separated, no header)::

        chrom  start  end  fwd_M  fwd_T  fwd_pct  rev_M  rev_T  rev_pct  M  T  pct

    Columns 4-6 are the forward-strand counts ``(N_meth, coverage, %)``;
    7-9 are reverse-strand; 10-12 are the strand-collapsed total used
    when the upstream pipeline merged the two Cs of each CpG dyad
    onto a single position. epykit consumes the combined triplet:

      * ``N_meth``     = col 10 (M)
      * ``coverage``   = col 11 (T)
      * ``N_unmeth``   = col 11 - col 10
      * ``methylation_percent`` = col 12

    This is the format emitted by a number of in-house WGBS pipelines
    that merge CpG dyads before downstream analysis (e.g. the ``.bed.gz``
    files distributed under GEO accession GSE263850).

    Parameters are identical to :func:`read_bismark`; ``path`` entries
    in the samplesheet point at ``.bed`` / ``.bed.gz`` files in the
    12-column layout.
    """
    return _read_methylation_samplesheet(
        samplesheet,
        pipeline="combined_strand_bed",
        source_format="combined_strand_bed",
        treatment_group=treatment_group,
        control_group=control_group,
        groups=groups,
        assembly=assembly,
        store_dir=store_dir,
        context=context,
        reference_fasta=reference_fasta,
    )


def load(path: str) -> MethylData:
    """Load a previously saved MethylData analysis directory."""
    return MethylData.load(path)


def _candidate_sample_ids_from_filename(name: str) -> list[str]:
    candidates = {name}
    suffixes = [
        ".deduplicated.bismark.cov.gz",
        ".bismark.cov.gz",
        ".deduplicated.cov.gz",
        ".cov.gz",
        ".deduplicated.bismark.cov",
        ".bismark.cov",
        ".deduplicated.cov",
        ".cov",
    ]
    for suf in suffixes:
        if name.endswith(suf):
            candidates.add(name[: -len(suf)])
    if "." in name:
        candidates.add(name.split(".")[0])
    return sorted(candidates, key=len, reverse=True)


def read_nfcore_methylseq(
    run_dir: str,
    treatment_group: str,
    control_group: str,
    assembly: str = "unknown",
    store_dir: str = "methyl_store",
    context: str = "CpG",
    samplesheet_name: str = "samplesheet.csv",
) -> MethylData:
    """Load methylation data directly from an nf-core/methylseq run directory."""
    run = Path(run_dir).resolve()
    cov_dir = run / "results" / "bismark" / "deduplicated"
    samplesheet_path = run / samplesheet_name

    if not cov_dir.exists():
        raise FileNotFoundError(
            f"Expected nf-core/methylseq bismark directory at: {cov_dir}"
        )
    if not samplesheet_path.exists():
        raise FileNotFoundError(
            f"Expected samplesheet at: {samplesheet_path}"
        )

    cov_files = sorted(cov_dir.glob("*.cov.gz")) + sorted(cov_dir.glob("*.cov"))
    if not cov_files:
        raise FileNotFoundError(f"No .cov/.cov.gz files found in {cov_dir}")

    sample_to_cov: dict[str, str] = {}
    for p in cov_files:
        for candidate in _candidate_sample_ids_from_filename(p.name):
            sample_to_cov.setdefault(candidate, str(p))

    with open(samplesheet_path) as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise ValueError(f"Samplesheet '{samplesheet_path}' is empty")

    if "sample_id" not in rows[0].keys() or "group" not in rows[0].keys():
        raise ValueError(
            "nf-core samplesheet requires at least 'sample_id' and 'group' columns"
        )

    obs_rows: list[dict] = []
    # Organize stores under a .cache subdirectory for a cleaner output layout
    analysis_root = Path(store_dir)
    cache_store_dir = str(analysis_root / ".cache" / "raw")

    Path(cache_store_dir).mkdir(parents=True, exist_ok=True)
    for row in rows:
        group = row["group"]
        if group not in (treatment_group, control_group):
            continue

        sample_id = row["sample_id"]
        cov_path = sample_to_cov.get(sample_id)
        if cov_path is None:
            raise FileNotFoundError(
                f"Could not match sample '{sample_id}' to .cov file in {cov_dir}. "
                f"Detected candidates: {sorted(sample_to_cov.keys())[:15]}"
            )

        obs_row = {
            "sample_id": sample_id,
            "group": group,
            "treatment": 1 if group == treatment_group else 0,
            "path": cov_path,
        }
        for key, value in row.items():
            if key not in {"sample_id", "group"}:
                obs_row[key] = value
        obs_rows.append(obs_row)

        logger.info("  %s (%s) <- %s", sample_id, group, cov_path)
        ensure_converted_sample(cov_path, sample_id, cache_store_dir, context=context)

    if not obs_rows:
        raise ValueError(
            "No samples matched treatment/control groups from nf-core samplesheet. "
            f"treatment_group={treatment_group}, control_group={control_group}"
        )

    n_sites_raw = _count_store_rows(cache_store_dir)
    uns = {
        "pipeline": "nf-core/methylseq",
        "nfcore_run": str(run),
        "samplesheet": str(samplesheet_path),
        "epykit_version": "0.1.0",
        "n_sites_raw": n_sites_raw,
        "_store_history": [{"step": "raw", "path": cache_store_dir, "n_sites": n_sites_raw}],
    }
    md = MethylData(
        obs=pl.DataFrame(obs_rows),
        store=cache_store_dir,
        assembly=assembly,
        context=context,
        uns=uns,
    )
    md._analysis_root = str(analysis_root)
    return md
