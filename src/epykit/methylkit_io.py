"""methylKit-compatible tabix export.

Writes tab-delimited methylation tables that methylKit's ``methylRawDB``
can ingest. Schema:

    chrBase  chr  base  strand  coverage  freqC  freqT

stored as bgzip + tabix-indexed text. Useful for sharing epykit output
with downstream tools that expect the methylKit input format.

Note on dependencies
--------------------
* Plain bgzip and tabix indexing requires ``pysam``. On Windows, pysam
  has no PyPI wheel; we still write the plain tab-delimited text and
  best-effort gzip the file, but the ``.tbi`` index step is silently
  skipped when pysam isn't importable.
"""

from __future__ import annotations

import gzip
import json
import logging
from pathlib import Path
from typing import Optional, Sequence

import polars as pl

logger = logging.getLogger(__name__)


_METHYLKIT_SCHEMA: tuple[str, ...] = (
    "chrBase", "chr", "base", "strand",
    "coverage", "freqC", "freqT",
)


def _write_one_sample(
    store_root: Path,
    sample: str,
    out_path: Path,
    *,
    chromosomes: Optional[Sequence[str]] = None,
) -> int:
    sample_dir = store_root / f"sample={sample}"
    if not sample_dir.exists():
        raise FileNotFoundError(
            f"sample={sample} not present under {store_root}"
        )
    written = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    chrom_dirs = sorted(sample_dir.glob("chrom=*"))
    if chromosomes:
        wanted = set(chromosomes)
        chrom_dirs = [
            d for d in chrom_dirs
            if d.name.removeprefix("chrom=") in wanted
        ]
    with gzip.open(out_path, "wt") as fh:
        fh.write("\t".join(_METHYLKIT_SCHEMA) + "\n")
        for chrom_dir in chrom_dirs:
            chrom = chrom_dir.name.removeprefix("chrom=")
            for part in sorted(chrom_dir.glob("part-*.parquet")):
                df = pl.read_parquet(
                    str(part),
                    columns=["pos", "strand", "N_meth", "coverage"],
                ).filter(pl.col("coverage") > 0)
                if len(df) == 0:
                    continue
                for row in df.iter_rows(named=True):
                    pos = int(row["pos"])
                    cov = int(row["coverage"])
                    n_meth = int(row["N_meth"])
                    freq_c = 100.0 * n_meth / max(cov, 1)
                    freq_t = 100.0 - freq_c
                    strand = row.get("strand") or "*"
                    fh.write(
                        f"{chrom}.{pos}\t{chrom}\t{pos}\t{strand}\t"
                        f"{cov}\t{freq_c:.6f}\t{freq_t:.6f}\n"
                    )
                    written += 1
    return written


def to_methylkit_tabix(
    md,
    output_dir: str,
    samples: Optional[list[str]] = None,
    *,
    chromosomes: Optional[Sequence[str]] = None,
) -> str:
    """Export per-sample methylKit-schema tab-delimited tables.

    Returns the output directory. A manifest file
    ``epykit_to_methylkit.json`` is written alongside the per-sample
    ``.txt.gz`` files mapping ``sample -> file -> treatment_assignment``
    (when ``md.obs["treatment"]`` exists).

    Tabix indexing is attempted via ``pysam`` and silently skipped on
    platforms where ``pysam`` isn't importable.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    store_root = Path(md.store)
    if samples is None:
        samples = md.obs.get_column("sample_id").to_list()

    manifest: dict = {"samples": [], "schema": list(_METHYLKIT_SCHEMA)}
    treatment_col = (
        md.obs.get_column("treatment").to_list()
        if "treatment" in md.obs.columns else None
    )
    obs_samples = md.obs.get_column("sample_id").to_list()
    for sample in samples:
        out_path = out / f"{sample}.methylraw.txt.gz"
        n_rows = _write_one_sample(
            store_root, sample, out_path, chromosomes=chromosomes,
        )
        # Best-effort tabix indexing.
        tbi_path: Optional[str] = None
        try:
            import pysam  # type: ignore
            tbi = pysam.tabix_index(
                str(out_path), preset=None, seq_col=1, start_col=2,
                end_col=2, force=True,
            )
            tbi_path = tbi
        except ImportError:
            logger.info(
                "pysam not installed; skipping tabix index for %s "
                "(install 'epykit[methylkit]' for tabix support).",
                sample,
            )
        except Exception as exc:
            logger.warning("tabix indexing failed for %s: %s", sample, exc)

        treatment = None
        if treatment_col is not None:
            try:
                treatment = int(treatment_col[obs_samples.index(sample)])
            except (ValueError, IndexError):
                treatment = None
        manifest["samples"].append({
            "sample_id": sample,
            "file": out_path.name,
            "tabix": (Path(tbi_path).name if tbi_path else None),
            "n_rows": n_rows,
            "treatment": treatment,
        })

    (out / "epykit_to_methylkit.json").write_text(
        json.dumps(manifest, indent=2),
    )
    return str(out)


__all__ = ["to_methylkit_tabix"]
