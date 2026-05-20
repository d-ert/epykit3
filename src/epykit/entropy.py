"""Methylation entropy: read-level Shannon entropy over CpG-window patterns.

For each genomic window of ``window_cpgs`` consecutive CpGs and each
sample, we collect every read that covers all CpGs in the window,
encode its methylation calls as a binary string (e.g. ``1010`` for a
4-CpG window), tally the 2^window_cpgs possible patterns, and compute
Shannon entropy.

The 'normalised_entropy' column is entropy divided by
``log2(2**window_cpgs) = window_cpgs`` so it sits in ``[0, 1]``. A
fully ordered window (one dominant pattern) has entropy ~0; a fully
disordered window (uniform over all patterns) has normalised entropy
~1.

Higher entropy reflects mixed methylation states across the reads at a
locus -- a signature of stochastic methylation, age-related drift, or
intra-tumour heterogeneity. Mean beta alone misses this; two reads with
patterns ``1100`` and ``0011`` have mean beta = 0.5 but max entropy.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Mapping, Optional

import numpy as np
import polars as pl

from .bam_io import read_methylation_calls

logger = logging.getLogger(__name__)


_ENTROPY_SCHEMA = {
    "sample_id":          pl.Utf8,
    "chrom":              pl.Utf8,
    "start":              pl.Int32,
    "end":                pl.Int32,
    "n_reads":            pl.Int32,
    "entropy":            pl.Float64,
    "normalised_entropy": pl.Float64,
}


def call_entropy(
    bam: Mapping[str, str | Path],
    *,
    window_cpgs: int = 4,
    min_reads: int = 10,
    chromosomes: Optional[list[str]] = None,
    caller: str = "bismark",
    min_baseq: int = 20,
    min_mapq: int = 10,
) -> pl.DataFrame:
    """Compute per-CpG-window methylation entropy across one or more samples.

    Parameters
    ----------
    bam
        ``{sample_id -> bam_path}``. BAMs must be coordinate-sorted and
        indexed.
    window_cpgs
        Number of consecutive CpGs in each entropy window. Must be in
        [2, 8]; 4 is the canonical choice (16 patterns).
    min_reads
        Skip windows with fewer than this many fully-covering reads.
    chromosomes
        Restrict analysis to these chromosomes; ``None`` for all.
    caller
        BAM dialect (``"bismark"`` or ``"methyldackel"``); see
        :mod:`epykit.bam_io`.

    Returns
    -------
    pl.DataFrame
        Long-form frame with columns
        ``(sample_id, chrom, start, end, n_reads, entropy,
        normalised_entropy)``.
    """
    if not (2 <= window_cpgs <= 8):
        raise ValueError(
            f"window_cpgs={window_cpgs}: must be in [2, 8] "
            "(beyond 8 the pattern space (256) gets too sparse for "
            "typical WGBS coverage)."
        )

    sample_frames: list[pl.DataFrame] = []
    log2_states = float(window_cpgs)  # log2(2**window_cpgs)

    for sample_id, bam_path in bam.items():
        logger.info("[entropy] %s: reading BAM", sample_id)
        meth_df = read_methylation_calls(
            bam_path, caller=caller, min_baseq=min_baseq, min_mapq=min_mapq,
        )
        if meth_df.height == 0:
            logger.warning("[entropy] %s: no usable reads", sample_id)
            continue
        if chromosomes is not None:
            meth_df = meth_df.filter(pl.col("chrom").is_in(chromosomes))
            if meth_df.height == 0:
                continue

        df = _entropy_one_sample(
            meth_df, window_cpgs=window_cpgs, min_reads=min_reads,
            log2_states=log2_states,
        )
        if df.height > 0:
            df = df.with_columns(pl.lit(sample_id).alias("sample_id"))
            sample_frames.append(df)

    if not sample_frames:
        return pl.DataFrame(schema=_ENTROPY_SCHEMA)

    return pl.concat(sample_frames, how="vertical_relaxed").sort(
        ["sample_id", "chrom", "start"]
    )


def _entropy_one_sample(
    meth_df: pl.DataFrame,
    *,
    window_cpgs: int,
    min_reads: int,
    log2_states: float,
) -> pl.DataFrame:
    """Per-sample entropy: walk each chrom's CpGs and score windows."""
    out_rows: list[dict[str, object]] = []
    for chrom_name, chrom_grp in meth_df.partition_by("chrom", as_dict=True).items():
        # Sort calls by (read_id, pos) so per-read patterns are stable.
        chrom_grp = chrom_grp.sort(["pos", "read_id"])
        # Unique sorted CpG positions on this chrom.
        cpg_positions = chrom_grp["pos"].unique().sort().to_numpy()
        if len(cpg_positions) < window_cpgs:
            continue

        # Index reads -> (pos, methylation) calls for quick lookup.
        read_calls: dict[str, dict[int, int]] = {}
        for row in chrom_grp.iter_rows(named=True):
            rid = row["read_id"]
            read_calls.setdefault(rid, {})[int(row["pos"])] = int(
                row["methylation_status"]
            )

        # Walk every window of window_cpgs consecutive CpGs.
        for i in range(len(cpg_positions) - window_cpgs + 1):
            window_pos = cpg_positions[i: i + window_cpgs]
            pattern_counts: dict[int, int] = {}
            n_reads = 0
            for calls in read_calls.values():
                # Reject reads that don't cover ALL CpGs in the window.
                pat = 0
                ok = True
                for j, p in enumerate(window_pos):
                    p_int = int(p)
                    if p_int not in calls:
                        ok = False
                        break
                    if calls[p_int] not in (0, 1):
                        ok = False
                        break
                    pat |= (calls[p_int] << j)
                if not ok:
                    continue
                pattern_counts[pat] = pattern_counts.get(pat, 0) + 1
                n_reads += 1

            if n_reads < min_reads:
                continue

            # Shannon entropy across the observed pattern frequencies.
            ent = 0.0
            for c in pattern_counts.values():
                p = c / n_reads
                if p > 0:
                    ent -= p * math.log2(p)
            out_rows.append({
                "chrom": chrom_name[0] if isinstance(chrom_name, tuple) else chrom_name,
                "start": int(window_pos[0]),
                "end": int(window_pos[-1]) + 1,
                "n_reads": int(n_reads),
                "entropy": float(ent),
                "normalised_entropy": float(ent / log2_states) if log2_states else 0.0,
            })

    if not out_rows:
        return pl.DataFrame(schema={k: v for k, v in _ENTROPY_SCHEMA.items()
                                    if k != "sample_id"})
    return pl.DataFrame(
        out_rows,
        schema={k: v for k, v in _ENTROPY_SCHEMA.items() if k != "sample_id"},
    )


def entropy(
    md,
    *,
    bam: Mapping[str, str | Path],
    window_cpgs: int = 4,
    min_reads: int = 10,
    chromosomes: Optional[list[str]] = None,
    caller: str = "bismark",
) -> None:
    """Run entropy and store results in ``md.varm["entropy"]``."""
    md_samples = set(md.obs.get_column("sample_id").to_list())
    missing = [s for s in bam if s not in md_samples]
    if missing:
        raise ValueError(f"bam keys not in md.obs.sample_id: {missing[:5]}")
    result = call_entropy(
        bam=bam, window_cpgs=window_cpgs, min_reads=min_reads,
        chromosomes=chromosomes, caller=caller,
    )
    md.varm["entropy"] = result
    md.uns["entropy"] = {
        "n_windows": int(result.height),
        "n_samples": len(bam),
        "window_cpgs": window_cpgs,
        "min_reads": min_reads,
    }


__all__ = ["call_entropy", "entropy"]
