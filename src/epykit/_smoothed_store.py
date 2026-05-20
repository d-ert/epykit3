"""Build a temporary methylstore of smoothed pseudo-counts.

Used by ``ep.tl.dmc(use_smoothed=True)`` (and indirectly by ``ep.tl.dmr``
when it inherits the flag) to feed BSmooth- or Gaussian-smoothed
methylation estimates into the existing DMC engine without rewriting
any of the eight test backends.

The transform is per (sample, chrom, CpG):

    coverage_smooth  = original coverage   (unchanged)
    N_meth_smooth    = round(beta_smooth * coverage_smooth)
    N_unmeth_smooth  = coverage_smooth - N_meth_smooth

Sites whose ``beta_smooth`` is NaN (i.e. the smoother couldn't fit a
local polynomial -- too few neighbours) fall back to the original
``N_meth`` / ``coverage`` from the raw store. This is intentional: we
never silently drop signal, and the smoother's own ``min_cpgs_for_smooth``
fallback to raw beta means most NaN cases are already handled upstream.

This module exists as a small focused unit so the DMC code path stays
clean -- ``process_chromosomes_dmc`` just takes a different
``methylstore_path`` and runs normally.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import polars as pl

logger = logging.getLogger(__name__)


def build_smoothed_pseudo_count_store(
    raw_store: Path,
    smooth_store: Path,
    samples: list[str],
    out_dir: Path,
) -> None:
    """Materialise a methylstore of smoothed pseudo-counts at ``out_dir``.

    Parameters
    ----------
    raw_store
        Path to the existing methylstore that ``tl.dmc`` would normally
        consume (``md.store``).
    smooth_store
        Path to the sidecar parquet store written by
        :func:`smooth_methylation_bsmooth` (or the Gaussian smoother).
        Hive-partitioned ``sample=*/chrom=*`` with columns
        ``chrom, pos, sample, beta_raw, beta_smooth``.
    samples
        Sample identifiers to materialise. Typically every sample on
        ``md.obs.sample_id``.
    out_dir
        Output root. Will be created if missing. Existing contents are
        not cleared -- callers are responsible for using a fresh
        :class:`tempfile.TemporaryDirectory`.
    """
    raw_store = Path(raw_store)
    smooth_store = Path(smooth_store)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        raw_sample_dir = raw_store / f"sample={sample}"
        smooth_sample_dir = smooth_store / f"sample={sample}"
        if not raw_sample_dir.exists():
            logger.warning(
                "[smoothed-store] sample %s missing from raw store; skipping",
                sample,
            )
            continue
        if not smooth_sample_dir.exists():
            raise FileNotFoundError(
                f"Smoothed sidecar missing sample={sample} at {smooth_sample_dir}. "
                "Did ep.pp.smooth(md) skip this sample?"
            )

        for chrom_dir in sorted(raw_sample_dir.glob("chrom=*")):
            chrom = chrom_dir.name.removeprefix("chrom=")
            raw_part = chrom_dir / "part-0.parquet"
            smooth_part = smooth_sample_dir / f"chrom={chrom}" / "part-0.parquet"
            if not raw_part.exists():
                continue

            raw_df = pl.read_parquet(str(raw_part))

            if not smooth_part.exists():
                # The smoother may have skipped a chrom with too few sites;
                # fall back to raw counts so downstream DMC still sees this
                # chrom (with the raw signal, not smoothed).
                logger.debug(
                    "[smoothed-store] %s/%s missing in smooth sidecar; using raw",
                    sample, chrom,
                )
                out_chunk = raw_df
            else:
                smooth_df = pl.read_parquet(
                    str(smooth_part), columns=["pos", "beta_smooth"]
                )
                # Left-join: every raw site keeps its strand / context /
                # coverage; beta_smooth comes from the sidecar.
                joined = raw_df.join(smooth_df, on="pos", how="left")
                pos = joined["pos"].to_numpy()
                cov = joined["coverage"].to_numpy().astype(np.int64)
                raw_n_meth = joined["N_meth"].to_numpy().astype(np.int64)
                beta = joined["beta_smooth"].to_numpy().astype(np.float64)

                # Where beta_smooth is NaN (rare -- smoother fell back to
                # raw beta for those sites anyway, but defensive), keep
                # the raw N_meth / N_unmeth.
                valid = np.isfinite(beta)
                pseudo_n_meth = np.where(
                    valid,
                    np.round(beta * cov).astype(np.int64),
                    raw_n_meth,
                )
                # Clamp into [0, coverage] in case rounding pushed beta
                # slightly outside [0, 1] (the smoother itself clips, but
                # belt-and-braces for the Gaussian path which clips at
                # the grid stage but can still drift on a float32 cast).
                pseudo_n_meth = np.clip(pseudo_n_meth, 0, cov)
                pseudo_n_unmeth = cov - pseudo_n_meth

                out_chunk = joined.with_columns([
                    pl.Series("N_meth", pseudo_n_meth.astype(np.int32)),
                    pl.Series("N_unmeth", pseudo_n_unmeth.astype(np.int32)),
                ]).drop("beta_smooth").select(raw_df.columns)

            out_part = out_dir / f"sample={sample}" / f"chrom={chrom}"
            out_part.mkdir(parents=True, exist_ok=True)
            out_chunk.write_parquet(str(out_part / "part-0.parquet"), compression="zstd")


__all__ = ["build_smoothed_pseudo_count_store"]
