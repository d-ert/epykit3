"""Backend dispatcher for per-chromosome compute.

The per-chromosome streaming pattern is identical across DMC, DMR (tile),
and DVC: iterate chromosomes, for each chrom build the canonical position
index, run a handler over it, stage the per-chrom result to a tempdir,
concatenate at the end. This module captures the iteration / dispatch
half of that pattern so the existing engines can opt into distributed
execution without rewriting their per-chrom math.

Default ``backend="sequential"`` is bit-identical to the prior in-line
loop -- the dispatcher just invokes the handler one chrom at a time on
the calling process. ``backend="dask"`` and ``backend="ray"`` are
optional (extras ``[distributed]`` / ``[ray]``) and submit one task per
chromosome to a worker pool.

The handler contract is intentionally minimal: a callable that takes a
single ``chrom: str`` argument and returns either a ``pl.DataFrame``
(the result for that chrom) or ``None`` (skip -- no rows, no warning).
Per-engine state (store path, sample lists, knobs) is captured by
closure or :func:`functools.partial`; both work under
``cloudpickle`` which Dask and Ray use by default.
"""

from __future__ import annotations

import logging
from typing import Callable, Iterable, Iterator, Optional

import polars as pl

logger = logging.getLogger(__name__)


ChromHandler = Callable[[str], Optional[pl.DataFrame]]


def run_chrom_pipeline(
    chromosomes: Iterable[str],
    handler: ChromHandler,
    *,
    backend: str = "sequential",
    n_workers: Optional[int] = None,
    label: str = "chrom",
) -> Iterator[tuple[str, pl.DataFrame]]:
    """Run one handler invocation per chromosome through the chosen backend.

    Parameters
    ----------
    chromosomes
        Iterable of chromosome names to process. Consumed eagerly so the
        order of results matches the order of inputs (important for the
        existing engines, which write per-chrom parquet files keyed by
        chromosome).
    handler
        Callable accepting a single chrom name. Returns the per-chrom
        result DataFrame, or None to skip.
    backend
        One of ``"sequential"`` (default), ``"dask"``, ``"ray"``.
    n_workers
        Number of workers for distributed backends. Ignored for
        sequential. ``None`` means the backend picks a default (Dask:
        cpu_count; Ray: cpu_count).
    label
        Short string used in progress logs (e.g. ``"DMC"``, ``"DVC"``).

    Yields
    ------
    (chrom, result_df) tuples in submission order. Chroms whose handler
    returned None or an empty frame are filtered out.
    """
    chromosomes = list(chromosomes)
    if not chromosomes:
        return

    backend = (backend or "sequential").lower()
    if backend == "sequential":
        yield from _run_sequential(chromosomes, handler, label)
    elif backend == "dask":
        yield from _run_dask(chromosomes, handler, n_workers, label)
    elif backend == "ray":
        yield from _run_ray(chromosomes, handler, n_workers, label)
    else:
        raise ValueError(
            f"Unknown compute backend {backend!r}. "
            "Use 'sequential' (default), 'dask', or 'ray'."
        )


def _emit(chrom: str, result: Optional[pl.DataFrame]) -> Optional[tuple[str, pl.DataFrame]]:
    """Filter empty / None results; return None to signal 'skip'."""
    if result is None:
        return None
    if len(result) == 0:
        return None
    return chrom, result


def _run_sequential(
    chromosomes: list[str],
    handler: ChromHandler,
    label: str,
) -> Iterator[tuple[str, pl.DataFrame]]:
    n = len(chromosomes)
    for i, chrom in enumerate(chromosomes):
        logger.info("[%s %d/%d] %s", label, i + 1, n, chrom)
        result = handler(chrom)
        emitted = _emit(chrom, result)
        if emitted is not None:
            yield emitted


def _run_dask(
    chromosomes: list[str],
    handler: ChromHandler,
    n_workers: Optional[int],
    label: str,
) -> Iterator[tuple[str, pl.DataFrame]]:
    try:
        from dask.distributed import Client, LocalCluster, get_client  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "Dask is required for backend='dask'. "
            "Install with: pip install 'epykit[distributed]'"
        ) from exc

    # Reuse an active client if the caller already opened one; otherwise
    # spin up a LocalCluster sized to n_workers.
    owned_cluster: Optional[LocalCluster] = None
    owned_client: Optional[Client] = None
    try:
        try:
            client = get_client()
            logger.info(
                "[%s] Using existing Dask client at %s for %d chrom(s)",
                label, client.scheduler.address, len(chromosomes),
            )
        except (ValueError, RuntimeError):
            owned_cluster = LocalCluster(n_workers=n_workers or None, processes=True)
            owned_client = Client(owned_cluster)
            client = owned_client
            logger.info(
                "[%s] Spun up LocalCluster (%d worker(s)) for %d chrom(s)",
                label, len(client.scheduler_info().get("workers", {})), len(chromosomes),
            )

        # Submit one future per chromosome. `pure=False` so Dask doesn't
        # de-duplicate by hash -- handler closures often look identical
        # to Dask's task hasher even though they target different chroms.
        futures = [
            client.submit(handler, chrom, pure=False, key=f"{label}-{chrom}")
            for chrom in chromosomes
        ]

        # Gather in submission order so downstream tempdir-staging logic
        # writes files in a deterministic order matching the sequential
        # path. (`as_completed` would be faster but would not match.)
        for i, (chrom, future) in enumerate(zip(chromosomes, futures)):
            logger.info("[%s %d/%d] %s (awaiting)", label, i + 1, len(chromosomes), chrom)
            result = future.result()
            emitted = _emit(chrom, result)
            if emitted is not None:
                yield emitted
    finally:
        if owned_client is not None:
            owned_client.close()
        if owned_cluster is not None:
            owned_cluster.close()


def _run_ray(
    chromosomes: list[str],
    handler: ChromHandler,
    n_workers: Optional[int],
    label: str,
) -> Iterator[tuple[str, pl.DataFrame]]:
    try:
        import ray  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "Ray is required for backend='ray'. "
            "Install with: pip install 'epykit[ray]'"
        ) from exc

    owned_init = False
    if not ray.is_initialized():
        ray.init(num_cpus=n_workers, ignore_reinit_error=True, log_to_driver=False)
        owned_init = True
        logger.info(
            "[%s] Initialised Ray (%s CPUs) for %d chrom(s)",
            label, n_workers or "auto", len(chromosomes),
        )
    else:
        logger.info(
            "[%s] Using existing Ray runtime for %d chrom(s)",
            label, len(chromosomes),
        )

    try:
        # ray.remote(handler) creates a remote function; ray gets the
        # cloudpickled closure. Submit one task per chromosome.
        remote_handler = ray.remote(handler)
        object_refs = [remote_handler.remote(chrom) for chrom in chromosomes]

        for i, (chrom, ref) in enumerate(zip(chromosomes, object_refs)):
            logger.info("[%s %d/%d] %s (awaiting)", label, i + 1, len(chromosomes), chrom)
            result = ray.get(ref)
            emitted = _emit(chrom, result)
            if emitted is not None:
                yield emitted
    finally:
        if owned_init:
            ray.shutdown()


__all__ = ["run_chrom_pipeline", "ChromHandler"]
