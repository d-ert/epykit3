"""Shared on-disk cache manifest helpers.

Two related but distinct cache layers live here:

  1. **Per-sample, per-step manifests** (``.epykit_<step>_manifest.json``
     inside each ``sample=*/`` dir). These are the originals from
     0.1-0.3: convert / filter / normalize fingerprint the upstream
     sample + params and skip the step if the fingerprint matches. Each
     step manages its own manifest.

  2. **Top-level pipeline manifest** (``.epykit_manifest.json`` at the
     analysis root). New in 0.4.0. Tracks completed *pipeline stages*
     across the whole analysis (raw, filtered, united, smoothed, dmc_lr,
     dmr_tile, ...) with their params + input signatures + sidecar paths.
     This is what powers the formal ``ep.pp.resume(md)`` /
     ``MethylData.resume_from(stage=...)`` API.

The two coexist: the per-sample manifests still gate per-sample work;
the pipeline manifest gates top-level stages and lets a long pipeline
restart from any completed stage without rerunning earlier work.
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Any, Optional


def file_signature(path: Path) -> dict[str, Any]:
    """Path + size + mtime fingerprint of a single file."""
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def load_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    with path.open() as handle:
        return json.load(handle)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def sample_dir(output_dir: Path, sample_name: str) -> Path:
    return output_dir / f"sample={sample_name}"


def expected_chrom_dirs(sd: Path) -> list[str]:
    return sorted(item.name for item in sd.glob("chrom=*") if item.is_dir())


def sample_is_complete(sd: Path, chroms: list[str]) -> bool:
    """True iff `sd` has exactly the expected chrom dirs and each has part-0."""
    if not sd.exists():
        return False
    if expected_chrom_dirs(sd) != sorted(chroms):
        return False
    return all((sd / chrom / "part-0.parquet").exists() for chrom in chroms)


# Top-level pipeline manifest (0.4.0 checkpoint/resume API)

_PIPELINE_MANIFEST_NAME = ".epykit_manifest.json"


def pipeline_manifest_path(analysis_root: Path | str) -> Path:
    """Return the canonical path to the pipeline-level manifest."""
    return Path(analysis_root) / _PIPELINE_MANIFEST_NAME


def manifest_read(analysis_root: Path | str) -> dict[str, Any]:
    """Read the pipeline manifest. Returns an empty skeleton if missing.

    Skeleton: ``{"epykit_version": <ver>, "stages": []}``. ``stages`` is
    a list of entries in order of completion.
    """
    mp = pipeline_manifest_path(analysis_root)
    if not mp.exists():
        from epykit import __version__  # local to avoid circular import at import time
        return {"epykit_version": __version__, "stages": []}
    with mp.open() as handle:
        return json.load(handle)


def manifest_write(analysis_root: Path | str, payload: dict[str, Any]) -> None:
    """Atomic-ish write of the pipeline manifest (write-temp + rename)."""
    mp = pipeline_manifest_path(analysis_root)
    mp.parent.mkdir(parents=True, exist_ok=True)
    tmp = mp.with_suffix(mp.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(mp)


def manifest_append(
    analysis_root: Path | str,
    stage: str,
    *,
    params: dict[str, Any],
    input_sig: str,
    output_path: str,
    extra: Optional[dict[str, Any]] = None,
) -> None:
    """Record one completed stage in the pipeline manifest.

    If a stage with the same ``name`` already exists, the prior entry is
    replaced (rerunning a stage updates the manifest in place). Order
    in the ``stages`` list reflects the *most recent* completion order.
    """
    payload = manifest_read(analysis_root)
    stages = [s for s in payload.get("stages", []) if s.get("name") != stage]
    entry: dict[str, Any] = {
        "name": stage,
        "params": params,
        "input_sig": input_sig,
        "output_path": output_path,
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    if extra:
        entry["extra"] = extra
    stages.append(entry)
    payload["stages"] = stages
    manifest_write(analysis_root, payload)


def manifest_find(
    analysis_root: Path | str,
    stage: str,
) -> Optional[dict[str, Any]]:
    """Return the manifest entry for ``stage`` if present, else None."""
    for entry in manifest_read(analysis_root).get("stages", []):
        if entry.get("name") == stage:
            return entry
    return None


def input_signature(*items: Any) -> str:
    """Stable hash of an arbitrary tuple of (paths, params, lists).

    Paths are fingerprinted by (resolved_path, size, mtime_ns) via
    :func:`file_signature` so that touching an upstream file invalidates
    downstream stages. Other items are JSON-serialised with sort_keys=True
    so dict-equivalent params hash to the same value.
    """
    parts: list[str] = []
    for item in items:
        if isinstance(item, (str, Path)) and Path(item).exists():
            sig = file_signature(Path(item))
            parts.append(json.dumps(sig, sort_keys=True))
        else:
            try:
                parts.append(json.dumps(item, sort_keys=True, default=str))
            except TypeError:
                parts.append(repr(item))
    raw = "\n".join(parts).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def upstream_sample_signature(input_sample_dir: Path) -> dict[str, Any]:
    """Fingerprint of an upstream sample directory in a partitioned store.

    Prefers any existing pipeline manifest inside the dir (raw / filtered /
    normalized) -- its content already captures the upstream lineage cheaply.
    Falls back to the on-disk chrom partition listing if no manifest exists.
    """
    for name in (
        ".epykit_raw_manifest.json",
        ".epykit_filter_manifest.json",
        ".epykit_normalize_manifest.json",
    ):
        mp = input_sample_dir / name
        if mp.exists():
            return {"manifest": name, "content": load_json(mp)}

    chroms = expected_chrom_dirs(input_sample_dir)
    return {
        "manifest": None,
        "chroms": chroms,
        "parts": {
            chrom: [
                file_signature(p)
                for p in sorted((input_sample_dir / chrom).glob("part-*.parquet"))
            ]
            for chrom in chroms
        },
    }
