"""Process-wide runtime configuration helpers.

Currently exposes one knob: :func:`set_tmp_dir`, which redirects every
internal ``tempfile.TemporaryDirectory`` call (DMC chrom staging, DMR
tile aggregation, DVC streaming, smoothed-DMC pseudo-count store, etc.)
to a user-chosen directory. Useful when the default OS tempdir (Windows:
``%TEMP%`` on ``C:\``) sits on a small or full drive and you want
transient work to go elsewhere without setting environment variables.
"""

from __future__ import annotations

import logging
import os
import tempfile
from pathlib import Path
from typing import Optional, Union

logger = logging.getLogger(__name__)


def set_tmp_dir(path: Optional[Union[str, Path]]) -> str:
    """Redirect all of epykit's transient temp files to ``path``.

    Sets :data:`tempfile.tempdir` to ``path`` (created if it doesn't
    exist) AND mirrors the value into the ``TMPDIR``/``TEMP``/``TMP``
    environment variables so any subprocesses epykit spawns (Dask
    workers, Ray actors) inherit the same setting.

    Pass ``None`` to restore the OS default behaviour.

    Examples
    --------
    Redirect all epykit transient work to a path on a larger drive::

        import epykit as ep
        ep.set_tmp_dir("D:/work/epykit_tmp")
        # ... every subsequent tl.dmc / tl.dmr / pp.smooth call now
        # stages its per-chrom parquet files under D:/work/epykit_tmp.

    Returns
    -------
    str
        The effective temp directory after the call (the OS default
        when ``path`` is None, otherwise the resolved absolute path).
    """
    if path is None:
        tempfile.tempdir = None
        for var in ("TMPDIR", "TEMP", "TMP"):
            os.environ.pop(var, None)
        effective = tempfile.gettempdir()
        logger.info("epykit tmp dir reset to OS default: %s", effective)
        return effective

    resolved = Path(path).expanduser().resolve()
    resolved.mkdir(parents=True, exist_ok=True)
    s = str(resolved)
    tempfile.tempdir = s
    for var in ("TMPDIR", "TEMP", "TMP"):
        os.environ[var] = s
    logger.info("epykit tmp dir set to: %s", s)
    return s


def get_tmp_dir() -> str:
    """Return the directory currently used for epykit's transient files."""
    return tempfile.tempdir or tempfile.gettempdir()


__all__ = ["set_tmp_dir", "get_tmp_dir"]
