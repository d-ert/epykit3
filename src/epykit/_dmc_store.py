"""Persistent, per-chromosome DMC result store.

DMC results scale linearly with genome size: ~22M CpG rows for a full
human WGBS dataset. The historical pattern of staging each chromosome
to a temp parquet and then ``pl.concat`` ing all 174 files into one
in-memory DataFrame breaks at the genome scale -- the eager concat
holds the list of per-chrom frames *and* the assembled output
simultaneously, and downstream consumers (BH correction, DMR sweep)
duplicate that table again.

``DMCStore`` is a thin handle around a directory of per-chromosome
parquet files plus a manifest. It lets BH and DMR stream the table
chromosome-by-chromosome (peak ~50 MB per chrom) and only materialises
the full DataFrame when a caller explicitly asks for one via
:meth:`to_dataframe`.

Layout::

    <store_dir>/
        .epykit_dmc_manifest.json
        chrom=chr1.parquet
        chrom=chr2.parquet
        ...

The manifest schema is defined in :data:`_MANIFEST_NAME` and written by
:func:`process_chromosomes_dmc` once all per-chrom files land. ``DMCStore``
can be re-opened on a populated directory via :meth:`DMCStore.open`.
"""

from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional

import polars as pl

from . import _cache

logger = logging.getLogger(__name__)


_MANIFEST_NAME = ".epykit_dmc_manifest.json"


def _chrom_filename(chrom: str) -> str:
    return f"chrom={chrom}.parquet"


@dataclass(frozen=True)
class DMCStore:
    """Handle to a persistent per-chromosome DMC result directory.

    Instances are immutable; the only mutation operation is
    :meth:`update_chrom`, which rewrites a single per-chrom parquet
    (used by streaming BH correction to attach qvalues).

    Attributes
    ----------
    path : Path
        Directory holding ``chrom=*.parquet`` files and the manifest.
    test : str
        Statistical test name (``"lr"``, ``"score"``, ...). Mirrors
        :func:`process_chromosomes_dmc`'s ``test`` argument; used for
        diagnostics and cache keys.
    """

    path: Path
    test: str
    _manifest: dict = field(default_factory=dict, repr=False, compare=False)

    @classmethod
    def open(cls, path: str | Path) -> "DMCStore":
        """Re-open an existing store directory by reading its manifest."""
        p = Path(path)
        manifest_path = p / _MANIFEST_NAME
        manifest = _cache.load_json(manifest_path)
        if manifest is None:
            raise FileNotFoundError(
                f"No DMC manifest at {manifest_path}. Either the directory "
                f"wasn't populated by process_chromosomes_dmc(), or the "
                f"manifest was lost. Re-run DMC to regenerate."
            )
        return cls(path=p, test=manifest.get("test", "unknown"), _manifest=manifest)

    # Manifest access

    @property
    def manifest(self) -> dict:
        return dict(self._manifest)

    def chroms(self) -> list[str]:
        """Chromosome names in manifest order (= submission order)."""
        return [c["name"] for c in self._manifest.get("chroms", [])]

    @property
    def total_sites(self) -> int:
        return int(self._manifest.get("total_sites", 0))

    @property
    def bh_applied(self) -> bool:
        return bool(self._manifest.get("bh_qvalues_applied", False))

    def __len__(self) -> int:
        return self.total_sites

    # Per-chrom IO

    def _chrom_path(self, chrom: str) -> Path:
        return self.path / _chrom_filename(chrom)

    def read_chrom(
        self,
        chrom: str,
        columns: Optional[list[str]] = None,
    ) -> pl.DataFrame:
        """Read one chromosome's parquet eagerly."""
        return pl.read_parquet(str(self._chrom_path(chrom)), columns=columns)

    def scan_chrom(self, chrom: str) -> pl.LazyFrame:
        """Lazy scan of one chromosome's parquet."""
        return pl.scan_parquet(str(self._chrom_path(chrom)))

    def iter_chroms(
        self,
        columns: Optional[list[str]] = None,
    ) -> Iterator[tuple[str, pl.DataFrame]]:
        """Yield ``(chrom, per-chrom DataFrame)`` in manifest order."""
        for entry in self._manifest.get("chroms", []):
            chrom = entry["name"]
            yield chrom, self.read_chrom(chrom, columns=columns)

    def to_dataframe(
        self,
        columns: Optional[list[str]] = None,
        *,
        preserve_enum_dtypes: bool = False,
    ) -> pl.DataFrame:
        """Eagerly assemble all chromosomes into one DataFrame.

        Use sparingly on genome-scale stores -- this is the operation
        that DMC's old assembly step performed and the one this class
        exists to avoid. Provided for back-compat with callers that
        truly need the full table (e.g. ``cli.py``'s
        ``results.write_parquet(...)``).

        ``chrom`` and ``strand`` are cast back to ``pl.Utf8`` by default
        so the returned DataFrame is drop-in compatible with code that
        joins on ``chrom`` (Polars rejects Enum-vs-Utf8 join keys). Pass
        ``preserve_enum_dtypes=True`` if you want the Enum dtype
        preserved -- saves ~10x memory on the chrom column but breaks
        joins with Utf8-keyed tables.
        """
        chroms = self.chroms()
        if not chroms:
            return pl.DataFrame()
        frames = []
        for entry in self._manifest.get("chroms", []):
            chrom = entry["name"]
            frames.append(self.read_chrom(chrom, columns=columns))
        df = pl.concat(frames, how="vertical")
        if not preserve_enum_dtypes:
            cast_exprs = []
            for col in ("chrom", "strand"):
                if col in df.columns and df.schema[col] != pl.Utf8:
                    cast_exprs.append(pl.col(col).cast(pl.Utf8))
            if cast_exprs:
                df = df.with_columns(cast_exprs)
        return df

    def update_chrom(self, chrom: str, df: pl.DataFrame) -> None:
        """Atomically overwrite one chromosome's parquet.

        Used by BH correction to attach qvalue / reject columns without
        materialising the full table. The atomic rename keeps the store
        consistent if the process is killed mid-write.
        """
        target = self._chrom_path(chrom)
        tmp = target.with_suffix(target.suffix + ".tmp")
        df.write_parquet(str(tmp))
        tmp.replace(target)

    def write_manifest(self, payload: dict) -> None:
        """Write the manifest atomically and update the cached copy.

        Mutates ``self._manifest`` in place because the dataclass is
        frozen at the attribute level (``object.__setattr__``-style),
        which is fine for the internal dict.
        """
        _cache.write_json(self.path / _MANIFEST_NAME, payload)
        # Replace contents of the cached dict without rebinding the attribute.
        self._manifest.clear()
        self._manifest.update(payload)

    def mark_bh_applied(self, qvalue_col: str = "qvalue", method: str = "fdr_bh") -> None:
        """Flip the manifest's ``bh_qvalues_applied`` flag and persist."""
        payload = dict(self._manifest)
        payload["bh_qvalues_applied"] = True
        payload["bh_qvalue_col"] = qvalue_col
        payload["bh_method"] = method
        self.write_manifest(payload)

    def cleanup(self) -> None:
        """Remove the store directory and all per-chrom files.

        For one-shot uses (permutation FDR, tile-DMR's internal DMC pass)
        that don't want to leave artifacts on disk.
        """
        if self.path.exists():
            shutil.rmtree(self.path, ignore_errors=True)


def has_qvalue(store_or_df: "DMCStore | pl.DataFrame") -> bool:
    """Return True if ``store_or_df`` carries a ``qvalue`` column."""
    if isinstance(store_or_df, DMCStore):
        for entry in store_or_df.manifest.get("chroms", []):
            schema = pl.read_parquet_schema(
                str(store_or_df.path / _chrom_filename(entry["name"]))
            )
            return "qvalue" in schema
        return False
    return "qvalue" in store_or_df.columns
