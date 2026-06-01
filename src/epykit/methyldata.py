from __future__ import annotations

from dataclasses import dataclass, field
import json
import logging
import os
import shutil
from pathlib import Path
from typing import Optional
import warnings

import polars as pl

logger = logging.getLogger(__name__)


@dataclass
class MethylData:
    """Central data object for WGBS methylation analysis.

    Preprocessing state (``filtered``, ``united``, ``smoothed``) is *derived*
    from ``uns`` and its ``_store_history`` log rather than stored as
    independent booleans -- see the ``state`` property and the ``_filtered``
    / ``_united`` / ``_smoothed`` aliases below. This means the flags can
    never drift from reality.
    """

    obs: pl.DataFrame
    store: str
    assembly: str = "unknown"
    context: str = "CpG"

    varm: dict[str, pl.DataFrame] = field(default_factory=dict)
    uns: dict = field(default_factory=dict)

    analysis_root: Optional[str] = field(default=None, repr=False)

    # --- Deprecated alias -------------------------------------------

    @property
    def _analysis_root(self) -> Optional[str]:
        """Deprecated; use ``analysis_root`` (the public name).

        Removed in 2.0.
        """
        warnings.warn(
            "MethylData._analysis_root is deprecated; use the public name "
            "MethylData.analysis_root. _analysis_root will be removed in 2.0.",
            DeprecationWarning, stacklevel=2,
        )
        return self.analysis_root

    @_analysis_root.setter
    def _analysis_root(self, value: Optional[str]) -> None:
        warnings.warn(
            "MethylData._analysis_root is deprecated; use the public name "
            "MethylData.analysis_root. _analysis_root will be removed in 2.0.",
            DeprecationWarning, stacklevel=2,
        )
        self.analysis_root = value

    # --- State (derived from uns) -------------------------------------

    @property
    def _filtered(self) -> bool:
        """True iff filter_coverage has been run (recorded in store history)."""
        history = self.uns.get("_store_history", [])
        return any(h.get("step") == "filtered" for h in history)

    @property
    def _united(self) -> bool:
        """True iff pp.unite has been called (records md.uns['unite'])."""
        return "unite" in self.uns

    @property
    def _smoothed(self) -> bool:
        """True iff pp.smooth has been called (records md.uns['smooth_path'])."""
        return "smooth_path" in self.uns

    @property
    def state(self) -> list[str]:
        """Ordered list of preprocessing steps applied to this object.

        Reads ``uns["_store_history"]`` for store-mutating steps (filtered,
        normalized) and appends ``united`` / ``smoothed`` if recorded in
        ``uns``. The result is suitable for ``__repr__`` / ``_repr_html_``.
        """
        history = self.uns.get("_store_history", [])
        steps = [h.get("step") for h in history if h.get("step")]
        if self._united and "united" not in steps:
            steps.append("united")
        if self._smoothed and "smoothed" not in steps:
            steps.append("smoothed")
        return steps

    @property
    def treatment_ids(self) -> list[str]:
        if "treatment" not in self.obs.columns:
            return []
        return (
            self.obs
            .filter(pl.col("treatment") == 1)
            .get_column("sample_id")
            .to_list()
        )

    @property
    def control_ids(self) -> list[str]:
        if "treatment" not in self.obs.columns:
            return []
        return (
            self.obs
            .filter(pl.col("treatment") == 0)
            .get_column("sample_id")
            .to_list()
        )

    @property
    def n_samples(self) -> int:
        return len(self.obs)

    @property
    def completed_stages(self) -> list[str]:
        """Names of stages recorded in the top-level pipeline manifest.

        Reads ``<analysis_root>/.epykit_manifest.json`` (the 0.4.0
        checkpoint/resume manifest). Always reflects what is on disk --
        unlike :py:attr:`state`, which is derived from ``uns`` and can
        diverge after a manual ``md.uns.pop()``.

        Returns an empty list when the analysis root has no manifest
        (e.g. an in-memory MethylData built by ``read_bismark`` without
        the ``store_dir`` argument).
        """
        from ._cache import manifest_read
        root = self.analysis_root or self.store
        if not root:
            return []
        try:
            return [s.get("name", "") for s in manifest_read(root).get("stages", [])]
        except (OSError, ValueError):
            return []

    def resume_from(self, stage: str) -> bool:
        """Re-hydrate ``uns`` / ``varm`` state for a completed stage.

        Looks up ``stage`` in the pipeline manifest. If present, restores
        any sidecar parquet results referenced by the manifest entry's
        ``output_path`` into the matching ``varm`` / ``uns`` slot and
        returns True. Returns False (without modifying anything) when
        no such stage was recorded.

        This is the read side of the 0.4.0 checkpoint/resume API; the
        write side is each ``pp.*`` / ``tl.*`` function appending to the
        manifest when invoked with ``resumable=True``.
        """
        from ._cache import manifest_find
        from pathlib import Path
        root = self.analysis_root or self.store
        if not root:
            return False
        entry = manifest_find(root, stage)
        if entry is None:
            return False
        out_path = entry.get("output_path")
        if not out_path:
            return False
        op = Path(out_path)
        if not op.is_absolute():
            op = Path(root) / op
        # varm/<key>.parquet -> re-load into varm
        if op.exists() and op.suffix == ".parquet" and stage.startswith(("dmc_", "dvc", "asm", "entropy")):
            self.varm[stage] = pl.read_parquet(str(op))
            return True
        # uns/<key>.parquet -> re-load into uns
        if op.exists() and op.suffix == ".parquet" and stage.startswith(("dmr_", "dvr", "pmd", "hmr", "lmr", "smooth")):
            self.uns[stage] = pl.read_parquet(str(op))
            return True
        # Stage references a directory (filtered/normalized stores etc.) --
        # just point md.store at it.
        if op.is_dir():
            self.store = str(op)
            return True
        return False

    def get_dmc(
        self,
        test: Optional[str] = None,
        annotated: bool = True,
    ) -> Optional[pl.DataFrame]:
        """Look up a DMC table by test name (explicit, recommended).

        Parameters
        ----------
        test : str, optional
            Test backend name (``"lr"``, ``"glm"``, ``"welch_t"``, ...). When
            ``None`` (default), returns the most-recently-written DMC table,
            as recorded by ``ep.tl.dmc`` in ``md.uns["dmc"]["last_key"]``.
        annotated : bool
            When True (default), prefer the ``*_annotated`` variant of the
            requested table if it exists (so plotting code that needs
            ``feature_type`` / ``cpg_context`` works out of the box).

        Returns
        -------
        pl.DataFrame or None
            The matching table, or None if no DMC has been run.
        """
        if test is None:
            key = self.uns.get("dmc", {}).get("last_key")
            if key is None:
                return None
        else:
            key = f"dmc_{test}"
        if annotated:
            ann_key = f"{key}_annotated"
            if ann_key in self.varm:
                return self.varm[ann_key]
        return self.varm.get(key)

    @property
    def dmc(self) -> Optional[pl.DataFrame]:
        """Most-recently-written DMC table (annotated if available).

        Equivalent to ``self.get_dmc(test=None, annotated=True)``. Use
        :meth:`get_dmc` with an explicit ``test=`` argument when running
        multiple tests in one session and you need a specific one. The
        legacy auto-pick-by-priority behaviour (glm > lr > score > ...) has
        been removed because it silently disagreed with the user's most
        recent call.
        """
        # Pointer-first resolution: ep.tl.dmc writes uns["dmc"]["last_key"]
        # on every run. If that's absent (older sessions), fall back to a
        # single existing key -- but never auto-prioritize, to avoid the
        # surprise documented in S5.
        last = self.get_dmc()
        if last is not None:
            return last
        # Fallback: if only one dmc_* table is present, return it.
        dmc_keys = [k for k in self.varm if k.startswith("dmc") and not k.endswith("_annotated")]
        if len(dmc_keys) == 1:
            key = dmc_keys[0]
            ann_key = f"{key}_annotated"
            return self.varm.get(ann_key, self.varm.get(key))
        return None

    @property
    def significant_dmcs(self) -> Optional[pl.DataFrame]:
        df = self.dmc
        if df is None:
            return None
        if "qvalue" in df.columns:
            return df.filter(pl.col("qvalue") < 0.05)
        if "pvalue" in df.columns:
            return df.filter(pl.col("pvalue") < 0.05)
        return df

    def save(self, path: str) -> None:
        """Persist obs/varm/uns + manifest to disk.

        Path interpretation:

        * If ``path`` contains any directory components (relative or
          absolute), the data is written there verbatim. ``load(path)``
          reads from the same place -- save and load are symmetric.
        * If ``path`` is a bare name (no separators) **and**
          ``analysis_root`` is set, the data is written under
          ``<analysis_root>/results/<path>``. This is the
          "analysis-project" convenience layout.
        * If ``path`` is a bare name with no ``analysis_root``, it's
          treated as a relative path in the current directory.

        The previous behaviour silently re-rooted every call (even
        absolute paths) under ``<analysis_root>/results/<basename>``,
        which broke save/load symmetry.
        """
        p = Path(path)
        has_components = p.is_absolute() or len(p.parts) > 1
        if self.analysis_root and not has_components:
            out = Path(self.analysis_root) / "results" / p.name
        else:
            out = p
        out.mkdir(parents=True, exist_ok=True)

        self.obs.write_parquet(str(out / "obs.parquet"))

        # Resolve the DMCStore reference (if any) so we can detect varm
        # entries that are *already* on disk in chrom-partitioned form
        # and just need to be linked / copied rather than re-encoded.
        # Encoding a 22M-row Polars DataFrame to a single parquet allocates
        # a second-copy buffer the same size as the table (~3 GB) -- that's
        # what makes naive `df.write_parquet(...)` OOM the host.
        dmc_meta       = self.uns.get("dmc", {}) if isinstance(self.uns.get("dmc"), dict) else {}
        dmc_store_path = dmc_meta.get("store_path")
        dmc_last_key   = dmc_meta.get("last_key")

        varm_format: dict[str, str] = {}
        for name, df in self.varm.items():
            target_single = out / f"varm_{name}.parquet"

            # Path 1: DMCStore-backed varm table -> copy per-chrom parquets
            # directly. Zero re-encoding, constant memory. The DMCStore
            # already carries BH-corrected q-values (apply_multiple_testing_
            # correction writes them back per chrom), so the on-disk files
            # are exactly the table held in self.varm[name].
            is_dmcstore_backed = (
                name == dmc_last_key
                and dmc_store_path is not None
                and (Path(dmc_store_path) / ".epykit_dmc_manifest.json").exists()
            )
            if is_dmcstore_backed:
                store_dir = Path(dmc_store_path)
                target_dir = out / f"varm_{name}"
                target_dir.mkdir(parents=True, exist_ok=True)
                # Hardlink where possible (no extra disk, instant); fall
                # back to a chunked file copy elsewhere (different drive /
                # filesystem that doesn't support hardlinks).
                for src in store_dir.glob("chrom=*.parquet"):
                    dst = target_dir / src.name
                    if dst.exists():
                        dst.unlink()
                    try:
                        os.link(src, dst)
                    except (OSError, NotImplementedError):
                        shutil.copyfile(src, dst)
                # Carry the manifest along so load() can verify integrity.
                manifest = store_dir / ".epykit_dmc_manifest.json"
                if manifest.exists():
                    shutil.copyfile(manifest, target_dir / manifest.name)
                varm_format[name] = "dmcstore"
                logger.info(
                    "save: %s linked from DMCStore at %s (no materialization)",
                    name, store_dir,
                )
                continue

            # Path 2: legacy / small varm -> single-file parquet.
            df.write_parquet(str(target_single))
            varm_format[name] = "parquet"

        serialisable_uns = self.uns.copy()
        for key, value in list(serialisable_uns.items()):
            if isinstance(value, pl.DataFrame):
                parquet_name = f"uns_{key}.parquet"
                value.write_parquet(str(out / parquet_name))
                serialisable_uns[key] = {"__parquet__": parquet_name}

        meta = {
            "store": self.store,
            "assembly": self.assembly,
            "context": self.context,
            # _filtered / _united / _smoothed are derived from uns; don't
            # persist them. They are recomputed on load() from the loaded
            # uns dict (which includes _store_history, unite, smooth_path).
            "_analysis_root": self.analysis_root,
            "varm_keys": list(self.varm.keys()),
            # Per-varm storage format: "parquet" (single file) or
            # "dmcstore" (chrom-partitioned dir). Older saves omit this
            # field and load() falls back to "parquet" for back-compat.
            "varm_format": varm_format,
            "uns": serialisable_uns,
        }
        (out / "methyldata.json").write_text(json.dumps(meta, indent=2, default=str))

    @classmethod
    def load(cls, path: str) -> "MethylData":
        out = Path(path)
        meta = json.loads((out / "methyldata.json").read_text())
        obs = pl.read_parquet(str(out / "obs.parquet"))
        varm_format = meta.get("varm_format", {})
        varm: dict[str, pl.DataFrame] = {}
        for key in meta.get("varm_keys", []):
            fmt = varm_format.get(key, "parquet")
            if fmt == "dmcstore":
                # Chrom-partitioned directory written by save() via direct
                # link-from-DMCStore. Streaming scan keeps the read on a
                # single pass and avoids a second-copy materialisation
                # buffer; we still collect into an eager frame for
                # back-compat with all the code that expects md.varm[key]
                # to be a pl.DataFrame.
                varm_dir = out / f"varm_{key}"
                varm[key] = (
                    pl.scan_parquet(str(varm_dir / "chrom=*.parquet"))
                    .collect()
                )
            else:
                varm[key] = pl.read_parquet(str(out / f"varm_{key}.parquet"))

        uns = meta.get("uns", {})
        for key, value in list(uns.items()):
            if isinstance(value, dict) and "__parquet__" in value:
                uns[key] = pl.read_parquet(str(out / value["__parquet__"]))

        md = cls(
            obs=obs,
            store=meta.get("store", ""),
            assembly=meta.get("assembly", "unknown"),
            context=meta.get("context", "CpG"),
            varm=varm,
            uns=uns,
            # _filtered / _united / _smoothed are properties derived from
            # uns -- nothing to pass through the constructor. Older saves
            # that include those keys in meta are silently ignored.
        )
        md.analysis_root = meta.get("_analysis_root")
        return md

    # --- Exports / reports --------------------------------------------

    def report(self, output: str, **kwargs) -> str:
        """Render a self-contained interactive HTML report.

        Thin wrapper for :func:`epykit.report.generate_report`. See its
        docstring for ``title``, ``gtf_path``, ``alpha`` etc.
        """
        from .report import generate_report
        return generate_report(self, output, **kwargs)

    def to_bedgraph(self, sample: str, output: str, *, value: str = "beta") -> str:
        from .export import to_bedgraph
        return to_bedgraph(self, sample, output, value=value)

    def to_bigwig(
        self, sample: str, output: str, *, value: str = "beta",
        chrom_sizes: Optional[dict] = None,
    ) -> str:
        from .export import to_bigwig
        return to_bigwig(self, sample, output, value=value, chrom_sizes=chrom_sizes)

    def dmcs_to_bed(
        self, output: str, *, alpha: float = 0.05,
        min_abs_diff: float = 0.0, test: Optional[str] = None,
    ) -> str:
        from .export import dmcs_to_bed
        return dmcs_to_bed(self, output, alpha=alpha,
                           min_abs_diff=min_abs_diff, test=test)

    def dmrs_to_bed(self, output: str) -> str:
        from .export import dmrs_to_bed
        return dmrs_to_bed(self, output)

    def to_anndata(self, **kwargs):
        """Return an AnnData of this MethylData (requires `pp.unite` first)."""
        from .anndata_io import to_anndata
        return to_anndata(self, **kwargs)

    def to_mudata(self, **kwargs):
        """Return a MuData with methylation as the ``'meth'`` modality."""
        from .anndata_io import to_mudata
        return to_mudata(self, **kwargs)

    def to_methylkit_tabix(self, output_dir: str, samples=None):
        """Export per-sample methylKit tabix-friendly tables."""
        from .methylkit_io import to_methylkit_tabix
        return to_methylkit_tabix(self, output_dir, samples=samples)

    def region_beta(
        self,
        chrom: str,
        start: int,
        end: int,
    ) -> pl.DataFrame:
        """Per-sample mean beta within ``chrom:start-end``.

        Returns columns: sample, mean_beta, n_cpgs, mean_coverage.
        """
        rows: list[dict] = []
        samples = self.obs.get_column("sample_id").to_list()
        store_root = Path(self.store)
        for s in samples:
            part = store_root / f"sample={s}" / f"chrom={chrom}" / "part-0.parquet"
            if not part.exists():
                rows.append({
                    "sample": s, "mean_beta": float("nan"),
                    "n_cpgs": 0, "mean_coverage": float("nan"),
                })
                continue
            df = (
                pl.read_parquet(str(part), columns=["pos", "N_meth", "coverage"])
                .filter((pl.col("pos") >= start) & (pl.col("pos") <= end))
            )
            if len(df) == 0:
                rows.append({
                    "sample": s, "mean_beta": float("nan"),
                    "n_cpgs": 0, "mean_coverage": float("nan"),
                })
                continue
            cov_sum = int(df.get_column("coverage").sum())
            meth_sum = int(df.get_column("N_meth").sum())
            rows.append({
                "sample": s,
                "mean_beta": float(meth_sum / max(cov_sum, 1)),
                "n_cpgs": int(len(df)),
                "mean_coverage": float(cov_sum / max(len(df), 1)),
            })
        return pl.DataFrame(rows)

    def __repr__(self) -> str:
        n_sites = self.uns.get("n_sites_filtered") or self.uns.get("n_sites_raw", "?")
        groups = "unknown"
        if "group" in self.obs.columns:
            grouped = self.obs.group_by("group").len().sort("group")
            groups = "  ".join(
                f"{r['group']} (n={r['len']})"
                for r in grouped.iter_rows(named=True)
            )

        status_str = ", ".join(self.state) if self.state else "raw"

        varm_str = ", ".join(self.varm.keys()) if self.varm else "none yet"
        uns_keys = ", ".join(sorted(self.uns.keys())) if self.uns else "none"

        n_sites_str = f"{n_sites:,}" if isinstance(n_sites, int) else str(n_sites)
        return (
            f"MethylData [{status_str}]\n"
            f"  assembly : {self.assembly}\n"
            f"  n_samples: {self.n_samples} ({groups})\n"
            f"  n_sites  : {n_sites_str}\n"
            f"  context  : {self.context}\n"
            f"  store    : {self.store}\n"
            f"  varm     : {varm_str}\n"
            f"  uns      : {uns_keys}\n"
        )

    def _repr_html_(self) -> str:
        """Render obs as a notebook-friendly table.

        Shows every column in ``self.obs`` rather than a hardcoded subset, so
        user-supplied covariates (sex, batch, age, ...) are visible. Floats are
        rounded to 4 significant figures; the ``treatment`` column, if present,
        renders as > (1) / o (0).
        """
        cols = list(self.obs.columns)

        def _fmt(value: object, col: str) -> str:
            if value is None:
                return "--"
            if col == "treatment":
                return ">" if value == 1 else "o" if value == 0 else str(value)
            if isinstance(value, float):
                if value != value:  # NaN
                    return "--"
                return f"{value:.4g}"
            return str(value)

        rows_html: list[str] = []
        for row in self.obs.iter_rows(named=True):
            cells = "".join(f"<td>{_fmt(row.get(c), c)}</td>" for c in cols)
            rows_html.append(f"<tr>{cells}</tr>")

        header_html = "".join(f"<th>{c}</th>" for c in cols)
        status = ", ".join(self.state) if hasattr(self, "state") and self.state else "raw"

        return f"""
        <div style="font-family:monospace;font-size:13px">
        <b>MethylData</b> [{status}] | assembly: {self.assembly} | context: {self.context}<br>
        <table border="1" style="border-collapse:collapse;margin:8px 0">
          <tr>{header_html}</tr>
          {''.join(rows_html)}
        </table>
        Results: {', '.join(self.varm.keys()) or 'none yet'}
        </div>
        """
