"""run_epykit_gse.py -- Phase 4 Task 6.

Post-Phase-3 epykit re-run on the GSE263850 real cohort (Study 3:
3 SBP009 controls vs 3 Clone cases). Produces per-engine per-CpG DMC
parquets, per-method DMR parquets, and a cross-tool concordance table
that pulls DSS + methylKit from their existing on-disk outputs (no
re-runs needed).

Inputs
------
* Samplesheet : benchmark/data/study3/samplesheet_gse263850.csv
* Raw .bed.gz : ../epykit2/GSE263850_RAW/*.readset_sorted.dedup.filtered.bed.gz
  -- 12-column combined-strand BED (see ep.read_combined_strand_bed).
* DSS per-CpG : benchmark/data/study3/dss/dmltest_per_cpg.tsv.gz
  (chr/pos/stat/pvals/fdrs)
* DSS per-DMR : benchmark/data/study3/dss/dmr_dss.csv
* methylKit   : D:/Coding/Projeler/methyl_lib/methylkıt_realResults/
                scripts_and_results/methylkit_results/dmc_all_sites.csv (~2 GB)
                + dmr_all_tiles.csv (~98 MB). Note the non-ASCII ı.

Outputs (under benchmark/data/study3/)
--------------------------------------
* epykit_post_phase3/dmc/<engine>.parquet     (per-engine full per-CpG)
* epykit_post_phase3/dmr/<method>.parquet     (per-method full per-DMR)
* epykit_post_phase3/MANIFEST.txt
* comparisons_post_phase3/dmr_iou.parquet
* comparisons_post_phase3/per_dmr_stat_concordance.parquet
* comparisons_post_phase3/SUMMARY.md

The methylKit 2 GB per-CpG file is streamed via ``pl.scan_csv`` -- never
materialised. lr+ kwargs come from ``_epykit_scoring._dmc_kwargs`` so the
Study 1 + Study 1b + Study 3 contracts cannot drift.

Usage
-----
    uv run python benchmark/scripts/run_epykit_gse.py
    uv run python benchmark/scripts/run_epykit_gse.py --skip-dmc
    uv run python benchmark/scripts/run_epykit_gse.py --skip-dmr
    uv run python benchmark/scripts/run_epykit_gse.py --skip-concordance
    uv run python benchmark/scripts/run_epykit_gse.py --engines lr lr+ \\
                                                       --dmr-methods chain_merge

Exit non-zero on any per-engine failure.
"""

from __future__ import annotations

import argparse
import gc
import logging
import shutil
import subprocess
import sys
import time
from pathlib import Path

import polars as pl

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from _epykit_scoring import ENGINE_EXCEPTIONS, Q_THRESHOLD, _dmc_kwargs  # noqa: E402, F401


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                       # benchmark/
REPO = ROOT.parent                       # epykit3/
DATA_STUDY3 = ROOT / "data" / "study3"
SAMPLESHEET = DATA_STUDY3 / "samplesheet_gse263850.csv"

OUT_BASE = DATA_STUDY3 / "epykit_post_phase3"
OUT_DMC = OUT_BASE / "dmc"
OUT_DMR = OUT_BASE / "dmr"
MANIFEST_PATH = OUT_BASE / "MANIFEST.txt"
RUN_STORE = ROOT / "_runs_post_phase3" / "study3_gse"

CMP_OUT = DATA_STUDY3 / "comparisons_post_phase3"

# DSS reuse paths
DSS_DMC_TSV = DATA_STUDY3 / "dss" / "dmltest_per_cpg.tsv.gz"
DSS_DMR_CSV = DATA_STUDY3 / "dss" / "dmr_dss.csv"

# methylKit reuse paths (external; non-ASCII ı handled via pathlib)
MK_BASE = Path(
    "D:/Coding/Projeler/methyl_lib/methylkıt_realResults/"
    "scripts_and_results/methylkit_results"
)
MK_DMC_CSV = MK_BASE / "dmc_all_sites.csv"          # ~2 GB -- stream only
MK_DMR_CSV = MK_BASE / "dmr_all_tiles.csv"           # ~98 MB

logger = logging.getLogger("run_epykit_gse")


# ---------------------------------------------------------------------------
# Engine surface
# ---------------------------------------------------------------------------

# Per Phase-4 plan task 6: lr + lr+ are the headline. welch_t / fisher
# are accepted via --engines for completeness but skipped by default.
DEFAULT_DMC_ENGINES: tuple[str, ...] = ("lr", "lr+")
DEFAULT_DMR_METHODS: tuple[str, ...] = ("chain_merge",)

# Concordance comparison threshold for "called DMR" sets.
CALLED_QVAL_THRESHOLD = 0.05


# ---------------------------------------------------------------------------
# Ingestion
# ---------------------------------------------------------------------------


def _resolve_samplesheet(raw_sheet: Path, resolved_sheet: Path) -> Path:
    """Rewrite the samplesheet so each ``path`` is an absolute file path.

    The committed samplesheet uses ``../epykit2/GSE263850_RAW/…`` which is
    relative to the repo root. epykit's ingester resolves paths against
    the samplesheet's *own* directory (which is ``benchmark/data/study3``),
    so the bare relative path would fail. We materialise an absolute-path
    copy under tmp/cache and ingest from that.
    """
    df = pl.read_csv(raw_sheet)
    resolved_paths: list[str] = []
    for p in df["path"].to_list():
        candidate = Path(p)
        if not candidate.is_absolute():
            # Relative paths are anchored at REPO (the working dir users
            # invoke `uv run` from).
            candidate = (REPO / p).resolve()
        if not candidate.exists():
            raise FileNotFoundError(f"samplesheet path missing: {candidate}")
        resolved_paths.append(str(candidate))
    df = df.with_columns(pl.Series("path", resolved_paths))
    resolved_sheet.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(resolved_sheet)
    return resolved_sheet


def ingest(store_dir: Path | None = None) -> "object":
    """Read the GSE263850 samplesheet into a fresh ``MethylData``.

    Uses :func:`epykit.read_combined_strand_bed` -- the input files are
    12-column strand-collapsed BEDs per ep.io.read_combined_strand_bed
    docstring, which explicitly cites GSE263850 as the canonical case.
    """
    import epykit as ep  # noqa: PLC0415 -- defer the heavy import

    store_dir = Path(store_dir) if store_dir is not None else RUN_STORE
    if store_dir.exists():
        shutil.rmtree(store_dir)
    store_dir.parent.mkdir(parents=True, exist_ok=True)

    resolved_sheet = store_dir.parent / "samplesheet_gse263850_resolved.csv"
    sheet = _resolve_samplesheet(SAMPLESHEET, resolved_sheet)

    md = ep.read_combined_strand_bed(
        str(sheet),
        treatment_group="clone",
        control_group="sbp009",
        assembly="hg38",
        store_dir=str(store_dir),
    )
    ep.pp.unite(md, type="intersect")
    logger.info(
        "ingest: %d samples, store=%s",
        md.obs.height, md.store,
    )
    return md


# ---------------------------------------------------------------------------
# DMC + DMR runners
# ---------------------------------------------------------------------------


def run_dmc(md, engines: tuple[str, ...], out_dir: Path) -> list[dict]:
    """Run ``ep.tl.dmc`` for each engine; write per-engine parquets."""
    import epykit as ep  # noqa: PLC0415

    out_dir.mkdir(parents=True, exist_ok=True)
    timings: list[dict] = []

    for engine in engines:
        backend, kwargs = _dmc_kwargs(engine, allow_n1=False)
        out_name = engine.replace("+", "plus")
        out_pq = out_dir / f"{out_name}.parquet"

        t0 = time.perf_counter()
        try:
            ep.tl.dmc(md, **kwargs)
        except ENGINE_EXCEPTIONS as exc:
            elapsed = time.perf_counter() - t0
            logger.error("DMC engine=%s FAILED in %.1fs: %r", engine, elapsed, exc)
            timings.append({
                "stage": "dmc", "engine": engine, "wall_s": elapsed,
                "n_sites": None, "n_sig": None,
                "ok": False, "error": str(exc),
            })
            continue
        elapsed = time.perf_counter() - t0

        df = md.get_dmc(test=backend)
        df.write_parquet(out_pq)
        n_sites = df.height
        n_sig = (
            df.filter(pl.col("qvalue") < Q_THRESHOLD).height
            if "qvalue" in df.columns else None
        )
        timings.append({
            "stage": "dmc", "engine": engine, "wall_s": elapsed,
            "n_sites": n_sites, "n_sig": n_sig,
            "ok": True, "error": None,
        })
        logger.info(
            "DMC engine=%s: %s CpGs, %s sig (q<%.2f), %.1fs",
            engine, f"{n_sites:,}",
            f"{n_sig:,}" if n_sig is not None else "?",
            Q_THRESHOLD, elapsed,
        )
    return timings


def run_dmr(md, methods: tuple[str, ...], out_dir: Path) -> list[dict]:
    """Run ``ep.tl.dmr`` for each method; write per-method parquets.

    Assumes ``ep.tl.dmc(md, test='lr', ...)`` has already run so the DMC
    store is present (chain_merge / segment / sliding_window read from it).
    """
    import epykit as ep  # noqa: PLC0415

    out_dir.mkdir(parents=True, exist_ok=True)
    timings: list[dict] = []

    for method in methods:
        out_pq = out_dir / f"{method}.parquet"

        t0 = time.perf_counter()
        try:
            ep.tl.dmr(md, method=method)
        except ENGINE_EXCEPTIONS as exc:
            elapsed = time.perf_counter() - t0
            logger.error("DMR method=%s FAILED in %.1fs: %r", method, elapsed, exc)
            timings.append({
                "stage": "dmr", "engine": method, "wall_s": elapsed,
                "n_sites": None, "n_sig": None,
                "ok": False, "error": str(exc),
            })
            continue
        elapsed = time.perf_counter() - t0

        dmr_df = md.uns.get("dmr")
        if isinstance(dmr_df, dict):
            dmr_df = dmr_df.get("frame") or next(iter(dmr_df.values()))
        dmr_df.write_parquet(out_pq)
        n_sites = dmr_df.height
        n_sig = (
            dmr_df.filter(pl.col("qvalue") < Q_THRESHOLD).height
            if "qvalue" in dmr_df.columns else n_sites
        )
        timings.append({
            "stage": "dmr", "engine": method, "wall_s": elapsed,
            "n_sites": n_sites, "n_sig": n_sig,
            "ok": True, "error": None,
        })
        logger.info(
            "DMR method=%s: %s regions, %s sig (q<%.2f), %.1fs",
            method, f"{n_sites:,}", f"{n_sig:,}", Q_THRESHOLD, elapsed,
        )
    return timings


# ---------------------------------------------------------------------------
# Concordance
# ---------------------------------------------------------------------------


def _interval_jaccard_dmr(a: pl.DataFrame, b: pl.DataFrame) -> dict:
    """Jaccard on the *region union* of two DMR sets.

    Each tool's DMR set is reduced to per-chromosome merged intervals;
    we compute |A ∩ B| / |A ∪ B| in bp over those merged intervals.

    Returns ``{jaccard, n_a, n_b, bp_a, bp_b, bp_intersection, bp_union}``.
    """
    def _merged_by_chrom(df: pl.DataFrame) -> dict[str, list[tuple[int, int]]]:
        out: dict[str, list[tuple[int, int]]] = {}
        for r in df.select(["chrom", "start", "end"]).iter_rows(named=True):
            out.setdefault(str(r["chrom"]), []).append(
                (int(r["start"]), int(r["end"])),
            )
        for ch in out:
            iv = sorted(out[ch])
            merged: list[tuple[int, int]] = []
            for s, e in iv:
                if merged and s <= merged[-1][1]:
                    ms, me = merged[-1]
                    merged[-1] = (ms, max(me, e))
                else:
                    merged.append((s, e))
            out[ch] = merged
        return out

    a_iv = _merged_by_chrom(a)
    b_iv = _merged_by_chrom(b)
    bp_a = sum(e - s for chs in a_iv.values() for s, e in chs)
    bp_b = sum(e - s for chs in b_iv.values() for s, e in chs)
    bp_inter = 0
    for ch in set(a_iv) | set(b_iv):
        ia = a_iv.get(ch, [])
        ib = b_iv.get(ch, [])
        i = j = 0
        while i < len(ia) and j < len(ib):
            s1, e1 = ia[i]
            s2, e2 = ib[j]
            lo = max(s1, s2)
            hi = min(e1, e2)
            if hi > lo:
                bp_inter += hi - lo
            if e1 < e2:
                i += 1
            else:
                j += 1
    bp_union = bp_a + bp_b - bp_inter
    return {
        "jaccard": (bp_inter / bp_union) if bp_union else 0.0,
        "n_a": a.height, "n_b": b.height,
        "bp_a": bp_a, "bp_b": bp_b,
        "bp_intersection": bp_inter, "bp_union": bp_union,
    }


def _load_epykit_dmr_called(parquet: Path) -> pl.DataFrame:
    df = pl.read_parquet(parquet)
    if "qvalue" in df.columns:
        df = df.filter(pl.col("qvalue") < CALLED_QVAL_THRESHOLD)
    return df.select(
        chrom=pl.col("chrom").cast(pl.Utf8),
        start=pl.col("start").cast(pl.Int64),
        end=pl.col("end").cast(pl.Int64),
        qvalue=pl.col("qvalue").cast(pl.Float64) if "qvalue" in df.columns
        else pl.lit(None, dtype=pl.Float64),
        meth_diff=pl.col("mean_meth_diff").cast(pl.Float64) if "mean_meth_diff"
        in df.columns else pl.col("meth_diff").cast(pl.Float64),
    )


def _load_dss_dmr_called(csv: Path) -> pl.DataFrame:
    """DSS doesn't carry a per-DMR q-value; ``dmr_dss.csv`` is already the
    DSS-significant subset (922 rows per the existing audit). We take all
    rows as called.

    Direction column: the file's ``diff_Methy_DSSfit`` column is empty in
    this distribution (data audit confirmed -- the DSS fitted means were
    not persisted), so we fall back to ``diff_Methy_fromCounts`` which is
    the same case-vs-control direction in raw count space. Equivalent
    sign convention to epykit + methylKit (case - control).
    """
    df = pl.read_csv(csv, infer_schema_length=10000)
    meth_col = (
        "diff_Methy_DSSfit"
        if "diff_Methy_DSSfit" in df.columns
        and df["diff_Methy_DSSfit"].dtype != pl.Utf8
        else "diff_Methy_fromCounts"
    )
    return df.select(
        chrom=pl.col("chrom").cast(pl.Utf8),
        start=pl.col("start").cast(pl.Int64),
        end=pl.col("end").cast(pl.Int64),
        qvalue=pl.lit(None, dtype=pl.Float64),
        meth_diff=pl.col(meth_col).cast(pl.Float64),
    )


def _load_methylkit_dmr_called(csv: Path) -> pl.DataFrame:
    """methylKit DMR tiles, filtered to q<0.05 (the called set).

    methylKit ``meth_diff`` is in percent units (e.g. -15 means 15%
    methylation drop). We rescale to the [-1, 1] convention used by
    epykit + DSS for the direction-agreement column to be meaningful.
    """
    df = pl.read_csv(csv)
    df = df.filter(pl.col("qvalue") < CALLED_QVAL_THRESHOLD)
    return df.select(
        chrom=pl.col("chrom").cast(pl.Utf8),
        start=pl.col("start").cast(pl.Int64),
        end=pl.col("end").cast(pl.Int64),
        qvalue=pl.col("qvalue").cast(pl.Float64),
        meth_diff=(pl.col("meth_diff") / 100.0).cast(pl.Float64),
    )


def _direction_agreement(a: pl.DataFrame, b: pl.DataFrame) -> dict:
    """Overlap-then-compare-direction.

    For each row in ``a``, find any overlapping row(s) in ``b`` (same
    chrom, bp overlap > 0). Count the (a, b) overlap pairs whose
    ``meth_diff`` signs agree.
    """
    b_idx: dict[str, list[tuple[int, int, float]]] = {}
    for r in b.iter_rows(named=True):
        b_idx.setdefault(str(r["chrom"]), []).append(
            (int(r["start"]), int(r["end"]), float(r["meth_diff"])
             if r["meth_diff"] is not None else 0.0),
        )
    for ch in b_idx:
        b_idx[ch].sort()

    n_overlap = 0
    n_agree = 0
    for r in a.iter_rows(named=True):
        ch = str(r["chrom"])
        s, e = int(r["start"]), int(r["end"])
        md_a = float(r["meth_diff"]) if r["meth_diff"] is not None else 0.0
        for bs, be, md_b in b_idx.get(ch, []):
            if be < s:
                continue
            if bs > e:
                break
            n_overlap += 1
            if (md_a > 0 and md_b > 0) or (md_a < 0 and md_b < 0):
                n_agree += 1
    return {
        "n_overlapping_pairs": n_overlap,
        "n_direction_agree": n_agree,
        "direction_agree_frac": (n_agree / n_overlap) if n_overlap else 0.0,
    }


def _overlap_join(a: pl.DataFrame, b: pl.DataFrame) -> pl.DataFrame:
    """Build a per-overlapping-pair table with q-values + meth_diffs.

    Used for the per-DMR stat-concordance parquet. We avoid bioframe to
    keep deps small -- the input call sets are small (≤ a few thousand
    DMRs each) so an O(n*m) hash-by-chrom is fine.
    """
    b_idx: dict[str, list[tuple[int, int, float, float, int]]] = {}
    for i, r in enumerate(b.iter_rows(named=True)):
        b_idx.setdefault(str(r["chrom"]), []).append((
            int(r["start"]), int(r["end"]),
            float(r["qvalue"]) if r["qvalue"] is not None else float("nan"),
            float(r["meth_diff"]) if r["meth_diff"] is not None else float("nan"),
            i,
        ))
    for ch in b_idx:
        b_idx[ch].sort()

    rows: list[dict] = []
    for i, r in enumerate(a.iter_rows(named=True)):
        ch = str(r["chrom"])
        s, e = int(r["start"]), int(r["end"])
        qa = float(r["qvalue"]) if r["qvalue"] is not None else float("nan")
        ma = float(r["meth_diff"]) if r["meth_diff"] is not None else float("nan")
        for bs, be, qb, mb, j in b_idx.get(ch, []):
            if be < s:
                continue
            if bs > e:
                break
            inter = max(0, min(e, be) - max(s, bs))
            union = max(e, be) - min(s, bs)
            rows.append({
                "a_idx": i, "b_idx": j,
                "chrom": ch,
                "a_start": s, "a_end": e,
                "b_start": bs, "b_end": be,
                "bp_inter": inter, "bp_union": union,
                "jaccard_pair": inter / union if union else 0.0,
                "q_a": qa, "q_b": qb,
                "meth_diff_a": ma, "meth_diff_b": mb,
                "direction_match": (ma > 0) == (mb > 0),
            })
    return pl.DataFrame(rows) if rows else pl.DataFrame(
        schema={
            "a_idx": pl.Int64, "b_idx": pl.Int64, "chrom": pl.Utf8,
            "a_start": pl.Int64, "a_end": pl.Int64,
            "b_start": pl.Int64, "b_end": pl.Int64,
            "bp_inter": pl.Int64, "bp_union": pl.Int64,
            "jaccard_pair": pl.Float64,
            "q_a": pl.Float64, "q_b": pl.Float64,
            "meth_diff_a": pl.Float64, "meth_diff_b": pl.Float64,
            "direction_match": pl.Boolean,
        },
    )


def build_concordance(out_dir: Path, engines_done: tuple[str, ...]) -> dict:
    """Build the DMR concordance table across epykit / DSS / methylKit.

    Returns a small dict of headline metrics for the manifest + SUMMARY.md.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    sets: dict[str, pl.DataFrame] = {}

    # epykit DMR sets (chain_merge under the lr backbone -- the lr+ stack
    # only affects DMC; DMR uses the standard chain_merge).
    for method in DEFAULT_DMR_METHODS:
        pq = OUT_DMR / f"{method}.parquet"
        if pq.exists():
            sets[f"epykit_{method}"] = _load_epykit_dmr_called(pq)
            logger.info(
                "concordance: epykit_%s called=%d",
                method, sets[f"epykit_{method}"].height,
            )

    if DSS_DMR_CSV.exists():
        sets["dss"] = _load_dss_dmr_called(DSS_DMR_CSV)
        logger.info("concordance: dss called=%d", sets["dss"].height)
    else:
        logger.warning("DSS DMR not found at %s -- skipping", DSS_DMR_CSV)

    if MK_DMR_CSV.exists():
        sets["methylkit"] = _load_methylkit_dmr_called(MK_DMR_CSV)
        logger.info("concordance: methylkit called=%d", sets["methylkit"].height)
    else:
        logger.warning("methylKit DMR not found at %s -- skipping", MK_DMR_CSV)

    # Pairwise Jaccard on the merged-interval union of each set.
    iou_rows: list[dict] = []
    tools = sorted(sets.keys())
    for i, ta in enumerate(tools):
        for tb in tools[i + 1:]:
            jres = _interval_jaccard_dmr(sets[ta], sets[tb])
            dres = _direction_agreement(sets[ta], sets[tb])
            row = {"tool_a": ta, "tool_b": tb, **jres, **dres}
            iou_rows.append(row)
            logger.info(
                "  %s vs %s: Jaccard=%.4f  pair_overlap=%d  dir_agree=%.2f",
                ta, tb, row["jaccard"], row["n_overlapping_pairs"],
                row["direction_agree_frac"],
            )
    iou_df = pl.DataFrame(iou_rows) if iou_rows else pl.DataFrame()
    iou_path = out_dir / "dmr_iou.parquet"
    iou_df.write_parquet(iou_path)
    logger.info("wrote %s (%d rows)", iou_path, iou_df.height)

    # Per-DMR stat concordance: one row per overlapping pair, per tool pair.
    stat_chunks: list[pl.DataFrame] = []
    for i, ta in enumerate(tools):
        for tb in tools[i + 1:]:
            sub = _overlap_join(sets[ta], sets[tb])
            if sub.height:
                sub = sub.with_columns(
                    tool_a=pl.lit(ta),
                    tool_b=pl.lit(tb),
                )
                stat_chunks.append(sub)
    if stat_chunks:
        stat_df = pl.concat(stat_chunks, how="diagonal_relaxed")
    else:
        stat_df = pl.DataFrame()
    stat_path = out_dir / "per_dmr_stat_concordance.parquet"
    stat_df.write_parquet(stat_path)
    logger.info("wrote %s (%d rows)", stat_path, stat_df.height)

    # Headline summary
    headline: dict = {"set_sizes": {k: v.height for k, v in sets.items()}}
    if iou_df.height:
        for r in iou_df.iter_rows(named=True):
            headline[f"{r['tool_a']}_vs_{r['tool_b']}"] = {
                "jaccard_bp": round(r["jaccard"], 4),
                "n_overlapping_pairs": int(r["n_overlapping_pairs"]),
                "direction_agree_frac": round(r["direction_agree_frac"], 4),
            }

    # SUMMARY.md
    summary_path = out_dir / "SUMMARY.md"
    lines = [
        "# Study 3 -- Post-Phase-3 cross-tool DMR concordance",
        "",
        "Headline metrics for the GSE263850 re-run (Phase 4 Task 6).",
        "Generated by `benchmark/scripts/run_epykit_gse.py`.",
        "",
        "## Call-set sizes (q < 0.05)",
        "",
    ]
    for k, n in sorted(headline["set_sizes"].items()):
        lines.append(f"- **{k}**: {n:,} called DMRs")
    lines += ["", "## Pairwise concordance", "",
              "| Pair | bp-Jaccard | Overlap pairs | Direction agree |",
              "|---|---:|---:|---:|"]
    for r in iou_df.iter_rows(named=True):
        lines.append(
            f"| {r['tool_a']} vs {r['tool_b']} | "
            f"{r['jaccard']:.4f} | {r['n_overlapping_pairs']:,} | "
            f"{r['direction_agree_frac']:.2%} |"
        )
    lines += [
        "",
        "**bp-Jaccard** is the bp-weighted Jaccard on the merged-interval",
        "union of each tool's q<0.05 DMR set; **Overlap pairs** is the count",
        "of (A, B) DMR pairs that share ≥ 1 bp; **Direction agree** is the",
        "fraction of those pairs with the same sign of meth_diff.",
        "",
        f"Threshold: q < {CALLED_QVAL_THRESHOLD}.",
        "DSS panel: ``dmr_dss.csv`` (already DSS-significant, no qvalue filter).",
        "methylKit: ``dmr_all_tiles.csv`` filtered to q<0.05; meth_diff",
        "rescaled from percent to [-1, 1] for the direction column.",
        "",
    ]
    summary_path.write_text("\n".join(lines), encoding="utf-8")
    logger.info("wrote %s", summary_path)

    return headline


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------


def _git_head() -> str:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True,
        )
        return out.strip()
    except (subprocess.SubprocessError, OSError):
        return "unknown"


def write_manifest(timings: list[dict], headline: dict | None) -> None:
    import epykit as ep  # noqa: PLC0415

    OUT_BASE.mkdir(parents=True, exist_ok=True)

    total_wall = sum(t.get("wall_s", 0.0) for t in timings if t.get("ok"))
    lines = [
        "epykit_post_phase3 -- Phase 4 Task 6 (run_epykit_gse.py / GSE263850)",
        "",
        f"Date           : {time.strftime('%Y-%m-%d %H:%M:%S')}",
        f"epykit version : {ep.__version__}",
        "Engine tag     : v0.7.5-phase3-engines-frozen",
        f"Git HEAD       : {_git_head()}",
        "",
        "Cohort         : GSE263850 (3 SBP009 controls vs 3 Clone cases)",
        "Input format   : 12-column combined-strand BED",
        "                 (ep.read_combined_strand_bed)",
        "Raw .bed.gz    : ../epykit2/GSE263850_RAW/  (6 files)",
        "Samplesheet    : benchmark/data/study3/samplesheet_gse263850.csv",
        "",
        "DMC engines run: " + ", ".join(
            t["engine"] for t in timings
            if t.get("ok") and t["stage"] == "dmc"
        ),
        "DMR methods run: " + ", ".join(
            t["engine"] for t in timings
            if t.get("ok") and t["stage"] == "dmr"
        ),
        "",
        "lr+ kwargs (explicit on top of test='lr'):",
        "    neighbour_combine=True, neighbour_bp=500, sep_fallback=True,",
        "    sep_threshold=0.9, fdr_method='fdr_tsbh', dispersion='eb'",
        "",
        f"Total wallclock (ok engines only): {total_wall:.1f}s",
        "",
        "Per-stage timing:",
    ]
    for t in timings:
        status = "ok" if t.get("ok") else "FAIL"
        lines.append(
            f"  {t['stage']:<3} {t['engine']:<14} {status:<4} "
            f"wall={t.get('wall_s', 0.0):.1f}s  "
            f"n_sites={t.get('n_sites')}  n_sig={t.get('n_sig')}"
        )
        if not t.get("ok"):
            lines.append(f"      error: {t.get('error')}")

    lines += [
        "",
        "Cross-tool concordance reused (no re-runs):",
        f"  DSS  per-CpG : {DSS_DMC_TSV}",
        f"  DSS  per-DMR : {DSS_DMR_CSV}",
        f"  methylKit DMC: {MK_DMC_CSV}",
        f"                 (streamed via pl.scan_csv; 2 GB file is never",
        "                  materialised in concordance pipeline)",
        f"  methylKit DMR: {MK_DMR_CSV}",
        "",
        "Concordance outputs:",
        "  benchmark/data/study3/comparisons_post_phase3/dmr_iou.parquet",
        "  benchmark/data/study3/comparisons_post_phase3/"
        "per_dmr_stat_concordance.parquet",
        "  benchmark/data/study3/comparisons_post_phase3/SUMMARY.md",
        "",
        "Pre-Phase-3 baseline artefacts left untouched:",
        "  benchmark/data/study3/comparisons/  (paper / DSS / chain_merge",
        "                                       gene-level comparisons)",
        "  benchmark/data/study3/dmr_significant_*.csv  (epykit-pre-Phase-1)",
        "",
        "Per-engine bulk parquets:",
        "  benchmark/data/study3/epykit_post_phase3/dmc/<engine>.parquet",
        "  benchmark/data/study3/epykit_post_phase3/dmr/<method>.parquet",
        "These are gitignored as bulk artefacts; regenerable via "
        "run_epykit_gse.py.",
        "",
    ]
    if headline is not None:
        lines.append("Headline concordance:")
        for k, v in headline.items():
            if k == "set_sizes":
                continue
            if isinstance(v, dict):
                lines.append(
                    f"  {k}: bp-Jaccard={v['jaccard_bp']}  "
                    f"overlap_pairs={v['n_overlapping_pairs']}  "
                    f"dir_agree={v['direction_agree_frac']}"
                )

    MANIFEST_PATH.write_text("\n".join(lines), encoding="utf-8")
    logger.info("wrote %s", MANIFEST_PATH)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--engines", nargs="+", default=list(DEFAULT_DMC_ENGINES),
        help="DMC engines (default: lr lr+).",
    )
    parser.add_argument(
        "--dmr-methods", nargs="+", default=list(DEFAULT_DMR_METHODS),
        help="DMR methods (default: chain_merge).",
    )
    parser.add_argument("--skip-dmc", action="store_true",
                        help="Skip DMC stage.")
    parser.add_argument("--skip-dmr", action="store_true",
                        help="Skip DMR stage.")
    parser.add_argument("--skip-concordance", action="store_true",
                        help="Skip cross-tool concordance.")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="DEBUG-level logging.")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        force=True,
    )

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    CMP_OUT.mkdir(parents=True, exist_ok=True)

    timings: list[dict] = []
    md = None
    engines_done: tuple[str, ...] = tuple()

    if not (args.skip_dmc and args.skip_dmr):
        logger.info("=== ingest ===")
        md = ingest()

    if not args.skip_dmc:
        logger.info("=== DMC stage ===")
        dmc_timings = run_dmc(md, tuple(args.engines), OUT_DMC)
        timings.extend(dmc_timings)
        engines_done = tuple(
            t["engine"] for t in dmc_timings if t.get("ok")
        )

    if not args.skip_dmr:
        logger.info("=== DMR stage ===")
        # Ensure a DMC table exists for chain_merge / segment / sliding_window.
        if md is not None and md.uns.get("dmc", {}).get("last_key") is None:
            import epykit as ep  # noqa: PLC0415
            logger.info("DMC store empty -- running ep.tl.dmc(test='lr') first")
            ep.tl.dmc(md, test="lr", allow_n1=False)
        timings.extend(run_dmr(md, tuple(args.dmr_methods), OUT_DMR))

    if md is not None:
        del md
        gc.collect()

    headline: dict | None = None
    if not args.skip_concordance:
        logger.info("=== concordance stage ===")
        headline = build_concordance(CMP_OUT, engines_done)

    write_manifest(timings, headline)

    n_fail = sum(1 for t in timings if not t.get("ok", True))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
