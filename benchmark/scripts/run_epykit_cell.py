"""run_epykit_cell.py -- run one epykit DMC engine on one simulator cell.

Standalone, single-cell, single-engine runner so the phi-sweep driver can
invoke epykit as a *monitored subprocess* (peak-RSS captured identically to the
R baselines via _resource_monitor). It reads the SAME ``bismark_cov/*.cov.gz``
files the R tools read (samples 1-3 = treat, 4-6 = ctrl, matching
run_methylkit_simulator.R / run_dss_simulator.R / write_samplesheet), so inputs
are byte-identical across tools.

Writes a per-CpG parquet: chrom, pos, pvalue, qvalue, meth_diff (+ the
*_combined columns when engine=lr+). The phi-sweep driver scores it with the
canonical score_dmc_parquet.

Usage:
    python run_epykit_cell.py --in-dir <cell>/bismark_cov --engine lr \\
        --out <cell>/epykit_lr.parquet --store-dir <scratch>/store
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import polars as pl

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from _epykit_scoring import _dmc_kwargs  # noqa: E402


def _build_samplesheet(cov_files: list[Path], n_per_group: int, sheet: Path) -> Path:
    """First ``n_per_group`` files -> treat, next -> ctrl (matches the R runners)."""
    rows = []
    for i, p in enumerate(cov_files):
        grp = "treat" if i < n_per_group else "ctrl"
        rows.append((f"{grp}_{(i % n_per_group) + 1}", grp, str(p)))
    sheet.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(rows, schema=["sample_id", "group", "path"], orient="row").write_csv(sheet)
    return sheet


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-dir", required=True, type=Path,
                    help="Directory of the 6 .cov.gz files (bismark_cov)")
    ap.add_argument("--out", required=True, type=Path, help="Output per-CpG parquet")
    ap.add_argument("--engine", required=True, help="lr | lr+ | welch_t | fisher")
    ap.add_argument("--assembly", default="sim")
    ap.add_argument("--store-dir", required=True, type=Path,
                    help="Scratch methylstore dir (wiped + recreated)")
    ap.add_argument("--n-per-group", type=int, default=3)
    ap.add_argument("--glob", default="*.cov.gz")
    args = ap.parse_args(argv)

    import epykit as ep

    cov_files = sorted(args.in_dir.glob(args.glob))
    if len(cov_files) != 2 * args.n_per_group:
        raise SystemExit(
            f"expected {2 * args.n_per_group} .cov.gz in {args.in_dir}, "
            f"found {len(cov_files)} (glob={args.glob!r})"
        )

    if args.store_dir.exists():
        shutil.rmtree(args.store_dir)
    args.store_dir.mkdir(parents=True, exist_ok=True)
    sheet = _build_samplesheet(cov_files, args.n_per_group, args.store_dir / "samplesheet.csv")

    md = ep.read_bismark(
        str(sheet), treatment_group="treat", control_group="ctrl",
        assembly=args.assembly, store_dir=str(args.store_dir),
    )
    ep.pp.unite(md, type="intersect")

    backend_test, kwargs = _dmc_kwargs(args.engine, allow_n1=False)
    ep.tl.dmc(md, **kwargs)
    df = md.get_dmc(test=backend_test)

    keep = [c for c in ("chrom", "pos", "pvalue", "qvalue", "meth_diff",
                        "pvalue_combined", "qvalue_combined") if c in df.columns]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.select(keep).write_parquet(args.out)
    print(f"wrote {args.out} ({df.height} CpGs, engine={args.engine}, cols={keep})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
