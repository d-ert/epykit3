"""run_tabix_vs_epykit.py -- focused head-to-head: methylKit tabix vs RAM vs epykit lr.

Answers one narrow question the phi-sweep leaves open: does epykit's out-of-core
engine actually beat methylKit's *on-disk tabix backend* on peak RAM (tabix's home
turf), and is it still faster -- at 1 core AND multi-core?

For ONE input cell (a dir of 2*n_per_group bismarkCoverage .cov.gz files) it runs,
each as a process-tree-monitored subprocess (peak RSS captured identically):
    - epykit_lr      (single-thread)
    - methylkit_ram  (in-memory backend)   at each --cores value
    - methylkit_tabix(on-disk tabix backend) at each --cores value

tabix changes only storage/IO, not statistics, so we compare wall_s + peak RSS
(rss and uss) + n_sites (n_sites must match RAM, a correctness check). Rows are
APPENDED to --out-csv (header written once); re-runnable across cells via a loop.

Usage:
    uv run python run_tabix_vs_epykit.py \
        --in-dir <cell>/bismark_cov --glob '*.cov.gz' --n-per-group 3 \
        --mincov 1 --assembly sim --dataset simulator --cell-id phi=0.0/seed=2026000 \
        --cores 1 8 --scratch /tmp/tve_scratch \
        --out-csv .../summaries/tabix_vs_epykit_per_run.csv
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from _resource_monitor import run_subprocess_monitored  # noqa: E402

ROW_FIELDS = [
    "dataset", "cell_id", "tool", "backend", "cores",
    "wall_s", "cpu_s", "rss_peak_mb", "rss_mean_mb", "uss_peak_mb",
    "cpu_percent_peak", "num_processes_peak", "n_sites", "returncode",
]


def _n_sites_tsv(path: Path) -> int:
    """Count data rows in a methylKit per-CpG TSV (minus header)."""
    if not path.exists():
        return -1
    with path.open() as fh:
        return max(sum(1 for _ in fh) - 1, 0)


def _n_sites_parquet(path: Path) -> int:
    if not path.exists():
        return -1
    import polars as pl
    return pl.read_parquet(path, columns=["pos"]).height


def _append_row(out_csv: Path, row: dict) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    new = not out_csv.exists()
    with out_csv.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ROW_FIELDS)
        if new:
            w.writeheader()
        w.writerow({k: row.get(k) for k in ROW_FIELDS})


def _record(out_csv, dataset, cell_id, tool, backend, cores, res, n_sites):
    row = dict(
        dataset=dataset, cell_id=cell_id, tool=tool, backend=backend, cores=cores,
        wall_s=round(res.get("wall_s") or 0, 2),
        cpu_s=round(res.get("cpu_s") or 0, 2),
        rss_peak_mb=res.get("rss_peak_mb"),
        rss_mean_mb=res.get("rss_mean_mb"),
        uss_peak_mb=res.get("uss_peak_mb"),
        cpu_percent_peak=res.get("cpu_percent_peak"),
        num_processes_peak=res.get("num_processes_peak"),
        n_sites=n_sites,
        returncode=res.get("returncode"),
    )
    _append_row(out_csv, row)
    flag = "OK" if res.get("returncode") == 0 else f"FAIL rc={res.get('returncode')}"
    print(f"  [{flag}] {tool:16s} cores={cores}  wall={row['wall_s']:>8}s  "
          f"rss_peak={row['rss_peak_mb']}MB  uss_peak={row['uss_peak_mb']}MB  "
          f"n_sites={n_sites}", flush=True)
    return row


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-dir", required=True, type=Path)
    ap.add_argument("--glob", default="*.cov.gz")
    ap.add_argument("--n-per-group", type=int, default=3)
    ap.add_argument("--mincov", type=int, default=1)
    ap.add_argument("--assembly", default="sim")
    ap.add_argument("--dataset", required=True, help="row label, e.g. simulator | gse263850")
    ap.add_argument("--cell-id", required=True, help="row label, e.g. phi=0.0/seed=2026000 | chr1")
    ap.add_argument("--cores", type=int, nargs="+", default=[1, 8])
    ap.add_argument("--scratch", required=True, type=Path, help="scratch dir (per-run subdirs wiped)")
    ap.add_argument("--out-csv", required=True, type=Path)
    ap.add_argument("--interval", type=float, default=0.2)
    ap.add_argument("--skip-ram", action="store_true", help="skip in-memory methylKit (e.g. OOM risk)")
    ap.add_argument("--skip-tabix", action="store_true", help="skip tabix methylKit")
    ap.add_argument("--skip-epykit", action="store_true")
    args = ap.parse_args(argv)

    in_dir = args.in_dir
    cov = sorted(in_dir.glob(args.glob))
    if len(cov) != 2 * args.n_per_group:
        raise SystemExit(f"expected {2*args.n_per_group} files matching {args.glob} in "
                         f"{in_dir}, found {len(cov)}")
    args.scratch.mkdir(parents=True, exist_ok=True)
    tag = f"{args.dataset}__{args.cell_id.replace('/', '_')}"
    print(f"=== {args.dataset} / {args.cell_id}  ({len(cov)} samples, mincov={args.mincov}) ===",
          flush=True)

    # ---- epykit lr (single-thread, out-of-core Parquet store) -------------
    if not args.skip_epykit:
        ep_out = args.scratch / f"{tag}__epykit_lr.parquet"
        ep_store = args.scratch / f"{tag}__epykit_store"
        cmd = [sys.executable, str(_HERE / "run_epykit_cell.py"),
               "--in-dir", str(in_dir), "--engine", "lr", "--out", str(ep_out),
               "--store-dir", str(ep_store), "--glob", args.glob,
               "--n-per-group", str(args.n_per_group), "--assembly", args.assembly]
        res = run_subprocess_monitored(cmd, interval=args.interval,
                                       stdout=None, stderr=None)
        _record(args.out_csv, args.dataset, args.cell_id, "epykit_lr", "parquet_ooc",
                1, res, _n_sites_parquet(ep_out))
        shutil.rmtree(ep_store, ignore_errors=True)

    # ---- methylKit, both backends, at each core count --------------------
    # Same unified runner for both -> identical code path except --dbtype.
    def _mk_cmd(dbtype, out, cores, dbdir=None):
        cmd = ["Rscript", str(_HERE / "run_methylkit_backend.R"),
               "--dbtype", dbtype, "--in-dir", str(in_dir), "--out", str(out),
               "--cores", str(cores), "--mincov", str(args.mincov),
               "--n-per-group", str(args.n_per_group), "--assembly", args.assembly]
        if dbdir is not None:
            cmd += ["--dbdir", str(dbdir)]
        return cmd

    for cores in args.cores:
        # in-memory (RAM) backend
        if not args.skip_ram:
            mk_out = args.scratch / f"{tag}__mk_ram_c{cores}.tsv"
            res = run_subprocess_monitored(_mk_cmd("memory", mk_out, cores),
                                           interval=args.interval, stdout=None, stderr=None)
            _record(args.out_csv, args.dataset, args.cell_id, "methylkit_ram",
                    "ram", cores, res, _n_sites_tsv(mk_out))

        # tabix (on-disk) backend
        if not args.skip_tabix:
            tb_out = args.scratch / f"{tag}__mk_tabix_c{cores}.tsv"
            tb_db = args.scratch / f"{tag}__tabixdb_c{cores}"
            res = run_subprocess_monitored(_mk_cmd("tabix", tb_out, cores, dbdir=tb_db),
                                           interval=args.interval, stdout=None, stderr=None)
            _record(args.out_csv, args.dataset, args.cell_id, "methylkit_tabix",
                    "tabix_ooc", cores, res, _n_sites_tsv(tb_out))
            shutil.rmtree(tb_db, ignore_errors=True)

    print(f"=== done {args.dataset}/{args.cell_id} -> {args.out_csv} ===", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
