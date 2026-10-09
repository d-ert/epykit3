"""Run the experimental spatial GLS backend on a read-only count store."""

import argparse
import json
import logging
import time
from pathlib import Path

import polars as pl

import epykit as ep


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samplesheet", default="samplesheet.csv")
    p.add_argument("--treatment-group", default="Het_AKAP11_KO")
    p.add_argument("--control-group", default="WT")
    p.add_argument("--input-store", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--chromosomes", nargs="+")
    p.add_argument("--n-perm", type=int, default=100)
    p.add_argument("--window-bp", type=int, default=500)
    p.add_argument("--perm-seed", type=int, default=42)
    p.add_argument("--workers", type=int, default=None)
    args = p.parse_args()
    out = args.output_dir.resolve()
    if out.exists():
        p.error(f"Output already exists: {out}")
    obs = pl.read_csv(args.samplesheet)
    if set(obs["group"].to_list()) != {args.treatment_group, args.control_group} or obs[
        "sample_id"
    ].n_unique() != len(obs):
        p.error("Samplesheet must contain unique samples in the two requested groups")
    obs = obs.with_columns(
        (pl.col("group") == args.treatment_group).cast(pl.Int32).alias("treatment")
    )
    md = ep.MethylData(
        obs=obs, store=str(args.input_store.resolve()), assembly="hg38", analysis_root=str(out)
    )
    started = time.monotonic()
    logging.basicConfig(level=logging.INFO)
    ep.tl.dmr(
        md,
        method="regional_gls",
        chromosomes=args.chromosomes,
        n_perm=args.n_perm,
        perm_seed=args.perm_seed,
        window_bp=args.window_bp,
        n_workers=args.workers,
        tsv=False,
    )
    out.mkdir(parents=True, exist_ok=True)
    # Preserve headers for a null result; generic exporter currently rejects
    # empty DMR tables, while this research runner must record empty outcomes.
    md.uns["dmr"].write_csv(out / "dmr.tsv", separator="\t")
    md.uns["dmr_candidates"].write_csv(out / "dmr_candidates.tsv", separator="\t")
    meta = dict(
        status="complete",
        experimental=True,
        parameters=md.uns["dmr_params"],
        seconds=time.monotonic() - started,
        input_store=md.store,
        samplesheet=str(Path(args.samplesheet).resolve()),
        discoveries=len(md.uns["dmr"]),
    )
    (out / "run.json").write_text(json.dumps(meta, indent=2))
    print(f"Saved {len(md.uns['dmr'])} experimental GLS regions to {out}", flush=True)


if __name__ == "__main__":
    main()
