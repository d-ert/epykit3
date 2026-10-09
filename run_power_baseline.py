"""Run a fresh intersection baseline or a separately labelled guarded union."""

import argparse
import json
from pathlib import Path

import polars as pl

import epykit as ep


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", default="samplesheet.csv")
    parser.add_argument("--treatment-group", default="Het_AKAP11_KO")
    parser.add_argument("--control-group", default="WT")
    parser.add_argument("--output-dir", type=Path, default=Path("results/baseline_intersection"))
    parser.add_argument("--input-store", type=Path, help="Read an existing filtered store without modifying it")
    parser.add_argument("--unite", choices=["intersect", "union"], default="intersect")
    parser.add_argument("--no-report", action="store_true")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        parser.error(f"Output already exists: {output}; choose a fresh directory")
    if args.input_store:
        source = args.input_store.resolve()
        obs = pl.read_csv(args.samplesheet)
        if set(obs["group"].to_list()) != {args.treatment_group, args.control_group}:
            parser.error("Samplesheet must contain exactly the requested case and control groups")
        if obs["sample_id"].n_unique() != len(obs):
            parser.error("Duplicate sample IDs")
        if not all((source / f"sample={s}").is_dir() for s in obs["sample_id"]):
            parser.error("Input store is missing a requested sample")
        obs = obs.with_columns((pl.col("group") == args.treatment_group).cast(pl.Int32).alias("treatment"))
        md = ep.MethylData(obs=obs, store=str(source), assembly="hg38", analysis_root=str(output / "methylstore"))
        md.uns["_store_history"] = [{"step": "filtered", "path": str(source)}]
    else:
        md = ep.read_bismark(samplesheet=args.samplesheet, treatment_group=args.treatment_group,
                             control_group=args.control_group, assembly="hg38",
                             store_dir=str(output / "methylstore"))
        ep.tl.qc(md)
        ep.pp.filter_coverage(md, lo_count=10, hi_perc=99.9)
    output.mkdir(parents=True, exist_ok=True)
    (output / "run.json").write_text(json.dumps({"status": "running", "unite": args.unite,
        "input_store": str(args.input_store.resolve()) if args.input_store else None,
        "min_samples_per_group": 2, "samplesheet": str(Path(args.samplesheet).resolve())}, indent=2))
    ep.pp.set_unite_type(md, type=args.unite)
    ep.tl.dmc(md, test="lr", min_samples_treatment=2, min_samples_control=2, materialize=False)
    ep.tl.dmr(md, method="chain_merge", preset="default")
    md.export_tables(str(output / "tables"), full=True)
    if not args.no_report:
        ep.generate_report(md, output=str(output / "report.html"), title="GSE263850 corrected baseline")
    metadata = {"status": "complete", "unite": args.unite, "dmc": md.uns["dmc"],
                "dmr": md.uns.get("dmr_params"), "input_store": md.store}
    (output / "run.json").write_text(json.dumps(metadata, indent=2, default=str))
    print(f"Saved corrected {args.unite} baseline: {output}")


if __name__ == "__main__":
    main()
