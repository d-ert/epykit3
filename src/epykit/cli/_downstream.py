"""Downstream commands: ``annotate``, ``qc-report``, ``smooth``, ``report``,
``aggregate-regions`` and ``export``, registered in that order."""

from __future__ import annotations

import argparse
from pathlib import Path

from ._common import _auto_tsv_path, _cli_tsv_opts, _write_table_local


def _cmd_annotate(args: argparse.Namespace):
    """Handler for 'annotate' subcommand."""
    import polars as pl

    sites = pl.read_parquet(args.input)

    if args.gtf:
        from ..annotate import annotate_features

        sites = annotate_features(
            sites,
            args.gtf,
            source="gtf",
            promoter_upstream_bp=args.promoter_upstream_bp,
            promoter_downstream_bp=args.promoter_downstream_bp,
        )
        print("Gene feature annotation complete.")

    if args.cpg_islands:
        from ..annotate import annotate_cpg_islands

        sites = annotate_cpg_islands(sites, cpg_island_bed=args.cpg_islands)
        print("CpG island annotation complete.")

    sites.write_parquet(args.output)
    print(f"Annotated results written to {args.output}")

    tsv_suppressed, tsv_path_opt, _, _ = _cli_tsv_opts(args)
    if not tsv_suppressed:
        tsv_path = tsv_path_opt or _auto_tsv_path(args.output)
        _write_table_local(sites, tsv_path)
        print(f"Annotated TSV: {tsv_path}")


def _cmd_qc_report(args: argparse.Namespace):
    """Handler for 'qc-report' subcommand."""
    import polars as pl

    from ..qc import coverage_uniformity, global_methylation_report

    samples = args.samples.split(",")
    tsv_suppressed = _cli_tsv_opts(args)[0]

    print("=== Global methylation report ===")
    meth_report = global_methylation_report(args.methylstore, samples)
    print(meth_report)
    if args.output_dir:
        out = Path(args.output_dir)
        out.mkdir(parents=True, exist_ok=True)
        meth_report.write_parquet(str(out / "global_methylation.parquet"))
        if not tsv_suppressed:
            _write_table_local(meth_report, str(out / "global_methylation.tsv"))

    print("\n=== Coverage uniformity report ===")
    cov_frames: list[pl.DataFrame] = []
    for sample in samples:
        try:
            cov_df = coverage_uniformity(args.methylstore, sample)
            print(f"\n{sample}")
            print(cov_df)
            cov_frames.append(cov_df)
        except ValueError as exc:
            print(f"  Warning: {exc}")

    if cov_frames and args.output_dir:
        combined = pl.concat(cov_frames)
        combined.write_parquet(str(Path(args.output_dir) / "coverage_uniformity.parquet"))
        if not tsv_suppressed:
            _write_table_local(
                combined,
                str(Path(args.output_dir) / "coverage_uniformity.tsv"),
            )
        print(f"\nQC reports written to {args.output_dir}")


def _cmd_smooth(args: argparse.Namespace):
    """Handler for 'smooth' subcommand (Gaussian-kernel smoothing)."""
    from ..dmr import smooth_methylation_gaussian

    samples = args.samples.split(",")
    smooth_path = args.output
    smooth_methylation_gaussian(
        args.methylstore,
        samples,
        bandwidth=args.bandwidth,
        output_path=smooth_path,
    )
    print(f"Smoothed betas written to {smooth_path}")


def _cmd_report(args: argparse.Namespace):
    """Handler for 'report' subcommand."""
    from ..methyldata import MethylData

    md = MethylData.load(args.md)
    kwargs: dict = {
        "alpha": args.alpha,
        "min_abs_diff": args.min_abs_diff,
        "self_contained": args.self_contained,
    }
    if args.title:
        kwargs["title"] = args.title
    if args.gtf:
        kwargs["gtf_path"] = args.gtf
    md.report(args.output, **kwargs)
    print(f"Report written: {args.output}")


def _cmd_aggregate_regions(args: argparse.Namespace):
    """Handler for 'aggregate-regions' subcommand."""
    from .. import pp as pp_mod
    from ..methyldata import MethylData

    md = MethylData.load(args.md)
    pp_mod.aggregate_regions(
        md,
        args.bed,
        region_id_col=args.region_id_col,
        output_store=args.output_store,
        min_cpgs_per_region=args.min_cpgs_per_region,
    )
    md.save(args.md)
    print(f"Region aggregation complete; methylstore = {md.store}")


def _cmd_export(args: argparse.Namespace):
    """Handler for 'export' subcommand."""
    from ..export import dmcs_to_bed, dmrs_to_bed, to_bedgraph, to_bigwig
    from ..methyldata import MethylData

    md = MethylData.load(args.md)
    fmt = args.export_cmd
    if fmt == "bedgraph":
        to_bedgraph(md, args.sample, args.output, value=args.value)
    elif fmt == "bigwig":
        to_bigwig(md, args.sample, args.output, value=args.value)
    elif fmt == "dmcs-bed":
        dmcs_to_bed(
            md, args.output, alpha=args.alpha, min_abs_diff=args.min_abs_diff, test=args.test
        )
    elif fmt == "dmrs-bed":
        dmrs_to_bed(md, args.output)
    else:
        raise SystemExit(f"unknown export format: {fmt}")
    print(f"Wrote {args.output}")


def register(sub: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    """Add the downstream subcommands to ``sub`` in help order."""
    # annotate
    p_ann = sub.add_parser("annotate", help="Annotate DMC/DMR results with genomic features")
    p_ann.add_argument(
        "--input", required=True, help="Parquet file from 'epykit dmc' or 'epykit dmr'"
    )
    p_ann.add_argument("--output", required=True)
    p_ann.add_argument(
        "--gtf",
        help="Ensembl/UCSC GTF/GFF3 for gene feature annotation",
    )
    p_ann.add_argument(
        "--cpg-islands",
        help="UCSC CpGIsland BED file for CpG context annotation",
    )
    p_ann.add_argument("--promoter-upstream-bp", type=int, default=2000)
    p_ann.add_argument("--promoter-downstream-bp", type=int, default=200)
    p_ann.add_argument(
        "--no-tsv",
        action="store_true",
        dest="no_tsv",
        default=False,
        help="Suppress the sibling .tsv auto-emit.",
    )
    p_ann.add_argument(
        "--tsv",
        dest="tsv_path",
        default=None,
        help="Override sibling table path. .csv suffix -> comma delim.",
    )
    # Deprecated csv* aliases (see dmc) -- honoured but warn; hidden from --help.
    p_ann.add_argument(
        "--no-csv",
        action="store_true",
        dest="no_csv",
        default=False,
        help=argparse.SUPPRESS,
    )
    p_ann.add_argument("--csv", dest="csv_path", default=None, help=argparse.SUPPRESS)
    p_ann.set_defaults(func=_cmd_annotate)

    # qc-report
    p_qc = sub.add_parser("qc-report", help="QC and coverage uniformity report")
    p_qc.add_argument("--methylstore", required=True)
    p_qc.add_argument(
        "--samples",
        required=True,
        help="Comma-separated list of sample IDs",
    )
    p_qc.add_argument(
        "--output-dir",
        help="Directory for Parquet QC output files (optional)",
    )
    p_qc.add_argument(
        "--no-tsv",
        action="store_true",
        dest="no_tsv",
        default=False,
        help="Suppress the sibling .tsv auto-emit alongside the parquets.",
    )
    # Deprecated csv alias (see dmc) -- honoured but warn; hidden from --help.
    p_qc.add_argument(
        "--no-csv",
        action="store_true",
        dest="no_csv",
        default=False,
        help=argparse.SUPPRESS,
    )
    p_qc.set_defaults(func=_cmd_qc_report)

    # smooth
    p_sm = sub.add_parser(
        "smooth",
        help="Gaussian-kernel methylation beta smoothing (approximates BSmooth)",
        description=(
            "Gaussian-kernel methylation beta smoothing. Approximates BSmooth "
            "via scipy.ndimage.gaussian_filter1d on a regular grid -- not a "
            "true local LOESS. ~500x faster than statsmodels LOESS."
        ),
    )
    p_sm.add_argument("--methylstore", required=True)
    p_sm.add_argument(
        "--samples",
        required=True,
        help="Comma-separated list of sample IDs",
    )
    p_sm.add_argument("--output", required=True, help="Output directory for smoothed beta chunks")
    p_sm.add_argument(
        "--bandwidth", type=int, default=1000, help="Smoothing bandwidth in bp (default 1000)"
    )
    p_sm.set_defaults(func=_cmd_smooth)

    # report
    p_rep = sub.add_parser(
        "report",
        help="Render an interactive HTML report from a saved MethylData",
    )
    p_rep.add_argument(
        "--md", required=True, help="Path to a directory previously written with md.save(...)"
    )
    p_rep.add_argument("--output", required=True, help="Output HTML file")
    p_rep.add_argument("--title", default=None)
    p_rep.add_argument("--gtf", default=None, help="Optional GTF for a TSS metaplot section")
    p_rep.add_argument("--alpha", type=float, default=0.05)
    p_rep.add_argument("--min-abs-diff", dest="min_abs_diff", type=float, default=0.1)
    p_rep.add_argument(
        "--self-contained",
        dest="self_contained",
        action="store_true",
        default=True,
        help="Embed Plotly inline so the HTML works offline (default)",
    )
    p_rep.add_argument(
        "--no-self-contained",
        dest="self_contained",
        action="store_false",
        help="Load Plotly from a CDN (smaller file, needs internet)",
    )
    p_rep.set_defaults(func=_cmd_report)

    # aggregate-regions
    p_agg = sub.add_parser(
        "aggregate-regions",
        help="Aggregate CpG counts to user-supplied BED regions",
    )
    p_agg.add_argument(
        "--md", required=True, help="Path to a directory previously written with md.save(...)"
    )
    p_agg.add_argument("--bed", required=True, help="BED file of regions to aggregate to")
    p_agg.add_argument("--output-store", default=None, help="Override output Parquet store path")
    p_agg.add_argument(
        "--region-id-col",
        default=None,
        help="BED column name to use as region_id (default: 'name' or chrom:start-end)",
    )
    p_agg.add_argument("--min-cpgs-per-region", type=int, default=1)
    p_agg.set_defaults(func=_cmd_aggregate_regions)

    # export
    p_exp = sub.add_parser("export", help="Export to BedGraph / BigWig / BED")
    exp_sub = p_exp.add_subparsers(dest="export_cmd", required=True)
    for name, help_text in (
        ("bedgraph", "Per-sample beta / coverage -> BedGraph"),
        ("bigwig", "Per-sample beta / coverage -> BigWig (pyBigWig)"),
    ):
        sp = exp_sub.add_parser(name, help=help_text)
        sp.add_argument("--md", required=True)
        sp.add_argument("--sample", required=True)
        sp.add_argument("--output", required=True)
        sp.add_argument("--value", default="beta", choices=["beta", "coverage", "N_meth"])
        sp.set_defaults(func=_cmd_export)
    sp_d = exp_sub.add_parser("dmcs-bed", help="DMCs -> 6-column BED")
    sp_d.add_argument("--md", required=True)
    sp_d.add_argument("--output", required=True)
    sp_d.add_argument("--alpha", type=float, default=0.05)
    sp_d.add_argument("--min-abs-diff", dest="min_abs_diff", type=float, default=0.0)
    sp_d.add_argument("--test", default=None)
    sp_d.set_defaults(func=_cmd_export)
    sp_r = exp_sub.add_parser("dmrs-bed", help="DMRs -> 6-column BED")
    sp_r.add_argument("--md", required=True)
    sp_r.add_argument("--output", required=True)
    sp_r.set_defaults(func=_cmd_export)
