"""Ingestion commands: ``convert``, ``filter`` and ``summary``."""

from __future__ import annotations

import argparse

from .. import filter
from ..convert import convert_sample


def _cmd_convert(args: argparse.Namespace):
    """Handler for 'convert' subcommand."""
    # --merge-cpg/--no-merge-cpg is a tri-state flag with a None sentinel
    # (neither flag given). Resolve None -> True so a bare ``epykit convert``
    # merges CpG dyads, matching the API default (convert_sample/read_bismark
    # default merge_strands=True). Explicit flags force True/False (D10).
    merge_strands = True if args.merge_cpg is None else args.merge_cpg
    convert_sample(
        args.input,
        args.sample_id,
        args.output_dir,
        context=args.context,
        reference_fasta=args.reference_fasta,
        merge_strands=merge_strands,  # CLI flag is --merge-cpg; param is merge_strands
        format=args.format,
        canonical_only=getattr(args, "canonical_only", False),
    )


def _cmd_filter(args: argparse.Namespace):
    """Handler for 'filter' subcommand."""
    filter.filter_sites(
        args.methylstore,
        args.output_dir,
        min_coverage=args.min_coverage,
        max_coverage_quantile=args.max_coverage_quantile,
        blacklist_bed=args.blacklist_bed,
        sample=args.sample,
    )


def _cmd_sample_summary(args: argparse.Namespace):
    """Handler for 'summary' subcommand."""
    df = filter.sample_summary(args.methylstore, args.sample, output_path=args.output)
    if not args.output:
        print(df)


def register(sub: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    """Add ``convert``, ``filter`` and ``summary`` to ``sub``, in that order."""
    # convert
    p_conv = sub.add_parser(
        "convert",
        help="Convert a Bismark .cov or MethylDackel .bedGraph file to Parquet",
    )
    p_conv.add_argument("--input", required=True)
    p_conv.add_argument("--sample-id", required=True)
    p_conv.add_argument("--output-dir", required=True)
    p_conv.add_argument(
        "--format",
        choices=["bismark", "methyldackel"],
        default="bismark",
        help=(
            "Source file format. 'bismark' (default) for .cov[.gz] files "
            "produced by bismark_methylation_extractor / bismark2bedGraph. "
            "'methyldackel' for .bedGraph[.gz] files produced by "
            "MethylDackel extract -- same 6-column layout, with the leading "
            "track header skipped automatically."
        ),
    )
    p_conv.add_argument(
        "--context",
        choices=["CpG", "CHG", "CHH"],
        default="CpG",
    )
    p_conv.add_argument("--reference-fasta")
    p_conv.add_argument(
        "--merge-cpg",
        dest="merge_cpg",
        action="store_true",
        help=(
            "Merge symmetric CpG dyads (+/- strand) into one record. This is "
            "the default (matching the Python API merge_strands=True); the flag "
            "is accepted for explicitness."
        ),
    )
    p_conv.add_argument(
        "--no-merge-cpg",
        dest="merge_cpg",
        action="store_false",
        help="Disable CpG-dyad merging; keep per-strand records.",
    )
    p_conv.add_argument(
        "--canonical-only",
        dest="canonical_only",
        action="store_true",
        default=False,
        help=(
            "Keep only the fixed human-style chromosome set (1-22, X, Y, M/MT, "
            "with or without a chr prefix) and drop every other contig before "
            "the partition write. The setting is part of the per-sample "
            "conversion cache. Default: keep every contig."
        ),
    )
    p_conv.set_defaults(merge_cpg=None, func=_cmd_convert)

    # filter
    p_filt = sub.add_parser("filter", help="Filter low-coverage CpGs")
    p_filt.add_argument("--methylstore", required=True)
    p_filt.add_argument("--output-dir", required=True)
    p_filt.add_argument("--min-coverage", type=int, default=10)
    p_filt.add_argument("--max-coverage-quantile", type=float, default=0.999)
    p_filt.add_argument("--blacklist-bed")
    p_filt.add_argument("--sample")
    p_filt.set_defaults(func=_cmd_filter)

    # summary
    p_sum = sub.add_parser("summary", help="Per-sample summary statistics")
    p_sum.add_argument("--methylstore", required=True)
    p_sum.add_argument("--sample", required=True)
    p_sum.add_argument("--output")
    p_sum.set_defaults(func=_cmd_sample_summary)
