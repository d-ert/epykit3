"""The ``dmr`` command: region calling from DMC results or read-pooled tiles."""

from __future__ import annotations

import argparse

from .._dmc_engines import PUBLIC_ENGINES
from ._common import (
    _add_min_samples_args,
    _auto_tsv_path,
    _cli_n1_and_footgun_checks,
    _cli_tsv_opts,
    _read_samplesheet_groups,
)


def _cmd_dmr(args: argparse.Namespace):
    """Handler for 'dmr' subcommand."""
    import polars as pl

    from ..dmr import (
        _DMR_DEFAULT_MIN_CPGS,
        _resolve_tile_chromosomes,
        apply_region_qfilter,
        call_dmr_chain_merge,
        call_dmr_sliding_window,
        call_dmr_tile_based,
        resolve_layer_min_cpgs,
    )

    # Mirror of tl.dmr: canonical_only is tile-only. The DMC-derived methods
    # inherit the chromosome universe of the DMC parquet, so the filter
    # belongs upstream; refuse rather than filter a finished q-value table.
    canonical_only = bool(getattr(args, "canonical_only", False))
    if canonical_only and args.method != "tile":
        raise SystemExit(
            f"error: --canonical-only applies to --method tile only; got "
            f"--method {args.method}. The DMC-derived methods use the chromosome "
            f"universe of the DMC parquet, so filter there: run epykit dmc "
            f"--canonical-only (or epykit convert --canonical-only), then call "
            f"epykit dmr without --canonical-only."
        )

    # M2 gate (mirror of tl.dmr): empirical_fdr is wired only for method=tile.
    # Pre-fix the CLI accepted --empirical-fdr against any method and silently
    # dropped it on chain_merge/sliding_window/segment, leaving users
    # thresholding an uncalibrated combined_qvalue as if it were FDR-controlled.
    if getattr(args, "empirical_fdr", False) and args.method != "tile":
        raise NotImplementedError(
            f"--empirical-fdr is currently implemented only for "
            f"--method=tile. Got --method={args.method!r}. Use --method=tile "
            f"or omit --empirical-fdr. (Follow-up: implement permutation FDR "
            f"for chain_merge/sliding_window/segment -- tracked in "
            f"docs/superpowers/plans/2026-06-07-epykit-audit-fixes.md "
            f"Batch-4 follow-up.)"
        )

    if args.method == "chain_merge":
        # --- DSS-style chain-merge path: takes a DMC parquet ---
        if not args.dmc_results:
            raise ValueError("method=chain_merge requires --dmc-results.")

        dmc_results = pl.read_parquet(args.dmc_results)
        # Resolve --min-cpgs at the layer (shared with tl.dmr) so bare CLI
        # chain_merge matches bare API: explicit --min-cpgs N wins; else the
        # active --preset's value; else 5 (the paper default). Passing the
        # concrete value overrides the engine's own preset resolution
        # deterministically. The resolver also validates --preset.
        cm_min_cpgs = resolve_layer_min_cpgs(args.min_cpgs, args.preset)
        dmr_results = call_dmr_chain_merge(
            dmc_results,
            preset=args.preset,
            alpha=args.alpha,
            min_abs_meth_diff=args.min_abs_meth_diff,
            dis_merge_bp=args.dis_merge_bp,
            min_cpgs=cm_min_cpgs,
            pct_sig=args.pct_sig,
            minlen_bp=args.minlen_bp,
            use_q_for_sig=args.use_q_for_sig,
        )
        # Region-level q-value post-filter, mirroring tl.dmr's chain_merge
        # branch exactly via the shared apply_region_qfilter helper. Pre-fix
        # the CLI skipped this, so CLI chain_merge output was less filtered
        # than the API (D11).
        dmr_results = apply_region_qfilter(dmr_results, getattr(args, "min_mean_qvalue", None))
    elif args.method == "tile":
        # --- tile-based path. Needs methylstore + samplesheet. ---
        if not args.methylstore or not args.samplesheet:
            raise ValueError(
                "method=tile requires --methylstore, --samplesheet, "
                "--treatment-group and --control-group."
            )

        treatment_samples, control_samples = _read_samplesheet_groups(
            args.samplesheet, args.treatment_group, args.control_group
        )
        args._samples = (treatment_samples, control_samples)
        _cli_n1_and_footgun_checks(args, unit="tiles")

        print(f"Treatment samples: {treatment_samples}")
        print(f"Control samples:   {control_samples}")
        print(f"Tile size:         {args.tile_size_bp} bp")
        print(f"Test:              {args.test}")
        if canonical_only:
            print("Chromosomes:       canonical set only (--canonical-only)")

        # Resolve the chromosome universe once, as tl.dmr does, and hand the
        # same explicit list to the observed run and every permutation.
        tile_chromosomes = _resolve_tile_chromosomes(
            args.methylstore, None, canonical_only=canonical_only
        )
        dmr_results = call_dmr_tile_based(
            methylstore_path=args.methylstore,
            samples_treatment=treatment_samples,
            samples_control=control_samples,
            tile_size_bp=args.tile_size_bp,
            test=args.test,
            chromosomes=tile_chromosomes,
            min_cpgs_per_tile=args.min_cpgs_per_tile,
            alpha=args.alpha,
            min_abs_meth_diff=args.min_abs_meth_diff,
            unite=args.unite,
            min_samples_treatment=args.min_samples_treatment,
            min_samples_control=args.min_samples_control,
        )
        # Optional stricter-than-alpha q-value post-filter on the per-tile
        # `qvalue` column, mirroring tl.dmr's tile branch via the shared
        # helper. Applied BEFORE empirical FDR so the permutation BH
        # correction sees the same region set as tl.dmr (which also filters
        # before empirical_fdr_for_dmr). Pre-fix the CLI tile branch skipped
        # this entirely, so `epykit dmr --method tile` diverged from
        # tl.dmr(method='tile') when --min-mean-qvalue != alpha (D11 gap).
        dmr_results = apply_region_qfilter(
            dmr_results,
            getattr(args, "min_mean_qvalue", None),
            candidate_cols=("qvalue",),
        )
        if getattr(args, "empirical_fdr", False) and len(dmr_results) > 0:
            from ..dmr import empirical_fdr_for_dmr

            dmr_results = empirical_fdr_for_dmr(
                methylstore_path=args.methylstore,
                samples_treatment=treatment_samples,
                samples_control=control_samples,
                observed_dmr=dmr_results,
                n_perm=args.n_perm,
                seed=args.perm_seed,
                tile_size_bp=args.tile_size_bp,
                test=args.test,
                chromosomes=tile_chromosomes,
                min_cpgs_per_tile=args.min_cpgs_per_tile,
                alpha=args.alpha,
                min_abs_meth_diff=args.min_abs_meth_diff,
                unite=args.unite,
                min_samples_treatment=args.min_samples_treatment,
                min_samples_control=args.min_samples_control,
            )
    elif args.method == "segment":
        # --- Rule-based segmentation path: takes a DMC parquet ---
        if not args.dmc_results:
            raise ValueError("method=segment requires --dmc-results.")
        from ..dmr_segment import call_dmr_rule_segment

        dmc_results = pl.read_parquet(args.dmc_results)
        dmr_results = call_dmr_rule_segment(
            dmc_results,
            min_cpgs=args.min_cpgs if args.min_cpgs is not None else _DMR_DEFAULT_MIN_CPGS,
            min_abs_meth_diff=args.min_abs_meth_diff,
            alpha=args.alpha,
        )
    else:
        # --- Legacy sliding-window path: takes a DMC parquet ---
        if not args.dmc_results:
            raise ValueError("method=sliding_window requires --dmc-results.")

        dmc_results = pl.read_parquet(args.dmc_results)
        dmr_results = call_dmr_sliding_window(
            dmc_results,
            window_bp=args.window_bp,
            step_bp=args.step_bp,
            min_cpgs=args.min_cpgs if args.min_cpgs is not None else _DMR_DEFAULT_MIN_CPGS,
            min_sites_significant=args.min_sites_significant,
            alpha=args.alpha,
            min_abs_meth_diff=args.min_abs_meth_diff,
        )
        # Same region-level q-value post-filter as tl.dmr's sliding_window
        # branch, via the shared apply_region_qfilter helper (BH-corrected
        # combined_qvalue, falling back to combined_pvalue).
        dmr_results = apply_region_qfilter(dmr_results, getattr(args, "min_mean_qvalue", None))

    dmr_results.write_parquet(args.output)
    print(f"DMR results written to {args.output}")
    print(f"Total DMRs called: {len(dmr_results):,}")
    if len(dmr_results) > 0 and "dmr_type" in dmr_results.columns:
        n_hyper = int((dmr_results["dmr_type"] == "hyper").sum())
        n_hypo = int((dmr_results["dmr_type"] == "hypo").sum())
        n_mixed = int((dmr_results["dmr_type"] == "mixed").sum())
        print(f"  Hyper: {n_hyper:,}  Hypo: {n_hypo:,}  Mixed: {n_mixed:,}")
        print(dmr_results.head(10))

    tsv_suppressed, tsv_path_opt, _, _ = _cli_tsv_opts(args)
    if not tsv_suppressed and len(dmr_results) > 0:
        from ..export import dmr_to_tsv
        from ..methyldata import MethylData

        md_tmp = MethylData(obs=pl.DataFrame({"sample_id": []}), store="")
        md_tmp.uns["dmr"] = dmr_results
        tsv_path = tsv_path_opt or _auto_tsv_path(args.output)
        dmr_to_tsv(md_tmp, tsv_path)
        print(f"DMR TSV: {tsv_path}")


def register(sub: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    """Add ``dmr`` to ``sub``."""
    # dmr
    p_dmr = sub.add_parser(
        "dmr", help="DMR calling (chain-merge, tile, segment, or sliding-window)"
    )
    p_dmr.add_argument(
        "--method",
        choices=["chain_merge", "tile", "sliding_window", "segment"],
        default="chain_merge",
        help=(
            "DMR algorithm. "
            "'chain_merge' (default) chains contiguous significant CpGs from "
            "a precomputed DMC parquet into DSS-style regions; tune via "
            "--preset. "
            "'tile' pools reads across CpGs within each fixed-size tile and "
            "runs one test per tile. "
            "'sliding_window' takes a precomputed DMC parquet and combines "
            "per-CpG p-values with signed Stouffer's Z (legacy). "
            "'segment' rule-based 3-state segmentation on meth_diff signal "
            "with Stouffer-combined per-segment p-values."
        ),
    )
    p_dmr.add_argument("--output", required=True)

    # Tile-method options
    p_dmr.add_argument("--methylstore", help="(tile only) Path to filtered Parquet methylstore.")
    p_dmr.add_argument("--samplesheet", help="(tile only) CSV: sample_id, group, path.")
    p_dmr.add_argument("--treatment-group", help="(tile only) Group label for treatment samples.")
    p_dmr.add_argument("--control-group", help="(tile only) Group label for control samples.")
    p_dmr.add_argument(
        "--tile-size-bp", type=int, default=1000, help="(tile only) Tile width in bp. Default 1000."
    )
    p_dmr.add_argument(
        "--min-cpgs-per-tile",
        type=int,
        default=5,
        help="(tile only) Minimum CpGs per tile per sample.",
    )
    p_dmr.add_argument(
        "--test",
        choices=list(PUBLIC_ENGINES),
        default="lr",
        help="(tile only) Statistical test applied to tile-level counts. "
        "Default 'lr': quasi-binomial LR with McCullagh-Nelder dispersion.",
    )
    p_dmr.add_argument(
        "--empirical-fdr",
        action="store_true",
        default=False,
        help="(tile only) permutation-based empirical FDR.",
    )
    p_dmr.add_argument(
        "--n-perm",
        type=int,
        default=100,
        help="(tile only) Number of permutations when --empirical-fdr is set.",
    )
    p_dmr.add_argument(
        "--perm-seed",
        type=int,
        default=42,
        help="(tile only) Seed for permutation RNG.",
    )
    p_dmr.add_argument(
        "--canonical-only",
        dest="canonical_only",
        action="store_true",
        default=False,
        help=(
            "(tile only) Test only the fixed human-style chromosome set (1-22, X, "
            "Y, M/MT, with or without a chr prefix) of the store's partitions; the "
            "same set is used by every --empirical-fdr permutation. The other "
            "methods inherit the DMC parquet's chromosomes and reject the flag: "
            "run epykit dmc --canonical-only instead. Default: every partition."
        ),
    )
    # Default is union, matching ep.tl.dmr (no pp.unite step). --unite forces
    # intersect (tiles covered in ALL samples); --no-unite explicitly selects
    # the default union, kept for backward compatibility.
    p_dmr.add_argument(
        "--unite",
        action="store_true",
        dest="unite",
        default=False,
        help="(tile only) Intersect: restrict to tiles covered in ALL samples "
        "(default: union, matching ep.tl.dmr).",
    )
    p_dmr.add_argument(
        "--no-unite",
        action="store_false",
        dest="unite",
        default=False,
        help="(tile only) Union: test tiles covered in at least one sample (the default).",
    )
    _add_min_samples_args(p_dmr, scope_help_prefix="(tile only) ")
    p_dmr.add_argument(
        "--allow-n1",
        action="store_true",
        default=False,
        help=(
            "(tile only) Permit n=1 per group. The pooled Fisher exact engine "
            "is used automatically in this case (it is anti-conservative -- do "
            "not trust borderline calls). Default is to refuse, since "
            "between-replicate variance is ignored under this fallback."
        ),
    )

    # Sliding-window-method options
    p_dmr.add_argument(
        "--dmc-results",
        help="(chain_merge, sliding_window, segment) Parquet file from 'epykit dmc'",
    )
    p_dmr.add_argument("--window-bp", type=int, default=500)
    p_dmr.add_argument("--step-bp", type=int, default=250)
    p_dmr.add_argument(
        "--min-cpgs",
        type=int,
        default=None,
        help="Minimum CpGs per DMR. Default when unset: the active --preset's "
        "value if a preset is given, otherwise 5 (chain_merge, "
        "sliding_window, and segment all share this default).",
    )
    p_dmr.add_argument("--min-sites-significant", type=int, default=3)

    # Chain-merge-method options (DSS callDMR semantics). Knob defaults match
    # call_dmr_chain_merge's signature so --preset bundles apply unless a knob
    # is overridden explicitly. min_cpgs comes from the preset / engine default.
    p_dmr.add_argument(
        "--preset",
        choices=["strict", "default", "permissive"],
        default=None,
        help="(chain_merge only) Parameter bundle from DMR_PRESETS. Explicit "
        "knob flags override the bundled value.",
    )
    p_dmr.add_argument(
        "--dis-merge-bp",
        type=int,
        default=500,
        help="(chain_merge only) Max bp gap between consecutive "
        "significant CpGs in a chain. Highest-leverage knob.",
    )
    p_dmr.add_argument(
        "--pct-sig",
        type=float,
        default=0.5,
        help="(chain_merge only) Min fraction of CpGs in a span that must be significant.",
    )
    p_dmr.add_argument(
        "--minlen-bp", type=int, default=50, help="(chain_merge only) Min DMR span length in bp."
    )
    p_dmr.add_argument(
        "--use-q-for-sig",
        action="store_true",
        default=False,
        help="(chain_merge only) Gate significance on qvalue "
        "instead of pvalue when a qvalue column is present.",
    )

    # Shared filters
    p_dmr.add_argument("--alpha", type=float, default=0.05)
    p_dmr.add_argument("--min-abs-meth-diff", type=float, default=0.1)
    p_dmr.add_argument(
        "--min-mean-qvalue",
        type=float,
        default=0.05,
        help="(chain_merge, sliding_window, tile) Region-level q-value cutoff "
        "applied as a post-filter, matching ep.tl.dmr (default 0.05). "
        "The filter is strict (q < cutoff), so set to a value above 1.0 "
        "(e.g. 1.1) to disable -- 1.0 still drops regions with q == 1.0.",
    )
    p_dmr.add_argument(
        "--no-tsv",
        action="store_true",
        dest="no_tsv",
        default=False,
        help="Suppress the sibling .tsv auto-emit.",
    )
    p_dmr.add_argument(
        "--tsv",
        dest="tsv_path",
        default=None,
        help="Override sibling table path. .csv suffix -> comma delim.",
    )
    # Deprecated csv* aliases (see dmc) -- honoured but warn; hidden from --help.
    p_dmr.add_argument(
        "--no-csv",
        action="store_true",
        dest="no_csv",
        default=False,
        help=argparse.SUPPRESS,
    )
    p_dmr.add_argument("--csv", dest="csv_path", default=None, help=argparse.SUPPRESS)
    p_dmr.set_defaults(func=_cmd_dmr)
