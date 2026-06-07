"""Command-line entry point for epykit.

Default DMC test is ``lr`` everywhere (CLI, Python API, docstrings) -- the
quasi-binomial likelihood-ratio chi-square with per-site McCullagh-Nelder
dispersion. Closed-form on streaming (S0_g, S1_g, Sigmam^2/n_g) accumulators,
recommended at n >= 2 replicates per group.

CLI surface:
* ``dmc`` -- per-CpG calling with ``--test {lr,glm,welch_t,fisher}``,
  ``--min-samples-treatment`` / ``--min-samples-control``
  filters, and ``--allow-n1`` to opt into the (anti-conservative) Fisher fallback when
  there are fewer than 2 replicates per group.
* ``dmr`` -- ``--method {chain_merge,tile,sliding_window,segment}`` (default
  ``chain_merge``). chain_merge / sliding_window / segment take a DMC parquet
  (``--dmc-results``); the tile path takes a methylstore + samplesheet and
  pools reads per tile.
"""

import argparse
import logging
import os
import sys
import warnings
from pathlib import Path
from .convert import convert_sample
from . import filter, dmc


def _auto_tsv_path(parquet_path: str, *, suffix: str = "") -> str:
    """Derive a sibling .tsv path from a --output parquet path.

    ``dmc.parquet`` -> ``dmc.significant.tsv`` (suffix=".significant")
    ``dmr.parquet`` -> ``dmr.tsv``             (suffix="")
    Strips a ``.parquet`` extension if present; otherwise appends.
    """
    p = Path(parquet_path)
    stem = p.stem if p.suffix.lower() == ".parquet" else p.name
    return str(p.with_name(f"{stem}{suffix}.tsv"))


def _cli_tsv_opts(args):
    """Resolve the --tsv* auto-emit options, honouring the deprecated --csv* aliases.

    epykit writes tab-delimited TSV by default, so ``--csv`` / ``--no-csv`` /
    ``--csv-full`` / ``--csv-alpha`` and ``EPYKIT_NO_AUTO_CSV`` were renamed to
    ``--tsv`` / ``--no-tsv`` / ``--tsv-full`` / ``--tsv-alpha`` and
    ``EPYKIT_NO_AUTO_TSV``. The old names still work (same code path) but emit a
    deprecation warning. Returns ``(suppressed, path, full, alpha)``.
    """
    env_csv = os.environ.get("EPYKIT_NO_AUTO_CSV") in ("1", "true", "True")
    used_old = (
        getattr(args, "no_csv", False)
        or getattr(args, "csv_path", None) is not None
        or getattr(args, "csv_full", False)
        or getattr(args, "csv_alpha", 0.05) != 0.05
        or env_csv
    )
    if used_old:
        logging.getLogger(__name__).warning(
            "The --csv / --no-csv / --csv-full / --csv-alpha flags and "
            "EPYKIT_NO_AUTO_CSV are deprecated aliases for --tsv / --no-tsv / "
            "--tsv-full / --tsv-alpha and EPYKIT_NO_AUTO_TSV (epykit writes "
            "tab-delimited TSV by default). The csv names will be removed in a "
            "future release."
        )
    suppressed = (
        getattr(args, "no_tsv", False)
        or getattr(args, "no_csv", False)
        or os.environ.get("EPYKIT_NO_AUTO_TSV") in ("1", "true", "True")
        or env_csv
    )
    path = getattr(args, "tsv_path", None) or getattr(args, "csv_path", None)
    full = getattr(args, "tsv_full", False) or getattr(args, "csv_full", False)
    tsv_alpha = getattr(args, "tsv_alpha", 0.05)
    csv_alpha = getattr(args, "csv_alpha", 0.05)
    alpha = tsv_alpha if tsv_alpha != 0.05 else csv_alpha
    return suppressed, path, full, alpha


def _write_table_local(df, path: str) -> str:
    """Write a raw Polars frame to ``path`` with suffix-derived delimiter.

    Mirror of ``export._write_table`` for handlers that hold a frame directly
    and don't need to wrap it in a stub ``MethylData``.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    sep = "," if str(path).lower().endswith(".csv") else "\t"
    df.write_csv(str(out), separator=sep)
    return str(out.resolve())


def _add_min_samples_args(p: argparse.ArgumentParser, scope_help_prefix: str = "") -> None:
    """Register ``--min-samples-treatment`` and ``--min-samples-control``."""
    p.add_argument(
        "--min-samples-treatment", type=int, default=0,
        dest="min_samples_treatment",
        help=(
            f"{scope_help_prefix}Per-site minimum number of treatment samples "
            f"with non-zero coverage. Sites failing the threshold are NaN'd "
            f"before FDR. Useful with --no-unite ."
        ),
    )
    p.add_argument(
        "--min-samples-control", type=int, default=0,
        dest="min_samples_control",
        help=f"{scope_help_prefix}Per-site minimum number of control samples "
             f"with non-zero coverage.",
    )


def _read_samplesheet_groups(samplesheet: str, treatment_group: str, control_group: str):
    import csv

    with open(samplesheet) as f:
        reader = csv.DictReader(f)
        samples_by_group: dict[str, list[str]] = {}
        for row in reader:
            group     = row["group"]
            sample_id = row["sample_id"]
            samples_by_group.setdefault(group, []).append(sample_id)

    treatment_samples = samples_by_group.get(treatment_group)
    control_samples   = samples_by_group.get(control_group)

    if not treatment_samples:
        raise ValueError(f"No samples found for group '{treatment_group}'")
    if not control_samples:
        raise ValueError(f"No samples found for group '{control_group}'")

    return treatment_samples, control_samples


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


def _cli_n1_and_footgun_checks(args, unit: str = "sites") -> None:
    """Mirror tl.* guards on the CLI side.

    May resolve ``args.test`` to 'fisher' when --allow-n1 is set at n<2.
    """
    treatment_samples, control_samples = args._samples  # set by caller
    n_min = min(len(treatment_samples), len(control_samples))
    if n_min < 2 and not args.allow_n1:
        raise SystemExit(
            f"error: at least 2 replicates per group required "
            f"(treatment={len(treatment_samples)}, control={len(control_samples)}). "
            f"Pass --allow-n1 to opt into the Fisher fallback."
        )
    # D12: --allow-n1 advertises a pooled-Fisher fallback, but the default
    # (and explicit) lr engine has no n=1 path, so the advertised fallback
    # never fired -- the n=1 run silently used lr. When the user is on the
    # lr/auto engine, resolve to fisher so the advertised behavior actually
    # happens. An explicit non-lr engine choice (glm/welch_t) is respected.
    if n_min < 2 and args.allow_n1 and getattr(args, "test", None) in (None, "lr", "auto"):
        warnings.warn(
            "n=1 per group with --allow-n1: resolving --test to 'fisher' "
            "(pooled Fisher exact) -- the lr engine has no n=1 fallback. "
            "Fisher is anti-conservative; do not trust borderline calls.",
            UserWarning, stacklevel=2,
        )
        args.test = "fisher"
    elif args.test == "fisher":
        warnings.warn(
            "test='fisher' is anti-conservative; prefer 'lr' at n >= 2.",
            UserWarning, stacklevel=2,
        )
    if (not args.unite) and args.min_samples_treatment == 0 and args.min_samples_control == 0:
        warnings.warn(
            f"Union mode (the default) + min_samples_*=0 will test {unit} "
            f"covered in only one sample per group. Recommended: "
            f"--min-samples-treatment 2 --min-samples-control 2, or pass "
            f"--unite to restrict to sites covered in all samples.",
            UserWarning, stacklevel=2,
        )


def _cmd_dmc(args: argparse.Namespace):
    """Handler for 'dmc' subcommand."""

    # formula/contrast path uses ALL samples from the
    # samplesheet rather than binary case/control. We build a tiny
    # MethylData on the fly so tl.dmc can resolve the contrast against
    # md.obs.
    if args.formula is not None or args.contrast is not None:
        from . import read_bismark, tl as _tl
        covariates = (
            [c.strip() for c in args.covariates.split(",")]
            if args.covariates else None
        )
        # All groups from the samplesheet
        import csv
        with open(args.samplesheet) as fh:
            groups = sorted({row["group"] for row in csv.DictReader(fh)})
        md = read_bismark(
            args.samplesheet,
            treatment_group=args.treatment_group,
            control_group=args.control_group,
            groups=groups,
            store_dir=str(args.methylstore),
        )
        _tl.dmc(
            md,
            test=args.test,
            formula=args.formula,
            contrast=args.contrast,
            covariates=covariates,
            min_samples_treatment=args.min_samples_treatment,
            min_samples_control=args.min_samples_control,
            dispersion=args.dispersion,
            reference=args.reference,
            fdr_method=args.fdr_method,
        )
        key = md.uns.get("dmc", {}).get("last_key", "dmc_glm_contrast")
        results = md.varm.get(key)
        if results is None:
            raise RuntimeError("dmc contrast path produced no results")
        results.write_parquet(args.output)
        print(f"DMC (contrast) results written to {args.output}")
        return

    treatment_samples, control_samples = _read_samplesheet_groups(
        args.samplesheet, args.treatment_group, args.control_group
    )
    args._samples = (treatment_samples, control_samples)
    _cli_n1_and_footgun_checks(args, unit="sites")

    print(f"Treatment samples: {treatment_samples}")
    print(f"Control samples:   {control_samples}")
    print(f"Test:              {args.test}")
    print(f"Unite mode:        {'intersect' if args.unite else 'union'}")
    if args.min_samples_treatment or args.min_samples_control:
        print(
            f"Per-site guards:   min_samples_treatment={args.min_samples_treatment}, "
            f"min_samples_control={args.min_samples_control}"
        )

    # Stream through a DMCStore so BH never holds the full DataFrame in
    # memory alongside its own pvalue/qvalue arrays. Materialise once
    # at the end purely for the parquet write.
    dmc_store = dmc.process_chromosomes_dmc(
        args.methylstore,
        treatment_samples,
        control_samples,
        test=args.test,
        unite=args.unite,
        min_samples_treatment=args.min_samples_treatment,
        min_samples_control=args.min_samples_control,
        dispersion=args.dispersion,
        reference=args.reference,
        return_store=True,
    )
    dmc_store = dmc.apply_multiple_testing_correction(dmc_store, method=args.fdr_method)
    results = dmc_store.to_dataframe()
    results.write_parquet(args.output)
    print(f"DMC results written to {args.output}")
    print(f"  DMC store at:      {dmc_store.path}")
    n_sig = int((results["qvalue"] < 0.05).sum()) if "qvalue" in results.columns else 0
    print(f"  Total sites tested: {len(results):,}")
    print(f"  Significant (q<0.05): {n_sig:,}")

    tsv_suppressed, tsv_path, tsv_full, tsv_alpha = _cli_tsv_opts(args)
    if not tsv_suppressed and len(results.columns) > 0:
        import polars as pl
        from .methyldata import MethylData
        from .export import dmc_to_tsv
        # Build a transient MethylData carrying just the dmc result so the
        # writer can re-use the same delegation path the API uses.
        obs = pl.DataFrame({"sample_id": treatment_samples + control_samples})
        md_tmp = MethylData(obs=obs, store=str(args.methylstore))
        md_tmp.varm["dmc_lr"] = results
        md_tmp.uns["dmc"] = {"last_key": "dmc_lr"}

        sig_path = tsv_path or _auto_tsv_path(
            args.output, suffix=".significant"
        )
        dmc_to_tsv(md_tmp, sig_path, alpha=tsv_alpha)
        print(f"  Significant TSV:    {sig_path}")
        if tsv_full:
            full_path = _auto_tsv_path(args.output)
            dmc_to_tsv(md_tmp, full_path, full=True)
            print(f"  Full TSV:           {full_path}")


def _cmd_dmr(args: argparse.Namespace):
    """Handler for 'dmr' subcommand."""
    import polars as pl
    from .dmr import (
        _DMR_DEFAULT_MIN_CPGS,
        apply_region_qfilter,
        call_dmr_chain_merge,
        call_dmr_sliding_window,
        call_dmr_tile_based,
        resolve_layer_min_cpgs,
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
        dmr_results = apply_region_qfilter(
            dmr_results, getattr(args, "min_mean_qvalue", None)
        )
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

        dmr_results = call_dmr_tile_based(
            methylstore_path=args.methylstore,
            samples_treatment=treatment_samples,
            samples_control=control_samples,
            tile_size_bp=args.tile_size_bp,
            test=args.test,
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
            from .dmr import empirical_fdr_for_dmr
            dmr_results = empirical_fdr_for_dmr(
                methylstore_path=args.methylstore,
                samples_treatment=treatment_samples,
                samples_control=control_samples,
                observed_dmr=dmr_results,
                n_perm=args.n_perm,
                seed=args.perm_seed,
                tile_size_bp=args.tile_size_bp,
                test=args.test,
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
        from .dmr_segment import call_dmr_rule_segment
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
        dmr_results = apply_region_qfilter(
            dmr_results, getattr(args, "min_mean_qvalue", None)
        )

    dmr_results.write_parquet(args.output)
    print(f"DMR results written to {args.output}")
    print(f"Total DMRs called: {len(dmr_results):,}")
    if len(dmr_results) > 0 and "dmr_type" in dmr_results.columns:
        n_hyper = int((dmr_results["dmr_type"] == "hyper").sum())
        n_hypo  = int((dmr_results["dmr_type"] == "hypo").sum())
        n_mixed = int((dmr_results["dmr_type"] == "mixed").sum())
        print(f"  Hyper: {n_hyper:,}  Hypo: {n_hypo:,}  Mixed: {n_mixed:,}")
        print(dmr_results.head(10))

    tsv_suppressed, tsv_path_opt, _, _ = _cli_tsv_opts(args)
    if not tsv_suppressed and len(dmr_results) > 0:
        from .methyldata import MethylData
        from .export import dmr_to_tsv
        md_tmp = MethylData(obs=pl.DataFrame({"sample_id": []}), store="")
        md_tmp.uns["dmr"] = dmr_results
        tsv_path = tsv_path_opt or _auto_tsv_path(args.output)
        dmr_to_tsv(md_tmp, tsv_path)
        print(f"DMR TSV: {tsv_path}")


def _cmd_annotate(args: argparse.Namespace):
    """Handler for 'annotate' subcommand."""
    import polars as pl

    sites = pl.read_parquet(args.input)

    if args.gtf:
        from .annotate import annotate_features
        sites = annotate_features(
            sites,
            args.gtf,
            source="gtf",
            promoter_upstream_bp=args.promoter_upstream_bp,
            promoter_downstream_bp=args.promoter_downstream_bp,
        )
        print("Gene feature annotation complete.")

    if args.cpg_islands:
        from .annotate import annotate_cpg_islands
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
    from .qc import global_methylation_report, coverage_uniformity

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
    from .dmr import smooth_methylation_gaussian

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
    from .methyldata import MethylData
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
    from .methyldata import MethylData
    from . import pp as pp_mod
    md = MethylData.load(args.md)
    pp_mod.aggregate_regions(
        md, args.bed,
        region_id_col=args.region_id_col,
        output_store=args.output_store,
        min_cpgs_per_region=args.min_cpgs_per_region,
    )
    md.save(args.md)
    print(f"Region aggregation complete; methylstore = {md.store}")


def _cmd_export(args: argparse.Namespace):
    """Handler for 'export' subcommand."""
    from .methyldata import MethylData
    from .export import to_bedgraph, to_bigwig, dmcs_to_bed, dmrs_to_bed
    md = MethylData.load(args.md)
    fmt = args.export_cmd
    if fmt == "bedgraph":
        to_bedgraph(md, args.sample, args.output, value=args.value)
    elif fmt == "bigwig":
        to_bigwig(md, args.sample, args.output, value=args.value)
    elif fmt == "dmcs-bed":
        dmcs_to_bed(md, args.output, alpha=args.alpha,
                    min_abs_diff=args.min_abs_diff, test=args.test)
    elif fmt == "dmrs-bed":
        dmrs_to_bed(md, args.output)
    else:
        raise SystemExit(f"unknown export format: {fmt}")
    print(f"Wrote {args.output}")


def _configure_logging(verbosity: int) -> None:
    """Configure logging only when running as a CLI.

    Library code never calls ``logging.basicConfig``; doing so at import time
    would override the host application's logging configuration. The CLI is
    allowed to configure logging because the user has explicitly invoked it.

    ``verbosity`` is the net of ``-v`` (count) minus ``-q`` (count):
      0  -> INFO (default)
      >=1 -> DEBUG
      <=-1 -> WARNING
    """
    if verbosity >= 1:
        level = logging.DEBUG
    elif verbosity <= -1:
        level = logging.WARNING
    else:
        level = logging.INFO
    # Guard against overriding handlers a host program (e.g. tests, notebooks)
    # may already have installed.
    if not logging.getLogger().handlers:
        logging.basicConfig(
            level=level,
            format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        )
    else:
        logging.getLogger().setLevel(level)


def build_parser() -> argparse.ArgumentParser:
    """Construct the epykit CLI argument parser.

    Extracted from ``main`` so tests can introspect flags/defaults without
    spawning a subprocess.
    """
    from . import __version__

    ap  = argparse.ArgumentParser(
        prog="epykit", description="Methylation Parquet store tools"
    )
    ap.add_argument(
        "--version", action="version", version=f"epykit {__version__}",
    )
    ap.add_argument("-v", "--verbose", action="count", default=0,
                    help="Increase logging verbosity (-v: DEBUG)")
    ap.add_argument("-q", "--quiet", action="count", default=0,
                    help="Decrease logging verbosity (-q: WARNING and above)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    # convert
    p_conv = sub.add_parser(
        "convert",
        help="Convert a Bismark .cov or MethylDackel .bedGraph file to Parquet",
    )
    p_conv.add_argument("--input",        required=True)
    p_conv.add_argument("--sample-id",    required=True)
    p_conv.add_argument("--output-dir",   required=True)
    p_conv.add_argument(
        "--format", choices=["bismark", "methyldackel"], default="bismark",
        help=(
            "Source file format. 'bismark' (default) for .cov[.gz] files "
            "produced by bismark_methylation_extractor / bismark2bedGraph. "
            "'methyldackel' for .bedGraph[.gz] files produced by "
            "MethylDackel extract -- same 6-column layout, with the leading "
            "track header skipped automatically."
        ),
    )
    p_conv.add_argument(
        "--context", choices=["CpG", "CHG", "CHH"], default="CpG",
    )
    p_conv.add_argument("--reference-fasta")
    p_conv.add_argument(
        "--merge-cpg", dest="merge_cpg", action="store_true",
        help=(
            "Merge symmetric CpG dyads (+/- strand) into one record. This is "
            "the default (matching the Python API merge_strands=True); the flag "
            "is accepted for explicitness."
        ),
    )
    p_conv.add_argument(
        "--no-merge-cpg", dest="merge_cpg", action="store_false",
        help="Disable CpG-dyad merging; keep per-strand records.",
    )
    p_conv.set_defaults(merge_cpg=None, func=_cmd_convert)

    # filter
    p_filt = sub.add_parser("filter", help="Filter low-coverage CpGs")
    p_filt.add_argument("--methylstore",           required=True)
    p_filt.add_argument("--output-dir",            required=True)
    p_filt.add_argument("--min-coverage",          type=int,   default=10)
    p_filt.add_argument("--max-coverage-quantile", type=float, default=0.999)
    p_filt.add_argument("--blacklist-bed")
    p_filt.add_argument("--sample")
    p_filt.set_defaults(func=_cmd_filter)

    # summary
    p_sum = sub.add_parser("summary", help="Per-sample summary statistics")
    p_sum.add_argument("--methylstore", required=True)
    p_sum.add_argument("--sample",      required=True)
    p_sum.add_argument("--output")
    p_sum.set_defaults(func=_cmd_sample_summary)

    # dmc
    p_dmc = sub.add_parser("dmc", help="Differential methylation calling (per-CpG)")
    p_dmc.add_argument("--methylstore",      required=True)
    p_dmc.add_argument("--samplesheet",      required=True,
                       help="CSV: sample_id, group, path")
    p_dmc.add_argument("--treatment-group",  required=True)
    p_dmc.add_argument("--control-group",    required=True)
    p_dmc.add_argument("--output",           required=True)
    p_dmc.add_argument(
        "--test",
        choices=[
            "lr", "glm", "welch_t",
            "fisher",
        ],
        default="lr",
        help=(
            "Statistical test (default: lr). "
            "lr -- Quasi-binomial likelihood-ratio chi-square with per-site "
            "McCullagh-Nelder dispersion. Closed-form on streaming "
            "accumulators; recommended default at n>=2. "
            "glm -- Binomial GLM with covariates (requires a design via "
            "--formula). "
            "welch_t -- Welch t on raw betas. "
            "fisher -- Fisher exact on reads pooled across replicates "
            "(anti-conservative, kept for backward compatibility; warns)."
        ),
    )
    p_dmc.add_argument(
        "--formula", default=None,
        help=(
            "patsy formula on md.obs columns (e.g. '~ group'). "
            "Triggers the GLM-contrast path; pair with --contrast."
        ),
    )
    p_dmc.add_argument(
        "--contrast", default=None,
        help=(
            "contrast specification. Either a single column "
            "name (continuous covariate primary effect), a factor name "
            "for a joint F-test (e.g. 'group'), or a patsy linear "
            "combination ('group[T.KO] - group[T.WT]')."
        ),
    )
    p_dmc.add_argument(
        "--covariates", default=None,
        help="Comma-separated list of nuisance covariate columns on md.obs.",
    )
    # Unite mode. Default is union (sites in at least one sample), matching
    # the ep.tl.dmc API with no prior pp.unite step. Pass --unite to restrict
    # to sites covered in ALL samples (intersect). --no-unite is kept as an
    # explicit, backward-compatible way to request the default (union).
    p_dmc.add_argument(
        "--unite", action="store_true", dest="unite", default=False,
        help="Intersect: restrict to sites covered in ALL samples "
             "(default: union -- sites in at least one sample, matching ep.tl.dmc).",
    )
    p_dmc.add_argument(
        "--no-unite", action="store_false", dest="unite", default=False,
        help="Union: include sites covered in at least one sample (the default).",
    )
    _add_min_samples_args(p_dmc)
    p_dmc.add_argument(
        "--allow-n1", action="store_true", default=False,
        help=(
            "Permit n=1 per group. The pooled Fisher exact engine is used "
            "automatically in this case (it is anti-conservative -- do not "
            "trust borderline calls). Default is to refuse, since "
            "between-replicate variance is ignored under this fallback."
        ),
    )
    # Parity with ep.tl.dmc defaults. Historically the CLI inherited
    # process_chromosomes_dmc's dispersion="site" default while the Python
    # API used "eb", so identical input produced different q-values. The
    # CLI now defaults to the same values tl.dmc uses (M-PKG2).
    p_dmc.add_argument(
        "--dispersion", choices=["site", "eb", "shrink", "chrom"], default="eb",
        help="Dispersion estimator for lr/glm (default: eb, matching ep.tl.dmc).",
    )
    p_dmc.add_argument(
        "--reference", choices=["adaptive", "F", "chi2"], default="adaptive",
        help="Reference distribution for the lr statistic (default: adaptive).",
    )
    p_dmc.add_argument(
        "--fdr-method", dest="fdr_method", default="fdr_bh",
        help="Multiple-testing correction method (default: fdr_bh).",
    )
    p_dmc.add_argument(
        "--no-tsv", action="store_true", dest="no_tsv", default=False,
        help="Suppress the sibling .significant.tsv auto-emit.",
    )
    p_dmc.add_argument(
        "--tsv", dest="tsv_path", default=None,
        help=(
            "Override sibling table path. Suffix .csv selects comma "
            "delimiter; otherwise tab. Implies the file is written."
        ),
    )
    p_dmc.add_argument(
        "--tsv-alpha", dest="tsv_alpha", type=float, default=0.05,
        help="qvalue threshold for the significant-only table. Default 0.05.",
    )
    p_dmc.add_argument(
        "--tsv-full", dest="tsv_full", action="store_true", default=False,
        help="Also write the full (unfiltered) table next to the parquet.",
    )
    # Deprecated csv* aliases (epykit writes TSV by default) -- still honoured,
    # but emit a deprecation warning. Hidden from --help to steer users to --tsv*.
    p_dmc.add_argument(
        "--no-csv", action="store_true", dest="no_csv", default=False,
        help=argparse.SUPPRESS,
    )
    p_dmc.add_argument("--csv", dest="csv_path", default=None, help=argparse.SUPPRESS)
    p_dmc.add_argument(
        "--csv-alpha", dest="csv_alpha", type=float, default=0.05,
        help=argparse.SUPPRESS,
    )
    p_dmc.add_argument(
        "--csv-full", dest="csv_full", action="store_true", default=False,
        help=argparse.SUPPRESS,
    )
    p_dmc.set_defaults(func=_cmd_dmc)

    # dmr
    p_dmr = sub.add_parser(
        "dmr", help="DMR calling (chain-merge, tile, segment, or sliding-window)")
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
    p_dmr.add_argument("--methylstore",
                       help="(tile only) Path to filtered Parquet methylstore.")
    p_dmr.add_argument("--samplesheet",
                       help="(tile only) CSV: sample_id, group, path.")
    p_dmr.add_argument("--treatment-group",
                       help="(tile only) Group label for treatment samples.")
    p_dmr.add_argument("--control-group",
                       help="(tile only) Group label for control samples.")
    p_dmr.add_argument("--tile-size-bp",       type=int,   default=1000,
                       help="(tile only) Tile width in bp. Default 1000.")
    p_dmr.add_argument("--min-cpgs-per-tile",  type=int,   default=5,
                       help="(tile only) Minimum CpGs per tile per sample.")
    p_dmr.add_argument(
        "--test",
        choices=[
            "lr", "glm", "welch_t",
            "fisher",
        ],
        default="lr",
        help="(tile only) Statistical test applied to tile-level counts. "
             "Default 'lr': quasi-binomial LR with McCullagh-Nelder dispersion.",
    )
    p_dmr.add_argument(
        "--empirical-fdr", action="store_true", default=False,
        help="(tile only) permutation-based empirical FDR.",
    )
    p_dmr.add_argument(
        "--n-perm", type=int, default=100,
        help="(tile only) Number of permutations when --empirical-fdr is set.",
    )
    p_dmr.add_argument(
        "--perm-seed", type=int, default=42,
        help="(tile only) Seed for permutation RNG.",
    )
    # Default is union, matching ep.tl.dmr (no pp.unite step). --unite forces
    # intersect (tiles covered in ALL samples); --no-unite explicitly selects
    # the default union, kept for backward compatibility.
    p_dmr.add_argument(
        "--unite", action="store_true", dest="unite", default=False,
        help="(tile only) Intersect: restrict to tiles covered in ALL samples "
             "(default: union, matching ep.tl.dmr).",
    )
    p_dmr.add_argument(
        "--no-unite", action="store_false", dest="unite", default=False,
        help="(tile only) Union: test tiles covered in at least one sample (the default).",
    )
    _add_min_samples_args(p_dmr, scope_help_prefix="(tile only) ")
    p_dmr.add_argument(
        "--allow-n1", action="store_true", default=False,
        help=(
            "(tile only) Permit n=1 per group. The pooled Fisher exact engine "
            "is used automatically in this case (it is anti-conservative -- do "
            "not trust borderline calls). Default is to refuse, since "
            "between-replicate variance is ignored under this fallback."
        ),
    )

    # Sliding-window-method options
    p_dmr.add_argument("--dmc-results",
                       help="(chain_merge, sliding_window, segment) "
                            "Parquet file from 'epykit dmc'")
    p_dmr.add_argument("--window-bp",            type=int,   default=500)
    p_dmr.add_argument("--step-bp",              type=int,   default=250)
    p_dmr.add_argument(
        "--min-cpgs", type=int, default=None,
        help="Minimum CpGs per DMR. Default when unset: the active --preset's "
             "value if a preset is given, otherwise 5 (chain_merge, "
             "sliding_window, and segment all share this default).",
    )
    p_dmr.add_argument("--min-sites-significant",type=int,   default=3)

    # Chain-merge-method options (DSS callDMR semantics). Knob defaults match
    # call_dmr_chain_merge's signature so --preset bundles apply unless a knob
    # is overridden explicitly. min_cpgs comes from the preset / engine default.
    p_dmr.add_argument(
        "--preset", choices=["strict", "default", "permissive"], default=None,
        help="(chain_merge only) Parameter bundle from DMR_PRESETS. Explicit "
             "knob flags override the bundled value.",
    )
    p_dmr.add_argument("--dis-merge-bp",  type=int,   default=500,
                       help="(chain_merge only) Max bp gap between consecutive "
                            "significant CpGs in a chain. Highest-leverage knob.")
    p_dmr.add_argument("--pct-sig",       type=float, default=0.5,
                       help="(chain_merge only) Min fraction of CpGs in a span "
                            "that must be significant.")
    p_dmr.add_argument("--minlen-bp",     type=int,   default=50,
                       help="(chain_merge only) Min DMR span length in bp.")
    p_dmr.add_argument("--use-q-for-sig", action="store_true", default=False,
                       help="(chain_merge only) Gate significance on qvalue "
                            "instead of pvalue when a qvalue column is present.")

    # Shared filters
    p_dmr.add_argument("--alpha",                type=float, default=0.05)
    p_dmr.add_argument("--min-abs-meth-diff",    type=float, default=0.1)
    p_dmr.add_argument(
        "--min-mean-qvalue", type=float, default=0.05,
        help="(chain_merge, sliding_window, tile) Region-level q-value cutoff "
             "applied as a post-filter, matching ep.tl.dmr (default 0.05). "
             "The filter is strict (q < cutoff), so set to a value above 1.0 "
             "(e.g. 1.1) to disable -- 1.0 still drops regions with q == 1.0.",
    )
    p_dmr.add_argument(
        "--no-tsv", action="store_true", dest="no_tsv", default=False,
        help="Suppress the sibling .tsv auto-emit.",
    )
    p_dmr.add_argument(
        "--tsv", dest="tsv_path", default=None,
        help="Override sibling table path. .csv suffix -> comma delim.",
    )
    # Deprecated csv* aliases (see dmc) -- honoured but warn; hidden from --help.
    p_dmr.add_argument(
        "--no-csv", action="store_true", dest="no_csv", default=False,
        help=argparse.SUPPRESS,
    )
    p_dmr.add_argument("--csv", dest="csv_path", default=None, help=argparse.SUPPRESS)
    p_dmr.set_defaults(func=_cmd_dmr)

    # annotate
    p_ann = sub.add_parser(
        "annotate", help="Annotate DMC/DMR results with genomic features"
    )
    p_ann.add_argument("--input",   required=True,
                       help="Parquet file from 'epykit dmc' or 'epykit dmr'")
    p_ann.add_argument("--output",  required=True)
    p_ann.add_argument(
        "--gtf",
        help="Ensembl/UCSC GTF/GFF3 for gene feature annotation",
    )
    p_ann.add_argument(
        "--cpg-islands",
        help="UCSC CpGIsland BED file for CpG context annotation",
    )
    p_ann.add_argument("--promoter-upstream-bp",   type=int, default=2000)
    p_ann.add_argument("--promoter-downstream-bp", type=int, default=200)
    p_ann.add_argument(
        "--no-tsv", action="store_true", dest="no_tsv", default=False,
        help="Suppress the sibling .tsv auto-emit.",
    )
    p_ann.add_argument(
        "--tsv", dest="tsv_path", default=None,
        help="Override sibling table path. .csv suffix -> comma delim.",
    )
    # Deprecated csv* aliases (see dmc) -- honoured but warn; hidden from --help.
    p_ann.add_argument(
        "--no-csv", action="store_true", dest="no_csv", default=False,
        help=argparse.SUPPRESS,
    )
    p_ann.add_argument("--csv", dest="csv_path", default=None, help=argparse.SUPPRESS)
    p_ann.set_defaults(func=_cmd_annotate)

    # qc-report
    p_qc = sub.add_parser("qc-report", help="QC and coverage uniformity report")
    p_qc.add_argument("--methylstore", required=True)
    p_qc.add_argument(
        "--samples", required=True,
        help="Comma-separated list of sample IDs",
    )
    p_qc.add_argument(
        "--output-dir",
        help="Directory for Parquet QC output files (optional)",
    )
    p_qc.add_argument(
        "--no-tsv", action="store_true", dest="no_tsv", default=False,
        help="Suppress the sibling .tsv auto-emit alongside the parquets.",
    )
    # Deprecated csv alias (see dmc) -- honoured but warn; hidden from --help.
    p_qc.add_argument(
        "--no-csv", action="store_true", dest="no_csv", default=False,
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
        "--samples", required=True,
        help="Comma-separated list of sample IDs",
    )
    p_sm.add_argument("--output",    required=True,
                      help="Output directory for smoothed beta chunks")
    p_sm.add_argument("--bandwidth", type=int, default=1000,
                      help="Smoothing bandwidth in bp (default 1000)")
    p_sm.set_defaults(func=_cmd_smooth)

    # report
    p_rep = sub.add_parser(
        "report",
        help="Render an interactive HTML report from a saved MethylData",
    )
    p_rep.add_argument("--md", required=True,
                       help="Path to a directory previously written with md.save(...)")
    p_rep.add_argument("--output", required=True, help="Output HTML file")
    p_rep.add_argument("--title", default=None)
    p_rep.add_argument("--gtf", default=None,
                       help="Optional GTF for a TSS metaplot section")
    p_rep.add_argument("--alpha", type=float, default=0.05)
    p_rep.add_argument("--min-abs-diff", dest="min_abs_diff",
                       type=float, default=0.1)
    p_rep.add_argument("--self-contained", dest="self_contained",
                       action="store_true", default=True,
                       help="Embed Plotly inline so the HTML works offline (default)")
    p_rep.add_argument("--no-self-contained", dest="self_contained",
                       action="store_false",
                       help="Load Plotly from a CDN (smaller file, needs internet)")
    p_rep.set_defaults(func=_cmd_report)

    # aggregate-regions
    p_agg = sub.add_parser(
        "aggregate-regions",
        help="Aggregate CpG counts to user-supplied BED regions",
    )
    p_agg.add_argument("--md", required=True,
                       help="Path to a directory previously written with md.save(...)")
    p_agg.add_argument("--bed", required=True, help="BED file of regions to aggregate to")
    p_agg.add_argument("--output-store", default=None,
                       help="Override output Parquet store path")
    p_agg.add_argument("--region-id-col", default=None,
                       help="BED column name to use as region_id (default: 'name' or chrom:start-end)")
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
        sp.add_argument("--value", default="beta",
                        choices=["beta", "coverage", "N_meth"])
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

    return ap


def main():
    # Help strings and log messages embed unicode (beta, ->, mu, ...). On
    # Windows the default console codec is cp1252 and argparse's
    # `--help` print crashes with UnicodeEncodeError before any
    # subcommand runs. Reconfigure both streams to UTF-8 with
    # replacement so we never crash on a glyph the terminal can't draw.
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            try:
                reconfigure(encoding="utf-8", errors="replace")
            except (OSError, ValueError):
                # Detached / non-text stream; nothing to do.
                pass

    ap = build_parser()
    args = ap.parse_args()
    _configure_logging(verbosity=args.verbose - args.quiet)
    try:
        args.func(args)
    except FileNotFoundError as exc:
        raise SystemExit(f"error: {exc}")


if __name__ == "__main__":
    main()