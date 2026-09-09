"""The ``dmc`` command: per-CpG differential methylation calling."""

from __future__ import annotations

import argparse

from .. import dmc
from .._dmc_engines import PUBLIC_ENGINES
from ._common import (
    _add_min_samples_args,
    _auto_tsv_path,
    _cli_n1_and_footgun_checks,
    _cli_tsv_opts,
    _read_samplesheet_groups,
)

# The engine that reads the DSS-style count smoothing; every other engine
# ignores the knobs, so the CLI refuses the combination instead.
_SMOOTHING_ENGINE = "lr"


def _smoothing_opts(args: argparse.Namespace) -> tuple[bool, int]:
    """``(smoothing, smoothing_span_bp)`` from the parsed flags.

    Exits with a usage error when smoothing is on with a non-positive span.
    The engine check runs separately, after ``--allow-n1`` may have
    resolved ``--test`` to ``fisher``.
    """
    smoothing = bool(getattr(args, "smoothing", False))
    span = int(getattr(args, "smoothing_span_bp", 500))
    if smoothing and span <= 0:
        raise SystemExit(
            f"error: --smoothing-span-bp must be a positive number of base pairs "
            f"when --smoothing is set (got {span})."
        )
    return smoothing, span


def _refuse_unconsumed_smoothing(args: argparse.Namespace, *, reason: str) -> None:
    """Exit when ``--smoothing`` was given to a path that does not read it."""
    if getattr(args, "smoothing", False):
        raise SystemExit(
            f"error: --smoothing is an option of the {_SMOOTHING_ENGINE} engine and "
            f"{reason}. Drop --smoothing, or run the binary treatment / control path "
            f"with --test {_SMOOTHING_ENGINE}."
        )


def _cmd_dmc(args: argparse.Namespace):
    """Handler for 'dmc' subcommand."""
    canonical_only = bool(getattr(args, "canonical_only", False))
    smoothing, smoothing_span_bp = _smoothing_opts(args)

    # formula/contrast path uses ALL samples from the
    # samplesheet rather than binary case/control. We build a tiny
    # MethylData on the fly so tl.dmc can resolve the contrast against
    # md.obs.
    if args.formula is not None or args.contrast is not None:
        # The GLM contrast engine never reads the smoothing knobs; refuse
        # before any store is opened rather than silently drop the flag.
        _refuse_unconsumed_smoothing(
            args, reason="the --formula / --contrast path does not consume it"
        )
        from .. import read_bismark
        from .. import tl as _tl

        covariates = [c.strip() for c in args.covariates.split(",")] if args.covariates else None
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
            canonical_only=canonical_only,
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
    # After the n=1 check: --allow-n1 may have resolved --test to fisher,
    # which does not read the smoothing knobs either.
    if args.test != _SMOOTHING_ENGINE:
        _refuse_unconsumed_smoothing(args, reason=f"--test {args.test} does not consume it")

    print(f"Treatment samples: {treatment_samples}")
    print(f"Control samples:   {control_samples}")
    print(f"Test:              {args.test}")
    print(f"Unite mode:        {'intersect' if args.unite else 'union'}")
    if args.min_samples_treatment or args.min_samples_control:
        print(
            f"Per-site guards:   min_samples_treatment={args.min_samples_treatment}, "
            f"min_samples_control={args.min_samples_control}"
        )
    if smoothing:
        print(f"Smoothing:         DSS-style box, span {smoothing_span_bp} bp")
    if canonical_only:
        print("Chromosomes:       canonical set only (--canonical-only)")

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
        smoothing=smoothing,
        smoothing_span_bp=smoothing_span_bp,
        canonical_only=canonical_only,
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

        from ..export import dmc_to_tsv
        from ..methyldata import MethylData

        # Build a transient MethylData carrying just the dmc result so the
        # writer can re-use the same delegation path the API uses.
        obs = pl.DataFrame({"sample_id": treatment_samples + control_samples})
        md_tmp = MethylData(obs=obs, store=str(args.methylstore))
        md_tmp.varm["dmc_lr"] = results
        md_tmp.uns["dmc"] = {"last_key": "dmc_lr"}

        sig_path = tsv_path or _auto_tsv_path(args.output, suffix=".significant")
        dmc_to_tsv(md_tmp, sig_path, alpha=tsv_alpha)
        print(f"  Significant TSV:    {sig_path}")
        if tsv_full:
            full_path = _auto_tsv_path(args.output)
            dmc_to_tsv(md_tmp, full_path, full=True)
            print(f"  Full TSV:           {full_path}")


def register(sub: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    """Add ``dmc`` to ``sub``."""
    # dmc
    p_dmc = sub.add_parser("dmc", help="Differential methylation calling (per-CpG)")
    p_dmc.add_argument("--methylstore", required=True)
    p_dmc.add_argument("--samplesheet", required=True, help="CSV: sample_id, group, path")
    p_dmc.add_argument("--treatment-group", required=True)
    p_dmc.add_argument("--control-group", required=True)
    p_dmc.add_argument("--output", required=True)
    p_dmc.add_argument(
        "--test",
        choices=list(PUBLIC_ENGINES),
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
        "--formula",
        default=None,
        help=(
            "patsy formula on md.obs columns (e.g. '~ group'). "
            "Triggers the GLM-contrast path; pair with --contrast."
        ),
    )
    p_dmc.add_argument(
        "--contrast",
        default=None,
        help=(
            "contrast specification. Either a single column "
            "name (continuous covariate primary effect), a factor name "
            "for a joint F-test (e.g. 'group'), or a patsy linear "
            "combination ('group[T.KO] - group[T.WT]')."
        ),
    )
    p_dmc.add_argument(
        "--covariates",
        default=None,
        help="Comma-separated list of nuisance covariate columns on md.obs.",
    )
    # Unite mode. Default is union (sites in at least one sample), matching
    # the ep.tl.dmc API with no prior pp.unite step. Pass --unite to restrict
    # to sites covered in ALL samples (intersect). --no-unite is kept as an
    # explicit, backward-compatible way to request the default (union).
    p_dmc.add_argument(
        "--unite",
        action="store_true",
        dest="unite",
        default=False,
        help="Intersect: restrict to sites covered in ALL samples "
        "(default: union -- sites in at least one sample, matching ep.tl.dmc).",
    )
    p_dmc.add_argument(
        "--no-unite",
        action="store_false",
        dest="unite",
        default=False,
        help="Union: include sites covered in at least one sample (the default).",
    )
    _add_min_samples_args(p_dmc)
    p_dmc.add_argument(
        "--allow-n1",
        action="store_true",
        default=False,
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
        "--dispersion",
        choices=["site", "eb", "shrink", "chrom"],
        default="eb",
        help="Dispersion estimator for lr/glm (default: eb, matching ep.tl.dmc).",
    )
    p_dmc.add_argument(
        "--reference",
        choices=["adaptive", "F", "chi2"],
        default="adaptive",
        help="Reference distribution for the lr statistic (default: adaptive).",
    )
    p_dmc.add_argument(
        "--fdr-method",
        dest="fdr_method",
        default="fdr_bh",
        help="Multiple-testing correction method (default: fdr_bh).",
    )
    # DSS-style count smoothing (an lr engine option, matching
    # ep.tl.dmc(smoothing=..., smoothing_span_bp=...)). Refused for the
    # other engines and for the --formula / --contrast path, which do not
    # read it, so the flag is never silently ignored.
    p_dmc.add_argument(
        "--smoothing",
        dest="smoothing",
        action="store_true",
        default=False,
        help=(
            "DSS-style per-sample count smoothing for --test lr: replace each "
            "sample's raw counts by a uniform-box average over the CpGs within "
            "+/- half the span before the test (DMLfit.multiFactor(smoothing=TRUE)). "
            "Default: off. Not accepted with the other engines or with --formula / "
            "--contrast."
        ),
    )
    p_dmc.add_argument(
        "--smoothing-span-bp",
        dest="smoothing_span_bp",
        type=int,
        default=500,
        help="Full smoothing window in bp for --smoothing (default: 500, the DSS default).",
    )
    p_dmc.add_argument(
        "--canonical-only",
        dest="canonical_only",
        action="store_true",
        default=False,
        help=(
            "Test only the fixed human-style chromosome set (1-22, X, Y, M/MT, "
            "with or without a chr prefix) of the store's partitions; other "
            "contigs are dropped before the test and the FDR correction. "
            "Applies to the binary and the --formula / --contrast path. "
            "Default: test every partition."
        ),
    )
    p_dmc.add_argument(
        "--no-tsv",
        action="store_true",
        dest="no_tsv",
        default=False,
        help="Suppress the sibling .significant.tsv auto-emit.",
    )
    p_dmc.add_argument(
        "--tsv",
        dest="tsv_path",
        default=None,
        help=(
            "Override sibling table path. Suffix .csv selects comma "
            "delimiter; otherwise tab. Implies the file is written."
        ),
    )
    p_dmc.add_argument(
        "--tsv-alpha",
        dest="tsv_alpha",
        type=float,
        default=0.05,
        help="qvalue threshold for the significant-only table. Default 0.05.",
    )
    p_dmc.add_argument(
        "--tsv-full",
        dest="tsv_full",
        action="store_true",
        default=False,
        help="Also write the full (unfiltered) table next to the parquet.",
    )
    # Deprecated csv* aliases (epykit writes TSV by default) -- still honoured,
    # but emit a deprecation warning. Hidden from --help to steer users to --tsv*.
    p_dmc.add_argument(
        "--no-csv",
        action="store_true",
        dest="no_csv",
        default=False,
        help=argparse.SUPPRESS,
    )
    p_dmc.add_argument("--csv", dest="csv_path", default=None, help=argparse.SUPPRESS)
    p_dmc.add_argument(
        "--csv-alpha",
        dest="csv_alpha",
        type=float,
        default=0.05,
        help=argparse.SUPPRESS,
    )
    p_dmc.add_argument(
        "--csv-full",
        dest="csv_full",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )
    p_dmc.set_defaults(func=_cmd_dmc)
