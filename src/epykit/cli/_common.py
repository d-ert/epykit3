"""Helpers shared by the epykit command modules.

Logging setup, the sibling-table (``--tsv``) options, the shared
``--min-samples-*`` arguments, samplesheet group reading and the n=1 and
footgun checks that ``dmc`` and ``dmr`` both run. Command modules import
from here; nothing here imports a command module.
"""

import argparse
import logging
import os
import warnings
from pathlib import Path


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
        # The CLI package logger: keeps the pre-split stderr line for this
        # warning identical to when the helper lived in epykit/cli.py.
        logging.getLogger("epykit.cli").warning(
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
        "--min-samples-treatment",
        type=int,
        default=0,
        dest="min_samples_treatment",
        help=(
            f"{scope_help_prefix}Per-site minimum number of treatment samples "
            f"with non-zero coverage. Sites failing the threshold are NaN'd "
            f"before FDR. Useful with --no-unite ."
        ),
    )
    p.add_argument(
        "--min-samples-control",
        type=int,
        default=0,
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
            group = row["group"]
            sample_id = row["sample_id"]
            samples_by_group.setdefault(group, []).append(sample_id)

    treatment_samples = samples_by_group.get(treatment_group)
    control_samples = samples_by_group.get(control_group)

    if not treatment_samples:
        raise ValueError(f"No samples found for group '{treatment_group}'")
    if not control_samples:
        raise ValueError(f"No samples found for group '{control_group}'")

    return treatment_samples, control_samples


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
            UserWarning,
            stacklevel=2,
        )
        args.test = "fisher"
    elif args.test == "fisher":
        warnings.warn(
            "test='fisher' is anti-conservative; prefer 'lr' at n >= 2.",
            UserWarning,
            stacklevel=2,
        )
    if (not args.unite) and args.min_samples_treatment == 0 and args.min_samples_control == 0:
        warnings.warn(
            f"Union mode (the default) + min_samples_*=0 will test {unit} "
            f"covered in only one sample per group. Recommended: "
            f"--min-samples-treatment 2 --min-samples-control 2, or pass "
            f"--unite to restrict to sites covered in all samples.",
            UserWarning,
            stacklevel=2,
        )


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
