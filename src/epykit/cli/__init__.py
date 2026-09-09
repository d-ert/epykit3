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

Package layout: ``_common`` holds the shared helpers. ``_ingest``, ``_dmc``,
``_dmr`` and ``_downstream`` each expose ``register(sub)`` and add their
subcommands; ``build_parser`` calls the four registrars in that order, which
fixes the top-level ``--help`` order.
"""

import argparse
import sys

from . import _dmc, _dmr, _downstream, _ingest
from ._common import _configure_logging


def build_parser() -> argparse.ArgumentParser:
    """Construct the epykit CLI argument parser.

    Extracted from ``main`` so tests can introspect flags/defaults without
    spawning a subprocess.
    """
    from .. import __version__

    ap = argparse.ArgumentParser(prog="epykit", description="Methylation Parquet store tools")
    ap.add_argument(
        "--version",
        action="version",
        version=f"epykit {__version__}",
    )
    ap.add_argument(
        "-v", "--verbose", action="count", default=0, help="Increase logging verbosity (-v: DEBUG)"
    )
    ap.add_argument(
        "-q",
        "--quiet",
        action="count",
        default=0,
        help="Decrease logging verbosity (-q: WARNING and above)",
    )
    sub = ap.add_subparsers(dest="cmd", required=True)
    _ingest.register(sub)
    _dmc.register(sub)
    _dmr.register(sub)
    _downstream.register(sub)
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
        raise SystemExit(f"error: {exc}") from exc
