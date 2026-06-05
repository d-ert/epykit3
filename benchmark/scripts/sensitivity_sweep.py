"""DMR-caller parameter sensitivity sweep.

Phase 1.3 of the GB resubmission plan asks for evidence that the
headline DMR numbers are robust to small perturbations of the locked
``DMR_PRESETS["default"]`` parameter set. This script provides that
evidence by running a one-at-a-time (OAT) sensitivity sweep:

  1. Load a DMC parquet (the output of an ``ep.tl.dmc`` run).
  2. Load a DMR truth table.
  3. For each of six tunables (``alpha``, ``min_abs_meth_diff``,
     ``dis_merge_bp``, ``min_cpgs``, ``pct_sig``, ``minlen_bp``)
     perturb that one parameter by -50% / 0 / +50% from the default
     while holding the others at default, call
     ``ep.dmr.call_dmr_chain_merge`` once per perturbation, and
     score against the truth via ``score_dmr_parquet``.
  4. Emit one parquet row per (parameter, value) cell with the
     headline metrics so a single table can show how each knob
     moves the result.

OAT keeps the design size at 1 + 6 * 2 = 13 cells instead of the
6^3 = 729 full factorial. The point is to *show robustness*, not
to characterise the full response surface -- 13 cells is enough for
a "does the headline survive these perturbations?" claim.

Usage:
    python sensitivity_sweep.py \\
        --dmc-parquet  benchmark/data/study1/.../dmc.parquet \\
        --truth-parquet benchmark/data/study1/ground_truth/dmr_truth.parquet \\
        --out          benchmark/data/study1/sensitivity_sweep.parquet \\
        --tool         epykit \\
        --scenario     cov10
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)

# Ensure the local scripts dir is importable for _epykit_scoring.
_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))


PERTURBATION_FACTORS = (0.5, 1.0, 1.5)

# Tunables and how to perturb each. ``cast`` coerces the final value
# back to the type the caller expects (ints stay ints).
TUNABLES = (
    ("alpha", float),
    ("min_abs_meth_diff", float),
    ("dis_merge_bp", int),
    ("min_cpgs", int),
    ("pct_sig", float),
    ("minlen_bp", int),
)


def _perturbed(default: dict, key: str, factor: float, cast) -> dict:
    """Return a copy of ``default`` with ``key`` scaled by ``factor``.

    Integer parameters are floored to ``>= 1`` after scaling so we
    never pass 0 or negative values to the engine.
    """
    new = dict(default)
    raw = default[key] * factor
    if cast is int:
        raw = max(1, int(round(raw)))
    new[key] = cast(raw)
    return new


def run_sweep(
    dmc_parquet: Path,
    truth_parquet: Path,
    *,
    tool: str,
    scenario: str,
    test: str = "lr",
    preset_name: str = "default",
) -> pl.DataFrame:
    """Run the OAT sweep and return one row per perturbation cell.

    Each row carries ``parameter`` (the knob being perturbed),
    ``factor`` (0.5 / 1.0 / 1.5), ``parameter_value`` (the actual
    value passed to the caller), ``n_dmrs``, plus the threshold-grid
    TPR/FPR/F1 columns from ``score_dmr_parquet``.
    """
    from epykit.dmr import DMR_PRESETS, call_dmr_chain_merge
    from _epykit_scoring import score_dmr_parquet

    if preset_name not in DMR_PRESETS:
        raise ValueError(
            f"unknown preset {preset_name!r}; available: "
            f"{sorted(DMR_PRESETS)}"
        )
    default = DMR_PRESETS[preset_name]

    truth = pl.read_parquet(str(truth_parquet))
    rows: list[dict] = []
    seen_baseline = False  # only emit factor=1.0 cell once

    for key, cast in TUNABLES:
        for factor in PERTURBATION_FACTORS:
            if factor == 1.0:
                if seen_baseline:
                    continue
                seen_baseline = True
            params = _perturbed(default, key, factor, cast)
            try:
                dmrs = call_dmr_chain_merge(
                    str(dmc_parquet), **params,
                )
            except Exception as exc:  # noqa: BLE001 -- log + record empty
                logger.warning(
                    "sweep cell failed: %s=%s (factor=%s): %s",
                    key, params[key], factor, exc,
                )
                rows.append({
                    "parameter": key, "factor": factor,
                    "parameter_value": params[key],
                    "n_dmrs": 0, "error": str(exc),
                })
                continue

            # Persist the per-cell DMR parquet to a temporary location
            # so score_dmr_parquet can read it back (it operates on
            # paths to keep the scoring contract uniform across tools).
            import tempfile
            with tempfile.NamedTemporaryFile(
                suffix=".parquet", delete=False,
            ) as tf:
                tmp = Path(tf.name)
            try:
                dmrs.write_parquet(str(tmp))
                scored = score_dmr_parquet(
                    tmp, truth,
                    tool=tool, scenario=scenario,
                    parameter=key, parameter_value=params[key],
                    test=test,
                )
            finally:
                tmp.unlink(missing_ok=True)

            # score_dmr_parquet emits one row per threshold; reduce to
            # the q<0.05 row + record the perturbation cell.
            scored_q = [r for r in scored
                        if r.get("threshold_kind") == "qvalue"]
            if not scored_q:
                rows.append({
                    "parameter": key, "factor": factor,
                    "parameter_value": params[key],
                    "n_dmrs": dmrs.height,
                    "error": "no qvalue row in scored output",
                })
                continue
            r = dict(scored_q[0])
            r.update({
                "parameter": key, "factor": factor,
                "parameter_value": params[key],
                "n_dmrs": dmrs.height,
                "preset": preset_name,
            })
            rows.append(r)
            logger.info(
                "sweep cell: %s=%s (factor=%s) -> n_dmrs=%d TPR=%.4f FPR=%.4f F1=%.4f",
                key, params[key], factor, dmrs.height,
                r.get("tpr", float("nan")),
                r.get("fpr", float("nan")),
                r.get("f1", float("nan")),
            )

    return pl.DataFrame(rows)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dmc-parquet", required=True, type=Path)
    parser.add_argument("--truth-parquet", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--tool", default="epykit",
                        help="Tool label for the scored rows")
    parser.add_argument("--scenario", default="default",
                        help="Scenario label for the scored rows")
    parser.add_argument("--test", default="lr",
                        help="DMC engine label for the scored rows")
    parser.add_argument("--preset", default="default",
                        choices=("strict", "default", "permissive"))
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    df = run_sweep(
        dmc_parquet=args.dmc_parquet,
        truth_parquet=args.truth_parquet,
        tool=args.tool,
        scenario=args.scenario,
        test=args.test,
        preset_name=args.preset,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.write_parquet(str(args.out))
    print(f"wrote {args.out} ({df.height} rows)")


if __name__ == "__main__":
    main()
