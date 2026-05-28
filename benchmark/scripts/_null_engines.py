"""Real-engine closures for run_null_calibration.py.

Each registry entry is a factory: ``factory(md) -> closure``, where
the closure has signature:
    closure(samples_treatment, samples_control, seed) -> np.ndarray of q-values

Factories capture the MethylData object once outside the shuffle loop
so the parquet store is loaded only once per (engine, scenario).

Usage in run_null_calibration.py main():
    from _null_engines import ENGINE_REGISTRY
    factory = ENGINE_REGISTRY["lr"]
    closure = factory(md)
    qvals = closure(samples_treatment=treat, samples_control=ctrl, seed=42)
"""
from __future__ import annotations

import copy
from typing import Callable

import numpy as np
import polars as pl

import epykit as ep


def _permute_md(md, samples_treatment: list[str], samples_control: list[str]):
    """Return a copy of md with group labels reassigned per the shuffle."""
    md_copy = copy.copy(md)
    label_map = {s: "treatment" for s in samples_treatment}
    label_map.update({s: "control" for s in samples_control})
    new_groups = [
        label_map.get(sid, grp)
        for sid, grp in zip(md.obs["sample_id"].to_list(), md.obs["group"].to_list())
    ]
    md_copy.obs = md.obs.with_columns(pl.Series("group", new_groups))
    return md_copy


def _dmc_engine(test_name: str, *, lr_plus: bool = False, glm: bool = False) -> Callable:
    """Factory for a DMC engine closure.

    ``glm=True`` routes through ``ep.tl.dmc(..., test="glm", formula="~ group",
    contrast="group")``. This is the contrast / multi-group path in
    :func:`epykit.tl.dmc`, which internally calls :func:`epykit._glm.build_design`
    to construct the full design (``~ group``) and the reduced design
    (intercept-only, equivalent to dropping the ``group`` term). The
    permuted ``group`` column on ``md_perm.obs`` is the regressor; under a
    calibrated null the resulting q-values should be approximately
    uniformly distributed on [0, 1].
    """
    def factory(md):
        def closure(samples_treatment, samples_control, seed):
            md_perm = _permute_md(md, samples_treatment, samples_control)
            if lr_plus:
                ep.tl.dmc(
                    md_perm, test="lr",
                    fdr_method="fdr_tsbh",
                    neighbour_combine=True,
                    sep_fallback=True,
                    dispersion="eb",
                )
            elif glm:
                # Full design ~ group (intercept + group[T.treatment]); reduced
                # design ~ 1 (intercept only). build_design is called inside
                # _run_dmc_contrast; we only need to pass the formula + contrast.
                ep.tl.dmc(
                    md_perm, test="glm",
                    formula="~ group",
                    contrast="group",
                )
            else:
                ep.tl.dmc(md_perm, test=test_name)
            df = md_perm.dmc
            if df is None or df.height == 0:
                return np.array([], dtype=np.float64)
            # Use qvalue if available, else pvalue.
            qcol = "qvalue" if "qvalue" in df.columns else "pvalue"
            return df[qcol].to_numpy().astype(np.float64)
        return closure
    return factory


ENGINE_REGISTRY: dict[str, Callable] = {
    # DMC engines.
    "lr":      _dmc_engine("lr"),
    "lr_plus": _dmc_engine("lr", lr_plus=True),
    "welch_t": _dmc_engine("welch_t"),
    "fisher":  _dmc_engine("fisher"),
    # glm: builds design matrices on every shuffle via ep.tl.dmc(formula=...)
    # which routes through _run_dmc_contrast -> _glm.build_design. The full
    # model is ~ group, the reduced model is the intercept-only ~ 1.
    "glm":     _dmc_engine("glm", glm=True),
}
