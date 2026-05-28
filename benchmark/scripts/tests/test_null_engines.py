"""Test the real-engine closures: each surviving DMC engine wraps
into the engine_fn(samples_treatment, samples_control, seed)->qvals
contract."""
from __future__ import annotations

import numpy as np
import pytest

from _null_engines import ENGINE_REGISTRY

SURVIVING_DMC = ["lr", "lr_plus", "welch_t", "fisher", "glm"]
# DMR engines are skipped here (they need a prior DMC run).
# The glm closure builds full + reduced design matrices on each shuffle via
# ep.tl.dmc(formula="~ group", contrast="group"), so no extra obs columns are
# required beyond the "group" column that the synth_md fixture already writes.


@pytest.mark.slow
@pytest.mark.parametrize("engine_name", SURVIVING_DMC)
def test_dmc_engine_closure_runs(synth_md_filtered, engine_name):
    """Each registered DMC engine: callable, returns 1D q-value array
    in [0, 1] (or NaN), deterministic across two seeded runs."""
    md = synth_md_filtered
    factory = ENGINE_REGISTRY[engine_name]
    closure = factory(md)
    treat = list(md.treatment_ids)
    ctrl = list(md.control_ids)
    out_a = closure(samples_treatment=treat, samples_control=ctrl, seed=42)
    out_b = closure(samples_treatment=treat, samples_control=ctrl, seed=42)
    assert isinstance(out_a, np.ndarray)
    assert out_a.ndim == 1
    np.testing.assert_array_equal(out_a, out_b, err_msg=f"{engine_name}: not deterministic")
    finite = out_a[np.isfinite(out_a)]
    if len(finite) > 0:
        assert ((finite >= 0) & (finite <= 1)).all(), (
            f"q-values out of [0, 1] for {engine_name}"
        )
