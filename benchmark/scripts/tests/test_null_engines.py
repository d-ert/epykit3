"""Test the real-engine closures: each surviving DMC engine wraps
into the engine_fn(samples_treatment, samples_control, seed)->qvals
contract."""
from __future__ import annotations

import numpy as np
import pytest

from _null_engines import ENGINE_REGISTRY, _permute_md

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


def test_permute_md_updates_treatment_column(synth_md):
    """Full swap of treatment <-> control must update obs.treatment (i64),
    not just obs.group. md.treatment_ids reads from obs.treatment, so
    failing to update it leaves the DMC engine running on the original
    label assignment -- the bug this regression test guards against.
    """
    md = synth_md
    orig_treat = list(md.treatment_ids)
    orig_ctrl = list(md.control_ids)
    assert orig_treat and orig_ctrl, "fixture sanity: both groups non-empty"

    # Full swap: original controls become treatment, and vice versa.
    md_copy = _permute_md(md, samples_treatment=orig_ctrl, samples_control=orig_treat)

    assert sorted(md_copy.treatment_ids) == sorted(orig_ctrl), (
        "treatment_ids should reflect the swapped labels"
    )
    assert sorted(md_copy.control_ids) == sorted(orig_treat), (
        "control_ids should reflect the swapped labels"
    )
    # obs.treatment must remain integer-valued (the property filters on == 1).
    treat_vals = set(md_copy.obs["treatment"].to_list())
    assert treat_vals == {0, 1}, (
        f"obs.treatment must be integer 0/1, got values {treat_vals}"
    )
    # obs.group should also be swapped (display consistency regression).
    group_by_sid = dict(
        zip(md_copy.obs["sample_id"].to_list(), md_copy.obs["group"].to_list())
    )
    for sid in orig_ctrl:
        assert group_by_sid[sid] == "treatment", (
            f"sample {sid} moved to treatment; obs.group must reflect that"
        )
    for sid in orig_treat:
        assert group_by_sid[sid] == "control", (
            f"sample {sid} moved to control; obs.group must reflect that"
        )
    # Original md must not be mutated.
    assert sorted(md.treatment_ids) == sorted(orig_treat)
    assert sorted(md.control_ids) == sorted(orig_ctrl)


def test_permute_md_partial_relabel(synth_md):
    """Mixed assignment: some original treatments stay treatment, some
    move to control; same for original controls. The treatment_ids
    property must reflect exactly the requested ``samples_treatment``
    list, regardless of original group. All samples must be covered by
    one of the two input lists (the contract used by the null
    calibration sweep).
    """
    md = synth_md
    orig_treat = list(md.treatment_ids)
    orig_ctrl = list(md.control_ids)
    assert len(orig_treat) >= 3 and len(orig_ctrl) >= 3, (
        "this test wants at least 3 samples per group"
    )

    # Build a mixed assignment that still covers every sample exactly once
    # (the in-practice contract). Use the first 3 originals from each group
    # to form a 3-vs-3 mix, and route any remaining samples deterministically.
    treat_first3 = orig_treat[:3]
    ctrl_first3 = orig_ctrl[:3]
    new_treat = [treat_first3[0], ctrl_first3[0], ctrl_first3[1]]
    new_ctrl = [treat_first3[1], treat_first3[2], ctrl_first3[2]]
    # Spill any extra samples into new_ctrl so every sample_id appears
    # exactly once across the two lists.
    leftovers = [s for s in orig_treat[3:] + orig_ctrl[3:]]
    new_ctrl = new_ctrl + leftovers

    md_copy = _permute_md(md, samples_treatment=new_treat, samples_control=new_ctrl)

    assert sorted(md_copy.treatment_ids) == sorted(new_treat)
    assert sorted(md_copy.control_ids) == sorted(new_ctrl)


@pytest.mark.slow
def test_closure_produces_different_qvalues_per_shuffle(synth_md_filtered):
    """Regression for the original bug: when the closure is called with
    two clearly different label assignments, the returned q-value arrays
    must differ. Before the fix, _permute_md only touched obs.group while
    md.treatment_ids reads obs.treatment, so every shuffle silently ran
    on the original assignment and returned the same q-values.

    A full label swap (treatment <-> control) is two-sidedly symmetric
    for the LR test and yields bit-identical q-values even after the
    fix, so we use an ASYMMETRIC mixed permutation instead: move one
    sample across the group boundary, leave the rest alone. The
    resulting test compares 3-vs-5 to the original 4-vs-4.

    ``ep.tl.dmc`` has a weak cache hit in ``src/epykit/dmc.py`` that can
    short-circuit recomputation when the per-chrom parquets already
    exist on disk; we wipe the DMC cache directory between the two
    closure calls so this test isolates the closure behaviour from the
    cache layer. The bug guarded against is independent of the cache:
    before the fix, the two calls produced bit-identical q-values even
    on the first invocation against a freshly-built methylstore (the
    real Phase 4 Task 7 symptom).
    """
    import shutil
    from pathlib import Path

    md = synth_md_filtered
    factory = ENGINE_REGISTRY["lr"]
    closure = factory(md)
    orig_treat = list(md.treatment_ids)
    orig_ctrl = list(md.control_ids)
    assert len(orig_treat) >= 2 and len(orig_ctrl) >= 2, (
        "fixture sanity: need at least 2 samples per group"
    )

    dmc_cache_dir = Path(md.store).parent / "dmc"

    # Assignment A: the original 4-vs-4.
    qvals_a = closure(samples_treatment=orig_treat, samples_control=orig_ctrl, seed=1)

    # Wipe the cache to force recompute under assignment B. This isolates
    # the regression -- we are testing the closure / permutation path,
    # not the DMC cache layer.
    if dmc_cache_dir.exists():
        shutil.rmtree(dmc_cache_dir)

    # Assignment B: move the last original treatment into control. Now
    # 3-vs-5 -- not a symmetric swap of A.
    mixed_treat = orig_treat[:-1]
    mixed_ctrl = orig_ctrl + [orig_treat[-1]]
    qvals_b = closure(samples_treatment=mixed_treat, samples_control=mixed_ctrl, seed=1)

    assert qvals_a.size > 0 and qvals_b.size > 0, "expected non-empty q-value arrays"
    assert not np.array_equal(qvals_a, qvals_b), (
        "q-values are bit-identical across two different label assignments; "
        "the closure is not actually permuting the labels."
    )
