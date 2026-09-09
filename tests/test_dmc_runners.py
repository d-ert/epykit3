"""The per-engine chromosome runners in ``epykit.dmc`` and the finalizer they feed.

One test per contract: the runner table names exactly the registry's
engines, a runner returns only per-site arrays (never a per-sample stack),
an empty chromosome short-circuits to the canonical schema, and an unknown
engine is refused before any runner runs.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import fields, is_dataclass
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from epykit._dmc_engines import ENGINES
from epykit._glm import build_design, resolve_contrast
from epykit.dmc import (
    _EMPTY_SCHEMA,
    _ENGINE_RUNNERS,
    EngineInput,
    _intersect_chrom,
    _process_one_chromosome,
)


def test_runner_table_matches_the_registry():
    assert set(_ENGINE_RUNNERS) == set(ENGINES)


def _engine_input(md, test: str, joint: bool = False) -> EngineInput:
    """One chromosome of the session fixture, with the design the GLM
    engines need. ``joint=True`` tests both design columns at once, which
    is the multi-group (k > 1) branch of ``glm_contrast``."""
    samples = md.treatment_ids + md.control_ids
    store = Path(md.store)
    design_full = design_reduced = contrast_matrix = None
    coef_idx = samples_all_ordered = group_labels = None
    if test in ("glm", "glm_contrast"):
        design_full, design_reduced, coef_idx, terms, _formula, info = build_design(
            md.obs, samples_ordered=samples, treatment_col="treatment", return_design_info=True
        )
    if test == "glm_contrast":
        samples_all_ordered = samples
        if joint:
            contrast_matrix = np.eye(design_full.shape[1])
            label_of = dict(zip(md.obs["sample_id"], md.obs["group"], strict=True))
            group_labels = [label_of[s] for s in samples]
        else:
            contrast_matrix, _label = resolve_contrast("treatment", terms, design_info=info)
    return EngineInput(
        methylstore_path=store,
        chrom="chr1",
        canonical_df=_intersect_chrom(store, "chr1", samples),
        samples_case=md.treatment_ids,
        samples_control=md.control_ids,
        test=test,
        min_samples_case=0,
        min_samples_control=0,
        dispersion="site",
        reference="adaptive",
        design_full=design_full,
        design_reduced=design_reduced,
        coef_idx=coef_idx,
        contrast_matrix=contrast_matrix,
        samples_all_ordered=samples_all_ordered,
        group_labels_per_sample=group_labels,
        glm_backend="cpu",
        smoothing=False,
        smoothing_span_bp=500,
        sep_fallback=False,
        sep_threshold=0.9,
    )


def _arrays(value, path: str = "") -> Iterator[tuple[str, np.ndarray]]:
    """Every array reachable from a result record, with its field path."""
    if isinstance(value, np.ndarray):
        yield path, value
    elif is_dataclass(value):
        for f in fields(value):
            yield from _arrays(getattr(value, f.name), f"{path}.{f.name}")
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from _arrays(item, f"{path}[{key!r}]")
    elif isinstance(value, tuple):
        for i, item in enumerate(value):
            yield from _arrays(item, f"{path}[{i}]")


@pytest.mark.parametrize(
    ("test", "joint"),
    [
        ("fisher", False),
        ("lr", False),
        ("welch_t", False),
        ("glm", False),
        ("glm_contrast", False),
        ("glm_contrast", True),
    ],
    ids=["fisher", "lr", "welch_t", "glm", "glm_contrast", "glm_contrast-joint"],
)
def test_runner_returns_only_per_site_arrays(synth_md_filtered, test, joint):
    """The streaming contract: every array a runner hands back has one entry
    per canonical site. Sample stacks and accumulators end with the runner."""
    inp = _engine_input(synth_md_filtered, test, joint=joint)
    res = _ENGINE_RUNNERS[test](inp)
    seen = dict(_arrays(res, "res"))
    assert seen, "runner returned no arrays"
    for path, arr in seen.items():
        assert arr.shape == (inp.n_sites,), f"{path} has shape {arr.shape}"
    assert (res.multigroup is not None) is joint
    if joint:
        assert set(res.multigroup.level_mean_beta) == {"treatment", "control"}


def test_empty_chromosome_returns_the_canonical_schema(tmp_path):
    empty = pl.DataFrame(
        {"pos": pl.Series([], dtype=pl.Int32), "strand": pl.Series([], dtype=pl.Utf8)}
    )
    out = _process_one_chromosome(tmp_path, "chr1", empty, ["a"], ["b"], "lr")
    assert out.height == 0
    assert dict(out.schema) == _EMPTY_SCHEMA


def test_unknown_engine_is_refused_before_any_runner(synth_md_filtered):
    md = synth_md_filtered
    inp = _engine_input(md, "lr")
    with pytest.raises(NotImplementedError, match="Test 'bogus' not implemented"):
        _process_one_chromosome(
            inp.methylstore_path,
            "chr1",
            inp.canonical_df,
            md.treatment_ids,
            md.control_ids,
            "bogus",
        )
