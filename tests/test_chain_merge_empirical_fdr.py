"""Permutation empirical FDR for chain_merge DMRs.

Two layers:

* ``empirical_fdr_for_chain_merge`` with the per-permutation engine
  ``_chain_merge_perm_survivors`` monkeypatched, so the harness and the
  count-ratio aggregation are exercised with hand-checkable expectations
  (mirrors ``test_dmr_region_fdr_mode.py`` for the tile path).
* The replay contract on the synth bundle: ``tl.dmr`` hands the observed DMC
  knobs, chromosome universe, multiple-testing correction and region filters
  to every permutation; rejections happen before the first permutation; the
  observed ``DMCStore`` is never touched; the seed makes the slow end-to-end
  run reproducible.
"""

from __future__ import annotations

import copy
import hashlib
from pathlib import Path

import numpy as np
import polars as pl
import pytest

import epykit as ep
from epykit import dmc as dmc_mod
from epykit import dmr as dmr_mod
from epykit._dmc_store import DMCStore

_SCHEMA = {
    "chrom": pl.Utf8,
    "start": pl.Int32,
    "end": pl.Int32,
    "n_cpgs": pl.Int32,
    "mean_meth_diff": pl.Float64,
    "combined_pvalue": pl.Float64,
}

_TREAT = ["t1", "t2", "t3"]
_CTRL = ["c1", "c2", "c3"]


def _observed(n=40, p=1e-6):
    return pl.DataFrame(
        {
            "chrom": ["chr1"] * n,
            "start": list(range(0, n * 1000, 1000)),
            "end": list(range(500, n * 1000 + 500, 1000)),
            "n_cpgs": [8] * n,
            "mean_meth_diff": [0.3] * n,
            "combined_pvalue": [p] * n,
        },
        schema=_SCHEMA,
    )


def _run(observed, **kwargs):
    base = dict(
        methylstore_path="/dev/null",
        samples_treatment=_TREAT,
        samples_control=_CTRL,
        observed_dmr=observed,
        chromosomes=["chr1"],
        dmc_kwargs={"test": "lr"},
        chain_merge_kwargs={"preset": "default"},
        min_mean_qvalue=0.05,
        n_perm=20,
        seed=0,
        n_jobs=1,
    )
    base.update(kwargs)
    return dmr_mod.empirical_fdr_for_chain_merge(**base)


# harness with the per-permutation engine stubbed


def test_chain_merge_region_mode_significant(monkeypatch):
    """40 observed regions at combined_pvalue=1e-6; each usable shuffle yields
    2 extreme null regions -> set FDR = 2/40 = 0.05, gradient <= that."""
    monkeypatch.setattr(
        dmr_mod,
        "_chain_merge_perm_survivors",
        lambda **kw: np.array([1e-6, 1e-6], dtype=np.float64),
    )
    out = _run(_observed(40, 1e-6), fdr_method="region")
    q = out.get_column("empirical_qvalue").to_numpy()
    assert "empirical_pvalue" in out.columns
    assert out.get_column("empirical_fdr_set")[0] == pytest.approx(0.05, abs=1e-6)
    assert np.all(q <= 0.05 + 1e-9)


def test_chain_merge_max_t_is_the_default(monkeypatch):
    monkeypatch.setattr(
        dmr_mod,
        "_chain_merge_perm_survivors",
        lambda **kw: np.array([1e-6, 1e-6], dtype=np.float64),
    )
    implicit = _run(_observed(40, 1e-6))
    explicit = _run(_observed(40, 1e-6), fdr_method="max_t")
    assert implicit.equals(explicit)
    # Every shuffle carries a null as extreme as the observed regions: the
    # min-P bar saturates and the set-level column is NaN in this mode.
    assert np.all(implicit.get_column("empirical_pvalue").to_numpy() > 0.5)
    assert np.all(np.isnan(implicit.get_column("empirical_fdr_set").to_numpy()))


def test_chain_merge_uses_combined_pvalue_not_pvalue(monkeypatch):
    """The harness reads the region statistic from ``combined_pvalue``."""
    monkeypatch.setattr(
        dmr_mod,
        "_chain_merge_perm_survivors",
        lambda **kw: np.array([0.5, 0.5], dtype=np.float64),
    )
    out = _run(_observed(10, 1e-8), fdr_method="region")
    assert np.nanmax(out.get_column("empirical_qvalue").to_numpy()) == 0.0


def test_chain_merge_small_n_warns_in_region_mode(monkeypatch):
    monkeypatch.setattr(
        dmr_mod, "_chain_merge_perm_survivors", lambda **kw: np.array([1e-6], dtype=np.float64)
    )
    with pytest.warns(UserWarning, match="underpowered"):
        _run(_observed(10, 1e-6), n_perm=5, fdr_method="region")


def test_chain_merge_empty_observed_returns_empty_cols(monkeypatch):
    def fake(**kw):
        raise AssertionError("no permutation must run on an empty observed table")

    monkeypatch.setattr(dmr_mod, "_chain_merge_perm_survivors", fake)
    out = _run(pl.DataFrame(schema=_SCHEMA), n_perm=5)
    assert out.height == 0
    for col in ("empirical_pvalue", "empirical_qvalue", "empirical_fdr_set"):
        assert out.schema[col] == pl.Float64


@pytest.mark.parametrize(
    ("kwargs", "exc", "pattern"),
    [
        ({"fdr_method": "bogus"}, ValueError, r"fdr_method"),
        ({"n_perm": 0}, ValueError, r"n_perm"),
        ({"chromosomes": []}, ValueError, r"chromosomes"),
        ({"dmc_kwargs": {"test": "glm"}}, NotImplementedError, r"two-group"),
        ({"dmc_kwargs": {"test": "lr", "design_full": 1}}, ValueError, r"dmc_kwargs"),
        ({"dmc_kwargs": {"test": "lr", "chromosomes": ["chr1"]}}, ValueError, r"dmc_kwargs"),
        ({"dmc_fdr_method": "region"}, ValueError, r"dmc_fdr_method"),
    ],
)
def test_chain_merge_rejects_before_any_permutation(monkeypatch, kwargs, exc, pattern):
    def fake(**kw):
        raise AssertionError("no permutation must run")

    monkeypatch.setattr(dmr_mod, "_chain_merge_perm_survivors", fake)
    with pytest.raises(exc, match=pattern):
        _run(_observed(5, 1e-6), **kwargs)


# per-permutation engine: replay contract


def _fake_engine(monkeypatch, *, chain_result, dmc_error=None):
    calls: dict[str, object] = {}

    def fake_process(**kw):
        calls["dmc"] = kw
        if dmc_error is not None:
            raise dmc_error
        return "STORE"

    def fake_correct(store, method="fdr_bh"):
        calls["correct"] = (store, method)
        return "CORRECTED"

    def fake_chain(store, **kw):
        calls["chain"] = (store, kw)
        return chain_result

    monkeypatch.setattr(dmc_mod, "process_chromosomes_dmc", fake_process)
    monkeypatch.setattr(dmc_mod, "apply_multiple_testing_correction", fake_correct)
    monkeypatch.setattr(dmr_mod, "call_dmr_chain_merge", fake_chain)
    return calls


def test_perm_engine_replays_knobs_correction_and_filter(monkeypatch):
    calls = _fake_engine(
        monkeypatch,
        chain_result=pl.DataFrame(
            {"combined_pvalue": [1e-4, 0.3, 2e-3], "combined_qvalue": [3e-4, 0.3, 4e-2]}
        ),
    )
    out = dmr_mod._chain_merge_perm_survivors(
        methylstore_path="/store",
        samples_treatment=["a", "d"],
        samples_control=["b", "c"],
        chromosomes=["chr2", "chr1"],
        dmc_kwargs={"test": "welch_t", "dispersion": "eb", "unite": False},
        dmc_fdr_method="fdr_tsbh",
        chain_merge_kwargs={"preset": "strict", "min_cpgs": 7},
        min_mean_qvalue=0.05,
    )
    dmc_call = calls["dmc"]
    assert dmc_call["methylstore_path"] == "/store"
    assert dmc_call["samples_treatment"] == ["a", "d"]
    assert dmc_call["samples_control"] == ["b", "c"]
    assert dmc_call["chromosomes"] == ["chr2", "chr1"]
    assert dmc_call["return_store"] is True
    assert {"test": "welch_t", "dispersion": "eb", "unite": False}.items() <= dmc_call.items()
    perm_dir = Path(dmc_call["out_dir"])
    assert perm_dir.name.startswith("epykit_dmr_perm_")
    assert not perm_dir.exists(), "the per-permutation store must be removed"
    assert calls["correct"] == ("STORE", "fdr_tsbh")
    assert calls["chain"] == ("CORRECTED", {"preset": "strict", "min_cpgs": 7})
    # The observed region filter (combined_qvalue < 0.05) drops the 0.3 row.
    np.testing.assert_allclose(out, [1e-4, 2e-3])


def test_perm_engine_failure_and_empty_are_distinct(monkeypatch):
    _fake_engine(monkeypatch, chain_result=None, dmc_error=RuntimeError("boom"))
    kwargs = dict(
        methylstore_path="/store",
        samples_treatment=["a"],
        samples_control=["b"],
        chromosomes=["chr1"],
        dmc_kwargs={"test": "lr"},
        dmc_fdr_method="fdr_bh",
        chain_merge_kwargs={},
        min_mean_qvalue=None,
    )
    assert dmr_mod._chain_merge_perm_survivors(**kwargs) is None
    _fake_engine(monkeypatch, chain_result=pl.DataFrame(schema=dmr_mod._DMR_EMPTY_SCHEMA))
    empty = dmr_mod._chain_merge_perm_survivors(**kwargs)
    assert empty is not None and empty.size == 0


# tl.dmr wiring on the synth bundle


@pytest.fixture(scope="module")
def dmc_md(synth_bundle, tmp_path_factory):
    """One lr DMC on the synth bundle for every test in this module."""
    store_dir = tmp_path_factory.mktemp("chain_merge_fdr")
    md = ep.read_bismark(
        synth_bundle.samplesheet,
        treatment_group="treatment",
        control_group="control",
        assembly="synth",
        store_dir=str(store_dir / "store"),
    )
    ep.pp.filter_coverage(md, lo_count=5, hi_perc=99.9)
    ep.pp.set_unite_type(md, type="intersect")
    ep.tl.dmc(md, test="lr")
    return md


@pytest.fixture
def md_lr(dmc_md):
    """Per-test view of the module DMC: private ``uns`` / ``obs`` / ``varm``."""
    md = copy.copy(dmc_md)
    md.uns = copy.deepcopy(dmc_md.uns)
    md.obs = dmc_md.obs.clone()
    md.varm = dict(dmc_md.varm)
    return md


def _capturing_engine(monkeypatch):
    captured: list[dict] = []

    def fake(**kw):
        captured.append(kw)
        return np.array([0.5], dtype=np.float64)

    monkeypatch.setattr(dmr_mod, "_chain_merge_perm_survivors", fake)
    return captured


def _refusing_engine(monkeypatch):
    def fake(**kw):
        raise AssertionError("no permutation must run")

    monkeypatch.setattr(dmr_mod, "_chain_merge_perm_survivors", fake)


def _store_fingerprint(store_path: str) -> dict[str, str]:
    root = Path(store_path)
    files = [root / ".epykit_dmc_manifest.json", *sorted(root.glob("chrom=*.parquet"))]
    return {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in files}


def test_tl_dmr_forwards_observed_knobs_universe_correction_and_filters(monkeypatch, md_lr):
    captured = _capturing_engine(monkeypatch)
    ep.tl.dmr(
        md_lr,
        method="chain_merge",
        preset="default",
        min_mean_qvalue=0.01,
        empirical_fdr=True,
        n_perm=3,
        perm_seed=0,
        fdr_method="region",
    )
    assert len(captured) == 3
    meta = md_lr.uns["dmc"]
    observed_store = DMCStore.open(meta["store_path"])
    for kw in captured:
        assert kw["methylstore_path"] == md_lr.store
        assert kw["chromosomes"] == observed_store.chroms()
        assert kw["dmc_kwargs"] == {
            "test": "lr",
            "unite": True,
            "min_samples_treatment": 0,
            "min_samples_control": 0,
            "dispersion": meta["dispersion"],
            "reference": meta["reference"],
            "smoothing": False,
            "sep_fallback": False,
            "sep_threshold": 0.9,
        }
        assert kw["dmc_fdr_method"] == meta["fdr_method"] == "fdr_bh"
        assert kw["chain_merge_kwargs"]["preset"] == "default"
        assert kw["chain_merge_kwargs"]["min_cpgs"] == md_lr.uns["dmr_params"]["min_cpgs"]
        assert kw["min_mean_qvalue"] == 0.01
    # Every permutation draws a genuine shuffle of the observed split.
    for kw in captured:
        assert sorted(kw["samples_treatment"] + kw["samples_control"]) == sorted(
            md_lr.treatment_ids + md_lr.control_ids
        )
        assert len(kw["samples_treatment"]) == len(md_lr.treatment_ids)
    params = md_lr.uns["dmr_params"]
    assert params["empirical_fdr"] is True
    assert params["n_perm"] == 3 and params["perm_seed"] == 0
    assert params["fdr_method"] == "region"
    assert isinstance(params["empirical_fdr_set"], float)
    dmr = md_lr.uns["dmr"]
    assert dmr.height > 0
    assert {"empirical_pvalue", "empirical_qvalue", "empirical_fdr_set"} <= set(dmr.columns)


def test_tl_dmr_records_no_permutation_keys_without_empirical_fdr(md_lr):
    ep.tl.dmr(md_lr, method="chain_merge")
    params = md_lr.uns["dmr_params"]
    assert params["empirical_fdr"] is False
    assert params["n_perm"] is None and params["fdr_method"] is None
    assert params["empirical_fdr_set"] is None
    assert "empirical_qvalue" not in md_lr.uns["dmr"].columns


def test_tl_dmr_matching_chromosome_restriction_is_accepted(monkeypatch, md_lr):
    captured = _capturing_engine(monkeypatch)
    universe = DMCStore.open(md_lr.uns["dmc"]["store_path"]).chroms()
    ep.tl.dmr(
        md_lr,
        method="chain_merge",
        chromosomes=list(reversed(universe)),
        empirical_fdr=True,
        n_perm=2,
        perm_seed=0,
    )
    assert len(captured) == 2
    assert captured[0]["chromosomes"] == universe


def test_tl_dmr_strata_column_restricts_shuffles(monkeypatch, md_lr):
    captured = _capturing_engine(monkeypatch)
    pairs = {s: f"pair{i}" for i, s in enumerate(md_lr.treatment_ids)}
    pairs.update({s: f"pair{i}" for i, s in enumerate(md_lr.control_ids)})
    md_lr.obs = md_lr.obs.with_columns(
        pl.col("sample_id").replace_strict(pairs, default=None).alias("pair")
    )
    ep.tl.dmr(
        md_lr,
        method="chain_merge",
        empirical_fdr=True,
        n_perm=4,
        perm_seed=1,
        empirical_strata="pair",
    )
    assert len(captured) == 4
    for kw in captured:
        # One sample of every pair lands in each group.
        assert sorted(pairs[s] for s in kw["samples_treatment"]) == sorted(set(pairs.values()))


@pytest.mark.parametrize(
    ("mutate", "exc", "pattern"),
    [
        (lambda md: md.uns["dmc"].__setitem__("test_used", "glm"), NotImplementedError, r"glm"),
        (
            lambda md: md.uns["dmc"].__setitem__("formula", "~ group + age"),
            NotImplementedError,
            r"formula",
        ),
        (
            lambda md: md.uns["dmc"].__setitem__("use_smoothed", True),
            NotImplementedError,
            r"use_smoothed",
        ),
        (lambda md: md.uns["dmc"].pop("dispersion"), ValueError, r"dispersion"),
        (lambda md: md.uns.pop("dmc"), ValueError, r"Run ep.tl.dmc"),
    ],
)
def test_tl_dmr_rejects_unsupported_dmc_before_any_permutation(
    monkeypatch, md_lr, mutate, exc, pattern
):
    _refusing_engine(monkeypatch)
    mutate(md_lr)
    with pytest.raises(exc, match=pattern):
        ep.tl.dmr(md_lr, method="chain_merge", empirical_fdr=True, n_perm=3)


def test_tl_dmr_rejects_missing_or_partial_strata_before_any_permutation(monkeypatch, md_lr):
    _refusing_engine(monkeypatch)
    with pytest.raises(ValueError, match=r"empirical_strata='nope' is not a column"):
        ep.tl.dmr(
            md_lr, method="chain_merge", empirical_fdr=True, n_perm=3, empirical_strata="nope"
        )
    first = md_lr.treatment_ids[0]
    md_lr.obs = md_lr.obs.with_columns(
        pl.when(pl.col("sample_id") == first).then(None).otherwise(pl.lit("s")).alias("stratum")
    )
    with pytest.raises(ValueError, match=r"does not cover"):
        ep.tl.dmr(
            md_lr, method="chain_merge", empirical_fdr=True, n_perm=3, empirical_strata="stratum"
        )


def test_tl_dmr_rejects_chromosome_restriction_outside_universe(monkeypatch, md_lr):
    _refusing_engine(monkeypatch)
    with pytest.raises(ValueError, match=r"rerun ep.tl.dmc"):
        ep.tl.dmr(md_lr, method="chain_merge", chromosomes=["chr1"], empirical_fdr=True, n_perm=3)


def test_tl_dmr_tile_missing_strata_raises_before_any_permutation(monkeypatch, md_lr):
    """The strict strata check applies to the tile harness as well."""

    def fake(**kw):
        raise AssertionError("no permutation must run")

    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based", fake)
    with pytest.raises(ValueError, match=r"empirical_strata='nope' is not a column"):
        ep.tl.dmr(md_lr, method="tile", empirical_fdr=True, n_perm=3, empirical_strata="nope")


@pytest.mark.slow
def test_chain_merge_empirical_fdr_end_to_end_is_reproducible(dmc_md):
    """Real permutations on the synth bundle: each shuffle recomputes the DMC,
    chain-merges and filters. The same seed yields the same table."""

    def run(seed):
        md = copy.copy(dmc_md)
        md.uns = copy.deepcopy(dmc_md.uns)
        md.obs = dmc_md.obs.clone()
        md.varm = dict(dmc_md.varm)
        ep.tl.dmr(
            md,
            method="chain_merge",
            empirical_fdr=True,
            n_perm=3,
            perm_seed=seed,
            fdr_method="region",
        )
        return md

    first = run(0)
    second = run(0)
    params = first.uns["dmr_params"]
    assert params["method"] == "chain_merge"
    assert params["empirical_fdr"] is True and params["fdr_method"] == "region"
    assert params["n_perm"] == 3 and params["perm_seed"] == 0
    dmr = first.uns["dmr"]
    assert dmr.height > 0, "the synth bundle seeds DMRs the default caller recovers"
    cols = ["empirical_pvalue", "empirical_qvalue", "empirical_fdr_set"]
    assert dmr.select(cols).equals(second.uns["dmr"].select(cols))
    assert params["empirical_fdr_set"] == second.uns["dmr_params"]["empirical_fdr_set"]
    assert 0.0 <= params["empirical_fdr_set"] <= 1.0
    q = dmr.get_column("empirical_qvalue").to_numpy()
    assert np.all((q >= 0) & (q <= 1))


@pytest.mark.slow
def test_observed_store_unchanged_after_serial_and_parallel_permutations(dmc_md):
    """Per-permutation DMC stores live in private temporary directories, so
    the observed manifest and parquet files are byte-identical afterwards,
    whether the permutations run serially or through joblib workers."""
    store_path = dmc_md.uns["dmc"]["store_path"]
    before = _store_fingerprint(store_path)
    chroms_before = DMCStore.open(store_path).chroms()

    def run(n_jobs):
        md = copy.copy(dmc_md)
        md.uns = copy.deepcopy(dmc_md.uns)
        md.obs = dmc_md.obs.clone()
        md.varm = dict(dmc_md.varm)
        ep.tl.dmr(
            md,
            method="chain_merge",
            empirical_fdr=True,
            n_perm=2,
            perm_seed=3,
            perm_n_jobs=n_jobs,
            fdr_method="region",
        )
        return md.uns["dmr"]

    serial = run(1)
    assert _store_fingerprint(store_path) == before
    parallel = run(2)
    assert _store_fingerprint(store_path) == before
    assert DMCStore.open(store_path).chroms() == chroms_before
    cols = ["empirical_pvalue", "empirical_qvalue", "empirical_fdr_set"]
    assert serial.select(cols).equals(parallel.select(cols))
