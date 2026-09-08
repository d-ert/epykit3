"""Integration tests for the ``fdr_method`` dispatch of ``empirical_fdr_for_dmr``.

Uses the monkeypatch pattern (the real ``call_dmr_tile_based`` is too heavy):
each "permutation" returns a scripted null survivor set, so the harness and
the estimator math are exercised end-to-end with hand-checkable expectations.
``max_t`` is the default and must reproduce the pre-``fdr_method`` numbers;
``region`` is the opt-in count-ratio construction.
"""

from __future__ import annotations

import warnings

import numpy as np
import polars as pl
import pytest

from epykit import dmr as dmr_mod

_SCHEMA = {
    "chrom": pl.Utf8,
    "start": pl.Int32,
    "end": pl.Int32,
    "n_cpgs": pl.Int32,
    "meth_diff": pl.Float64,
    "pvalue": pl.Float64,
}

_TREAT = ["t1", "t2", "t3"]
_CTRL = ["c1", "c2", "c3"]


def _observed(n=50, p=1e-6, pvalues=None):
    pvals = list(pvalues) if pvalues is not None else [p] * n
    n = len(pvals)
    return pl.DataFrame(
        {
            "chrom": ["chr1"] * n,
            "start": list(range(0, n * 1000, 1000)),
            "end": list(range(1000, n * 1000 + 1000, 1000)),
            "n_cpgs": [10] * n,
            "meth_diff": [0.3] * n,
            "pvalue": pvals,
        },
        schema=_SCHEMA,
    )


def _run(observed, *, n_perm=20, fake=None, **kwargs):
    return dmr_mod.empirical_fdr_for_dmr(
        methylstore_path="/dev/null",
        samples_treatment=_TREAT,
        samples_control=_CTRL,
        observed_dmr=observed,
        n_perm=n_perm,
        seed=0,
        n_jobs=1,
        **kwargs,
    )


def _scripted_engine(script):
    """Per-permutation engine: each entry is an exception to raise, ``None``
    for a clean zero-survivor run, or the null p-values to return."""
    calls = {"i": 0}

    def fake(**kwargs):
        entry = script[calls["i"] % len(script)]
        calls["i"] += 1
        if isinstance(entry, Exception):
            raise entry
        if entry is None:
            return pl.DataFrame({"pvalue": []}, schema={"pvalue": pl.Float64})
        return pl.DataFrame({"pvalue": list(entry)})

    return fake


# max_t: default, unchanged numbers


def test_max_t_is_the_default(monkeypatch):
    """Omitting ``fdr_method`` and passing ``"max_t"`` give identical tables,
    and the new ``empirical_fdr_set`` column is NaN in that mode."""
    monkeypatch.setattr(
        dmr_mod, "call_dmr_tile_based", lambda **kw: pl.DataFrame({"pvalue": [1e-6, 1e-6]})
    )
    implicit = _run(_observed(50, 1e-6))
    explicit = _run(_observed(50, 1e-6), fdr_method="max_t")
    assert implicit.equals(explicit)
    assert set(implicit.columns) >= {"empirical_pvalue", "empirical_qvalue", "empirical_fdr_set"}
    assert np.all(np.isnan(implicit.get_column("empirical_fdr_set").to_numpy()))


def test_max_t_with_failed_and_empty_permutations_matches_main(monkeypatch):
    """One permutation raises, one yields zero survivors, three succeed.
    Main's max_t excluded both kinds from the denominator, so
    emp_p = (#successful perms with min null p <= p_obs + 1) / (3 + 1)."""
    script = [
        RuntimeError("engine blew up"),
        None,
        [0.5, 0.9],  # min 0.5
        [0.05, 0.7],  # min 0.05
        [0.001],  # min 0.001
    ]
    observed = _observed(pvalues=[0.01, 0.3, 1e-5])
    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based", _scripted_engine(script))
    implicit = _run(observed, n_perm=5)
    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based", _scripted_engine(script))
    explicit = _run(observed, n_perm=5, fdr_method="max_t")
    assert implicit.equals(explicit)

    emp_p = implicit.get_column("empirical_pvalue").to_numpy()
    # p_obs=0.01: one min null (0.001) <= it -> (1+1)/4; p_obs=0.3: two
    # (0.05, 0.001) -> (2+1)/4; p_obs=1e-5: none -> (0+1)/4.
    assert emp_p == pytest.approx([2 / 4, 3 / 4, 1 / 4])
    assert np.all(np.isnan(implicit.get_column("empirical_fdr_set").to_numpy()))


def test_max_t_does_not_exclude_self_or_mirror(monkeypatch):
    """The legacy construction counts every permutation that produced a
    null, including assignments equal to the observed split or its mirror."""
    seen: list[tuple[str, ...]] = []

    def fake(**kw):
        seen.append(tuple(sorted(kw["samples_treatment"])))
        return pl.DataFrame({"pvalue": [0.5]})

    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based", fake)
    out = _run(_observed(pvalues=[0.01]), n_perm=30, fdr_method="max_t")
    self_or_mirror = {tuple(sorted(_TREAT)), tuple(sorted(_CTRL))}
    n_excludable = sum(1 for s in seen if s in self_or_mirror)
    assert n_excludable > 0, "seed 0 draws no self/mirror assignment; adjust n_perm"
    # All 30 perms contribute a min null p of 0.5 > 0.01: emp_p = 1/(30+1).
    assert out.get_column("empirical_pvalue")[0] == pytest.approx(1 / 31)


def test_max_t_mode_still_saturates(monkeypatch):
    """Two extreme nulls in every shuffle: the min-P bar makes every observed
    tile look unremarkable."""
    monkeypatch.setattr(
        dmr_mod, "call_dmr_tile_based", lambda **kw: pl.DataFrame({"pvalue": [1e-6, 1e-6]})
    )
    out = _run(_observed(50, 1e-6), fdr_method="max_t")
    emp_p = out.get_column("empirical_pvalue").to_numpy()
    assert np.all(emp_p > 0.5), f"max-T should saturate; got max emp_p={emp_p.max()}"


# region: opt-in count-ratio


def test_region_mode_significant_where_maxt_saturates(monkeypatch):
    """50 real tiles at p=1e-6; every usable shuffle yields only 2 extreme
    nulls. Count-ratio: V(1e-6)=2, R=50 -> q=0.04 (significant)."""
    monkeypatch.setattr(
        dmr_mod, "call_dmr_tile_based", lambda **kw: pl.DataFrame({"pvalue": [1e-6, 1e-6]})
    )
    out = _run(_observed(50, 1e-6), fdr_method="region")
    q = out.get_column("empirical_qvalue").to_numpy()
    assert np.all(q < 0.05), f"region q not significant: max={q.max()}"
    assert out.get_column("empirical_fdr_set")[0] == pytest.approx(0.04, abs=1e-6)


def test_region_distinguishes_failed_empty_and_excluded(monkeypatch):
    """Failed runs and self/mirror assignments leave the null; a clean
    zero-survivor run stays as a zero contribution."""
    seen: list[bool] = []

    def fake(**kw):
        is_self = set(kw["samples_treatment"]) in (set(_TREAT), set(_CTRL))
        seen.append(is_self)
        if is_self:
            # Would poison the null if counted: a decoy that beats every target.
            return pl.DataFrame({"pvalue": [1e-12] * 100})
        idx = len(seen) - 1
        if idx % 3 == 1:
            raise RuntimeError("engine failure")
        if idx % 3 == 2:
            return pl.DataFrame({"pvalue": []}, schema={"pvalue": pl.Float64})
        return pl.DataFrame({"pvalue": [0.5]})

    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based", fake)
    out = _run(_observed(pvalues=[0.01, 0.02]), n_perm=30, fdr_method="region")
    assert any(seen), "seed 0 draws no self/mirror assignment; adjust n_perm"
    non_self = [i for i, s in enumerate(seen) if not s]
    n_used = sum(1 for i in non_self if i % 3 != 1)
    n_with_null = sum(1 for i in non_self if i % 3 == 0)
    # Every usable null (0.5) is above both observed p-values: q = 0 and the
    # set-level FDR is the mean survivor count over the usable runs only.
    assert out.get_column("empirical_qvalue").to_numpy() == pytest.approx([0.0, 0.0])
    assert out.get_column("empirical_fdr_set")[0] == pytest.approx(n_with_null / n_used / 2)


def test_region_zero_usable_assignments_returns_nan_with_warning(monkeypatch):
    """Every permutation fails: NaN estimates and a UserWarning, never 0."""

    def fake(**kw):
        raise RuntimeError("no engine")

    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based", fake)
    with pytest.warns(UserWarning, match="undefined"):
        out = _run(_observed(pvalues=[0.01, 0.02]), n_perm=4, fdr_method="region")
    assert np.all(np.isnan(out.get_column("empirical_qvalue").to_numpy()))
    assert np.all(np.isnan(out.get_column("empirical_pvalue").to_numpy()))
    assert np.all(np.isnan(out.get_column("empirical_fdr_set").to_numpy()))


def test_region_preserves_nonfinite_observed_as_nan(monkeypatch):
    monkeypatch.setattr(
        dmr_mod, "call_dmr_tile_based", lambda **kw: pl.DataFrame({"pvalue": [0.5]})
    )
    out = _run(_observed(pvalues=[0.01, float("nan"), 0.02]), n_perm=6, fdr_method="region")
    q = out.get_column("empirical_qvalue").to_numpy()
    p = out.get_column("empirical_pvalue").to_numpy()
    assert np.isnan(q[1]) and np.isnan(p[1])
    assert np.isfinite(q[[0, 2]]).all() and np.isfinite(p[[0, 2]]).all()


def test_region_small_n_warns(monkeypatch):
    """Permutation inference at <4/group is underpowered -> UserWarning."""
    monkeypatch.setattr(
        dmr_mod, "call_dmr_tile_based", lambda **kw: pl.DataFrame({"pvalue": [1e-6]})
    )
    with pytest.warns(UserWarning, match="underpowered"):
        _run(_observed(10, 1e-6), n_perm=10, fdr_method="region")


def test_max_t_small_n_does_not_warn(monkeypatch):
    monkeypatch.setattr(
        dmr_mod, "call_dmr_tile_based", lambda **kw: pl.DataFrame({"pvalue": [1e-6]})
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _run(_observed(10, 1e-6), n_perm=10, fdr_method="max_t")


# validation and edge cases


def test_empty_observed_returns_three_null_columns():
    out = dmr_mod.empirical_fdr_for_dmr(
        methylstore_path="/dev/null",
        samples_treatment=_TREAT,
        samples_control=_CTRL,
        observed_dmr=pl.DataFrame(schema=_SCHEMA),
        n_perm=5,
    )
    assert out.height == 0
    for col in ("empirical_pvalue", "empirical_qvalue", "empirical_fdr_set"):
        assert out.schema[col] == pl.Float64


@pytest.mark.parametrize(
    "kwargs",
    [{"fdr_method": "bogus"}, {"n_perm": 0}, {"n_perm": -3}],
)
def test_invalid_request_raises_before_any_permutation(monkeypatch, kwargs):
    def fake(**kw):
        raise AssertionError("no permutation must run")

    monkeypatch.setattr(dmr_mod, "call_dmr_tile_based", fake)
    with pytest.raises(ValueError, match=r"fdr_method|n_perm"):
        _run(
            _observed(5, 1e-6),
            n_perm=kwargs.get("n_perm", 5),
            **{k: v for k, v in kwargs.items() if k != "n_perm"},
        )
