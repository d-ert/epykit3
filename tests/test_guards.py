"""Layer 5: footgun-guard tests.

Verifies the safety nets that prevent users from silently getting wrong
answers:

* n=1 per group -> ``ValueError`` unless ``allow_n1=True`` (B6).
* ``allow_n1=True`` runs Fisher but emits a ``UserWarning`` (B6).
* ``unite='union'`` + ``min_samples_*=0`` warns once (B8).
* Explicit ``test='fisher'`` emits a one-shot ``UserWarning`` per session
  (B6 second half).
"""

from __future__ import annotations

import warnings

import polars as pl
import pytest

from tests.fixtures.synth import SimConfig, generate



# n=1 per group guard


@pytest.fixture
def synth_md_n1(tmp_path):
    """A MethylData with only 1 replicate per group -- exercises the n=1 path."""
    import epykit as ep
    cfg = SimConfig(
        n_per_group=1,
        chromosomes=("chr1",),
        cpgs_per_chrom=200,
        n_scattered_dmcs=20,
        n_dmrs=2,
        dmr_size_cpgs=5,
        seed=7,
    )
    result = generate(cfg, tmp_path / "n1")
    md = ep.read_bismark(
        result["samplesheet"],
        treatment_group="treatment",
        control_group="control",
        store_dir=str(tmp_path / "store"),
    )
    ep.pp.filter_coverage(md, lo_count=2, hi_perc=99.9)
    ep.pp.unite(md, type="intersect")
    return md


def test_n1_dmc_raises_without_allow_n1(synth_md_n1):
    """``tl.dmc(md)`` with n=1 must refuse statistical inference."""
    import epykit as ep
    with pytest.raises(ValueError, match=r"(?i)at least 2 replicates"):
        ep.tl.dmc(synth_md_n1, test="auto")


def test_n1_dmc_explicit_test_also_raises(synth_md_n1):
    """Even with an explicit test name, n=1 should still be refused."""
    import epykit as ep
    with pytest.raises(ValueError, match=r"(?i)at least 2 replicates"):
        ep.tl.dmc(synth_md_n1, test="lr")


def test_n1_dmc_with_allow_n1_runs_and_warns(synth_md_n1):
    """``allow_n1=True`` opts in: function runs but emits UserWarning."""
    import epykit as ep
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ep.tl.dmc(synth_md_n1, test="auto", allow_n1=True)
    fired = [
        w for w in caught
        if issubclass(w.category, UserWarning)
        and ("fisher" in str(w.message).lower() or "n<2" in str(w.message).lower())
    ]
    assert fired, "expected a UserWarning explaining the n=1 fallback"
    # Auto-resolved test should be fisher.
    assert synth_md_n1.uns["dmc"]["test_used"] == "fisher"


def test_n1_dmr_tile_raises_without_allow_n1(synth_md_n1):
    """The tile-based DMR path also enforces n>=2."""
    import epykit as ep
    with pytest.raises(ValueError, match=r"(?i)at least 2 replicates"):
        ep.tl.dmr(synth_md_n1, method="tile", tile_size_bp=500, min_cpgs_per_tile=2)



# union + zero-zero min_samples footgun


def test_union_with_zero_min_samples_warns(synth_md, tmp_path):
    """``pp.unite('union') + tl.dmc(min_samples_*=0)`` must warn -- testing
    sites covered in only one sample per group is the textbook footgun."""
    import epykit as ep
    ep.pp.filter_coverage(synth_md, lo_count=5, hi_perc=99.9)
    ep.pp.unite(synth_md, type="union")

    with pytest.warns(UserWarning, match=r"unite='union'"):
        ep.tl.dmc(synth_md, test="lr",
                  min_samples_treatment=0, min_samples_control=0)


def test_union_with_explicit_min_samples_does_not_warn(synth_md):
    """Providing sensible min_samples_* with union mode suppresses the warning."""
    import epykit as ep
    ep.pp.filter_coverage(synth_md, lo_count=5, hi_perc=99.9)
    ep.pp.unite(synth_md, type="union")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ep.tl.dmc(synth_md, test="lr",
                  min_samples_treatment=2, min_samples_control=2)
    union_warns = [
        w for w in caught
        if issubclass(w.category, UserWarning) and "unite='union'" in str(w.message)
    ]
    assert union_warns == [], "should not fire union warning when min_samples_*>=2"


def test_intersect_mode_never_warns_about_union(synth_md_filtered):
    """unite='intersect' should not trigger the union footgun warning."""
    import epykit as ep
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ep.tl.dmc(synth_md_filtered, test="lr",
                  min_samples_treatment=0, min_samples_control=0)
    union_warns = [
        w for w in caught
        if issubclass(w.category, UserWarning) and "unite='union'" in str(w.message)
    ]
    assert union_warns == []


# Removed min_samples_case kwarg


def test_min_samples_case_kwarg_rejected_on_dmc(synth_md_filtered):
    """``tl.dmc(min_samples_case=...)`` was removed; raises TypeError."""
    import epykit as ep
    with pytest.raises(TypeError):
        ep.tl.dmc(synth_md_filtered, test="lr",
                  min_samples_case=2, min_samples_control=2)



# Explicit fisher warning (one-shot per session)


def test_explicit_fisher_emits_user_warning(synth_md_filtered):
    """Picking ``test='fisher'`` outside of the auto fallback should warn
    so the user knows what they signed up for."""
    import epykit as ep
    # Reset the one-shot gate so a previous test doesn't suppress this run.
    import epykit.tl as _tl_module
    _tl_module._FISHER_WARNED = False

    with pytest.warns(UserWarning, match="fisher"):
        ep.tl.dmc(synth_md_filtered, test="fisher")


def test_fisher_warning_fires_only_once_per_session(synth_md_filtered):
    """The *tl-level* fisher warning should fire only once per session.

    Two distinct fisher warnings exist after Cleanup 2:

    * tl-level (``tl._warn_fisher_once``): one-shot per session, gated
      by ``_FISHER_WARNED``. Message: "test='fisher' ignores
      between-replicate variance; ...". This test asserts that gate works.
    * dmc-level (``_validate_sample_size_and_warn``): once per
      ``process_chromosomes_dmc`` call. Message: "test='fisher' pools
      reads across replicates; ...". *Intentionally* not session-gated --
      direct-API users get warned every call.

    Both messages contain "Prefer test='lr'", so we identify the
    tl-level one by its unique phrase "ignores between-replicate".
    """
    import epykit as ep
    import epykit.tl as _tl_module

    _tl_module._FISHER_WARNED = False
    ep.tl.dmc(synth_md_filtered, test="fisher")  # primes the gate

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ep.tl.dmc(synth_md_filtered, test="fisher")
    tl_level_warns = [
        w for w in caught
        if issubclass(w.category, UserWarning)
        and "ignores between-replicate" in str(w.message)
    ]
    assert tl_level_warns == [], (
        "tl-level fisher warning should be one-shot per session "
        f"(saw {len(tl_level_warns)} on second call)"
    )



# pp ordering checks (no normalize before filter, no smooth before filter)


def test_normalize_coverage_before_filter_raises(synth_md):
    """``pp.normalize_coverage`` must refuse to run before ``filter_coverage``."""
    import epykit as ep
    with pytest.raises(ValueError, match=r"filter_coverage"):
        ep.pp.normalize_coverage(synth_md, method="median")


def test_smooth_before_filter_raises(synth_md):
    """``pp.smooth`` must refuse to run before ``filter_coverage`` ()."""
    import epykit as ep
    with pytest.raises(ValueError, match=r"filter_coverage"):
        ep.pp.smooth(synth_md, bandwidth=500)
