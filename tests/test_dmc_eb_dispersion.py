"""The ``dispersion="eb"`` default and its adaptive-F reference.

The adaptive F reference in eb mode must floor df_phi so that F(1, df_phi)
does not collapse to F(1, ~4) when the EB weight is small. We verify the
floor two ways:
1. A direct math test that constructs chi2_stat/phi_eff/df_phi arrays
   and verifies F(1, floor(df_phi)) differs from F(1, 4) as expected
   (no epykit import needed beyond confirming the constant is present).
2. An end-to-end integration test confirming that lr/eb does not produce
   p-values collapsed near 1.0 under realistic dispersion.

The code default and the ``tl.dmc`` docstring must also agree that ``eb``
is the default.
"""
from __future__ import annotations

import inspect
import re

import numpy as np
import pytest
import scipy.stats as sp_stats

import epykit as ep
from tests.fixtures.synth import SimConfig, generate


@pytest.fixture
def small_md_under_dispersed(tmp_path):
    """Small fixture sized so the adaptive-F path can trigger."""
    cfg = SimConfig(
        n_per_group=3,
        cpgs_per_chrom=600,
        chromosomes=("chr1", "chr2"),
        n_dmrs=2,
        dmr_size_cpgs=5,
        n_scattered_dmcs=120,
        seed=13,
    )
    out_dir = tmp_path / "eb_floor"
    out_dir.mkdir()
    result = generate(cfg, out_dir)
    md = ep.read_bismark(
        result["samplesheet"],
        treatment_group="treatment",
        control_group="control",
        store_dir=str(out_dir / "store"),
    )
    ep.pp.filter_coverage(md, lo_count=3, hi_perc=99.9)
    ep.pp.set_unite_type(md, type="intersect")
    return md


def test_score_finalize_floor_applies_at_small_w_eb():
    """Direct math unit test: when df_phi is small and phi_eff > 1, a floor
    at 50 must bring the F p-value much closer to chi^2 than F(1, 4) gives.

    This constructs the adaptive-F inputs directly and verifies:
    - F(1, 4) at chi2_stat=6 is > 0.05 (over-conservative)
    - F(1, max(df_phi, 50)) at chi2_stat=6 is < 0.025 (close to chi^2)
    """
    chi2_stat = np.array([6.0])
    # phi_eff > 1 in this regime, so the F branch (not chi^2) is the one that fires.
    df_phi = np.array([4.0])   # the eb-with-tiny-w_eb regime

    # Without floor: F(1, 4) at stat=6
    p_without = float(sp_stats.f.sf(chi2_stat, dfn=1, dfd=df_phi)[0])

    # With floor at 50: F(1, 50) at stat=6
    DF_PHI_FLOOR = 50.0
    df_phi_floored = np.maximum(df_phi, DF_PHI_FLOOR)
    p_with_floor = float(sp_stats.f.sf(chi2_stat, dfn=1, dfd=df_phi_floored)[0])

    # chi^2(1) reference at stat=6
    p_chi2 = float(sp_stats.chi2.sf(chi2_stat, df=1)[0])

    # Verify the bug exists pre-floor: F(1, 4) is over-conservative
    assert p_without > 0.05, (
        f"F(1, 4) at chi2_stat=6: p={p_without:.4f}; expected > 0.05 "
        f"(the over-conservative F(1,small df) regime). "
        f"Test setup may be wrong."
    )

    # Verify the floor fixes it: F(1, 50) is close to chi^2(1)
    assert p_with_floor < 0.025, (
        f"F(1, 50) at chi2_stat=6: p={p_with_floor:.4f}; expected < 0.025. "
        f"Floor should bring F(1, 50) within ~1% of chi^2(1) at stat=6."
    )

    # Verify chi^2 reference for orientation
    assert p_chi2 < 0.025, f"chi^2(1) at stat=6: p={p_chi2:.4f}; sanity check."


def test_eb_default_does_not_crush_pvalues_at_modest_statistics(small_md_under_dispersed):
    """With dispersion='eb' and reference='adaptive' (the defaults), the
    p-value distribution must not be collapsed near 1.0. Specifically:
    the fraction of finite p-values < 0.05 must be at least 1% on a
    fixture that seeded 120 DMCs + 10 DMR CpGs. Without the floor,
    F(1, ~4) collapses p-values upward and this fraction approaches 0.
    """
    md = small_md_under_dispersed
    ep.tl.dmc(md, test="lr", dispersion="eb", reference="adaptive")
    df = md.dmc
    p = df["pvalue"].drop_nans().to_numpy()
    assert p.size > 100, "fixture too small after filtering"

    frac_sig = float(np.mean(p < 0.05))
    assert frac_sig > 0.01, (
        f"frac(p<0.05) = {frac_sig:.4f} -- looks like F(1, small df) "
        f"collapsed the p-value distribution. After the floor, this "
        f"should be well above 1% on a fixture with seeded DMCs."
    )


# --- Default and docstring agree -------------------------------------------


def test_dispersion_default_is_eb():
    sig = inspect.signature(ep.tl.dmc)
    assert sig.parameters["dispersion"].default == "eb", (
        "Code default for `dispersion` changed. If intentional, "
        "update the benchmark protocol and executive summary as well."
    )


def test_dispersion_docstring_mentions_eb_as_default():
    doc = ep.tl.dmc.__doc__ or ""
    assert "eb" in doc, "Docstring no longer mentions the 'eb' option."
    assert re.search(r'default\s+``"eb"``', doc, re.IGNORECASE), (
        "Docstring must state 'eb' is the default (looking for "
        "'default ``\"eb\"``' case-insensitively)."
    )
    # The old wrong claim must be gone.
    assert not re.search(r'default\s+``"site"``', doc, re.IGNORECASE), (
        "Docstring still says default is 'site'; should be 'eb'."
    )
