"""Regression guard: ``power_stack="lr+"`` must produce output identical to
setting its four bundled knobs explicitly.

The lr+ bundle (CLAUDE.md / the dispatcher in ``tl.py``) expands to::

    neighbour_combine=True, sep_fallback=True, fdr_method="fdr_tsbh",
    dispersion="eb"           # dispersion="eb" is already the tl.dmc default

This test pins the equivalence so the bundle cannot silently drift from its
documented expansion. It complements the dispatch-level checks in
``test_lr_improvements.py`` (which verify the knobs *flip*, not that the full
per-CpG output is bit-identical to the explicit form).
"""
from __future__ import annotations

from polars.testing import assert_frame_equal

import epykit as ep
from tests.fixtures.synth import SimConfig, generate


def _make_md(tmp_subdir, n_per_group: int = 3):
    """Small two-group MethylData; fixed seed so two calls are identical."""
    cfg = SimConfig(
        n_per_group=n_per_group,
        chromosomes=("chr1",),
        cpgs_per_chrom=500,
        n_scattered_dmcs=50,
        n_dmrs=2,
        dmr_size_cpgs=5,
        seed=123,
    )
    result = generate(cfg, tmp_subdir / "sim")
    md = ep.read_bismark(
        result["samplesheet"],
        treatment_group="treatment",
        control_group="control",
        store_dir=str(tmp_subdir / "store"),
    )
    ep.pp.filter_coverage(md, lo_count=2, hi_perc=99.9)
    ep.pp.set_unite_type(md, type="intersect")
    return md


def test_power_stack_lr_plus_equals_explicit_knobs(tmp_path):
    md_bundle = _make_md(tmp_path / "bundle")
    md_explicit = _make_md(tmp_path / "explicit")

    # The bundle.
    ep.tl.dmc(md_bundle, test="lr", power_stack="lr+")
    # The documented expansion, with the dispatcher disabled.
    ep.tl.dmc(
        md_explicit,
        test="lr",
        power_stack="off",
        neighbour_combine=True,
        sep_fallback=True,
        fdr_method="fdr_tsbh",
        dispersion="eb",
    )

    a = md_bundle.varm["dmc_lr"].sort(["chrom", "pos"])
    b = md_explicit.varm["dmc_lr"].sort(["chrom", "pos"])

    assert a.columns == b.columns, (
        f"power_stack='lr+' columns {a.columns} != explicit knobs {b.columns}"
    )
    assert_frame_equal(a, b, check_exact=False, rel_tol=1e-9, abs_tol=1e-12)
