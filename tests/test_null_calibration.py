"""Real-engine null-calibration test (M-PKG3).

The single most important property of a differential-methylation caller is
that, under the null (no true difference between groups), its p-values are
*not anti-conservative* -- i.e. the false-positive rate at nominal alpha does
not exceed alpha (beyond sampling tolerance). The benchmark paper claims the
stronger property that epykit's FPR is *tighter* than the R baselines, so a
conservative (sub-uniform) distribution is expected and acceptable; an
*anti-conservative* one is a genuine calibration bug.

Before this test, ``tests/test_calibration.py`` only validated the K-S/FDR
scaffolding against *mock* uniform/Beta engines -- the real ``ep.tl.dmc``
engines had no executing null-calibration check anywhere in ``tests/``. This
fills that gap for ``lr``, ``glm`` (group contrast), and ``welch_t``.

Marked ``slow`` (builds a store + runs three engines); the slow CI job runs it.

NOTE on the assertion direction: we bound the FPR from *above* only. A
conservative engine (FPR below nominal) is valid and is what the paper claims,
so a hard lower bound / two-sided K-S-uniformity assertion would wrongly fail
it. The upper bound is the meaningful guard -- it fails an anti-conservative
engine, which is the real bug. If a future change makes an engine
anti-conservative here, do NOT loosen this bound; that is a finding to fix.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

import epykit as ep
from tests.fixtures.synth import SimConfig, generate

pytestmark = pytest.mark.slow


# Per-threshold upper bounds on the null rejection rate. Generous enough to
# absorb sampling noise on ~10k sites, tight enough to catch a genuinely
# anti-conservative engine (which would roughly double these).
_FPR_UPPER = {0.01: 0.025, 0.05: 0.075, 0.10: 0.130}


def _null_pvalues(md, **dmc_kwargs):
    ep.tl.dmc(md, min_samples_treatment=0, min_samples_control=0, **dmc_kwargs)
    p = md.dmc["pvalue"].to_numpy()
    return p[np.isfinite(p)]


def _build_null_md(tmp_path, store_name):
    """A true-null cohort: every site has zero true effect, so both groups are
    drawn from the same baseline + replicate-noise + binomial model."""
    cfg = SimConfig(
        n_per_group=4,
        n_scattered_dmcs=0,
        n_dmrs=0,
        dmc_effect=0.0,
        dmr_effect=0.0,
        seed=20260606,
    )
    res = generate(cfg, tmp_path / store_name)
    md = ep.read_bismark(
        res["samplesheet"], treatment_group="treatment", control_group="control",
        store_dir=str(tmp_path / f"{store_name}_store"),
    )
    ep.pp.filter_coverage(md, lo_count=5, hi_perc=99.9)
    ep.pp.set_unite_type(md, type="intersect")
    return md


@pytest.mark.parametrize("name,kwargs", [
    ("lr", dict(test="lr")),
    ("welch_t", dict(test="welch_t")),
    ("glm", dict(formula="~ group", contrast="group")),
])
def test_engine_not_anticonservative_under_null(tmp_path, name, kwargs):
    md = _build_null_md(tmp_path, f"null_{name}")
    p = _null_pvalues(md, **kwargs)

    assert p.size > 1000, f"{name}: too few finite p-values ({p.size}) to assess FPR"

    ks = stats.kstest(p, "uniform")
    report = {a: float((p < a).mean()) for a in _FPR_UPPER}
    msg = (
        f"{name}: n={p.size}, FPR@0.01={report[0.01]:.4f}, "
        f"FPR@0.05={report[0.05]:.4f}, FPR@0.10={report[0.10]:.4f}, "
        f"KS_D={ks.statistic:.4f}, KS_p={ks.pvalue:.3g}"
    )
    print(msg)

    for alpha, upper in _FPR_UPPER.items():
        fpr = report[alpha]
        assert fpr <= upper, (
            f"{name} is ANTI-CONSERVATIVE under the null: empirical FPR at "
            f"nominal alpha={alpha} is {fpr:.4f} (> {upper}). {msg}"
        )
