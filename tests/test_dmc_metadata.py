"""The contract for ``md.uns["dmc"]`` and the per-engine DMC column set.

``ep.tl.dmc`` records one metadata dict per run. Readers are
``MethylData.get_dmc`` / ``.dmc_store`` / ``.save`` (``last_key``,
``store_path``), ``report.py`` (``test_used``, ``test_requested``,
``fdr_method``), ``multiqc_export.py`` and ``pl.dashboard`` (``test_used``,
``n_sites``, ``unite``), ``export.export_tables`` and the CLI
(``last_key``), and ``test_resume.py`` (``resumed``). Every code path that
writes the record must produce the same key set; fields that do not apply
to a path are present with ``None``.

Drift recorded on the unmodified tree (commit 1 of the refactor series):

* the main path wrote 24 keys and lacked ``resumed`` and the five
  contrast fields ``formula``, ``contrast``, ``design_terms``,
  ``covariates``, ``treatment_col``;
* the ``resumable=True`` cache-hit path wrote 19 keys: it had ``resumed``
  but lacked ``materialized``, ``use_smoothed``, ``smooth_method``,
  ``smoothing``, ``smoothing_span_bp``, ``store_path`` and the five
  contrast fields;
* the ``formula`` / ``contrast`` path wrote 16 keys: it had the five
  contrast fields but lacked ``materialized``, ``empirical_fdr``,
  ``n_perm``, ``perm_seed``, ``power_stack``, ``sep_fallback``,
  ``sep_threshold``, ``neighbour_combine``, ``neighbour_bp``,
  ``use_smoothed``, ``smooth_method``, ``smoothing``,
  ``smoothing_span_bp`` and ``resumed``.

The ``uns["dmc"]`` key-set tests therefore failed on all three paths before
commit 2 unified the writer.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

import epykit as ep
from epykit.methyldata import MethylData
from tests.fixtures.synth import SimConfig, generate

# Main-path keys as written on the unmodified tree, in their original
# order, followed by the path-specific fields that the unified writer
# carries on every path.
CANONICAL_UNS_KEYS = [
    "test_requested",
    "test_used",
    "n_sites",
    "materialized",
    "unite",
    "min_samples_treatment",
    "min_samples_control",
    "dispersion",
    "reference",
    "empirical_fdr",
    "n_perm",
    "perm_seed",
    "power_stack",
    "sep_fallback",
    "sep_threshold",
    "neighbour_combine",
    "neighbour_bp",
    "fdr_method",
    "last_key",
    "use_smoothed",
    "smooth_method",
    "smoothing",
    "smoothing_span_bp",
    "store_path",
    # resume-hit path
    "resumed",
    # formula / contrast path
    "formula",
    "contrast",
    "design_terms",
    "covariates",
    "treatment_col",
]


def _assert_canonical_keys(record: dict, path: str) -> None:
    missing = sorted(set(CANONICAL_UNS_KEYS) - set(record))
    extra = sorted(set(record) - set(CANONICAL_UNS_KEYS))
    assert not missing and not extra, f"{path} path: missing={missing} extra={extra}"


def _assert_store_path_is_live(record: dict) -> None:
    store_dir = Path(record["store_path"])
    assert (store_dir / ".epykit_dmc_manifest.json").exists(), store_dir


@pytest.fixture(scope="module")
def samplesheet(tmp_path_factory) -> str:
    cfg = SimConfig(
        n_per_group=3,
        chromosomes=("chr1",),
        cpgs_per_chrom=200,
        n_dmrs=1,
        n_scattered_dmcs=40,
    )
    return generate(cfg, tmp_path_factory.mktemp("dmc_meta_synth"))["samplesheet"]


@pytest.fixture
def md(samplesheet, tmp_path) -> MethylData:
    """Fresh filtered + united MethylData per test, with its own analysis root
    so ``resumable=True`` starts from an empty pipeline manifest."""
    md = ep.read_bismark(
        samplesheet,
        treatment_group="treatment",
        control_group="control",
        store_dir=str(tmp_path / "store"),
    )
    ep.pp.filter_coverage(md, lo_count=3, hi_perc=99.9)
    ep.pp.set_unite_type(md, type="intersect")
    return md


@pytest.fixture
def md_n1(tmp_path) -> MethylData:
    """One sample per group: the only cohort on which ``test="fisher"`` runs
    (behind ``allow_n1=True``)."""
    cfg = SimConfig(
        n_per_group=1,
        chromosomes=("chr1",),
        cpgs_per_chrom=200,
        n_dmrs=1,
        n_scattered_dmcs=40,
    )
    res = generate(cfg, tmp_path / "synth_n1")
    md = ep.read_bismark(
        res["samplesheet"],
        treatment_group="treatment",
        control_group="control",
        store_dir=str(tmp_path / "store_n1"),
    )
    ep.pp.filter_coverage(md, lo_count=3, hi_perc=99.9)
    ep.pp.set_unite_type(md, type="intersect")
    return md


# uns["dmc"] key set


def test_main_path_writes_canonical_record(md):
    ep.tl.dmc(md, test="lr", tsv=False)
    rec = md.uns["dmc"]
    _assert_canonical_keys(rec, "main")

    assert rec["last_key"] == "dmc_lr"
    assert "dmc_lr" in md.varm
    _assert_store_path_is_live(rec)
    assert md.dmc_store.path == Path(rec["store_path"])
    assert rec["test_requested"] == "lr"
    assert rec["test_used"] == "lr"
    assert rec["n_sites"] == len(md.varm["dmc_lr"]) > 0
    assert rec["materialized"] is True
    assert rec["unite"] is True
    assert rec["resumed"] is False
    for key in ("formula", "contrast", "design_terms", "covariates", "treatment_col"):
        assert rec[key] is None, key


def test_resume_hit_path_writes_canonical_record(md):
    ep.tl.dmc(md, test="lr", resumable=True, tsv=False)
    computed = md.uns["dmc"]

    ep.tl.dmc(md, test="lr", resumable=True, tsv=False)
    rec = md.uns["dmc"]
    _assert_canonical_keys(rec, "resume-hit")

    assert computed["resumed"] is False  # first call computed
    assert rec["resumed"] is True
    assert rec["last_key"] == "dmc_lr"
    assert rec["n_sites"] == computed["n_sites"] == len(md.varm["dmc_lr"])
    assert rec["materialized"] is True  # the sidecar was loaded onto md.varm
    # The cache hit never opens a DMCStore, so no store path is recorded and
    # md.dmc_store stays None on this path.
    assert rec["store_path"] is None
    assert md.dmc_store is None
    # Everything the resume signature covers is carried over unchanged.
    for key in (
        "test_requested",
        "test_used",
        "unite",
        "min_samples_treatment",
        "min_samples_control",
        "dispersion",
        "reference",
        "empirical_fdr",
        "n_perm",
        "perm_seed",
        "power_stack",
        "sep_fallback",
        "sep_threshold",
        "neighbour_combine",
        "neighbour_bp",
        "fdr_method",
    ):
        assert rec[key] == computed[key], key


def test_contrast_path_writes_canonical_record(md):
    ep.tl.dmc(md, formula="~ group", contrast="group", tsv=False)
    rec = md.uns["dmc"]
    _assert_canonical_keys(rec, "contrast")

    assert rec["last_key"] == "dmc_glm_contrast"
    assert "dmc_glm_contrast" in md.varm
    _assert_store_path_is_live(rec)
    assert md.dmc_store.path == Path(rec["store_path"])
    assert rec["test_requested"] == "auto"
    assert rec["test_used"] == "glm_contrast"
    assert rec["n_sites"] == len(md.varm["dmc_glm_contrast"]) > 0
    assert rec["materialized"] is True
    assert rec["resumed"] is False
    assert rec["formula"] == "~ group"
    assert rec["contrast"] == "group"
    assert rec["design_terms"] == ["Intercept", "group[T.treatment]"]
    assert rec["covariates"] is None
    assert rec["treatment_col"] == "treatment"
    assert rec["fdr_method"] == "fdr_bh"
    # Knobs the contrast path does not consume are recorded as not applicable.
    for key in (
        "empirical_fdr",
        "n_perm",
        "perm_seed",
        "power_stack",
        "sep_fallback",
        "sep_threshold",
        "neighbour_combine",
        "neighbour_bp",
        "use_smoothed",
        "smooth_method",
        "smoothing",
        "smoothing_span_bp",
    ):
        assert rec[key] is None, key


def test_contrast_path_rejects_invalid_power_stack(md):
    """An unknown ``power_stack`` raises the same ValueError as on the binary
    path instead of being ignored, and nothing is written."""
    with pytest.raises(ValueError, match="power_stack must be one of"):
        ep.tl.dmc(md, formula="~ group", contrast="group", power_stack="bogus", tsv=False)
    assert "dmc" not in md.uns
    assert "dmc_glm_contrast" not in md.varm


def test_contrast_path_rejects_materialize_false(md):
    """The contrast path always assembles the full result onto ``md.varm``,
    so ``materialize=False`` is refused instead of being silently ignored."""
    with pytest.raises(ValueError, match="materialize=False is not supported on the formula"):
        ep.tl.dmc(md, formula="~ group", contrast="group", materialize=False, tsv=False)
    assert "dmc" not in md.uns
    assert "dmc_glm_contrast" not in md.varm


def test_contrast_path_ignores_power_stack_with_a_notice(md, caplog):
    """``power_stack="lr+"`` has no GLM knob to switch on: the result table
    and the record equal the default run's, and one INFO line says so."""
    ep.tl.dmc(md, formula="~ group", contrast="group", tsv=False)
    baseline = md.varm["dmc_glm_contrast"]
    baseline_rec = dict(md.uns["dmc"])

    with caplog.at_level(logging.INFO, logger="epykit._dmc_stages"):
        ep.tl.dmc(md, formula="~ group", contrast="group", power_stack="lr+", tsv=False)

    notices = [r.getMessage() for r in caplog.records if "power_stack" in r.getMessage()]
    assert notices == [
        "power_stack='lr+' is ignored on the formula / contrast path: the GLM has no lr+ knobs."
    ]
    assert md.varm["dmc_glm_contrast"].equals(baseline)
    assert md.uns["dmc"] == baseline_rec
    assert md.uns["dmc"]["power_stack"] is None


# Per-engine column set


def _expected_columns(engine: str) -> list[str]:
    """Canonical DMC schema as emitted today, in order.

    ``strand`` is part of the site key. The effect column is
    ``log2_odds_ratio_pooled`` for pooled-count engines and
    ``coef_treatment_log2`` for the GLM (P1-11); ``log2_odds_ratio`` is the
    NaN-filled transitional column slated for removal in 1.2. ``qvalue``
    and ``reject`` are appended by the FDR step. GLM adds the
    ``coef_treatment`` / ``coef_se`` extras.
    """
    effect = "coef_treatment_log2" if engine == "glm" else "log2_odds_ratio_pooled"
    cols = [
        "chrom",
        "pos",
        "strand",
        "n_case",
        "n_control",
        "mean_beta_case",
        "mean_beta_control",
        "pvalue",
        effect,
        "log2_odds_ratio",
        "meth_diff",
        "meth_diff_ci_lo",
        "meth_diff_ci_hi",
    ]
    if engine == "glm":
        cols += ["coef_treatment", "coef_se"]
    return [*cols, "qvalue", "reject"]


@pytest.mark.parametrize(
    "engine,kwargs",
    [
        ("lr", {}),
        ("welch_t", {}),
        ("glm", {"formula": "~ treatment"}),
    ],
)
def test_engine_emits_canonical_columns(md, engine, kwargs):
    ep.tl.dmc(md, test=engine, tsv=False, **kwargs)
    df = md.varm[md.uns["dmc"]["last_key"]]
    assert df.columns == _expected_columns(engine)


def test_fisher_emits_canonical_columns(md_n1):
    ep.tl.dmc(md_n1, test="fisher", allow_n1=True, tsv=False)
    assert md_n1.uns["dmc"]["last_key"] == "dmc_fisher"
    df = md_n1.varm["dmc_fisher"]
    assert df.columns == _expected_columns("fisher")
