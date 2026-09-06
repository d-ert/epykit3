"""The stage split of ``tl.dmc`` (``epykit._dmc_stages``).

One test per behaviour the split has to keep or fix:

* ``open_input_store`` builds the smoothed pseudo-count store in a temp
  directory and removes it on both normal and exceptional exit;
* ``materialize=False`` never assembles the full table on the binary path;
* the ``resumable=True`` cache hit emits the transitional-column warning
  like every other path (it was silent before the split);
* every warning ``tl.dmc`` raises names the caller's file, whether it comes
  from a stage, from the context manager or from a shared ``tl`` helper.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import polars as pl
import pytest

import epykit as ep
import epykit.tl as tl_mod
from epykit._dmc_config import DMCConfig
from epykit._dmc_stages import DMCPlan, open_input_store
from epykit._dmc_store import DMCStore
from epykit.methyldata import MethylData
from tests.fixtures.synth import SimConfig, generate

HERE = Path(__file__).resolve()


# Smoothed input store


def _write_part(root: Path, sample: str, chrom: str, df: pl.DataFrame) -> None:
    part_dir = root / f"sample={sample}" / f"chrom={chrom}"
    part_dir.mkdir(parents=True)
    df.write_parquet(part_dir / "part-0.parquet")


def _raw_chunk(chrom: str, pos: list[int], n_meth: list[int], coverage: list[int]) -> pl.DataFrame:
    n = len(pos)
    return pl.DataFrame(
        {
            "chrom": [chrom] * n,
            "pos": pos,
            "strand": ["+"] * n,
            "context": ["CpG"] * n,
            "N_meth": n_meth,
            "N_unmeth": [c - m for c, m in zip(coverage, n_meth, strict=True)],
            "coverage": coverage,
            "sample": ["S"] * n,
        }
    )


@pytest.fixture
def smoothed_md(tmp_path: Path) -> MethylData:
    """A one-sample raw store and its smoothed sidecar, hand-built so every
    branch of the pseudo-count transform has one site.

    chr1: ``beta_smooth`` 0.26 (rounds up), 1.2 (clips to coverage), -0.1
    (clips to zero) and NaN (falls back to the raw counts). chr2 is absent
    from the sidecar, so the whole chromosome falls back to raw counts.
    """
    raw = tmp_path / "raw"
    smooth = tmp_path / "smooth"
    _write_part(raw, "S", "chr1", _raw_chunk("chr1", [10, 20, 30, 40], [1, 1, 1, 7], [10] * 4))
    _write_part(raw, "S", "chr2", _raw_chunk("chr2", [5], [4], [8]))
    _write_part(
        smooth,
        "S",
        "chr1",
        pl.DataFrame(
            {
                "chrom": ["chr1"] * 4,
                "pos": [10, 20, 30, 40],
                "sample": ["S"] * 4,
                "beta_raw": [0.1, 0.1, 0.1, 0.7],
                "beta_smooth": [0.26, 1.2, -0.1, float("nan")],
            }
        ),
    )
    md = MethylData(obs=pl.DataFrame({"sample_id": ["S"], "treatment": [1]}), store=str(raw))
    md.uns["smooth_path"] = str(smooth)
    return md


def _binary_plan(cfg: DMCConfig) -> DMCPlan:
    return DMCPlan(
        cfg=cfg,
        mode="binary",
        selected_test="lr",
        unite=False,
        smooth_method="bsmooth" if cfg.use_smoothed else None,
        key="dmc_lr_smoothed" if cfg.use_smoothed else "dmc_lr",
        tsv=None,
    )


def test_open_input_store_builds_and_removes_the_smoothed_store(smoothed_md):
    plan = _binary_plan(DMCConfig(use_smoothed=True))
    raw_columns = pl.read_parquet(
        Path(smoothed_md.store) / "sample=S" / "chrom=chr1" / "part-0.parquet"
    ).columns

    with (
        pytest.warns(DeprecationWarning, match="smoothing=True"),
        open_input_store(smoothed_md, plan) as store_path,
    ):
        tmp = Path(store_path)
        assert tmp.is_dir() and tmp != Path(smoothed_md.store)
        chr1 = pl.read_parquet(tmp / "sample=S" / "chrom=chr1" / "part-0.parquet")
        assert chr1.columns == raw_columns
        assert chr1["coverage"].to_list() == [10, 10, 10, 10]
        # round(0.26 * 10) = 3; 1.2 * 10 clips to 10; -0.1 * 10 clips to 0;
        # NaN keeps the raw N_meth of 7. N_unmeth is coverage - N_meth.
        assert chr1["N_meth"].to_list() == [3, 10, 0, 7]
        assert chr1["N_unmeth"].to_list() == [7, 0, 10, 3]
        # chr2 has no smoothed sidecar: raw counts pass through unchanged.
        chr2 = pl.read_parquet(tmp / "sample=S" / "chrom=chr2" / "part-0.parquet")
        assert chr2["N_meth"].to_list() == [4]
        assert chr2["N_unmeth"].to_list() == [4]
    assert not tmp.exists()

    # An exception inside the with block removes the temp store too.
    with (
        pytest.raises(RuntimeError, match="engine failed"),
        pytest.warns(DeprecationWarning),
        open_input_store(smoothed_md, plan) as store_path,
    ):
        tmp = Path(store_path)
        assert tmp.is_dir()
        raise RuntimeError("engine failed")
    assert not tmp.exists()


def test_open_input_store_yields_the_raw_store_by_default(smoothed_md):
    plan = _binary_plan(DMCConfig())
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with open_input_store(smoothed_md, plan) as store_path:
            assert store_path == smoothed_md.store


def test_open_input_store_requires_a_smoothed_sidecar(smoothed_md):
    del smoothed_md.uns["smooth_path"]
    plan = _binary_plan(DMCConfig(use_smoothed=True))
    with (
        pytest.warns(DeprecationWarning),
        pytest.raises(ValueError, match=r"use_smoothed=True requires ep\.pp\.smooth"),
        open_input_store(smoothed_md, plan),
    ):
        pass


# End-to-end paths on a small cohort


@pytest.fixture(scope="module")
def samplesheet(tmp_path_factory) -> str:
    cfg = SimConfig(
        n_per_group=3,
        chromosomes=("chr1",),
        cpgs_per_chrom=200,
        n_dmrs=1,
        n_scattered_dmcs=40,
    )
    return generate(cfg, tmp_path_factory.mktemp("dmc_stages_synth"))["samplesheet"]


def _read(samplesheet: str, store_dir: Path) -> MethylData:
    md = ep.read_bismark(
        samplesheet,
        treatment_group="treatment",
        control_group="control",
        store_dir=str(store_dir),
    )
    ep.pp.filter_coverage(md, lo_count=3, hi_perc=99.9)
    ep.pp.set_unite_type(md, type="intersect")
    return md


@pytest.fixture
def md(samplesheet, tmp_path) -> MethylData:
    """Fresh filtered + united MethylData per test, with its own analysis
    root so ``resumable=True`` starts from an empty pipeline manifest."""
    return _read(samplesheet, tmp_path / "store")


def test_materialize_false_never_assembles_the_table(md, monkeypatch):
    """The binary path with ``materialize=False`` keeps the streaming
    DMCStore as the source of truth: ``DMCStore.to_dataframe`` is never
    called, no ``dmc_*`` table lands on ``md.varm``, and the metadata still
    records the site count from the store manifest."""
    calls: list[str] = []
    original = DMCStore.to_dataframe

    def spy(self, *args, **kwargs):
        calls.append(str(self.path))
        return original(self, *args, **kwargs)

    monkeypatch.setattr(DMCStore, "to_dataframe", spy)

    ep.tl.dmc(md, test="lr", materialize=False)

    assert calls == []
    assert not any(k.startswith("dmc") for k in md.varm), dict(md.varm)
    rec = md.uns["dmc"]
    assert rec["materialized"] is False
    assert rec["last_key"] == "dmc_lr"
    assert rec["n_sites"] == md.dmc_store.total_sites > 0


def test_resume_hit_emits_the_transitional_column_warning(md):
    """Issue 13: the ``log2_odds_ratio`` FutureWarning is emitted on every
    path, including the ``resumable=True`` cache hit that was silent before
    the split, and it still names the caller."""
    with pytest.warns(FutureWarning, match="log2_odds_ratio"):
        ep.tl.dmc(md, test="lr", resumable=True, tsv=False)
    assert md.uns["dmc"]["resumed"] is False

    with pytest.warns(FutureWarning, match="log2_odds_ratio") as record:
        ep.tl.dmc(md, test="lr", resumable=True, tsv=False)
    assert md.uns["dmc"]["resumed"] is True
    hits = [w for w in record if "log2_odds_ratio" in str(w.message)]
    assert len(hits) == 1
    assert Path(hits[0].filename).resolve() == HERE


def _warnings_from(fn) -> list[warnings.WarningMessage]:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fn()
    return list(caught)


def _assert_named_here(caught: list[warnings.WarningMessage], needle: str) -> None:
    hits = [w for w in caught if needle in str(w.message)]
    assert hits, f"no warning containing {needle!r}; got {[str(w.message)[:60] for w in caught]}"
    files = {Path(w.filename).resolve() for w in hits}
    assert files == {HERE}, f"{needle!r} points at {files}, not the caller"


def test_dmc_warnings_point_at_the_caller(md):
    """Warnings raised inside the stages keep their ``stacklevel`` aimed at
    the caller of ``tl.dmc``: the csv alias and the union footgun come from
    shared ``tl`` helpers called by ``plan_run``, ``use_smoothed`` from the
    ``open_input_store`` context manager, the transitional-column notice
    from ``finish``, and the explicit-Fisher notice from the shared one-shot
    helper. Each one is checked against this file, not a module in epykit."""
    ep.pp.smooth(md, method="gaussian", bandwidth=1000)
    ep.pp.set_unite_type(md, type="union")
    tl_mod._FISHER_WARNED = False

    caught = _warnings_from(
        lambda: ep.tl.dmc(md, test="lr", use_smoothed=True, tsv=False, csv_alpha=0.01)
    )
    _assert_named_here(caught, "`csv` / `csv_full` / `csv_alpha` arguments are deprecated")
    _assert_named_here(caught, "unite='union' with min_samples_treatment")
    _assert_named_here(caught, "use_smoothed=True")
    _assert_named_here(caught, "'log2_odds_ratio' column is deprecated")

    caught = _warnings_from(lambda: ep.tl.dmc(md, test="fisher", tsv=False))
    _assert_named_here(caught, "ignores between-replicate variance")


def test_n1_fallback_warning_points_at_the_caller(tmp_path):
    """The n<2 Fisher fallback notice comes from ``_auto_test_simple``, two
    helpers below ``tl.dmc`` before the split, where it named ``tl.py``. It
    now names the caller like every other ``tl.dmc`` warning."""
    cfg = SimConfig(
        n_per_group=1,
        chromosomes=("chr1",),
        cpgs_per_chrom=200,
        n_dmrs=1,
        n_scattered_dmcs=40,
    )
    sheet = generate(cfg, tmp_path / "synth_n1")["samplesheet"]
    md = _read(sheet, tmp_path / "store_n1")
    tl_mod._FISHER_WARNED = False

    caught = _warnings_from(lambda: ep.tl.dmc(md, test="auto", allow_n1=True, tsv=False))
    _assert_named_here(caught, "n<2 per group: falling back to Fisher exact")
    _assert_named_here(caught, "ignores between-replicate variance")
    assert md.uns["dmc"]["test_used"] == "fisher"
