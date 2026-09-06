"""D-minor: the CLI ``epykit dmc`` contrast path must forward
``--dispersion`` / ``--reference`` / ``--fdr-method`` to ``tl.dmc``.

Pre-fix ``_cmd_dmc``'s formula/contrast branch called ``_tl.dmc(...)`` with
only test/formula/contrast/covariates/min_samples_*, silently dropping the
three knobs the binary path forwards (cli.py binary path passes
``args.dispersion`` / ``args.reference`` / ``args.fdr_method``). So
``epykit dmc --formula '~ group' --contrast group --dispersion site
--fdr-method fdr_by`` produced different q-values from the equivalent API
call.
"""

from __future__ import annotations

import argparse
import types

import polars as pl


def _fake_md():
    """A minimal stand-in for MethylData carrying just uns/varm dicts."""
    return types.SimpleNamespace(uns={}, varm={})


def test_cli_contrast_forwards_dispersion_reference_fdr(tmp_path, monkeypatch):
    import epykit
    import epykit.tl as ep_tl
    from epykit.cli._dmc import _cmd_dmc

    # A samplesheet with a 'group' column is all the contrast handler reads
    # directly (it derives the group set before delegating to tl.dmc).
    samplesheet = tmp_path / "samples.csv"
    samplesheet.write_text(
        "sample_id,group,path\n"
        "s1,treatment,a.cov\n"
        "s2,control,b.cov\n"
    )

    captured: dict = {}

    def _fake_read_bismark(*_args, **_kwargs):
        return _fake_md()

    def _fake_dmc(md, **kwargs):
        captured.update(kwargs)
        md.uns["dmc"] = {"last_key": "dmc_glm_contrast"}
        md.varm["dmc_glm_contrast"] = pl.DataFrame({"chrom": ["chr1"], "pos": [1]})

    monkeypatch.setattr(epykit, "read_bismark", _fake_read_bismark)
    monkeypatch.setattr(ep_tl, "dmc", _fake_dmc)

    out = tmp_path / "out.parquet"
    args = argparse.Namespace(
        formula="~ group",
        contrast="group",
        covariates=None,
        samplesheet=str(samplesheet),
        treatment_group="treatment",
        control_group="control",
        methylstore=str(tmp_path / "store"),
        test="lr",
        min_samples_treatment=0,
        min_samples_control=0,
        dispersion="site",       # non-default (default is "eb")
        reference="F",           # non-default (default is "adaptive")
        fdr_method="fdr_by",     # non-default (default is "fdr_bh")
        output=str(out),
    )

    _cmd_dmc(args)

    assert captured.get("dispersion") == "site", "--dispersion not forwarded"
    assert captured.get("reference") == "F", "--reference not forwarded"
    assert captured.get("fdr_method") == "fdr_by", "--fdr-method not forwarded"
