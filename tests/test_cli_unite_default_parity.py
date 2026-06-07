"""D-minor: bare ``epykit dmc`` / ``epykit dmr --method tile`` must default
to the SAME unite mode as the bare API.

The API (``ep.tl.dmc`` / ``ep.tl.dmr``) derives unite from
``md.uns["unite"]``: with no prior ``ep.pp.unite(...)`` step the key is unset
and unite resolves to **False (union)**. Pre-fix the CLI defaulted
``--no-unite`` to ``dest="unite", default=True`` -- so bare CLI dmc/dmr
*intersected* while bare API *unioned*. This pins the parity fix: bare CLI now
resolves unite to False (union), while ``--unite`` still forces intersect.
"""

from __future__ import annotations

from epykit.cli import build_parser

_DMC_REQUIRED = [
    "--methylstore", "store",
    "--samplesheet", "s.csv",
    "--treatment-group", "t",
    "--control-group", "c",
    "--output", "out.parquet",
]

_DMR_REQUIRED = [
    "--methylstore", "store",
    "--samplesheet", "s.csv",
    "--treatment-group", "t",
    "--control-group", "c",
    "--output", "out.parquet",
    "--method", "tile",
]


def _parse(argv):
    return build_parser().parse_args(argv)


def test_dmc_bare_defaults_to_union():
    """Bare ``epykit dmc`` resolves unite to False (union), matching the API."""
    args = _parse(["dmc", *_DMC_REQUIRED])
    assert args.unite is False, "bare CLI dmc must default to union (unite=False)"


def test_dmc_unite_forces_intersect():
    args = _parse(["dmc", *_DMC_REQUIRED, "--unite"])
    assert args.unite is True, "--unite must force intersect"


def test_dmc_no_unite_is_union():
    args = _parse(["dmc", *_DMC_REQUIRED, "--no-unite"])
    assert args.unite is False, "--no-unite must select union"


def test_dmr_tile_bare_defaults_to_union():
    args = _parse(["dmr", *_DMR_REQUIRED])
    assert args.unite is False, "bare CLI dmr must default to union (unite=False)"


def test_dmr_unite_forces_intersect():
    args = _parse(["dmr", *_DMR_REQUIRED, "--unite"])
    assert args.unite is True, "--unite must force intersect on dmr"


def test_cli_api_bare_unite_parity():
    """The CLI-bare unite default must equal the API-bare resolution.

    API resolution (no pp.unite step): ``md.uns.get("unite") is None`` ->
    ``unite=False``. So both must be False.
    """
    cli_unite = _parse(["dmc", *_DMC_REQUIRED]).unite
    # Mirror tl.dmc's resolution with an unset md.uns["unite"].
    uns: dict = {}
    unite_info = uns.get("unite")
    api_unite = (unite_info is not None) and (unite_info.get("type") == "intersect")
    assert cli_unite == api_unite is False, (
        f"CLI-bare ({cli_unite}) and API-bare ({api_unite}) unite must both "
        f"be False (union)"
    )
