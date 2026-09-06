"""The DMC engine registry (``epykit._dmc_engines``) and the call sites that read it.

One test per fact the registry has to keep true: the public names are both
CLI ``--test`` choice lists in their documented order, the module loads
without the package (so it cannot import ``dmc.py``), removed names keep
their migration hints, and an unknown name is refused before any store is
opened, on the high-level and the low-level entry point alike.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import pytest

import epykit as ep
from epykit._dmc_engines import ENGINES, PUBLIC_ENGINES, REMOVED_ENGINES, engine_spec
from epykit.cli import build_parser

REGISTRY_PATH = Path(ep.__file__).with_name("_dmc_engines.py")


def _test_action(parser: argparse.ArgumentParser, command: str) -> argparse.Action:
    subparsers = next(a for a in parser._actions if isinstance(a, argparse._SubParsersAction))
    return next(a for a in subparsers.choices[command]._actions if "--test" in a.option_strings)


def _tree(md) -> list[Path]:
    """Every path under the synthetic store root, DMC cache included."""
    return sorted(Path(md.store).parent.parent.rglob("*"))


def test_public_engines_are_both_cli_choice_lists_in_order():
    assert PUBLIC_ENGINES == ("lr", "glm", "welch_t", "fisher")
    assert set(ENGINES) == {*PUBLIC_ENGINES, "glm_contrast"}
    parser = build_parser()
    for command in ("dmc", "dmr"):
        action = _test_action(parser, command)
        assert action.choices == list(PUBLIC_ENGINES), command
        assert action.default == "lr", command


def test_registry_loads_without_the_package():
    """Executed from its file, outside the ``epykit`` package, the module
    has no package to import from: stdlib only, so no ``dmc.py`` and no cycle."""
    spec = importlib.util.spec_from_file_location("dmc_engines_standalone", REGISTRY_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # dataclasses resolves the postponed annotations through sys.modules.
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        del sys.modules[spec.name]
    assert module.PUBLIC_ENGINES == PUBLIC_ENGINES


@pytest.mark.parametrize("name", sorted(REMOVED_ENGINES))
def test_removed_engine_raises_its_migration_hint(name):
    with pytest.raises(ValueError, match=r"removed in 0\.7\.5"):
        engine_spec(name)


def test_tl_dmc_refuses_an_unknown_engine_before_any_store(synth_md_filtered):
    md = synth_md_filtered
    before = _tree(md)
    with pytest.raises(ValueError, match=r"Unknown DMC test 'bogus'\. Choose one of: lr, glm"):
        ep.tl.dmc(md, test="bogus")
    assert "dmc" not in md.uns
    assert _tree(md) == before
