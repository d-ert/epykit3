"""D10: ``epykit convert`` CpG-merge default parity with the Python API.

The API (``convert_sample`` / ``read_bismark``) defaults ``merge_strands=True``.
The CLI ``--merge-cpg`` / ``--no-merge-cpg`` is a tri-state argparse flag with a
``None`` sentinel default. Pre-fix, ``_cmd_convert`` forwarded that ``None``
straight through, so a bare ``epykit convert`` did NOT merge CpG dyads while the
equivalent API call did -- a silent CLI/API divergence.

The fix resolves the ``None`` sentinel to ``True`` in ``_cmd_convert`` so the CLI
matches the API. Explicit ``--merge-cpg`` / ``--no-merge-cpg`` still force
True / False respectively.
"""

from __future__ import annotations

import argparse


def _capture_merge_strands(monkeypatch, *, merge_cpg) -> object:
    """Run ``_cmd_convert`` with a captured ``convert_sample`` and return the
    ``merge_strands`` value it was called with."""
    import epykit.cli as cli

    captured: dict = {}

    def _fake_convert_sample(*args, **kwargs):
        captured.update(kwargs)
        captured["_positional"] = args
        return None

    monkeypatch.setattr(cli, "convert_sample", _fake_convert_sample)

    args = argparse.Namespace(
        input="in.cov",
        sample_id="s1",
        output_dir="out",
        context="CpG",
        reference_fasta=None,
        merge_cpg=merge_cpg,
        format="bismark",
    )
    cli._cmd_convert(args)
    return captured["merge_strands"]


def test_bare_convert_defaults_to_merge(monkeypatch):
    """Bare ``epykit convert`` (no merge flag -> argparse None) must resolve to
    ``merge_strands=True``, matching the API default."""
    assert _capture_merge_strands(monkeypatch, merge_cpg=None) is True


def test_no_merge_cpg_disables_merge(monkeypatch):
    """``--no-merge-cpg`` (argparse False) must forward ``merge_strands=False``."""
    assert _capture_merge_strands(monkeypatch, merge_cpg=False) is False


def test_explicit_merge_cpg_enables_merge(monkeypatch):
    """``--merge-cpg`` (argparse True) must forward ``merge_strands=True``."""
    assert _capture_merge_strands(monkeypatch, merge_cpg=True) is True
