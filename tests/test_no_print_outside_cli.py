"""Enforce the load-bearing logging convention from CLAUDE.md.

Library code under ``epykit.*`` emits progress through the stdlib
``logging`` module and never calls :func:`print`. The CLI entry point
(``epykit.cli``) is the only module allowed to use :func:`print`,
where it produces final user-facing result lines on stdout.

This invariant is what lets host applications and notebooks consume
epykit without stdout pollution. A regression here would silently
break any caller that captures stdout.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src" / "epykit"
ALLOWED = {"cli.py"}  # the only module permitted to call print()


def _iter_source_files():
    for path in sorted(SRC_ROOT.rglob("*.py")):
        if path.name in ALLOWED:
            continue
        # __pycache__ etc. excluded by rglob("*.py")
        yield path


def _print_call_lines(path: Path) -> list[int]:
    """Return the line numbers of any top-level print() calls in `path`.

    Uses ast.walk so a print() nested anywhere (function body, method,
    conditional, lambda) is caught.
    """
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    except SyntaxError as exc:  # pragma: no cover - hard fail if any file is broken
        pytest.fail(f"Could not parse {path}: {exc}")
    hits = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            func = node.func
            if isinstance(func, ast.Name) and func.id == "print":
                hits.append(node.lineno)
    return hits


def test_no_print_calls_outside_cli():
    """Every library module under epykit/ must use logging, not print()."""
    offenders = {}
    for path in _iter_source_files():
        hits = _print_call_lines(path)
        if hits:
            offenders[str(path.relative_to(SRC_ROOT.parent.parent))] = hits

    if offenders:
        lines = ["Found print() calls in library code (use `logger.info` instead):"]
        for rel, line_nos in offenders.items():
            lines.append(f"  {rel}: lines {line_nos}")
        lines.append(
            "If a print() is genuinely needed (e.g., a new CLI subcommand), "
            "move the call site into src/epykit/cli.py or add the module to "
            "the ALLOWED set in this test with a CLAUDE.md justification."
        )
        pytest.fail("\n".join(lines))


def test_cli_module_is_present():
    """Sanity-check the ALLOWED list still matches the real CLI module."""
    assert (SRC_ROOT / "cli.py").is_file(), (
        "epykit.cli is the canonical print()-bearing module. If it has been "
        "renamed or split, update ALLOWED in this test."
    )
