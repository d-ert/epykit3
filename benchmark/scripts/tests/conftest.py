"""Pytest configuration for the benchmark-script test suite.

These tests are independent of the main epykit test suite — run them via
`uv run pytest benchmark/scripts/tests/` from the repo root. They exercise
the simulator, null-calibration runner, and CI helpers without touching
the epykit package internals.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Make the scripts directory importable as a flat package so tests can do
# `from simulate_piao import simulate_dmc` rather than messing with PYTHONPATH.
SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
