"""Dependency-floor guards.

These tests pin API features the codebase actually depends on, so a wrong
lower-bound in pyproject.toml (which imports fine but breaks at call time)
cannot pass silently. See M-PKG1 in the pre-submission review.
"""

from __future__ import annotations

import polars as pl


def test_polars_pivot_on_keyword_available():
    """epykit uses ``DataFrame.pivot(on=...)`` (qc.py, pl/_compute.py, tl.py).

    The ``on=`` keyword replaced ``columns=`` in polars 1.0.0; a <1.0 install
    raises ``TypeError`` here. This is the floor declared in pyproject.toml
    (``polars>=1.0``).
    """
    df = pl.DataFrame({"i": [1, 1], "k": ["a", "b"], "v": [1.0, 2.0]})
    wide = df.pivot(values="v", index="i", on="k")
    assert "a" in wide.columns and "b" in wide.columns
    assert wide.height == 1
