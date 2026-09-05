"""Dependency-floor guards.

These tests pin API features the codebase actually depends on, so a wrong
lower-bound in pyproject.toml (which imports fine but breaks at call time)
cannot pass silently. See M-PKG1 in the pre-submission review.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
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


def test_pandas_nullable_int_map_to_numpy_na_value():
    """annotate.py maps gene ids through a nullable ``Int64`` TSS lookup and
    materialises it with ``to_numpy(dtype=float64, na_value=nan)``.

    pandas is imported directly (annotate.py, filter.py, pl/composer.py) and
    declared as ``pandas>=2.0`` in pyproject.toml; this is the call pattern
    that floor stands behind.
    """
    tss = pd.Series([100, 200], index=["g1", "g2"], dtype="Int64")
    out = (
        pd.Series(["g1", "absent", "g2"])
        .map(tss)
        .to_numpy(dtype=np.float64, na_value=np.nan)
    )
    assert out[0] == 100.0 and np.isnan(out[1]) and out[2] == 200.0
