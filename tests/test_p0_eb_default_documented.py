"""P0-3: docstring and code default for `dispersion` must agree."""
from __future__ import annotations

import inspect
import epykit as ep


def test_dispersion_default_is_eb():
    sig = inspect.signature(ep.tl.dmc)
    assert sig.parameters["dispersion"].default == "eb", (
        "Code default for `dispersion` changed. If intentional, "
        "update PROTOCOL.md and the EXECUTIVE_SUMMARY as well."
    )


def test_dispersion_docstring_mentions_eb_as_default():
    import re
    doc = ep.tl.dmc.__doc__ or ""
    assert "eb" in doc, "Docstring no longer mentions the 'eb' option."
    assert re.search(r'default\s+``"eb"``', doc, re.IGNORECASE), (
        "Docstring must state 'eb' is the default (looking for "
        "'default ``\"eb\"``' case-insensitively)."
    )
    # The old wrong claim must be gone.
    assert not re.search(r'default\s+``"site"``', doc, re.IGNORECASE), (
        "Docstring still says default is 'site'; should be 'eb'."
    )
