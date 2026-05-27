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
    doc = ep.tl.dmc.__doc__ or ""
    # Both: 'eb' appears in the choices set AND default points to it.
    assert "eb" in doc, "Docstring no longer mentions the 'eb' option."
    assert 'Default ``"eb"' in doc or 'default ``"eb"' in doc, (
        "Docstring must state 'eb' is the default."
    )
    # The old wrong claim must be gone.
    assert 'Default ``"site"' not in doc, (
        "Docstring still says default is 'site'; should be 'eb'."
    )
