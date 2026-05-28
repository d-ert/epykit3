"""Phase 3 cleanup: dropped engines must raise ValueError with a
migration hint pointing the user at a surviving engine."""
from __future__ import annotations

import pytest

import epykit as ep


@pytest.mark.parametrize("engine,hint_substring", [
    ("logit_t", "welch_t"),
    ("bb_lr",   "lr"),
    ("score",   "lr"),
    ("cmh",     "formula='~ group + batch'"),
])
def test_dropped_engine_raises_with_migration_hint(
    synth_md_filtered, engine, hint_substring,
):
    """Each dropped engine raises ValueError; the message includes
    text pointing at the recommended replacement."""
    with pytest.raises(ValueError) as exc:
        ep.tl.dmc(synth_md_filtered, test=engine)
    msg = str(exc.value)
    assert "removed in 0.7.5" in msg, f"missing version note in: {msg}"
    assert hint_substring in msg, (
        f"missing migration hint '{hint_substring}' in: {msg}"
    )
