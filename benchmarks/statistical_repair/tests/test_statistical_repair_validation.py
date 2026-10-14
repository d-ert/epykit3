"""Independent all-null readouts must include every declared completed run."""

import pytest


def test_null_readout_rejects_incomplete_or_nonnull_evidence():
    from scripts.run_statistical_repair_validation import null_validation_summary

    row = dict(config_id="a", track="empirical", status="complete", all_null=True, gls_calls=0, dmc_calls=0)
    with pytest.raises(ValueError, match="complete"):
        null_validation_summary([row], expected=2)
    with pytest.raises(ValueError, match="null"):
        null_validation_summary([{**row, "all_null": False}], expected=1)


def test_null_readout_uses_any_discovery_not_average_call_count():
    from scripts.run_statistical_repair_validation import null_validation_summary

    rows = [dict(config_id="a", track="empirical", status="complete", all_null=True, gls_calls=n, dmc_calls=0)
            for n in [0, 4, 0, 0]]
    result = null_validation_summary(rows, expected=4)
    assert result["a/empirical"]["gls"]["probability_any_discovery"] == .25
    assert result["a/empirical"]["dmc"]["probability_any_discovery"] == 0
