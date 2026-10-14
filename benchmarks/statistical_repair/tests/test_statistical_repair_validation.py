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


def test_null_final_verification_failure_prevents_successful_summary(tmp_path):
    import importlib.util
    import json
    from pathlib import Path
    from types import SimpleNamespace
    from unittest.mock import patch

    pl = pytest.importorskip('polars')
    ep = pytest.importorskip('epykit')
    from epykit import _regional_gls
    # Repair the archived harness while running processes retain executed root code.
    workspace = Path(__file__).resolve().parents[1]
    path = workspace / 'development/epykit_statistical_repair/benchmarks/statistical_repair/scripts/run_statistical_repair_validation.py'
    spec = importlib.util.spec_from_file_location('isolated_null_validator', path)
    validator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(validator)
    data, output = tmp_path / 'data', tmp_path / 'out'
    name = 'empirical_null_seed1'
    (data / name).mkdir(parents=True)
    (data / 'suite.json').write_text(json.dumps({'runs':[dict(dataset=name, track='empirical', seed=1,
                                                            all_null=True, validated=True)]}))
    (data / name / 'manifest.json').write_text('{}')
    sheet = output / 'derived_counts' / name / 'samplesheet.csv'
    sheet.parent.mkdir(parents=True)
    sheet.write_text('fixture')
    fake_md = SimpleNamespace(obs=None, store=str(tmp_path / 'store'), dmc=pl.DataFrame({'qvalue':[.9]}))
    candidate = pl.DataFrame({'combined_qvalue':[.9], 'maxT_pvalue':[.9]})
    source = dict(source_signature='fixture', evaluation_signature='fixture', commit='fixture')
    with (
        patch.object(validator, 'ROOT', workspace),
        patch.object(validator, 'snapshot_source', return_value=source),
        patch.object(validator, 'store_identity', return_value='fixture'),
        patch.object(validator, 'verify_evaluation', side_effect=[None, ValueError('Live evaluation dependency changed: fixture')]),
        patch.object(ep, 'read_bismark', return_value=fake_md),
        patch.object(ep, 'MethylData', return_value=fake_md),
        patch.object(ep, 'set_tmp_dir'),
        patch.object(ep.pp, 'filter_coverage'),
        patch.object(ep.pp, 'set_unite_type'),
        patch.object(ep.tl, 'dmc'),
        patch.object(_regional_gls, 'run_regional_gls', return_value=(candidate, candidate, {'null_regions':0})),
    ):
        with pytest.raises(ValueError, match='complete'):
            validator.main(['--dataset-root',str(data),'--output',str(output),'--configs','cov10_boundNone'])
    manifest = output / 'runs' / name / 'cov10_boundNone' / 'run.json'
    assert json.loads(manifest.read_text())['status'] == 'failed'
    assert json.loads((output / 'per_run.json').read_text())[0]['status'] == 'failed'
    assert not (output / 'summary.json').exists()
