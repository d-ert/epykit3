"""Delivery verification must prove complete immutable result sets."""

import importlib.util
import json
from pathlib import Path

import pytest


def delivery_module():
    path = Path(__file__).resolve().parents[1] / 'development/repair_delivery/verify_results.py'
    spec = importlib.util.spec_from_file_location('repair_delivery_verifier', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_delivery_rejects_incomplete_or_changed_result_set(tmp_path):
    module = delivery_module()
    run = tmp_path / 'runs/pilot/dataset/config'
    run.mkdir(parents=True)
    value = run / 'metrics.json'
    value.write_text('{}\n')
    manifest = run / 'run.json'
    identity = dict(source_signature='one', evaluation_signature='two')
    state = dict(status='complete', identity=identity, outputs={'metrics.json':module.sha256(value)})
    manifest.write_text(json.dumps(state))
    assert module.verify_run_set(tmp_path, 1, 'one', 'two')['complete_runs'] == 1
    with pytest.raises(ValueError, match='count'):
        module.verify_run_set(tmp_path, 2, 'one', 'two')
    manifest.write_text(json.dumps({**state, 'status':'failed'}))
    with pytest.raises(ValueError, match='complete'):
        module.verify_run_set(tmp_path, 1, 'one', 'two')
    manifest.write_text(json.dumps(state))
    value.write_text('changed\n')
    with pytest.raises(ValueError, match='hash'):
        module.verify_run_set(tmp_path, 1, 'one', 'two')


def test_delivery_rejects_changed_frozen_source(tmp_path):
    module = delivery_module()
    package = tmp_path / 'source/epykit'
    package.mkdir(parents=True)
    source = package / '__init__.py'
    source.write_text('version=1\n')
    evaluator = tmp_path / 'evaluation_source/scripts/runner.py'
    evaluator.parent.mkdir(parents=True)
    evaluator.write_text('original\n')
    info = {'source_sha256':{'__init__.py':module.sha256(source)},
            'evaluation_sha256_by_file':{'scripts/runner.py':module.sha256(evaluator)}}
    assert module.verify_frozen_files(tmp_path, info)['source_files'] == 1
    evaluator.write_text('altered\n')
    with pytest.raises(ValueError, match='hash'):
        module.verify_frozen_files(tmp_path, info)
