#!/usr/bin/env python3
"""Verify final repair source, complete runs, tables and all recorded hashes."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20), b''):
            digest.update(block)
    return digest.hexdigest()


def check_hash(path, expected):
    if not path.is_file() or sha256(path) != expected:
        raise ValueError(f'File hash mismatch: {path}')


def verify_frozen_files(directory, info):
    for name, expected in info['source_sha256'].items():
        check_hash(directory / 'source/epykit' / name, expected)
    for name, expected in info['evaluation_sha256_by_file'].items():
        check_hash(directory / 'evaluation_source' / name, expected)
    return dict(source_files=len(info['source_sha256']), evaluation_files=len(info['evaluation_sha256_by_file']))


def verify_run_set(directory, expected_count, source_signature, evaluation_signature):
    manifests = sorted(directory.glob('runs/**/run.json'))
    if len(manifests) != expected_count:
        raise ValueError(f'Run count differs: {directory}: {len(manifests)} != {expected_count}')
    outputs = 0
    for path in manifests:
        state = json.loads(path.read_text())
        if state['status'] != 'complete':
            raise ValueError(f'Run is not complete: {path}')
        identity = state['identity']
        if identity['source_signature'] != source_signature or identity['evaluation_signature'] != evaluation_signature:
            raise ValueError(f'Run identity differs: {path}')
        if not state.get('outputs'):
            raise ValueError(f'No verified outputs: {path}')
        for name, expected in state['outputs'].items():
            check_hash(path.parent / name, expected)
            outputs += 1
    return dict(complete_runs=len(manifests), verified_outputs=outputs)


def verify_count_store(manifest):
    info = json.loads(manifest.read_text())
    signature = hashlib.sha256(json.dumps(info['sha256'], sort_keys=True).encode()).hexdigest()
    if signature != info['signature']:
        raise ValueError(f'Count-store signature differs: {manifest}')
    for name, expected in info['sha256'].items():
        check_hash(Path(info['store']) / name, expected)
    return len(info['sha256'])


def verify_tables(directory):
    import polars as pl
    total = dict(candidates=0, selected=0, null_datasets=0)
    for path in sorted(directory.glob('runs/**/candidates.parquet')):
        candidates = pl.read_parquet(path)
        selected = pl.read_parquet(path.with_name('selected.parquet'))
        expected = candidates.filter(pl.col('combined_qvalue') <= .05)
        if not expected.equals(selected):
            raise ValueError(f'Selected table disagrees with recorded inference: {path}')
        params = json.loads(path.with_name('params.json').read_text())
        for col in ['statistic','standard_error','combined_pvalue','combined_qvalue','maxT_pvalue']:
            if candidates[col].null_count() or not candidates[col].is_finite().all():
                raise ValueError(f'Invalid finite {col}: {path}')
        for col in ['combined_pvalue','combined_qvalue','maxT_pvalue']:
            if len(candidates) and not candidates[col].is_between(0,1).all():
                raise ValueError(f'Invalid probability {col}: {path}')
        if len(candidates) and candidates['standard_error'].min() <= 0:
            raise ValueError(f'Nonpositive standard error: {path}')
        if not params.get('experimental') or not params.get('maxT_available'):
            raise ValueError(f'Inference metadata missing required limitations/diagnostics: {path}')
        lower = params['minimum_full_orbit_two_sided_pvalue']
        if len(candidates) and candidates['maxT_pvalue'].min()+1e-14 < lower:
            raise ValueError(f'Maximum-score diagnostic exceeds exact permutation resolution: {path}')
        metrics = json.loads(path.with_name('metrics.json').read_text())
        n = metrics.get('n_calls', metrics.get('n_dmr_calls', metrics.get('gls_calls')))
        if n != len(selected):
            raise ValueError(f'Metric call count differs: {path}')
        if metrics.get('status','complete') != 'complete':
            raise ValueError(f'Metric failed: {path}')
        total['candidates'] += len(candidates)
        total['selected'] += len(selected)
        total['null_datasets'] += int(metrics.get('all_null', False))
    return total


def verify_final(workspace=ROOT):
    results = {}
    source_signature = None
    for name, expected in [('final_full',72), ('final_real',4), ('final_null_holdout',200)]:
        directory = workspace / f'results/statistical_repair_{name}'
        info = json.loads((directory / 'source_provenance.json').read_text())
        if source_signature is not None and source_signature != info['source_signature']:
            raise ValueError('Final result sets use different inference source')
        source_signature = info['source_signature']
        frozen = verify_frozen_files(directory, info)
        for rel, expected_hash in info['source_sha256'].items():
            check_hash(workspace / 'development/epykit_statistical_repair/src/epykit' / rel, expected_hash)
        runs = verify_run_set(directory, expected, source_signature, info['evaluation_signature'])
        tables = verify_tables(directory)
        count_files = sum(verify_count_store(p) for p in directory.glob('inputs/*.json'))
        count_files += sum(verify_count_store(p) for p in directory.glob('runs/**/input_store.json'))
        for path in directory.glob('runs/**/run.json'):
            identity = json.loads(path.read_text())['identity']
            if 'truth_sha256' in identity:
                for rel, expected_hash in identity['truth_sha256'].items():
                    check_hash(workspace / 'data/pilot' / identity['dataset'] / rel, expected_hash)
            if 'suite_sha256' in identity:
                dataset = workspace / 'data/statistical_repair_null_holdout' / identity['dataset']
                check_hash(dataset / 'manifest.json', identity['manifest_sha256'])
                check_hash(dataset.parent / 'suite.json', identity['suite_sha256'])
        results[name] = dict(commit=info['commit'], **frozen, **runs, **tables, count_files=count_files)
    return dict(status='verified', source_signature=source_signature, result_sets=results)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify_final(args.workspace)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
