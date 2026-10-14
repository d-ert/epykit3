#!/usr/bin/env python3
"""Package committed repair code and existing evidence; never execute inference."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20), b''):
            h.update(block)
    return h.hexdigest()


def copy_file(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    workspace, repo, output = args.workspace.resolve(), args.repo.resolve(), args.output.resolve()
    if output.exists():
        raise ValueError('Preserve existing bundles; choose a new directory')
    if subprocess.check_output(['git','-C',str(repo),'status','--porcelain'], text=True):
        raise ValueError('Commit the handoff before packaging')
    commit = subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip()
    output.mkdir(parents=True)
    tracked = subprocess.check_output(['git','-C',str(repo),'ls-files','-z'],text=True).split('\0')
    for name in filter(None,tracked):
        if name.startswith(('src/','tests/','benchmarks/statistical_repair/','docs/reviews/','docs/superpowers/')) or name in ['docs/STATISTICAL_REPAIR_HANDOFF.md','pyproject.toml','uv.lock','README.md']:
            copy_file(repo/name, output/'repair_source'/name)
    for folder in ['scripts','wgbs_sim','configs']:
        for path in (workspace/folder).rglob('*'):
            if path.is_file() and path.suffix in {'.py','.json','.tsv'} and '__pycache__' not in path.parts:
                copy_file(path, output/'benchmark_workspace'/path.relative_to(workspace))
    copy_file(workspace/'pyproject.toml', output/'benchmark_workspace/pyproject.toml')
    copy_file(workspace/'argon_chat.txt', output/'reference/argon_chat.txt')
    copy_file(repo/'docs/STATISTICAL_REPAIR_HANDOFF.md', output/'FINAL_REPORT.md')
    result_sets=['candidate_full','candidate_real','candidate_null_development','final_full','final_real','final_null_holdout']
    for label in result_sets:
        source=workspace/'results'/('statistical_repair_'+label)
        target=output/'results'/label
        for name in ['source_provenance.json','aggregate.json','metrics.json','per_run.tsv','per_run.json','summary.json']:
            if (source/name).is_file():copy_file(source/name,target/name)
        for manifest in source.glob('runs/**/run.json'):
            state=json.loads(manifest.read_text())
            names=['run.json','metrics.json']
            if state['status']=='complete':
                names+=['selected.parquet','candidates.parquet','dmc_selected.parquet','gene_comparison.tsv','scenario_recovery.json']
            for name in names:
                path=manifest.with_name(name)
                if path.is_file():copy_file(path,target/path.relative_to(source))
        # Preserve the actual executed code identity, including the pre-fix validator.
        for path in (source/'evaluation_source').rglob('*.py'):
            copy_file(path,target/path.relative_to(source))
    with (output/'epykit_statistical_repair.patch').open('wb') as stream:
        subprocess.run(['git','-C',str(repo),'diff','--binary','0bc083e','HEAD'],stdout=stream,check=True)
    subprocess.run(['git','-C',str(repo),'bundle','create',str(output/'epykit_statistical_repair.bundle'),'stats-repair-20261014'],check=True)
    subprocess.run(['git','-C',str(repo),'bundle','verify',str(output/'epykit_statistical_repair.bundle')],check=True,capture_output=True)
    (output/'START_HERE.md').write_text(f'''# Epykit repair handoff

Committed branch: `stats-repair-20261014`, HEAD `{commit}`.
Last inference change: `d91378f`; later commits preserve tooling and evidence.
Read FINAL_REPORT.md, then repair_source/docs/reviews/ and the tests.

No benchmark should be restarted automatically. All three benchmark jobs were
stopped by the user. Candidate runs completed, final-source sweeps are partial.
No nominal5% FDR or universally superior recovery claim is supported.

The source tree and tests are in repair_source/. The exact frozen executed
harness dependencies are under results/*/evaluation_source where available;
benchmark_workspace/ contains the current fixed harness. Successful preserved
outputs were integrity-checked independently of the old validator failure path.
Raw count stores, full DMC tables, truth input files and large original raw data
remain in the original workspace; source manifests retain their hashes/paths.
Results here include candidate/selected region tables, run status and summaries.

To recover the full committed branch into a new directory without changing your
existing checkout:

```sh
git clone -b stats-repair-20261014 epykit_statistical_repair.bundle epykit_repaired
```

Alternatively review/apply the supplied patch against the matching starting
source. The external ../epykit3 checkout was not modified or pushed.

Suggested external review: audit the observed-cell GLS covariance/design and
mirror symmetry, depth-relative EB reference, adaptive multiscale discovery,
complete-scan FDP estimation and maxT assumptions, plus remaining adjacent-tile
heuristics. Separate software correctness from inferential calibration. Do not
use paper concordance as ground truth or partial null runs as completed holdout
validation. Argon advice was followed in intent; invalid formulas were corrected.
''')
    files={str(p.relative_to(output)):digest(p) for p in sorted(output.rglob('*')) if p.is_file()}
    (output/'inventory.json').write_text(json.dumps(dict(commit=commit,sha256=files),indent=2)+'\n')
    for name, expected in files.items():
        if digest(output/name)!=expected:raise ValueError(f'Packaged hash differs: {name}')
    print(json.dumps(dict(output=str(output),commit=commit,verified_files=len(files))))


if __name__=='__main__':
    main()
