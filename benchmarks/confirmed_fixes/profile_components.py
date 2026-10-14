"""Time nested discovery/fit components while running the frozen workspace harness."""

import argparse
import cProfile
import importlib
import json
import runpy
import sys
import time
from functools import wraps
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('harness_args', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    sys.path.insert(0, str(args.repo.resolve()/'src'))
    sys.path.insert(1, str(args.workspace.resolve()))
    stages = {}
    counts = {'complete_regions': 0, 'masked_regions': 0}

    def instrument(module_name, function_name):
        module = importlib.import_module(module_name)
        original = getattr(module, function_name)

        @wraps(original)
        def timed(*a, **kw):
            if function_name == 'fit_regions':
                mask = kw.get('obs_mask')
                for lo, hi in a[4]:
                    key = 'complete_regions' if mask is None or mask[:, lo:hi].all() else 'masked_regions'
                    counts[key] += 1
            start = time.perf_counter()
            try:
                return original(*a, **kw)
            finally:
                stages[function_name] = stages.get(function_name, 0.)+time.perf_counter()-start
        setattr(module, function_name, timed)

    for module, name in [('epykit._regional_gls', '_prepare_label_scan'),
                         ('epykit._regional_gls_native', 'fit_regions'),
                         ('epykit._regional_gls_discovery', '_smooth_tricube'),
                         ('epykit._regional_gls', 'run_regional_gls')]:
        instrument(module, name)
    args.output.mkdir(parents=True, exist_ok=True)
    extra = args.harness_args[1:] if args.harness_args[:1] == ['--'] else args.harness_args
    sys.argv = ['run_statistical_repair_benchmarks.py', '--repo', str(args.repo.resolve()),
                '--output', str(args.output/'runs'), *extra]
    profiler = cProfile.Profile()
    profiler.enable()
    try:
        runpy.run_path(str(args.workspace/'scripts/run_statistical_repair_benchmarks.py'), run_name='__main__')
    except SystemExit as exc:
        if exc.code:
            raise
    finally:
        profiler.disable()
        profiler.dump_stats(str(args.output/'python.prof'))
        (args.output/'stages.json').write_text(json.dumps(dict(seconds=stages, counts=counts,
            semantics='Nested inclusive times; cold compilation included. Warm repeats required for optimization decisions.'), indent=2)+'\n')


if __name__ == '__main__':
    main()
