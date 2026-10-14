"""Opt-in, subprocess-isolated timing of fixed-seed confirmed-fix fixtures."""

import argparse
import hashlib
import json
import os
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np


def scalar_tss(chroms, mids, index):
    names = np.full(len(mids), '', object)
    distances = np.full(len(mids), np.iinfo(np.int32).min, np.int64)
    for i, chrom in enumerate(chroms):
        pos, genes, strands = index[chrom]
        center = int(mids[i])
        j = int(np.searchsorted(pos, center, side='left'))
        candidates = ([j-1] if j else []) + ([j] if j < len(pos) else [])
        best = min(candidates, key=lambda k: abs(int(pos[k])-center))
        names[i] = str(genes[best])
        distances[i] = (-1 if strands[best] == '-' else 1)*(center-int(pos[best]))
    return names, distances


def fixture(component, size, samples, k, dropout, invalid):
    rng = np.random.default_rng(84012)
    if component == 'glm':
        from epykit._glm import irls_binomial_batch
        x = np.c_[np.ones(10), np.repeat([0, 1], 5), np.linspace(-1, 1, 10)]
        n = np.full((size, 10), 50)
        m = rng.binomial(n, .4)
        # One nonconvergent site keeps baseline full-batch IRLS running.
        m[-1] = 0
        if invalid:
            n[0, 5:] = m[0, 5:] = 0
        return lambda: irls_binomial_batch(m, n, x, return_cov=True), [x, m, n]
    if component == 'tss':
        import epykit.annotate as annotation
        chroms = np.array(['chr1', 'chr2', 'chr3', 'chr4'])[np.arange(size) % 4]
        mids = rng.integers(0, 100_000_000, size).astype(float)+.5
        index = {chrom: (np.sort(rng.integers(0, 100_000_000, 5000)),
                         np.array([f'g{i}' for i in range(5000)]),
                         np.where(np.arange(5000) % 2, '-', '+')) for chrom in np.unique(chroms)}
        function = getattr(annotation, '_nearest_tss_arrays', scalar_tss)
        return lambda: function(chroms, mids, index), [chroms, mids, *[v for row in index.values() for v in row]]
    from epykit._regional_gls_masked import _spectral_prepare
    p = np.cumsum(rng.integers(5, 50, k)).astype(float)
    y = rng.normal(size=(samples, k))
    h = np.sqrt(rng.integers(5, 50, k)).astype(float)
    g = np.r_[np.zeros(samples//2), np.ones(samples-samples//2)]
    mask = rng.random((samples, k)) >= dropout
    mask[[0, samples//2], :] = True
    return lambda: _spectral_prepare(p, y, h, g, mask, 120.), [p, y, h, g, mask]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--component', choices=['glm', 'tss', 'masked'], required=True)
    parser.add_argument('--size', type=int, default=100000)
    parser.add_argument('--samples', type=int, default=6)
    parser.add_argument('--k', type=int, default=100)
    parser.add_argument('--dropout', type=float, default=.2)
    parser.add_argument('--invalid', action='store_true')
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--inner', type=int, default=1)
    args = parser.parse_args()
    sys.path.insert(0, str(args.repo.resolve()/'src'))
    import numba

    import epykit
    numba.set_num_threads(args.workers)
    function, arrays = fixture(args.component, args.size, args.samples, args.k, args.dropout, args.invalid)
    fingerprint = hashlib.sha256()
    for array in arrays:
        fingerprint.update(str(array.shape).encode())
        fingerprint.update(array.tobytes())
    start = time.perf_counter()
    function()
    cold = time.perf_counter()-start
    elapsed = []
    for _ in range(args.repeats):
        start = time.perf_counter()
        for _ in range(args.inner):
            function()
        elapsed.append((time.perf_counter()-start)/args.inner)
    package = args.repo.resolve()/'src/epykit'
    source = hashlib.sha256(b''.join(p.read_bytes() for p in sorted(package.glob('*.py')))).hexdigest()
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result = dict(component=args.component, dimensions=dict(size=args.size, samples=args.samples, k=args.k),
                  dropout=args.dropout, invalid=args.invalid, workers=args.workers,
                  commit=subprocess.check_output(['git', '-C', str(args.repo), 'rev-parse', 'HEAD'], text=True).strip(),
                  source_hash=source, input_hash=fingerprint.hexdigest(), import_path=epykit.__file__,
                  platform=platform.platform(), processor=platform.processor(), python=sys.version,
                  numpy=np.__version__, numba=numba.__version__, blas_threads={key: os.getenv(key) for key in
                  ['OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'OMP_NUM_THREADS']},
                  cold_seconds=cold, warm_seconds=elapsed, median_seconds=float(np.median(elapsed)),
                  peak_rss_bytes=rss if sys.platform == 'darwin' else rss*1024)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
