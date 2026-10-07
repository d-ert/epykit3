"""Multiplicity for the complete geometry-defined chain-merge search.

The family contains every contiguous CpG span meeting the count and length
constraints, without crossing a gap larger than the chain limit. It is defined
before inspecting effects or p-values. Unreported spans receive p=1.

BY accommodates dependence between overlapping interval hypotheses and the
data-dependent selection of reported spans, provided each fixed interval's
raw p-value is valid. It cannot repair miscalibrated CpG p-values or within-span
correlation ignored by signed Stouffer combination.
"""

import numpy as np
from scipy.special import digamma


def count_intervals(positions, min_cpgs, minlen_bp, max_gap):
    """Count admissible inclusive spans without enumerating them (O(n log n))."""
    pos = np.asarray(positions, dtype=np.int64)
    if not len(pos):
        return 0
    # Each gap-connected segment has its own allowed right boundary.
    cuts = np.r_[0, np.flatnonzero(np.diff(pos) > max_gap) + 1, len(pos)]
    segment_ends = np.repeat(cuts[1:], np.diff(cuts))
    first_end = np.maximum(
        np.arange(len(pos)) + min_cpgs - 1,
        np.searchsorted(pos, pos + minlen_bp - 1),
    )
    return int(np.maximum(segment_ends - first_end, 0).sum(dtype=np.int64))


def adjust_selected(pvalues, family_size):
    """BY adjusted values as if all unreported family members had p=1.

    Only selected spans need materialization. The harmonic factor uses the
    complete family size; computing it with digamma avoids allocating a vector
    proportional to the number of possible intervals.
    """
    p = np.asarray(pvalues, dtype=np.float64)
    if family_size < len(p):
        raise ValueError("Interval family is smaller than the reported candidate set")
    q = np.full_like(p, np.nan)
    valid = np.isfinite(p) & (p >= 0) & (p <= 1)
    if not np.any(valid):
        return q
    indices = np.flatnonzero(valid)
    order = np.argsort(p[valid], kind="stable")
    harmonic = float(digamma(float(family_size) + 1) + np.euler_gamma)
    ranked = p[indices[order]] * (float(family_size) * harmonic)
    ranked /= np.arange(1, len(order) + 1)
    q[indices[order]] = np.minimum(1., np.minimum.accumulate(ranked[::-1])[::-1])
    return q
