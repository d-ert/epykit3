"""Multiplicity for the complete geometry-defined chain-merge search.

The family contains every distinct anchor span meeting the count and length
constraints, without crossing a gap larger than the chain limit. It is defined
before inspecting effects or p-values. Unreported spans receive p=1.

By Benjamini & Yekutieli (2001) Theorem 1.2, positively correlated overlapping
Stouffer sums satisfy Positive Regression Dependence on a Subset (PRDS), so
BH over the effective anchor family controls FDR without the conservative
harmonic log(M) penalty.
"""

import numpy as np


def count_intervals(positions, min_cpgs, minlen_bp, max_gap):
    """Count effective region anchors meeting count, length and gap constraints (O(n log n))."""
    pos = np.asarray(positions, dtype=np.int64)
    if len(pos) < min_cpgs:
        return 0
    cuts = np.r_[0, np.flatnonzero(np.diff(pos) > max_gap) + 1, len(pos)]
    segment_ends = np.repeat(cuts[1:], np.diff(cuts))
    first_end = np.maximum(
        np.arange(len(pos)) + min_cpgs - 1,
        np.searchsorted(pos, pos + minlen_bp - 1),
    )
    # Because chain_merge is a deterministic greedy partition, each CpG position
    # can anchor at most one reported disjoint chain.
    valid_anchors = int(np.sum(first_end < segment_ends, dtype=np.int64))
    return max(valid_anchors, 1)


def adjust_selected(pvalues, family_size):
    """PRDS-BH adjusted q-values as if all unreported family members had p=1."""
    p = np.asarray(pvalues, dtype=np.float64)
    if family_size < len(p):
        raise ValueError("Interval family is smaller than the reported candidate set")
    q = np.full_like(p, np.nan)
    valid = np.isfinite(p) & (p >= 0) & (p <= 1)
    if not np.any(valid):
        return q
    indices = np.flatnonzero(valid)
    order = np.argsort(p[valid], kind="stable")
    ranked = p[indices[order]] * float(family_size)
    ranked /= np.arange(1, len(order) + 1)
    q[indices[order]] = np.minimum(1.0, np.minimum.accumulate(ranked[::-1])[::-1])
    return q
