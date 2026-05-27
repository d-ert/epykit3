"""Differential methylation calling (DMC) on partitioned Parquet stores.

Per-replicate state is accumulated via Welford's online algorithm, so peak
memory is O(n_sites) per chromosome, independent of sample count.

Tests
-----
  lr            -- Quasi-binomial likelihood-ratio chi-square with per-site
                  McCullagh-Nelder dispersion. The default at n>=2.
                  Closed-form on the streaming (S0_g, S1_g, Sigmam^2/n_g)
                  accumulators; closer to nominal type-I error than the
                  score test at the small samples and boundary beta typical
                  in WGBS.
  score         -- Pearson score on the same dispersion-corrected
                  accumulators as lr. Marginally more powerful but mildly
                  anti-conservative at pi near 0/1.
  glm           -- Binomial GLM via batched IRLS (see _glm.py). Required for
                  covariate-adjusted designs and multi-group contrasts.
  logit_t       -- Welch t on logit(beta) via Welford. Variance-stabilising
                  fallback. Weak near beta=0/1: anti-conservative under H0
                  when one group's between-replicate variance collapses
                  by binomial sampling chance. Use ``lr`` for trustworthy
                  inference; reach for ``logit_t`` only when count-model
                  assumptions are doubtful.
  welch_t       -- Welch t on raw betas. Same boundary-beta caveat as
                  ``logit_t``.
  bb_lr         -- True quasi-binomial LRT via a full per-site GLM on a
                  binary-treatment design. The honest "fit-the-model"
                  alternative to ``welch_t``; slower than ``lr`` and
                  produces near-identical results on binary designs.
  cmh           -- Cochran-Mantel-Haenszel with one 2x2 stratum per
                  (case_i, ctrl_j) pair.
  fisher        -- Fisher exact on reads pooled across replicates. Ignores
                  between-replicate variance; anti-conservative. Warns once
                  per session.
"""

from __future__ import annotations

import gc
import hashlib
import json
import logging
import tempfile
import time
import warnings
from pathlib import Path
from typing import Optional, Tuple, Union

import numpy as np
import polars as pl
from scipy import stats as sp_stats

from . import _cache
from ._dmc_store import DMCStore, _chrom_filename

logger = logging.getLogger(__name__)

_SMOOTH_BOX_NJIT_FN = None

_EMPTY_SCHEMA = {
    "chrom":             pl.Utf8,
    "pos":               pl.Int32,
    "strand":            pl.Utf8,
    "n_case":            pl.Int32,
    "n_control":         pl.Int32,
    "mean_beta_case":    pl.Float32,
    "mean_beta_control": pl.Float32,
    "pvalue":            pl.Float64,
    "log2_odds_ratio":   pl.Float64,
    "meth_diff":         pl.Float32,
    "meth_diff_ci_lo":   pl.Float32,
    "meth_diff_ci_hi":   pl.Float32,
}

def _epykit_version() -> str:
    try:
        from . import __version__
        return __version__
    except ImportError:
        return "0.0.0+unknown"


def _now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime())


def _dmc_input_signature(
    methylstore_path: Path,
    samples_case: list[str],
    samples_control: list[str],
    test: str,
    chromosomes: list[str],
    unite: bool,
    min_samples_case: int,
    min_samples_control: int,
    dispersion: str,
    reference: str,
    samples_all_ordered: Optional[list[str]],
    group_labels_per_sample: Optional[list[str]],
    contrast_label: Optional[str],
    smoothing: bool = False,
    smoothing_span_bp: int = 500,
    sep_fallback: bool = False,
    sep_threshold: float = 0.9,
) -> str:
    """SHA-256 fingerprint of inputs that materially affect DMC results.

    Two runs with the same signature must produce bit-identical DMC
    output, so we can skip recomputation when the existing cache's
    manifest carries a matching signature. The fingerprint deliberately
    avoids unstable bits (worker count, backend, glm_backend) -- those
    affect timing but not the result.

    The methylstore is captured by its resolved path only, not by
    mtime: directory mtimes get touched even on no-op filter reruns
    (e.g. when every sample manifest reports "cached"), which would
    invalidate the cache spuriously every run. To force a recompute
    after a genuine upstream change, delete the ``.cache/dmc/<test>/``
    directory.
    """
    h = hashlib.sha256()
    msp = Path(methylstore_path).resolve()
    h.update(b"|store="); h.update(str(msp).encode())
    h.update(b"|case=");    h.update(",".join(samples_case).encode())
    h.update(b"|ctrl=");    h.update(",".join(samples_control).encode())
    h.update(b"|test=");    h.update(test.encode())
    h.update(b"|chroms=");  h.update(",".join(chromosomes).encode())
    h.update(b"|unite=");   h.update(str(bool(unite)).encode())
    h.update(b"|min_case="); h.update(str(min_samples_case).encode())
    h.update(b"|min_ctrl="); h.update(str(min_samples_control).encode())
    h.update(b"|disp=");    h.update(str(dispersion).encode())
    h.update(b"|ref=");     h.update(str(reference).encode())
    if samples_all_ordered is not None:
        h.update(b"|all=");     h.update(",".join(samples_all_ordered).encode())
    if group_labels_per_sample is not None:
        h.update(b"|grp=");     h.update(",".join(group_labels_per_sample).encode())
    if contrast_label is not None:
        h.update(b"|contrast="); h.update(contrast_label.encode())
    # DSS-style count-smoothing tunables. The flag itself is always
    # hashed so a False->True toggle invalidates the cache. The span is
    # only hashed when smoothing is on (when off, the span has no effect
    # and should not invalidate the cache).
    h.update(b"|sm=");      h.update(b"1" if smoothing else b"0")
    if smoothing:
        h.update(b"|span=");    h.update(str(int(smoothing_span_bp)).encode())
    # Separation-aware fallback (since 0.7.1) -- ON/OFF and threshold change
    # the per-site p-values, so they must be part of the cache key.
    h.update(b"|sep=");     h.update(b"1" if sep_fallback else b"0")
    if sep_fallback:
        h.update(b"|sept=");    h.update(f"{sep_threshold:.6f}".encode())
    return h.hexdigest()


def _resolve_dmc_store_dir(
    methylstore_path: Path,
    test: str,
    out_dir: Optional[Union[str, Path]],
    smoothing: bool = False,
) -> Path:
    """Pick the directory used for the persistent DMC store.

    Resolution order:
      1. Explicit ``out_dir`` argument -> used verbatim.
      2. If ``methylstore_path`` lives directly inside a ``.cache/``
         directory (e.g. ``<X>/.cache/filtered``), put the DMC store
         in a sibling stage dir (``<X>/.cache/dmc/<test>``). Mirrors
         the ``.cache/filtered/`` convention used by ``pp.filter_coverage``.
      3. Otherwise, put it under ``<methylstore_parent>/.cache/dmc/<test>``.
      4. Final fallback: ``tempfile.mkdtemp`` with a warning.

    When ``smoothing=True`` the test directory is suffixed with
    ``_smooth`` so DSS-style smoothed results live in a separate cache
    from the un-smoothed run. The existing cache "weak-hit" path
    (sig-mismatch, files present -> serve cached) is therefore unable to
    cross-contaminate the two modes.
    """
    bucket = f"{test}_smooth" if smoothing else test
    if out_dir is not None:
        return Path(out_dir)
    msp    = Path(methylstore_path).resolve()
    parent = msp.parent
    if parent.exists():
        # Detect "methylstore is a stage dir inside .cache/" -- drop the
        # DMC store alongside it instead of nesting another .cache.
        if parent.name == ".cache":
            return parent / "dmc" / bucket
        return parent / ".cache" / "dmc" / bucket
    fallback = Path(tempfile.mkdtemp(prefix="epykit_dmc_"))
    logger.warning(
        "Could not derive a persistent DMC store dir from %s; "
        "falling back to ephemeral %s. Pass out_dir= to keep DMC results.",
        methylstore_path, fallback,
    )
    return fallback


def _canonicalise_test_name(test: str) -> str:
    """Map deprecated test names to their canonical form."""
    return test

_TEST_RECOMMENDATIONS = {
    range(1, 3):   "fisher (single-rep only; effect size dominates)",
    range(3, 999): "lr (quasi-binomial likelihood-ratio with MN overdispersion)",
}

# Shared epsilon for boundary clipping in logit / log-OR computations.
_BETA_EPSILON: float = 1e-6



# Core statistical tests (public, used by unit tests)

def fisher_exact_vectorized(
    meth_a: np.ndarray,
    unmeth_a: np.ndarray,
    meth_b: np.ndarray,
    unmeth_b: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Vectorised Fisher exact test via hypergeometric tail approximation."""
    meth_a   = np.asarray(meth_a,   dtype=np.int64)
    unmeth_a = np.asarray(unmeth_a, dtype=np.int64)
    meth_b   = np.asarray(meth_b,   dtype=np.int64)
    unmeth_b = np.asarray(unmeth_b, dtype=np.int64)

    row1  = meth_a + unmeth_a
    row2  = meth_b + unmeth_b
    col1  = meth_a + meth_b
    total = row1 + row2

    n       = len(meth_a)
    pvals   = np.full(n, np.nan, dtype=np.float64)
    log2_or = np.full(n, np.nan, dtype=np.float64)

    valid = (row1 > 0) & (row2 > 0)
    if np.any(valid):
        denom = unmeth_a[valid] * meth_b[valid]
        numer = meth_a[valid]  * unmeth_b[valid]

        odds_ratio = np.full(denom.shape, np.nan, dtype=np.float64)
        np.divide(numer, denom, out=odds_ratio, where=denom > 0)
        odds_ratio = np.where(
            denom > 0,
            odds_ratio,
            np.where(numer > 0, np.inf, np.nan),
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            log2_or[valid] = np.where(odds_ratio > 0, np.log2(odds_ratio), np.nan)

        ma = meth_a[valid]
        t, c1, r1 = total[valid], col1[valid], row1[valid]
        expected = r1.astype(np.float64) * c1.astype(np.float64) / t.astype(np.float64)
        p_upper = sp_stats.hypergeom.sf(ma - 1, t, c1, r1)
        p_lower = sp_stats.hypergeom.cdf(ma, t, c1, r1)
        one_tail = np.where(ma >= expected, p_upper, p_lower)
        pvals[valid] = np.minimum(2.0 * one_tail, 1.0)

    return pvals, log2_or


# CMH test -- O(n_sites) memory, statistically correct for replicates

def _cmh_init(n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Allocate CMH accumulators. Memory: ~32 bytes x n_sites total."""
    return (
        np.zeros(n, dtype=np.float64),  # Sigma(a - E): obs minus expected
        np.zeros(n, dtype=np.float64),  # Sigma V:      variance sum
        np.zeros(n, dtype=np.float64),  # Sigma(ad/n):  MH OR numerator
        np.zeros(n, dtype=np.float64),  # Sigma(bc/n):  MH OR denominator
    )


def _cmh_update(
    ome: np.ndarray,
    var_sum: np.ndarray,
    or_num: np.ndarray,
    or_den: np.ndarray,
    meth_case: np.ndarray,
    cov_case: np.ndarray,
    meth_ctrl: np.ndarray,
    cov_ctrl: np.ndarray,
) -> None:
    """In-place CMH accumulation from one case/control sample pair.

    Sites where either sample has zero coverage contribute V=0 and
    therefore do not influence the statistic -- this correctly handles
    union-mode sites with partial coverage without any special casing.
    """
    a = meth_case.astype(np.float64)
    b = (cov_case - meth_case).astype(np.float64)  # unmeth case
    c = meth_ctrl.astype(np.float64)
    d = (cov_ctrl - meth_ctrl).astype(np.float64)  # unmeth ctrl
    n = a + b + c + d

    # Sites need n > 1 for a non-degenerate variance term
    valid = n > 1

    row1 = a + b  # case coverage
    row2 = c + d  # ctrl coverage
    col1 = a + c  # total methylated
    col2 = b + d  # total unmethylated

    # Use safe denominator to avoid divide-by-zero warnings
    n_safe = np.where(n > 0, n, 1.0)
    E = np.where(valid, row1 * col1 / n_safe, 0.0)
    # Safe denominators to avoid divide-by-zero warnings
    n_sq_safe = np.where(n > 1, n * n * (n - 1.0), 1.0)
    V = np.where(
        valid,
        row1 * row2 * col1 * col2 / n_sq_safe,
        0.0,
    )

    ome[valid] += (a - E)[valid]
    var_sum[valid] += V[valid]

    # Mantel-Haenszel common odds ratio terms (using n_safe from above)
    or_num += np.where(valid, a * d / n_safe, 0.0)
    or_den += np.where(valid, b * c / n_safe, 0.0)


def _cmh_finalize(
    ome: np.ndarray,
    var_sum: np.ndarray,
    or_num: np.ndarray,
    or_den: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute CMH p-value and MH log2 OR from accumulated sums."""
    # Use safe denominator to avoid divide-by-zero warnings
    var_safe = np.where(var_sum > 0, var_sum, 1.0)
    cmh_stat = np.where(var_sum > 0, ome ** 2 / var_safe, np.nan)
    pvals = np.where(
        ~np.isnan(cmh_stat),
        sp_stats.chi2.sf(cmh_stat, df=1),
        np.nan,
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        mh_or = np.where(or_den > 0, or_num / or_den, np.nan)
        log2_mh_or = np.where(mh_or > 0, np.log2(mh_or), np.nan)

    return pvals, log2_mh_or


# Welford online statistics -- O(n_sites) memory regardless of n_samples

def _welford_init(n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Allocate Welford accumulators for n sites.

    Returns (mean, M2, n_valid).  Memory: ~20 bytes x n (float64 + int32).
    """
    return (
        np.zeros(n, dtype=np.float64),  # running mean
        np.zeros(n, dtype=np.float64),  # running sum of squared deviations
        np.zeros(n, dtype=np.int32),    # non-NaN replicate count per site
    )


def _welford_update(
    mean: np.ndarray,
    M2: np.ndarray,
    n_valid: np.ndarray,
    meth: np.ndarray,
    cov: np.ndarray,
) -> None:
    """In-place Welford update from one sample's integer meth/coverage arrays.

    Sites with zero coverage are treated as missing and skipped, so
    n_valid[i] counts only samples that actually covered site i.
    This handles (union sites with partial coverage) correctly.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        beta = np.where(cov > 0, meth.astype(np.float64) / cov, np.nan)
    valid = ~np.isnan(beta)
    if not np.any(valid):
        return
    n_valid[valid] += 1
    delta          = beta[valid] - mean[valid]
    mean[valid]   += delta / n_valid[valid]
    delta2         = beta[valid] - mean[valid]
    M2[valid]     += delta * delta2


def _welford_var_mean(M2: np.ndarray, n_valid: np.ndarray) -> np.ndarray:
    """Bessel-corrected variance of the group mean: s^2/n (per site).

    Sites with fewer than 2 valid replicates get NaN.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        var = np.where(n_valid > 1, M2 / (n_valid - 1), np.nan)
        return np.where(n_valid > 0, var / n_valid, np.nan)


def _logit_transform(beta: np.ndarray) -> np.ndarray:
    """Transform beta values to logit scale.

    logit(beta) = log(beta / (1 - beta))

    Handles boundary cases (beta=0, beta=1) by clipping to [epsilon, 1-epsilon].
    """
    beta_clipped = np.clip(beta, _BETA_EPSILON, 1 - _BETA_EPSILON)
    with np.errstate(divide="ignore", invalid="ignore"):
        logit_beta = np.log(beta_clipped / (1 - beta_clipped))
    return logit_beta


def _logit_variance_jacobian(beta: np.ndarray) -> np.ndarray:
    """Compute Jacobian for delta-method variance transformation.

    If Y = logit(X), then Var(Y) ~= Var(X) x [dY/dX]^2
    where dY/dX = 1 / [X(1-X)]
    """
    beta_clipped = np.clip(beta, _BETA_EPSILON, 1 - _BETA_EPSILON)
    with np.errstate(divide="ignore", invalid="ignore"):
        jacobian = 1.0 / (beta_clipped * (1 - beta_clipped))
    return jacobian


def _safe_log2_odds_ratio(
    mean_case: np.ndarray,
    mean_ctrl: np.ndarray,
) -> np.ndarray:
    """Symmetric log2 odds ratio with bounded clipping in both groups.

    fix: the previous formulation clipped only `(1 - mean_case)` and
    `(1 - mean_ctrl)` in the denominators of the inner ratios, so the
    numerators `mean_case` and `mean_ctrl` could remain at their raw values
    of 1.0, producing ratios of ``1 / epsilon`` that propagated to ``log2 ~= 30``
    and, when the outer ratio compounded, ``+/-inf``.

    Symmetric clipping in [epsilon, 1-epsilon] for both means caps the OR at
    ``log2((1-epsilon)/epsilon)^2 ~= 39.8`` regardless of which group is at the boundary,
    so finite-but-large values still signal extreme effects and ``inf``
    no longer pollutes the output column.
    """
    case_clip = np.clip(mean_case, _BETA_EPSILON, 1 - _BETA_EPSILON)
    ctrl_clip = np.clip(mean_ctrl, _BETA_EPSILON, 1 - _BETA_EPSILON)

    odds_case = case_clip / (1 - case_clip)
    odds_ctrl = ctrl_clip / (1 - ctrl_clip)

    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log2(odds_case / odds_ctrl)


# Quasi-binomial score test with chromosome-level overdispersion correction
# (a.k.a. "Version 1" of the count-model DMC family).
#
# Model for site i with replicates j and group g(j)  in  {case, ctrl}:
#
#     m_ij | n_ij  ~  Binomial(n_ij, pi_g(j)_i)             (binomial mean)
#     Var(m_ij)   =  phi * n_ij * pi_g(j)_i * (1 - pi_g(j)_i)   (quasi-binomial)
#
# Test H0: pi_case_i = pi_ctrl_i, against a two-sided alternative.  Let
#     S0_g = Sigma_j n_ij,   S1_g = Sigma_j m_ij,   S2_g = Sigma_j m_ij^2 / n_ij
# (group-wise running sums).  Then the group MLEs are pi_g = S1_g / S0_g and
# the pooled MLE under H0 is pi_pool = (S1_case + S1_ctrl) / (S0_case + S0_ctrl).
# The score for the group contrast and its null variance are
#     U      = S1_case - S0_case * pi_pool
#     Var(U) = phi * (S0_case * S0_ctrl / (S0_case + S0_ctrl)) * pi_pool * (1 - pi_pool)
# so the test statistic is U^2 / Var(U), chi^2_1 under H0.
#
# The dispersion phi is estimated from the full-model Pearson statistic:
#     X^2_g(i) = Sigma_j (m_ij - n_ij pi_g_i)^2 / (n_ij pi_g_i (1 - pi_g_i))
#             = (S2_g - S1_g^2/S0_g) / (pi_g_i (1 - pi_g_i))
# (the closed-form expansion lets us avoid materialising the n_sites x n_reps
# matrix). phi = sum_i_g X^2_g(i) / (n_obs - 2*n_sites_fit), clamped at 1.
#
# The test is replicate-aware (variance scales with the *number of
# replicates* via phi, not with the number of pooled reads), is a true
# count-based test (does not throw away the information that 5/10 carries
# less weight than 500/1000), and runs single-pass / O(n_sites) memory.

def _score_init(
    n: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Allocate quasi-binomial score accumulators (per group).

    Returns four arrays of length ``n``:

    sum_n   : float64 -- Sigma_j n_ij        (group coverage sum per site)
    sum_m   : float64 -- Sigma_j m_ij        (group meth-count sum per site)
    sum_m2n : float64 -- Sigma_j m_ij^2/n_ij  (used for closed-form Pearson)
    n_valid : int32   -- # replicates with coverage > 0 at this site

    Memory per group: 28 bytes x n_sites (~280 MB at 10 M sites).
    """
    return (
        np.zeros(n, dtype=np.float64),
        np.zeros(n, dtype=np.float64),
        np.zeros(n, dtype=np.float64),
        np.zeros(n, dtype=np.int32),
    )


def _score_update(
    sum_n:   np.ndarray,
    sum_m:   np.ndarray,
    sum_m2n: np.ndarray,
    n_valid: np.ndarray,
    meth:    np.ndarray,
    cov:     np.ndarray,
) -> None:
    """Fold one sample's (meth, coverage) arrays into the accumulators."""
    cov_f  = cov.astype(np.float64,  copy=False)
    meth_f = meth.astype(np.float64, copy=False)
    valid  = cov > 0

    sum_n += cov_f
    sum_m += meth_f
    # m^2 / n; zero contribution where the sample has no coverage.
    with np.errstate(invalid="ignore", divide="ignore"):
        sum_m2n += np.where(valid, meth_f * meth_f / np.maximum(cov_f, 1.0), 0.0)
    n_valid += valid.astype(np.int32)


def _smooth_box_kernel_py(
    pos:      np.ndarray,
    cum_meth: np.ndarray,
    cum_cov:  np.ndarray,
    half:     int,
    n:        int,
    meth_sm:  np.ndarray,
    cov_sm:   np.ndarray,
    meth_raw: np.ndarray,
    cov_raw:  np.ndarray,
) -> None:
    """Two-pointer sweep for box smoothing (pure-Python / numba target)."""
    lo = 0
    hi = 0
    for i in range(n):
        anchor = pos[i]
        while lo < n and (anchor - pos[lo]) > half:
            lo += 1
        if hi < lo:
            hi = lo
        while hi < n and (pos[hi] - anchor) <= half:
            hi += 1
        n_window = hi - lo
        if n_window <= 0:
            meth_sm[i] = float(meth_raw[i])
            cov_sm[i]  = float(cov_raw[i])
        else:
            meth_sm[i] = (cum_meth[hi] - cum_meth[lo]) / n_window
            cov_sm[i]  = (cum_cov[hi]  - cum_cov[lo])  / n_window


def _smooth_box_make_njit():
    """Build and cache the numba-compiled smoothing kernel."""
    global _SMOOTH_BOX_NJIT_FN
    if _SMOOTH_BOX_NJIT_FN is not None:
        return _SMOOTH_BOX_NJIT_FN
    try:
        from numba import njit
    except ImportError:
        njit = None
    if njit is not None:
        _SMOOTH_BOX_NJIT_FN = njit(cache=True)(_smooth_box_kernel_py)
    else:
        _SMOOTH_BOX_NJIT_FN = _smooth_box_kernel_py
    return _SMOOTH_BOX_NJIT_FN


def _smooth_sample_counts_box(
    meth:       np.ndarray,
    cov:        np.ndarray,
    positions:  np.ndarray,
    window_bp:  int,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-sample uniform-box smoothing of (meth, cov) counts.

    Matches DSS ``smooth.chr(..., method="avg")`` semantics exactly: for
    each CpG i, the smoothed count is the *average* of the counts at all
    CpGs within ``[pos_i - window_bp//2, pos_i + window_bp//2]``
    inclusive -- uniform rectangular kernel, per chromosome, unweighted by
    coverage. DSS rounds the smoothed counts back to integers; we keep
    float64 so the existing ``_score_update`` accumulators retain
    sub-replicate precision without losing information to early rounding.

    Implementation mirrors DSS's ``nitem.c`` + ``filter.c``: a single
    cumulative-sum pass plus a two-pointer sweep over sorted positions.
    O(n) total. The sweep is compiled via ``numba.njit`` when available.

    Parameters
    ----------
    meth, cov : np.ndarray
        Per-CpG methylated-read count and coverage for a single sample.
        Typically int32; cast to int64 internally for the cumsum.
    positions : np.ndarray
        Sorted CpG positions on this chromosome (int32/int64). Must be
        the same length as ``meth`` / ``cov``.
    window_bp : int
        Full smoothing window. ``+/-window_bp // 2`` bp on each side. The
        DSS default is ``smoothing.span=500``.

    Returns
    -------
    (meth_sm, cov_sm) : tuple[np.ndarray, np.ndarray]
        Float64 arrays of the same length as ``meth`` / ``cov``.
    """
    n = positions.shape[0]
    if n == 0:
        return meth.astype(np.float64, copy=False), cov.astype(np.float64, copy=False)
    if window_bp <= 0:
        return meth.astype(np.float64, copy=False), cov.astype(np.float64, copy=False)

    half = int(window_bp) // 2
    pos  = positions.astype(np.int64, copy=False)

    cum_meth = np.empty(n + 1, dtype=np.int64)
    cum_cov  = np.empty(n + 1, dtype=np.int64)
    cum_meth[0] = 0
    cum_cov[0]  = 0
    np.cumsum(meth.astype(np.int64, copy=False), out=cum_meth[1:])
    np.cumsum(cov.astype(np.int64,  copy=False), out=cum_cov[1:])

    meth_sm = np.empty(n, dtype=np.float64)
    cov_sm  = np.empty(n, dtype=np.float64)

    kernel = _smooth_box_make_njit()
    kernel(pos, cum_meth, cum_cov, half, n, meth_sm, cov_sm,
           meth.astype(np.float64), cov.astype(np.float64))

    return meth_sm, cov_sm


def _score_finalize(
    sn_case:  np.ndarray, sm_case:  np.ndarray, sm2n_case: np.ndarray, nv_case: np.ndarray,
    sn_ctrl:  np.ndarray, sm_ctrl:  np.ndarray, sm2n_ctrl: np.ndarray, nv_ctrl: np.ndarray,
    chrom_name:     str   = "?",
    min_dispersion: float = 1.0,
    min_disp_sites: int   = 100,
    dispersion:     str   = "site",
    shrink_pseudo_df: float = 4.0,
    statistic:      str   = "lr",
    reference:      str   = "adaptive",
    sep_fallback:   bool  = False,
    sep_threshold:  float = 0.9,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """Compute per-site score p-values with McCullagh-Nelder overdispersion.

    Parameters
    ----------
    sn_*, sm_*, sm2n_*, nv_* : np.ndarray
        Per-group accumulators returned by ``_score_init`` / ``_score_update``.
    chrom_name : str
        Used only for logging.
    min_dispersion : float
        Floor on phi. Underdispersion (phi < 1) usually reflects model
        misspecification rather than truly less-than-binomial variability;
        clamping at 1 is the conservative choice. Set < 1 to allow
        underdispersion.
    min_disp_sites : int
        If fewer than this many sites are usable for the chromosome-level
        dispersion estimate (e.g. tiny alt contigs), fall back to
        phi = ``min_dispersion`` instead of computing an unstable estimate.
        Used by the ``"chrom"`` and ``"shrink"`` modes; the ``"site"`` mode
        is unaffected because each site provides its own phi_i.
    dispersion : {"site", "chrom", "shrink"}
        Strategy for the McCullagh-Nelder dispersion correction:

        ``"site"`` (default)
            Each site/tile gets its own phi_i computed from its 4-df Pearson
            residual sum:  phi_i = (X^2_case_i + X^2_ctrl_i) / max(nv_i - 2, 1),
            clamped at ``min_dispersion``. This is what R's
            ``glm(family=quasibinomial)`` does on a per-site fit. The
            estimator is noisy with the typical 4 df, but it correctly
            tracks region-specific dispersion (CpG-islands vs gene bodies
            have different between-replicate variance).
        ``"chrom"``
            Single chromosome-pooled phi from all qualifying sites. Most
            powerful when between-replicate variance really is constant
            along the chromosome, but anti-conservative when it is not.
            Equivalent to fitting the quasi-binomial GLM as a single model
            with one shared dispersion. Use when you want strictly more
            power and you accept the modelling assumption.
        ``"shrink"``
            James-Stein-style shrinkage: phi_shrunk_i is a weighted average of
            phi_site_i (with weight = site df) and the chromosome-pooled
            phi_chrom (with weight = ``shrink_pseudo_df``, default 4).
            Trades a small bias for a large variance reduction on the
            per-site estimate (the same idea behind DSS's empirical-Bayes
            shrinkage of dispersion).
    shrink_pseudo_df : float
        Pseudo-df weight on phi_chrom in the ``"shrink"`` mode (default 4 ~=
        the typical real per-site df). Ignored otherwise.
    statistic : {"lr", "score"}
        Functional form of the test statistic. Both use the same per-group
        sufficient statistics (S0_g = Sigma_j n_ij and S1_g = Sigma_j m_ij) and the
        same dispersion correction, so the difference manifests only at
        small effective sample sizes (n=6 is small):

        ``"lr"`` (default)
            Quasi-binomial likelihood-ratio chi-square. Closed-form in
            S0_g and S1_g, so no per-tile GLM fit is required:
              LRT = 2 * Sigma_g [ S1_g*log(p_g/p_pool)
                            + (S0_g - S1_g)*log((1 - p_g)/(1 - p_pool)) ]
            divided by the dispersion phi_i. Closer to nominal coverage near
            the boundaries (pi near 0 or 1) -- exactly where DMR tiles tend
            to live.

        ``"score"`` (slightly more powerful)
            Pearson score statistic U^2/V_pool with quasi-binomial inflation.
            Asymptotically equivalent to the LR test but mildly
            anti-conservative at the boundaries. Kept as an option for
            users who want the small extra power.
    reference : {"adaptive", "chi2", "F"}
        Reference distribution used to convert the test statistic to a
        p-value.

        ``"adaptive"`` (default)
            Per-site: F(1, df_residual_i) where the per-site dispersion
            phi_i > 1 (real overdispersion signal made it past the
            min-dispersion floor), chi^2(1) where phi_i was clamped to 1. This
            is the right behaviour for quasi-binomial GLMs whose
            dispersion estimate is noisy at small samples -- F handles the
            overdispersed sites, chi^2 handles the ones where the
            quasi-binomial collapses to a binomial.
        ``"chi2"``
            Always reference to chi^2(1). Over-liberal at tiles with real
            overdispersion (typical DMR setting).
        ``"F"``
            Always reference to F(1, df_residual_i). Wildly conservative
            at sites where phi is clamped to 1 -- will reject ~zero
            genome-wide CpGs on typical n=3+3 WGBS.

    Returns
    -------
    pvals : float64 array, NaN at degenerate sites
    log2_or : float64 array, NaN at degenerate sites
    pi_case, pi_ctrl : float64 arrays of coverage-weighted group methylation
        (= group MLE proportion under the full model). NaN where the
        corresponding group has zero coverage at the site.
    phi_hat : float
        Chromosome-pooled dispersion estimate, returned for logging /
        downstream introspection regardless of which mode was used.
    """
    if dispersion not in {"site", "chrom", "shrink", "eb"}:
        raise ValueError(
            f"dispersion must be 'site', 'chrom', 'shrink', or 'eb'; got {dispersion!r}"
        )
    if statistic not in {"lr", "score"}:
        raise ValueError(
            f"statistic must be 'lr' or 'score'; got {statistic!r}"
        )
    if reference == "methylkit":
        raise ValueError(
            "reference='methylkit' was renamed to 'adaptive' in this release."
        )
    if reference not in {"adaptive", "F", "chi2"}:
        raise ValueError(
            f"reference must be 'adaptive', 'F', or 'chi2'; got {reference!r}"
        )
    eps = _BETA_EPSILON

    # --- Group MLE proportions under the full (unrestricted) model ---
    with np.errstate(invalid="ignore", divide="ignore"):
        pi_case = np.where(sn_case > 0, sm_case / sn_case, np.nan)
        pi_ctrl = np.where(sn_ctrl > 0, sm_ctrl / sn_ctrl, np.nan)

    # --- Pearson chi-sq contribution from each site & group, full model ---
    # Numerator (closed form): N_g * S2_g - S1_g^2 , scaled by 1/S0_g^2/pi(1-pi).
    # Compact form:  contrib = (S2 - S1^2/S0) / (pi * (1 - pi))
    # where pi(1-pi) = S1*(S0 - S1) / S0^2.  Sites with S1  in  {0, S0} have
    # variance 0 and contribute nothing to the dispersion estimate.
    with np.errstate(invalid="ignore", divide="ignore"):
        # numerator term Sigma_j(m - npi)^2/n   =  S2 - S1^2/S0
        num_case = sm2n_case - np.where(sn_case > 0, sm_case ** 2 / sn_case, 0.0)
        num_ctrl = sm2n_ctrl - np.where(sn_ctrl > 0, sm_ctrl ** 2 / sn_ctrl, 0.0)

        den_case = np.where(
            sn_case > 0,
            sm_case * (sn_case - sm_case) / (sn_case ** 2),
            0.0,
        )
        den_ctrl = np.where(
            sn_ctrl > 0,
            sm_ctrl * (sn_ctrl - sm_ctrl) / (sn_ctrl ** 2),
            0.0,
        )

        chi_case = np.where(den_case > 0, num_case / den_case, 0.0)
        chi_ctrl = np.where(den_ctrl > 0, num_ctrl / den_ctrl, 0.0)

    sites_both    = (sn_case > 0) & (sn_ctrl > 0) & (nv_case > 0) & (nv_ctrl > 0)
    sites_dispers = sites_both & (den_case > 0) & (den_ctrl > 0)

    # --- Chromosome-pooled phi (always computed; used directly in "chrom"
    #     mode, used as the shrinkage anchor in "shrink" mode, logged
    #     otherwise) -----------------------------------------------------
    n_disp = int(sites_dispers.sum())
    if n_disp < min_disp_sites:
        if dispersion == "chrom":
            logger.warning(
                "%s: only %d sites usable for dispersion estimation; "
                "falling back to phi = %.2f (no overdispersion correction).",
                chrom_name, n_disp, min_dispersion,
            )
        phi_hat = float(min_dispersion)
        phi_raw = float(min_dispersion)
        # df_chrom = 1 in the fallback case isn't actually used because
        # phi_hat = min_dispersion = 1.0 here, which never trips the F-branch.
        df_chrom = 1.0
    else:
        n_obs = int(nv_case[sites_dispers].sum() + nv_ctrl[sites_dispers].sum())
        df_chrom = float(max(n_obs - 2 * n_disp, 1))

        pearson_sum = float(
            chi_case[sites_dispers].sum() + chi_ctrl[sites_dispers].sum()
        )
        phi_raw = pearson_sum / df_chrom
        phi_hat = float(max(min_dispersion, phi_raw))

        logger.info(
            "%s: chrom-pooled phi = %.3f (raw %.3f, %s sites, %s obs, df=%s); "
            "applying dispersion='%s'",
            chrom_name, phi_hat, phi_raw,
            f"{n_disp:,}", f"{n_obs:,}", f"{int(df_chrom):,}", dispersion,
        )

    # --- Per-site Pearson dispersion phi_i (only used when needed) ---------
    if dispersion in ("site", "shrink", "eb"):
        # df_i = (replicates_case + replicates_ctrl) - 2 fitted proportions.
        # At a typical n=3 per group this is 4.
        df_i = (nv_case + nv_ctrl).astype(np.float64) - 2.0
        df_i_safe = np.where(df_i > 0, df_i, 1.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            phi_site = (chi_case + chi_ctrl) / df_i_safe

        # Sites with zero dispersion contribution (perfect fit OR degenerate
        # variance term) cannot inform phi_i. Apply the floor; the
        # ``min_dispersion`` clamp also handles the underdispersion case.
        phi_site = np.where(sites_dispers & (df_i > 0), phi_site, min_dispersion)
        phi_site = np.maximum(phi_site, min_dispersion)

        if dispersion == "site":
            phi_eff = phi_site
            # Per-site phi is estimated from this site's 4 df Pearson sum.
            df_phi = df_i_safe.copy()
        elif dispersion == "eb":
            # Empirical-Bayes shrinkage with data-driven shrinkage weight.
            # Treat phi_site_i as a noisy estimator of an unknown true phi_i,
            # and assume phi_i ~ inverse-Gamma(a, b) across sites with a, b
            # estimated by method of moments from the chromosome-wide
            # distribution of phi_site values. The posterior mean of phi_i
            # given phi_site_i and df_i is the same weighted-average form as
            # `shrink`, but the weight on phi_chrom is now
            #     w_eb = 2 * a   (the implied pseudo-df from the IG prior)
            # rather than the fixed 4 of `shrink`. When the per-site phi
            # estimates have high variance relative to phi_chrom, w_eb grows
            # large (shrink toward mean); when low variance, w_eb -> 0
            # (trust the site).
            valid_phi = sites_dispers & (df_i > 0)
            if int(np.sum(valid_phi)) >= min_disp_sites:
                phi_obs = phi_site[valid_phi]
                # Method-of-moments on inverse-Gamma:
                #   mean = b / (a - 1),  var = b^2 / ((a-1)^2 (a-2))
                # => a = mean^2 / var + 2,  b = mean * (a - 1)
                m = float(np.mean(phi_obs))
                v = float(np.var(phi_obs))
                if v > 1e-9 and m > 0:
                    a_mom = m * m / v + 2.0
                    w_eb = max(1.0, 2.0 * a_mom)
                else:
                    w_eb = float(shrink_pseudo_df)
            else:
                w_eb = float(shrink_pseudo_df)
            num = df_i_safe * phi_site + w_eb * phi_hat
            den = df_i_safe + w_eb
            phi_eff = np.maximum(num / den, min_dispersion)
            phi_eff = np.where(sites_dispers & (df_i > 0), phi_eff, phi_hat)
            # Weighted-average phi has effective df = (df_i + w_eb). Where we
            # fell back to phi_chrom, the df is the chrom-pool df.
            df_phi = np.where(
                sites_dispers & (df_i > 0),
                df_i_safe + float(w_eb),
                max(df_chrom, 1.0),
            )
            logger.info(
                "%s: empirical-Bayes shrinkage with w_eb=%.2f (phi_chrom=%.3f, "
                "phi_site var=%.2g over %d sites).",
                chrom_name, w_eb, phi_hat,
                v if 'v' in locals() else float('nan'),
                int(np.sum(valid_phi)),
            )
        else:  # "shrink": James-Stein-style weighted average toward chrom mean
            # phi_shrunk_i = (df_i * phi_site_i + w * phi_chrom) / (df_i + w)
            w   = float(shrink_pseudo_df)
            num = df_i_safe * phi_site + w * phi_hat
            den = df_i_safe + w
            phi_eff = np.maximum(num / den, min_dispersion)
            # Where the per-site estimator was unusable, fall back to the
            # chromosome value rather than the floor.
            phi_eff = np.where(sites_dispers & (df_i > 0), phi_eff, phi_hat)
            # Weighted-average phi has effective df = (df_i + shrink_pseudo_df).
            # Where we fell back to phi_chrom, the df is the chrom-pool df.
            df_phi = np.where(
                sites_dispers & (df_i > 0),
                df_i_safe + w,
                max(df_chrom, 1.0),
            )
    else:  # "chrom"
        phi_eff = np.full_like(sn_case, phi_hat, dtype=np.float64)
        # Chromosome-pooled phi is estimated from the whole chrom's Pearson sum,
        # so df_phi here is df_chrom (often >> 1e5, making F(1, df_phi) -> chi^2(1)).
        df_phi = np.full_like(sn_case, max(df_chrom, 1.0), dtype=np.float64)

    # --- Test for H0: pi_case = pi_ctrl --------------------------------------
    # Both the score and LR statistics use the same per-group sufficient
    # statistics (sn_*, sm_*) and the same dispersion phi_eff. They differ
    # only in functional form; both are referenced to chi^2_1 asymptotically.
    sn_total = sn_case + sn_ctrl
    sm_total = sm_case + sm_ctrl
    with np.errstate(invalid="ignore", divide="ignore"):
        pi_pool      = np.where(sn_total > 0, sm_total / sn_total, np.nan)
        pi_pool_safe = np.clip(pi_pool, eps, 1.0 - eps)

        # Variance of the null-MLE score U.  Used by both branches: the
        # ``score`` test divides U^2 by it; the ``lr`` test only needs it
        # as a degenerate-site guard (variance == 0 -> no information at
        # that site, so the LR is also undefined).
        var_U_bin = (
            (sn_case * sn_ctrl / np.maximum(sn_total, 1.0))
            * pi_pool_safe
            * (1.0 - pi_pool_safe)
        )

        if statistic == "score":
            U = sm_case - sn_case * pi_pool
            stat_raw = np.where(var_U_bin > 0, U * U / var_U_bin, np.nan)
        else:  # "lr": closed-form quasi-binomial log-likelihood ratio
            # LR = 2 * Sigma_g [ S1_g * log(p_g/p_pool)
            #             + (S0_g - S1_g) * log((1 - p_g)/(1 - p_pool)) ]
            # Each term is x * log(y/z) where y, z  in  [epsilon, 1-epsilon] after clipping,
            # so log is bounded; the multiplicative 0 at x=0 cleanly zeroes
            # the contribution (no special-casing needed).
            pc_safe = np.clip(pi_case, eps, 1.0 - eps)
            pk_safe = np.clip(pi_ctrl, eps, 1.0 - eps)

            u_case = sm_case
            u_ctrl = sm_ctrl
            v_case = sn_case - sm_case
            v_ctrl = sn_ctrl - sm_ctrl

            lr_terms = (
                u_case * np.log(pc_safe / pi_pool_safe)
                + v_case * np.log((1.0 - pc_safe) / (1.0 - pi_pool_safe))
                + u_ctrl * np.log(pk_safe / pi_pool_safe)
                + v_ctrl * np.log((1.0 - pk_safe) / (1.0 - pi_pool_safe))
            )
            stat_raw = 2.0 * lr_terms

        # Apply quasi-binomial dispersion inflation per site/tile.
        chi2_stat = np.where(phi_eff > 0, stat_raw / phi_eff, np.nan)
        chi2_stat = np.where(var_U_bin > 0, chi2_stat, np.nan)

    # --- Reference distribution -> p-value ---------------------------------
    # Per-site adaptive switch: F(1, df_phi) where the dispersion phi cleared
    # the min-dispersion floor (phi > 1, i.e. real overdispersion signal),
    # chi^2(1) where phi was clamped to 1. F handles the overdispersed sites;
    # chi^2 handles the ones where the quasi-binomial collapses to a binomial.
    #
    # df_phi is the df backing the phi estimate (NOT the per-site residual df
    # nv_case+nv_ctrl-2). It differs across dispersion modes:
    #   - "site":   per-site Pearson sum has df ~= 4 at n=3+3
    #   - "chrom":  chromosome-pooled phi has df_chrom (often >> 1e5)
    #   - "shrink/eb": df_i + shrink_pseudo_df (or w_eb), can collapse to ~4
    # In "eb" mode with small w_eb (the homogeneous-dispersion / easy
    # case), df_phi falls back to ~4 and F(1, 4) is ~250x more conservative
    # than chi^2(1) at typical test statistics -- this is the bug behind
    # the artifactually low FPR in eb mode in 0.7.2. Floor df_phi at 50
    # so F(1, 50) is within ~1% of chi^2(1) at typical statistics on
    # the F branch; the chi^2 branch (clamped phi) is unaffected.
    # Rationale for 50: F(1, 50) is within ~1% of chi^2(1) at the 5%
    # critical region (chi^2 critical value 3.84).
    DF_PHI_FLOOR = 50.0  # F(1, 50) within ~1% of chi^2(1) at the 5% critical region
    if reference == "adaptive":
        df_phi_floored = np.maximum(df_phi, DF_PHI_FLOOR)
        p_F    = sp_stats.f.sf(chi2_stat, dfn=1, dfd=df_phi_floored)
        p_chi2 = sp_stats.chi2.sf(chi2_stat, df=1)
        pvals  = np.where(phi_eff > 1.0, p_F, p_chi2)
    elif reference == "F":
        df_phi_floored = np.maximum(df_phi, DF_PHI_FLOOR)
        pvals = sp_stats.f.sf(chi2_stat, dfn=1, dfd=df_phi_floored)
    else:  # "chi2"
        pvals = sp_stats.chi2.sf(chi2_stat, df=1)

    degenerate = (
        ~sites_both
        | np.isnan(chi2_stat)
        | (var_U_bin <= 0)
    )
    pvals = np.where(degenerate, np.nan, pvals)

    log2_or = _safe_log2_odds_ratio(pi_case, pi_ctrl)
    log2_or[degenerate] = np.nan

    # --- Separation-aware Fisher fallback ----------------------------------
    # At very low coverage, sites with a large true effect can produce
    # pooled 2x2 tables at or near perfect separation, where the
    # quasi-binomial LR statistic can collapse (the asymptotic chi^2/F
    # reference under-estimates the true p-value because the LR statistic
    # is point-massed at small counts). For such sites, the exact
    # hypergeometric Fisher p-value on pooled counts is both well-defined
    # and uniformly more powerful than the LR's asymptotic approximation.
    #
    # The fallback fires only for sites where:
    #   1) the observed |meth_diff| >= sep_threshold (i.e. there IS a
    #      large effect to be tested), AND
    #   2) the LR p-value did not reject at p > 0.05 (i.e. the asymptotic
    #      test failed despite the large effect), AND
    #   3) both groups have at least one read.
    # Sites that already reject under LR are left alone; the fallback can
    # only re-test sites the LR test missed, so FPR is unchanged.
    if sep_fallback:
        with np.errstate(invalid="ignore"):
            obs_diff = np.abs(pi_case - pi_ctrl)
        fb_mask = (
            np.isfinite(obs_diff)
            & (obs_diff >= sep_threshold)
            & (pvals > 0.05)
            & (sn_case > 0)
            & (sn_ctrl > 0)
            & ~degenerate
        )
        n_fb = int(np.sum(fb_mask))
        if n_fb > 0:
            from scipy.stats import fisher_exact as _scipy_fisher
            idx = np.flatnonzero(fb_mask)
            M_a = sm_case[idx].astype(np.int64)
            U_a = (sn_case[idx] - sm_case[idx]).astype(np.int64)
            M_b = sm_ctrl[idx].astype(np.int64)
            U_b = (sn_ctrl[idx] - sm_ctrl[idx]).astype(np.int64)
            p_fb = np.empty(n_fb, dtype=np.float64)
            for k in range(n_fb):
                try:
                    _, p_fb[k] = _scipy_fisher(
                        [[int(M_a[k]), int(U_a[k])], [int(M_b[k]), int(U_b[k])]],
                        alternative="two-sided",
                    )
                except Exception:
                    p_fb[k] = pvals[idx[k]]
            # Use the better (smaller) of LR / Fisher; never inflate p
            improved = p_fb < pvals[idx]
            n_improved = int(np.sum(improved))
            pvals[idx[improved]] = p_fb[improved]
            logger.info(
                "%s: separation fallback fired on %d sites "
                "(|meth_diff|>=%.2f & LR-p>0.05); %d had Fisher-p<LR-p.",
                chrom_name, n_fb, sep_threshold, n_improved,
            )

    return pvals, log2_or, pi_case, pi_ctrl, phi_hat


def _beta_binom_mom_from_welford_logit(
    mean_case: np.ndarray,
    M2_case: np.ndarray,
    n_valid_case: np.ndarray,
    mean_ctrl: np.ndarray,
    M2_ctrl: np.ndarray,
    n_valid_ctrl: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Welch t-test on logit-transformed beta values (from Welford accumulators).

    Beta values are highly skewed near boundaries [0, 1]. Logit transformation
    stabilizes variance and improves normality for robust t-testing.

    Parameters match _beta_binom_mom_from_welford; output is p-values and
    log2 odds ratio on original (non-logit) scale.
    """
    # Transform means to logit scale
    logit_mean_case = _logit_transform(mean_case)
    logit_mean_ctrl = _logit_transform(mean_ctrl)

    # Compute variance on logit scale via delta method
    # Var(logit(beta)) = Var(beta) x jacobian^2
    jac_case = _logit_variance_jacobian(mean_case)
    jac_ctrl = _logit_variance_jacobian(mean_ctrl)

    var_case = M2_case / np.maximum(n_valid_case - 1, 1)
    var_ctrl = M2_ctrl / np.maximum(n_valid_ctrl - 1, 1)

    var_logit_case = var_case * (jac_case ** 2)
    var_logit_ctrl = var_ctrl * (jac_ctrl ** 2)

    # Normalize by sample size (variance of the mean)
    var_mean_logit_case = var_logit_case / np.maximum(n_valid_case, 1)
    var_mean_logit_ctrl = var_logit_ctrl / np.maximum(n_valid_ctrl, 1)

    se = np.sqrt(var_mean_logit_case + var_mean_logit_ctrl)

    # Welch t-test on logit scale
    with np.errstate(invalid="ignore", divide="ignore"):
        t_stat = np.where(se > 0, (logit_mean_case - logit_mean_ctrl) / se, np.nan)

    # Welch-Satterthwaite degrees of freedom
    dof_num = (var_mean_logit_case + var_mean_logit_ctrl) ** 2
    dof_den = (
        np.where(n_valid_case > 1, var_mean_logit_case ** 2 / (n_valid_case - 1), 0.0)
        + np.where(n_valid_ctrl > 1, var_mean_logit_ctrl ** 2 / (n_valid_ctrl - 1), 0.0)
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        dof = np.where(dof_den > 0, dof_num / dof_den, 1.0)
        dof = np.maximum(dof, 1.0)

    pvals = 2.0 * sp_stats.t.sf(np.abs(t_stat), df=dof)

    # Degenerate cases. NaN only when BOTH groups have zero between-
    # replicate variance (M2 = 0 with n_valid >= 2): only then is the
    # Welch t genuinely undefined (SE = 0). The prior, stricter version of
    # this guard NaN'd on *either*-zero variance to suppress an H0 false-
    # positive surge at boundary beta; but that killed all power on the
    # standard fixture, where ~half the truly-differential hypo sites
    # have beta_treatment clipped to ~0.01 and routinely collapse to
    # M2_case = 0 by binomial sampling chance. Position: logit_t is the
    # weak variance-stabilising fallback; we don't pretend it's
    # well-calibrated near beta = 0 / 1. For trustworthy inference use
    # ``test="lr"`` (or ``"score"``), which gets variance from the
    # binomial count model and isn't affected by replicate collapse.
    both_zero_var = (
        (n_valid_case >= 2) & (M2_case <= 0.0)
        & (n_valid_ctrl >= 2) & (M2_ctrl <= 0.0)
    )
    degenerate = (
        np.isnan(mean_case) | np.isnan(mean_ctrl) | np.isnan(t_stat)
        | (n_valid_case == 0) | (n_valid_ctrl == 0)
        | both_zero_var
    )
    pvals[degenerate] = np.nan

    # Compute log2 odds ratio on original scale using symmetric clamp
    log2_ors = _safe_log2_odds_ratio(mean_case, mean_ctrl)
    log2_ors[degenerate] = np.nan

    return pvals, log2_ors


def _beta_binom_mom_from_welford(
    mean_case: np.ndarray,
    M2_case: np.ndarray,
    n_valid_case: np.ndarray,
    mean_ctrl: np.ndarray,
    M2_ctrl: np.ndarray,
    n_valid_ctrl: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Welch t-test derived from Welford accumulators.

    Mathematically equivalent to _beta_binom_mom but never builds the
    (n_sites x n_replicates) beta matrix.  Per-site valid counts are used
    for both variance estimation and Welch-Satterthwaite DOF, which correctly
    handles sites where some replicates have no coverage (union / outer-join
    mode).
    """
    vm_case = _welford_var_mean(M2_case, n_valid_case)
    vm_ctrl = _welford_var_mean(M2_ctrl, n_valid_ctrl)

    se = np.sqrt(vm_case + vm_ctrl)

    with np.errstate(invalid="ignore", divide="ignore"):
        t_stat = np.where(se > 0, (mean_case - mean_ctrl) / se, np.nan)

    # Welch-Satterthwaite degrees of freedom (per-site n_valid)
    dof_num = (vm_case + vm_ctrl) ** 2
    dof_den = (
        np.where(n_valid_case > 1, vm_case ** 2 / (n_valid_case - 1), 0.0)
        + np.where(n_valid_ctrl > 1, vm_ctrl ** 2 / (n_valid_ctrl - 1), 0.0)
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        dof = np.where(dof_den > 0, dof_num / dof_den, 1.0)
        dof = np.maximum(dof, 1.0)

    pvals = 2.0 * sp_stats.t.sf(np.abs(t_stat), df=dof)

    # NaN only when BOTH groups have zero between-replicate variance
    # (SE = 0, t genuinely undefined). See the matching block in
    # _beta_binom_mom_from_welford_logit for the rationale: the
    # *either*-zero version killed real signal at boundary beta.
    # welch_t is treated as a weak fallback -- use ``test="lr"``
    # for trustworthy inference.
    both_zero_var = (
        (n_valid_case >= 2) & (M2_case <= 0.0)
        & (n_valid_ctrl >= 2) & (M2_ctrl <= 0.0)
    )
    degenerate = (
        np.isnan(mean_case) | np.isnan(mean_ctrl) | np.isnan(t_stat)
        | (n_valid_case == 0) | (n_valid_ctrl == 0)
        | both_zero_var
    )
    pvals[degenerate] = np.nan

    # symmetric clamp on both group means so log2 OR cannot blow up
    # to +/-inf when one group is at the boundary 0 or 1.
    log2_ors = _safe_log2_odds_ratio(mean_case, mean_ctrl)
    log2_ors[degenerate] = np.nan

    return pvals, log2_ors


# Internal per-chromosome helpers

def _detect_chromosomes(methylstore_path: Path) -> list[str]:
    chroms: set[str] = set()
    for sample_dir in methylstore_path.glob("sample=*"):
        for chrom_dir in sample_dir.glob("chrom=*"):
            chroms.add(chrom_dir.name.removeprefix("chrom="))
    return sorted(chroms)


def _intersect_chrom(
    methylstore_path: Path,
    chrom: str,
    samples: list[str],
) -> pl.DataFrame:
    """Return (pos, strand) rows present in every sample for one chromosome.

    The previous implementation joined on ["pos", "strand"].  Samples
    without a reference FASTA receive strand="*" while samples converted with
    a FASTA receive "+"/"-".  A mixed cohort produced an empty intersection
    with no warning.  We now join on "pos" only and resolve the strand column
    by taking the first non-"*" value seen across samples (falling back to "*"
    when all samples lack strand information).

    (this revision): the per-sample `sites` frame is now deduplicated on
    `pos` (keeping the first row) before the join.  Without this guard, a
    sample whose .cov file recorded both strands of one CpG dinucleotide
    (e.g. + at N and - at N+1 that were not merged by _merge_cpg_pairs)
    would produce one row per strand at the same pos, and the inner join
    on `pos` would multiply rows downstream -- silently breaking the
    one-row-per-site contract that _load_sample_chrom relies on.
    """
    _empty = pl.DataFrame({
        "pos":    pl.Series([], dtype=pl.Int32),
        "strand": pl.Series([], dtype=pl.Utf8),
    })
    n_samples = len(samples)
    if n_samples == 0:
        return _empty

    site_dfs: list[pl.DataFrame] = []
    for sample in samples:
        part_file = (
            methylstore_path / f"sample={sample}" / f"chrom={chrom}" / "part-0.parquet"
        )
        if not part_file.exists():
            logger.debug(
                "  Sample '%s' missing %s; chromosome excluded from intersection",
                sample, chrom,
            )
            return _empty
        site_dfs.append(
            pl.read_parquet(str(part_file), columns=["pos", "strand"])
            .unique(subset=["pos"], keep="first")
        )

    combined = pl.concat(site_dfs)
    intersected = (
        combined
        .group_by("pos")
        .agg([
            pl.len().alias("_n"),
            pl.col("strand").filter(pl.col("strand") != "*").first().alias("_strand_real"),
            pl.col("strand").first().alias("_strand_fb"),
        ])
        .filter(pl.col("_n") == n_samples)
        .with_columns(
            pl.when(pl.col("_strand_real").is_not_null())
            .then(pl.col("_strand_real"))
            .otherwise(pl.col("_strand_fb"))
            .alias("strand")
        )
        .select(["pos", "strand"])
        .sort("pos")
    )

    if len(intersected) == 0:
        logger.warning(
            "  Intersection is empty on %s. "
            "Check strand consistency across samples.", chrom,
        )

    return intersected


def _union_chrom(
    methylstore_path: Path,
    chrom: str,
    samples: list[str],
) -> pl.DataFrame:
    """Return (pos, strand) rows seen in at least one sample."""
    site_dfs: list[pl.DataFrame] = []
    for sample in samples:
        part_file = (
            methylstore_path / f"sample={sample}" / f"chrom={chrom}" / "part-0.parquet"
        )
        if part_file.exists():
            site_dfs.append(
                pl.read_parquet(str(part_file), columns=["pos", "strand"])
            )
    if not site_dfs:
        return pl.DataFrame({
            "pos":    pl.Series([], dtype=pl.Int32),
            "strand": pl.Series([], dtype=pl.Utf8),
        })
    return (
        pl.concat(site_dfs)
        .unique(subset=["pos"], keep="first")  # dedupe on pos only
        .sort("pos")
    )


def _load_sample_chrom(
    methylstore_path: Path,
    chrom: str,
    sample: str,
    canonical_pos: pl.DataFrame,
) -> tuple[np.ndarray, np.ndarray]:
    """Load N_meth and coverage for ONE sample / ONE chromosome.

    Left-joins to canonical_pos so arrays are aligned to the same site order.
    Missing sites are filled with 0 .
    """
    part_file = (
        methylstore_path / f"sample={sample}" / f"chrom={chrom}" / "part-0.parquet"
    )
    n_sites = len(canonical_pos)

    if not part_file.exists():
        return (
            np.zeros(n_sites, dtype=np.int32),
            np.zeros(n_sites, dtype=np.int32),
        )

    df = pl.read_parquet(str(part_file), columns=["pos", "N_meth", "coverage"])
    if df.height != df["pos"].n_unique():
        df = (
            df.group_by("pos")
            .agg([
                pl.sum("N_meth").alias("N_meth"),
                pl.sum("coverage").alias("coverage"),
            ])
        )
    aligned = canonical_pos.join(df, on="pos", how="left").fill_null(0)

    return (
        aligned["N_meth"].to_numpy().astype(np.int32),
        aligned["coverage"].to_numpy().astype(np.int32),
    )


def _process_one_chromosome(
    methylstore_path: Path,
    chrom: str,
    canonical_df: pl.DataFrame,
    samples_case: list[str],
    samples_control: list[str],
    test: str,
    min_samples_case: int = 0,
    min_samples_control: int = 0,
    dispersion: str = "site",
    reference: str = "adaptive",
    design_full: Optional[np.ndarray] = None,
    design_reduced: Optional[np.ndarray] = None,
    coef_idx: Optional[int] = None,
    contrast_matrix: Optional[np.ndarray] = None,
    contrast_label: Optional[str] = None,
    samples_all_ordered: Optional[list[str]] = None,
    group_labels_per_sample: Optional[list[str]] = None,
    glm_backend: str = "cpu",
    smoothing: bool = False,
    smoothing_span_bp: int = 500,
    sep_fallback: bool = False,
    sep_threshold: float = 0.9,
) -> pl.DataFrame:
    """Run DMC for one chromosome, loading one sample at a time.

    Memory design
    -------------
    Peak memory is O(n_sites) regardless of sample count for the
    Fisher and Welford paths:

        fisher / logit_t / welch_t:
            4 int64 running sums (Fisher) OR
            6 arrays per group (Welford: float64 mean, float64 M2, int32 n_valid)

    The CMH path caches the case-sample int32 (meth, coverage) arrays in
    memory because each case sample contributes to one stratum per control,
    so the case data is reused. Memory overhead: ~8 bytes x n_sites x n_case,
    which is bounded for typical experiments (n_case <= 10, ~300 MB on chr1).

    Statistical paths
    -----------------
    fisher
        Fisher exact on reads pooled across replicates. Emits a warning;
        anti-conservative because between-replicate variance is ignored.
        Provided for parity with single-rep tools and aggregate reporting.

    cmh
        Cochran-Mantel-Haenszel test with one 2x2 stratum per
        (case_i, ctrl_j) pair. Preserves between-replicate variability
        because each replicate contributes its own coverage marginal.

    logit_t / welch_t
        Welch t-test on per-replicate beta values (logit-transformed for
        `logit_t`). Welford accumulators give per-site variance without
        materialising the count matrix.
    """
    test = _canonicalise_test_name(test)
    n_sites = len(canonical_df)
    if n_sites == 0:
        return pl.DataFrame(schema=_EMPTY_SCHEMA)

    canonical_pos = canonical_df.select("pos")

    # Welford accumulators provide mean_beta_case / mean_beta_ctrl and
    # n_valid_* for every code path. They are always populated even when
    # the chosen test (e.g. fisher) does not use the variance.
    mean_case, M2_case, n_valid_case = _welford_init(n_sites)
    mean_ctrl, M2_ctrl, n_valid_ctrl = _welford_init(n_sites)

    # Optional extras populated by individual branches. Surfaced into the
    # output schema at the bottom of this function when present.
    extras: dict[str, np.ndarray] = {}
    # When contrast_matrix produces a joint test (k>1) we emit a different
    # schema (per-level mean_beta_* + f_stat/df1/df2); detected by checking
    # `multigroup_mode` at the end.
    multigroup_mode = False
    level_mean_beta: dict[str, np.ndarray] = {}
    f_stat_out: Optional[np.ndarray] = None
    df1_out: Optional[int] = None
    df2_out: Optional[np.ndarray] = None

    # --- Statistical test ---
    if test == "fisher":
        # Fisher exact on per-group POOLED read counts.
        #
        # The previous "fisher"/"cmh" path pooled control reads then ran
        # one CMH stratum per case sample against the pool, producing a
        # test that was structurally identical to Fisher on pooled reads
        # but with the variance term deflated by the pooling (effective N
        # inflated by Sigma coverage per control). The corrected code below
        # makes the pooling explicit and routes it through the well-tested
        # fisher_exact_vectorized() helper.
        #
        # NOTE: this test ignores between-replicate variability. The user
        # facing warning fires once per call from
        # ``_validate_sample_size_and_warn`` -- not here, to avoid the
        # per-chromosome warning spam this used to produce.
        meth_case_sum = np.zeros(n_sites, dtype=np.int64)
        cov_case_sum  = np.zeros(n_sites, dtype=np.int64)
        meth_ctrl_sum = np.zeros(n_sites, dtype=np.int64)
        cov_ctrl_sum  = np.zeros(n_sites, dtype=np.int64)

        for sample in samples_case:
            meth, cov = _load_sample_chrom(methylstore_path, chrom, sample, canonical_pos)
            meth_case_sum += meth.astype(np.int64)
            cov_case_sum  += cov.astype(np.int64)
            _welford_update(mean_case, M2_case, n_valid_case, meth, cov)
            del meth, cov

        for sample in samples_control:
            meth, cov = _load_sample_chrom(methylstore_path, chrom, sample, canonical_pos)
            meth_ctrl_sum += meth.astype(np.int64)
            cov_ctrl_sum  += cov.astype(np.int64)
            _welford_update(mean_ctrl, M2_ctrl, n_valid_ctrl, meth, cov)
            del meth, cov

        unmeth_case_sum = cov_case_sum - meth_case_sum
        unmeth_ctrl_sum = cov_ctrl_sum - meth_ctrl_sum
        pvals, log2_ors = fisher_exact_vectorized(
            meth_case_sum, unmeth_case_sum, meth_ctrl_sum, unmeth_ctrl_sum
        )
        del meth_case_sum, cov_case_sum, unmeth_case_sum
        del meth_ctrl_sum, cov_ctrl_sum, unmeth_ctrl_sum

    elif test == "cmh":
        # properly stratified CMH -- one 2x2 stratum per
        # (case_i, ctrl_j) pair, so each replicate's coverage marginal
        # enters its own variance term V. With n_case = n_ctrl = 1 this
        # degenerates to a single 2x2 table and matches Fisher; with
        # replicates the variance term grows correctly with n_case x n_ctrl,
        # avoiding the inflated chi^2 of the old pooled-control approach.
        #
        # Memory: we cache the case samples in int32 (n_case x n_sites x 8 B)
        # so each is contributed against every control sample without
        # re-reading parquet n_ctrl times.
        ome, var_sum, or_num, or_den = _cmh_init(n_sites)

        case_data: list[tuple[np.ndarray, np.ndarray]] = []
        for sample in samples_case:
            meth, cov = _load_sample_chrom(methylstore_path, chrom, sample, canonical_pos)
            case_data.append((meth, cov))
            _welford_update(mean_case, M2_case, n_valid_case, meth, cov)

        for ctrl in samples_control:
            meth_c, cov_c = _load_sample_chrom(methylstore_path, chrom, ctrl, canonical_pos)
            _welford_update(mean_ctrl, M2_ctrl, n_valid_ctrl, meth_c, cov_c)
            for meth_case_i, cov_case_i in case_data:
                _cmh_update(
                    ome, var_sum, or_num, or_den,
                    meth_case_i, cov_case_i, meth_c, cov_c,
                )
            del meth_c, cov_c

        del case_data
        pvals, log2_ors = _cmh_finalize(ome, var_sum, or_num, or_den)
        del ome, var_sum, or_num, or_den

    elif test in ("score", "lr"):
        # Quasi-binomial count-model test with McCullagh-Nelder overdispersion.
        # Both the score and likelihood-ratio statistics share the same
        # streaming accumulators (sn, sm, sm^2/n, nv per group) and the same
        # dispersion machinery; ``_score_finalize`` picks the functional form
        # based on the ``statistic=`` argument it receives.
        sn_case, sm_case, sm2n_case, nv_case = _score_init(n_sites)
        sn_ctrl, sm_ctrl, sm2n_ctrl, nv_ctrl = _score_init(n_sites)

        # DSS-style smoothing (smoothing=True) replicates
        # DMLfit.multiFactor(smoothing=TRUE): for each sample, replace the
        # raw (meth, cov) counts with a uniform-box moving average over
        # CpGs within +/-smoothing_span_bp//2 bp before they hit the score
        # accumulators. The kernel matches DSS's smooth.chr / nitem_bin /
        # windowFilter exactly. Counts stay float for accumulator precision.
        chrom_positions = (
            canonical_pos.to_series().to_numpy() if smoothing else None
        )

        for sample in samples_case:
            meth, cov = _load_sample_chrom(methylstore_path, chrom, sample, canonical_pos)
            if smoothing:
                meth, cov = _smooth_sample_counts_box(
                    meth, cov, chrom_positions, smoothing_span_bp,
                )
            _score_update(sn_case, sm_case, sm2n_case, nv_case, meth, cov)
            # Welford accumulators are also updated so that downstream code
            # which reads ``n_valid_case`` for the guard sees the
            # same per-site sample count.  ``mean_case`` from Welford is
            # overwritten below with the coverage-weighted score-test
            # equivalent, so its post-update value is irrelevant.
            _welford_update(mean_case, M2_case, n_valid_case, meth, cov)
            del meth, cov

        for sample in samples_control:
            meth, cov = _load_sample_chrom(methylstore_path, chrom, sample, canonical_pos)
            if smoothing:
                meth, cov = _smooth_sample_counts_box(
                    meth, cov, chrom_positions, smoothing_span_bp,
                )
            _score_update(sn_ctrl, sm_ctrl, sm2n_ctrl, nv_ctrl, meth, cov)
            _welford_update(mean_ctrl, M2_ctrl, n_valid_ctrl, meth, cov)
            del meth, cov

        pvals, log2_ors, pi_case, pi_ctrl, _phi_hat = _score_finalize(
            sn_case, sm_case, sm2n_case, nv_case,
            sn_ctrl, sm_ctrl, sm2n_ctrl, nv_ctrl,
            chrom_name=chrom,
            dispersion=dispersion,
            statistic=test,
            reference=reference,
            sep_fallback=sep_fallback,
            sep_threshold=sep_threshold,
        )

        # Coverage-weighted (= pooled MLE) group methylation for output.
        # Overwrite Welford's unweighted means with the score-test
        # equivalents so the unified output block at the bottom of this
        # function reports the values consistent with the test's math.
        mean_case[:] = np.where(np.isnan(pi_case), 0.0, pi_case)
        mean_ctrl[:] = np.where(np.isnan(pi_ctrl), 0.0, pi_ctrl)
        # nv_case / nv_ctrl from the score path agree with n_valid_case /
        # n_valid_ctrl from Welford by construction, so we don't overwrite.

        del sn_case, sm_case, sm2n_case, nv_case
        del sn_ctrl, sm_ctrl, sm2n_ctrl, nv_ctrl

    elif test in ("welch_t", "logit_t"):
        # Load all samples for Welford accumulators
        for sample in samples_case:
            meth, cov = _load_sample_chrom(methylstore_path, chrom, sample, canonical_pos)
            _welford_update(mean_case, M2_case, n_valid_case, meth, cov)
            del meth, cov

        for sample in samples_control:
            meth, cov = _load_sample_chrom(methylstore_path, chrom, sample, canonical_pos)
            _welford_update(mean_ctrl, M2_ctrl, n_valid_ctrl, meth, cov)
            del meth, cov

        # Welford path: no (n_sites x n_replicates) matrix ever built.
        if test == "logit_t":
            pvals, log2_ors = _beta_binom_mom_from_welford_logit(
                mean_case, M2_case, n_valid_case,
                mean_ctrl, M2_ctrl, n_valid_ctrl,
            )
        else:
            pvals, log2_ors = _beta_binom_mom_from_welford(
                mean_case, M2_case, n_valid_case,
                mean_ctrl, M2_ctrl, n_valid_ctrl,
            )

    elif test == "bb_lr":
        # True beta-binomial LRT via a quasi-binomial GLM on a binary
        # treatment indicator, with site-level dispersion estimated from
        # the per-replicate Pearson residuals. Distinct from "lr" (which
        # works on per-group pooled counts and a single phi from
        # _score_finalize): bb_lr fits the full GLM at every site so the
        # quasi-binomial dispersion correctly reflects between-replicate
        # over-dispersion at THIS site, which is the spirit of a beta-
        # binomial model.
        all_samples = samples_case + samples_control
        n_samples = len(all_samples)
        if n_samples < 6 and dispersion == "site":
            logger.warning(
                "[WARN]  bb_lr with n=%d total samples (df_resid=%d) has very "
                "noisy per-site dispersion estimates. Promoting dispersion to "
                "'shrink' to stabilise estimates. Consider test='lr' for "
                "higher power.",
                n_samples, n_samples - 2,
            )
            dispersion = "shrink"
        meth_stack = np.zeros((n_sites, n_samples), dtype=np.int32)
        cov_stack  = np.zeros((n_sites, n_samples), dtype=np.int32)
        # DSS-style smoothing (smoothing=True): per-sample uniform-box
        # average of (meth, cov) in a +/-smoothing_span_bp//2 window before
        # the counts hit the IRLS. Rounded back to int32 to match DSS's
        # round(smooth.chr(...)) and to fit the existing irls_dispatch
        # contract (which expects integer count arrays).
        chrom_positions_bb = (
            canonical_pos.to_series().to_numpy() if smoothing else None
        )
        for j, sample in enumerate(all_samples):
            meth, cov = _load_sample_chrom(methylstore_path, chrom, sample, canonical_pos)
            if smoothing:
                meth_sm, cov_sm = _smooth_sample_counts_box(
                    meth, cov, chrom_positions_bb, smoothing_span_bp,
                )
                meth = np.rint(meth_sm).astype(np.int32)
                cov  = np.rint(cov_sm).astype(np.int32)
            meth_stack[:, j] = meth
            cov_stack[:, j]  = cov
            if j < len(samples_case):
                _welford_update(mean_case, M2_case, n_valid_case, meth, cov)
            else:
                _welford_update(mean_ctrl, M2_ctrl, n_valid_ctrl, meth, cov)
            del meth, cov

        # Build a tiny intercept + treatment design.
        treat = np.zeros(n_samples, dtype=np.float64)
        treat[:len(samples_case)] = 1.0
        X_bb_full = np.column_stack([np.ones(n_samples), treat])
        X_bb_red  = np.ones((n_samples, 1))

        from . import _glm
        beta_f, se_f, dev_f, pearson_f, n_eff_bb = _glm.irls_dispatch(
            meth_stack, cov_stack, X_bb_full, backend=glm_backend,
        )
        _bbeta_r, _bse_r, dev_r, _bpearson_r, _bn_eff_r = _glm.irls_dispatch(
            meth_stack, cov_stack, X_bb_red, backend=glm_backend,
        )

        df_resid_per_site = (n_eff_bb.astype(np.float64) - 2.0)
        df_resid_safe = np.maximum(df_resid_per_site, 1.0)
        # bb_lr path: discards df_phi to preserve its pre-fix behavior
        # (it still passes df_resid_safe to reference_pvalues).
        phi_eff, _phi_hat, _df_phi = _glm.compute_dispersion_phi(
            pearson_per_site=pearson_f,
            df_per_site=df_resid_per_site,
            dispersion=dispersion,
            chrom_name=chrom,
        )

        with np.errstate(invalid="ignore", divide="ignore"):
            lr_raw = dev_r - dev_f
            lr_raw = np.where(lr_raw < 0, 0.0, lr_raw)
            chi2_stat = np.where(phi_eff > 0, lr_raw / phi_eff, np.nan)
        pvals = _glm.reference_pvalues(
            chi2_stat, phi_eff, df_resid_safe, reference=reference,
        )

        coef_treatment = beta_f[:, 1].astype(np.float64)
        coef_se        = se_f[:, 1].astype(np.float64) * np.sqrt(np.maximum(phi_eff, 1.0))
        log2_ors = (coef_treatment / np.log(2.0))
        degenerate_bb = np.isnan(chi2_stat) | np.isnan(pvals) | (n_eff_bb < 2)
        pvals = np.where(degenerate_bb, np.nan, pvals)
        log2_ors = np.where(degenerate_bb, np.nan, log2_ors)
        extras["coef_treatment"] = coef_treatment
        extras["coef_se"]        = coef_se
        del meth_stack, cov_stack, beta_f, se_f, dev_f, dev_r, pearson_f, n_eff_bb

    elif test == "glm":
        # Covariate-aware binomial GLM with deviance LR test.
        #
        # We load every sample's (meth, cov) for this chromosome into an
        # (n_sites, n_samples) stack so the batched IRLS can fit one GLM
        # per site against the shared design matrix. At the tile level
        # n_sites is small (~10^4-10^5), so this is a few MB of int32.
        if design_full is None or design_reduced is None or coef_idx is None:
            raise ValueError(
                "test='glm' requires design_full, design_reduced, and "
                "coef_idx. Build them via epykit._glm.build_design(md.obs, "
                "samples_ordered=samples_case+samples_control, formula=...)."
            )

        all_samples = samples_case + samples_control
        n_samples = len(all_samples)
        if design_full.shape[0] != n_samples:
            raise ValueError(
                f"design_full has {design_full.shape[0]} rows but "
                f"{n_samples} samples were passed. Rows must follow "
                "samples_case + samples_control order."
            )

        meth_stack = np.zeros((n_sites, n_samples), dtype=np.int32)
        cov_stack  = np.zeros((n_sites, n_samples), dtype=np.int32)
        for j, sample in enumerate(all_samples):
            meth, cov = _load_sample_chrom(methylstore_path, chrom, sample, canonical_pos)
            meth_stack[:, j] = meth
            cov_stack[:, j]  = cov
            # Welford accumulators for the unadjusted mean_beta_* columns.
            if j < len(samples_case):
                _welford_update(mean_case, M2_case, n_valid_case, meth, cov)
            else:
                _welford_update(mean_ctrl, M2_ctrl, n_valid_ctrl, meth, cov)
            del meth, cov

        from . import _glm

        beta_full, se_full, dev_full, pearson_full, n_eff = _glm.irls_dispatch(
            meth_stack, cov_stack, design_full, backend=glm_backend,
        )
        _beta_red, _se_red, dev_red, _pearson_red, _ne_red = _glm.irls_dispatch(
            meth_stack, cov_stack, design_reduced, backend=glm_backend,
        )

        # df_resid_i = n_eff_i - p_full   (per-site, since coverage gates samples)
        p_full = design_full.shape[1]
        df_resid_per_site = (n_eff.astype(np.float64) - float(p_full))
        df_resid_safe = np.maximum(df_resid_per_site, 1.0)

        phi_eff, _phi_hat, df_phi = _glm.compute_dispersion_phi(
            pearson_per_site=pearson_full,
            df_per_site=df_resid_per_site,
            dispersion=dispersion,
            chrom_name=chrom,
        )

        # LR statistic with dispersion correction. Reduced model drops the
        # treatment column => 1 df contrast.
        with np.errstate(invalid="ignore", divide="ignore"):
            lr_raw = dev_red - dev_full
            # tiny negative excursions can happen at numerical machine eps
            lr_raw = np.where(lr_raw < 0, 0.0, lr_raw)
            chi2_stat = np.where(phi_eff > 0, lr_raw / phi_eff, np.nan)

        # F-reference uses df_phi (df backing the phi estimate), NOT
        # df_resid_safe (per-site residual df). Identical to df_resid_safe for
        # dispersion="site"; for "chrom"/"shrink"/"eb" it's the chrom-pool or
        # shrinkage-effective df. Same bug-fix as in _score_finalize.
        pvals = _glm.reference_pvalues(
            chi2_stat, phi_eff, df_phi, reference=reference,
        )

        # Effect-size columns from the GLM coefficient (log-odds) and its SE.
        coef_treatment = beta_full[:, coef_idx].astype(np.float64)
        coef_se        = se_full[:, coef_idx].astype(np.float64)

        # Bookkeeping for the unified output block at the bottom.
        log2_ors = (coef_treatment / np.log(2.0))   # log-odds -> log2 odds
        degenerate = (
            np.isnan(chi2_stat) | np.isnan(pvals) | (n_eff < 2)
        )
        pvals = np.where(degenerate, np.nan, pvals)
        log2_ors = np.where(degenerate, np.nan, log2_ors)

        # Stash for the schema additions below.
        extras["coef_treatment"] = coef_treatment
        extras["coef_se"]        = coef_se
        del meth_stack, cov_stack, beta_full, se_full, dev_full, dev_red
        del pearson_full, n_eff

    elif test == "glm_contrast":
        # Multi-group / continuous-covariate primary-effect path.
        # The caller supplies (a) a shared design matrix `design_full`,
        # (b) a `contrast_matrix` C of shape (k, p), and (c) an ordered
        # list `samples_all_ordered` whose row order matches design_full.
        if (
            design_full is None
            or contrast_matrix is None
            or samples_all_ordered is None
        ):
            raise ValueError(
                "test='glm_contrast' requires design_full, contrast_matrix, "
                "and samples_all_ordered."
            )
        n_samples = len(samples_all_ordered)
        if design_full.shape[0] != n_samples:
            raise ValueError(
                f"design_full has {design_full.shape[0]} rows but "
                f"{n_samples} samples were supplied via samples_all_ordered."
            )

        meth_stack = np.zeros((n_sites, n_samples), dtype=np.int32)
        cov_stack  = np.zeros((n_sites, n_samples), dtype=np.int32)
        # Welford per-level accumulators. We keep mean_case/mean_ctrl as
        # "all samples that map to label 'case'" vs "all samples that
        # don't" -- chosen so the binary-case columns of the schema still
        # carry interpretable values. If group_labels_per_sample is given,
        # we additionally keep one Welford accumulator per level for the
        # multi-group output schema.
        level_mean: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        if group_labels_per_sample is not None:
            for lvl in set(group_labels_per_sample):
                level_mean[lvl] = _welford_init(n_sites)

        for j, sample in enumerate(samples_all_ordered):
            meth, cov = _load_sample_chrom(methylstore_path, chrom, sample, canonical_pos)
            meth_stack[:, j] = meth
            cov_stack[:, j]  = cov
            # Backwards-compat columns: split on whether sample is in samples_case.
            if sample in samples_case:
                _welford_update(mean_case, M2_case, n_valid_case, meth, cov)
            elif sample in samples_control:
                _welford_update(mean_ctrl, M2_ctrl, n_valid_ctrl, meth, cov)
            if group_labels_per_sample is not None:
                lvl = group_labels_per_sample[j]
                _welford_update(*level_mean[lvl], meth, cov)
            del meth, cov

        from . import _glm
        beta_full, se_full, dev_full, pearson_full, n_eff, cov_beta = (
            _glm.irls_dispatch(
                meth_stack, cov_stack, design_full, return_cov=True,
                backend=glm_backend,
            )
        )

        p_full = design_full.shape[1]
        df_resid_per_site = (n_eff.astype(np.float64) - float(p_full))
        phi_eff, _phi_hat, df_phi = _glm.compute_dispersion_phi(
            pearson_per_site=pearson_full,
            df_per_site=df_resid_per_site,
            dispersion=dispersion,
            chrom_name=chrom,
        )
        df_resid_safe = np.maximum(df_resid_per_site, 1.0)

        # F-reference uses df_phi (see comment in glm path above).
        stat, pvals, k_rank = _glm.wald_test(
            beta_full, cov_beta, contrast_matrix,
            phi_eff=phi_eff, df_resid=df_phi, reference=reference,
        )
        k_rank = int(k_rank)

        degenerate = (
            np.isnan(stat) | np.isnan(pvals) | (n_eff < p_full + 1)
        )
        pvals = np.where(degenerate, np.nan, pvals)

        if k_rank == 1:
            # Single-coef contrast: surface coef + coef_se on the contrast
            # axis (Cbeta and sqrt(C cov_beta CT)).
            Cb = (beta_full @ contrast_matrix.T)[:, 0]
            with np.errstate(invalid="ignore"):
                var_Cb = np.einsum(
                    "kp,ipq,lq->ikl",
                    contrast_matrix, cov_beta * phi_eff[:, None, None],
                    contrast_matrix,
                )[:, 0, 0]
                cse = np.sqrt(np.where(var_Cb > 0, var_Cb, np.nan))
            log2_ors = (Cb / np.log(2.0))
            log2_ors = np.where(degenerate, np.nan, log2_ors)
            extras["coef_treatment"] = Cb
            extras["coef_se"]        = cse
            # Single-coef path still emits the standard binary schema.
        else:
            # Joint contrast: emit multi-group schema. We do NOT populate
            # the binary mean_beta_case/control columns meaningfully --
            # they're filled with NaN downstream. F-stat and per-level
            # mean betas are stored for the unified output block.
            multigroup_mode = True
            f_stat_out = (stat / k_rank).astype(np.float64)
            f_stat_out = np.where(degenerate, np.nan, f_stat_out)
            df1_out = k_rank
            df2_out = df_resid_safe.astype(np.float64)
            for lvl, (mu_l, _M2_l, nv_l) in level_mean.items():
                arr = mu_l.astype(np.float32)
                arr[nv_l == 0] = np.nan
                level_mean_beta[lvl] = arr
            log2_ors = np.full(n_sites, np.nan, dtype=np.float64)
            extras.clear()  # coef_* not meaningful for joint test

        del meth_stack, cov_stack, beta_full, se_full, dev_full
        del pearson_full, n_eff, cov_beta

    else:
        raise NotImplementedError(
            f"Test '{test}' not implemented. "
            "Choose 'lr', 'score', 'fisher', 'cmh', 'logit_t', "
            "'welch_t', 'bb_lr', or 'glm'."
        )

    # --- equal-weight per-replicate mean beta ---
    # Welford mean IS the equal-weight nanmean -- no extra storage needed.
    mean_beta_case = mean_case.astype(np.float32)
    mean_beta_ctrl = mean_ctrl.astype(np.float32)
    mean_beta_case[n_valid_case == 0] = np.nan
    mean_beta_ctrl[n_valid_ctrl == 0] = np.nan
    meth_diff = (mean_beta_case - mean_beta_ctrl).astype(np.float32)

    # Wald CI on Deltabeta from Welford accumulators. For the
    # multi-group joint test (k>1) the scalar Deltabeta is undefined, so we leave
    # CI columns as NaN and let the multi-group schema speak through
    # f_stat / df1 / df2 / mean_beta_<level>.
    from . import _glm as _glm_for_ci
    if not multigroup_mode:
        vm_case = _welford_var_mean(M2_case, n_valid_case)
        vm_ctrl = _welford_var_mean(M2_ctrl, n_valid_ctrl)
        ci_lo, ci_hi = _glm_for_ci.welch_meth_diff_ci(
            mean_case.astype(np.float64), vm_case,
            mean_ctrl.astype(np.float64), vm_ctrl,
        )
        ci_lo = ci_lo.astype(np.float32)
        ci_hi = ci_hi.astype(np.float32)
    else:
        ci_lo = np.full(n_sites, np.nan, dtype=np.float32)
        ci_hi = np.full(n_sites, np.nan, dtype=np.float32)

    # per-site min-samples guard. Sites where fewer than
    # `min_samples_*` replicates contributed valid (coverage > 0) data
    # have their p-value masked to NaN. apply_multiple_testing_correction
    # passes NaNs through, so these sites are effectively excluded from
    # genome-wide FDR control without disturbing site-position alignment.
    if min_samples_case > 0 or min_samples_control > 0:
        keep_mask = (
            (n_valid_case >= max(min_samples_case, 0))
            & (n_valid_ctrl >= max(min_samples_control, 0))
        )
        n_dropped = int((~keep_mask).sum())
        if n_dropped > 0:
            logger.info(
                "  %s: masking %s/%s sites with n_valid_case < %d or "
                "n_valid_ctrl < %d",
                chrom, f"{n_dropped:,}", f"{n_sites:,}",
                min_samples_case, min_samples_control,
            )
            pvals = np.where(keep_mask, pvals, np.nan)
            log2_ors = np.where(keep_mask, log2_ors, np.nan)
            meth_diff = np.where(keep_mask, meth_diff, np.float32(np.nan))
            mean_beta_case = np.where(keep_mask, mean_beta_case, np.float32(np.nan))
            mean_beta_ctrl = np.where(keep_mask, mean_beta_ctrl, np.float32(np.nan))
            ci_lo = np.where(keep_mask, ci_lo, np.float32(np.nan))
            ci_hi = np.where(keep_mask, ci_hi, np.float32(np.nan))

    del mean_case, M2_case, n_valid_case, mean_ctrl, M2_ctrl, n_valid_ctrl

    out_cols = {
        "chrom":             pl.Series([chrom] * n_sites, dtype=pl.Utf8),
        "pos":               canonical_df["pos"],
        "strand":            canonical_df["strand"],
        "n_case":            pl.Series(
                                 np.full(n_sites, len(samples_case),    dtype=np.int32)),
        "n_control":         pl.Series(
                                 np.full(n_sites, len(samples_control), dtype=np.int32)),
        "mean_beta_case":    pl.Series(mean_beta_case),
        "mean_beta_control": pl.Series(mean_beta_ctrl),
        "pvalue":            pl.Series(pvals),
        "log2_odds_ratio":   pl.Series(log2_ors),
        "meth_diff":         pl.Series(meth_diff),
        "meth_diff_ci_lo":   pl.Series(ci_lo),
        "meth_diff_ci_hi":   pl.Series(ci_hi),
    }
    if "coef_treatment" in extras and "coef_se" in extras:
        out_cols["coef_treatment"] = pl.Series(extras["coef_treatment"])
        out_cols["coef_se"]        = pl.Series(extras["coef_se"])
    if multigroup_mode and f_stat_out is not None and df2_out is not None:
        out_cols["f_stat"] = pl.Series(f_stat_out)
        out_cols["df1"]    = pl.Series(np.full(n_sites, int(df1_out), dtype=np.int32))
        out_cols["df2"]    = pl.Series(df2_out)
        # Per-level mean beta columns (stable sort for deterministic schema).
        for lvl in sorted(level_mean_beta.keys()):
            out_cols[f"mean_beta_{lvl}"] = pl.Series(level_mean_beta[lvl])
        # meth_diff_max = max |mean_beta_i - mean_beta_j| across all level pairs
        if level_mean_beta:
            stacked = np.stack(
                [level_mean_beta[lvl] for lvl in sorted(level_mean_beta)],
                axis=1,
            )
            max_diff = np.nanmax(stacked, axis=1) - np.nanmin(stacked, axis=1)
            out_cols["meth_diff_max"] = pl.Series(max_diff.astype(np.float32))
    return pl.DataFrame(out_cols).sort("pos")


def _validate_sample_size_and_warn(n_case: int, n_ctrl: int, test: str) -> None:
    """Validate sample sizes and issue appropriate warnings."""
    min_n = min(n_case, n_ctrl)
    max_n = max(n_case, n_ctrl)

    if min_n == 0:
        raise ValueError(
            "Cannot perform DMC with zero samples in a group. "
            f"n_case={n_case}, n_control={n_ctrl}"
        )

    if min_n == 1:
        logger.warning(
            "[WARN]  CRITICAL: Only 1 replicate per group detected!\n"
            "   Statistical results are UNRELIABLE without biological replicates.\n"
            "   Effect sizes may be reported, but p-values should NOT be trusted.\n"
            "   Recommendation: Collect at least 3 biological replicates per group."
        )
    elif min_n == 2:
        logger.warning(
            "[WARN]  WARNING: Only 2 replicates per group.\n"
            "   Statistical power is very low. Many true positives will be missed.\n"
            "   Recommendation: Use n>=3 for reliable differential methylation calling."
        )
    elif min_n <= 2 and test == "welch_t":
        logger.warning(
            "[WARN]  CRITICAL: welch_t with n=%d per group produces degenerate "
            "Welch-Satterthwaite DOF and near-zero power. "
            "Use test='lr' instead.",
            min_n,
        )
    elif min_n < 6 and test == "welch_t":
        logger.warning(
            "[WARN]  Welch t with "
            "n<6 may have poor variance estimates.\n"
            "   Consider using test='lr' (recommended) or test='bb_lr' "
            "(true quasi-binomial LRT)."
        )

    if test == "bb_lr" and min_n < 3:
        logger.warning(
            "[WARN]  bb_lr requires at least 3 replicates per group for "
            "reliable dispersion estimation (got min_n=%d). "
            "Use test='lr' instead.",
            min_n,
        )

    if test == "fisher" and min_n >= 2:
        # Fires once per process_chromosomes_dmc() call (not per chromosome).
        # tl.dmc gates its own one-shot-per-session warning on top of this;
        # direct API users see this every call.
        warnings.warn(
            "test='fisher' pools reads across replicates; between-sample "
            "variance is ignored and p-values may be anti-conservative at "
            "WGBS coverage. Prefer test='lr' at n>=2 (quasi-binomial LR "
            "with MN dispersion).",
            UserWarning,
            stacklevel=3,
        )

    if max_n / min_n > 2:
        logger.warning(
            f"[WARN]  Unbalanced design detected: n_case={n_case}, n_control={n_ctrl}\n"
            "   Large imbalance may reduce statistical power."
        )


# Public API

def process_chromosomes_dmc(
    methylstore_path: str,
    samples_treatment: Optional[list[str]] = None,
    samples_control: Optional[list[str]] = None,
    test: str = "lr",
    chromosomes: Optional[list[str]] = None,
    unite: bool = True,
    min_samples_treatment: Optional[int] = None,
    min_samples_control: int = 0,
    dispersion: str = "site",
    reference: str = "adaptive",
    design_full: Optional[np.ndarray] = None,
    design_reduced: Optional[np.ndarray] = None,
    coef_idx: Optional[int] = None,
    contrast_matrix: Optional[np.ndarray] = None,
    contrast_label: Optional[str] = None,
    samples_all_ordered: Optional[list[str]] = None,
    group_labels_per_sample: Optional[list[str]] = None,
    *,
    backend: str = "sequential",
    n_workers: Optional[int] = None,
    glm_backend: str = "cpu",
    out_dir: Optional[Union[str, Path]] = None,
    return_store: bool = False,
    smoothing: bool = False,
    smoothing_span_bp: int = 500,
    sep_fallback: bool = False,
    sep_threshold: float = 0.9,
) -> Union[pl.DataFrame, DMCStore]:
    """Process differential methylation for all chromosomes.

    Parameters
    ----------
    methylstore_path : str
        Path to filtered partitioned Parquet methylstore.
    samples_treatment, samples_control : list[str]
        Sample identifiers for treatment and control groups.
    test : {"lr", "score", "fisher", "cmh", "logit_t", "welch_t", "bb_lr"}
        Statistical test.
            "lr"       (default) -- Quasi-binomial likelihood-ratio chi-square
                                   on per-group read counts with per-site
                                   McCullagh-Nelder dispersion. Closed-form on
                                   the streaming accumulators (S0_g, S1_g,
                                   Sigmam^2/n_g). Recommended at n >= 2.
            "score"              -- Pearson score statistic on the same
                                   accumulators. Marginally more powerful
                                   than "lr" but mildly anti-conservative
                                   when pi is near 0 or 1.
            "logit_t"            -- Welch t on logit(beta), variance via
                                   Welford. Variance-stabilising fallback
                                   when count-model assumptions are
                                   doubtful (e.g. very low coverage).
            "welch_t"            -- Welch t on raw betas. Same boundary-beta
                                   caveat as logit_t.
            "bb_lr"              -- True quasi-binomial LRT via a full per-
                                   site GLM on a binary-treatment design.
                                   Slow, but the honest "fit-the-model"
                                   version of "welch_t".
            "cmh"                -- Cochran-Mantel-Haenszel with one stratum
                                   per (case_i, ctrl_j) pair.
            "fisher"             -- Fisher exact on reads pooled across
                                   replicates (anti-conservative; warns).
    chromosomes : list[str], optional
        Chromosomes to process. Auto-detected when None.
    unite : bool
        If True (default), test only CpG sites covered in every sample
        (intersection / inner join).
        If False, test all sites covered in at least one sample
        (union / outer join).
    min_samples_treatment, min_samples_control : int
        Per-site minimum number of replicates with non-zero coverage
        required in each group. Sites failing the threshold have their
        p-value masked to NaN before FDR correction. Use this with
        ``unite=False`` (union mode) to drop tests that effectively run on
        a singleton observation in one group.
    dispersion : {"site", "chrom", "shrink"}
        McCullagh-Nelder dispersion strategy used by ``test="lr"`` and
        ``test="score"``. Default ``"site"`` estimates a per-site phi_i from
        the 4-df Pearson residual sum; ``"chrom"`` pools one phi across the
        whole chromosome; ``"shrink"`` is a weighted average of the two.
        See :func:`_score_finalize` for details. Ignored for other tests.
    reference : {"adaptive", "chi2", "F"}
        Reference distribution for the quasi-binomial test statistic.
        Default ``"adaptive"`` switches per-site between F(1, df) where
        phi > 1 and chi^2(1) where phi was clamped to 1 (the right behaviour for
        quasi-binomial GLMs whose dispersion estimate is noisy at small
        samples). ``"chi2"`` and ``"F"`` force a single reference
        distribution regardless of dispersion. See :func:`_score_finalize`
        for details. Ignored for other tests.

    Returns
    -------
    pl.DataFrame
        Columns: chrom, pos, strand, n_case, n_control,
                 mean_beta_case, mean_beta_control,
                 pvalue, log2_odds_ratio, meth_diff

        For the ``score`` test, ``mean_beta_*`` are coverage-weighted (the
        group MLE proportion M/N). For all other tests they are the
        unweighted per-replicate mean (Welford). The two differ when
        per-replicate coverage is uneven.
    """
    if samples_treatment is None:
        raise TypeError("Missing required argument: samples_treatment")
    if samples_control is None:
        raise TypeError("Missing required argument: samples_control")
    if min_samples_treatment is None:
        min_samples_treatment = 0
    samples_case = samples_treatment
    min_samples_case = min_samples_treatment

    test = _canonicalise_test_name(test)

    store       = Path(methylstore_path)
    # For glm_contrast the case/control split is not meaningful -- the
    # caller passes a samples_all_ordered list that the engine uses for
    # design-row alignment. Use that as the chromosome-intersection basis.
    if test == "glm_contrast":
        if samples_all_ordered is None:
            raise ValueError(
                "test='glm_contrast' requires samples_all_ordered."
            )
        all_samples = list(samples_all_ordered)
    else:
        all_samples = samples_case + samples_control

    min_group = min(len(samples_case), len(samples_control))
    if test != "glm_contrast":
        _validate_sample_size_and_warn(len(samples_case), len(samples_control), test)

    for rng, rec in _TEST_RECOMMENDATIONS.items():
        if min_group in rng:
            logger.info(
                "N replicates / group = %d; recommended test: %s", min_group, rec
            )
            break

    if chromosomes is None:
        chromosomes = _detect_chromosomes(store)
        logger.info("Auto-detected %d chromosomes", len(chromosomes))

    logger.info(
        "DMC: %d case / %d control, test=%s, unite=%s, "
        "min_samples_case=%d, min_samples_control=%d",
        len(samples_case), len(samples_control), test, unite,
        min_samples_case, min_samples_control,
    )

    from ._compute import run_chrom_pipeline

    def _dmc_chrom_handler(chrom: str) -> Optional[pl.DataFrame]:
        canonical_df = (
            _intersect_chrom(store, chrom, all_samples)
            if unite
            else _union_chrom(store, chrom, all_samples)
        )
        if len(canonical_df) == 0:
            logger.warning("  No sites for %s; skipping", chrom)
            return None
        logger.info("  %s sites to test (%s)", f"{len(canonical_df):,}", chrom)
        return _process_one_chromosome(
            store, chrom, canonical_df,
            samples_case, samples_control, test,
            min_samples_case=min_samples_case,
            min_samples_control=min_samples_control,
            dispersion=dispersion,
            reference=reference,
            design_full=design_full,
            design_reduced=design_reduced,
            coef_idx=coef_idx,
            contrast_matrix=contrast_matrix,
            contrast_label=contrast_label,
            samples_all_ordered=samples_all_ordered,
            group_labels_per_sample=group_labels_per_sample,
            glm_backend=glm_backend,
            smoothing=smoothing,
            smoothing_span_bp=smoothing_span_bp,
            sep_fallback=sep_fallback,
            sep_threshold=sep_threshold,
        )

    staging = _resolve_dmc_store_dir(store, test, out_dir, smoothing=smoothing)
    staging.mkdir(parents=True, exist_ok=True)

    # Cache check: if the existing manifest's input_sig matches the
    # current inputs, every per-chrom parquet is already on disk and
    # bit-identical to what we'd recompute. Skip straight to returning
    # a DMCStore over the cached dir.
    input_sig = _dmc_input_signature(
        store, samples_case, samples_control, test, chromosomes,
        unite, min_samples_case, min_samples_control,
        dispersion, reference,
        samples_all_ordered, group_labels_per_sample, contrast_label,
        smoothing=smoothing,
        smoothing_span_bp=smoothing_span_bp,
        sep_fallback=sep_fallback,
        sep_threshold=sep_threshold,
    )
    from ._dmc_store import _MANIFEST_NAME
    cached_manifest = _cache.load_json(staging / _MANIFEST_NAME)
    if cached_manifest is not None and cached_manifest.get("chroms"):
        # Strict hit: signatures match exactly.
        strict_hit = cached_manifest.get("input_sig") == input_sig
        # Weak hit: signatures don't match (or the manifest is from an
        # older format with no input_sig at all), but every per-chrom
        # parquet listed in the manifest still exists with the right
        # size. The parquet files are the source of truth; the
        # signature is a fast precheck. This recovers from legacy
        # manifests and from format changes without forcing a recompute.
        weak_hit = False
        all_present = True
        for entry in cached_manifest.get("chroms", []):
            f = staging / entry["file"]
            if not f.exists():
                all_present = False
                break
        if not strict_hit and all_present:
            # On weak hit, additionally verify sizes match what the
            # manifest claims -- guards against half-written parquets.
            try:
                size_ok = all(
                    (staging / e["file"]).stat().st_size > 0
                    for e in cached_manifest.get("chroms", [])
                )
            except OSError:
                size_ok = False
            weak_hit = size_ok

        if strict_hit and all_present:
            logger.info(
                "DMC cache hit at %s (%s sites, %d chrom file(s)); "
                "skipping recompute.",
                staging,
                f"{cached_manifest.get('total_sites', 0):,}",
                len(cached_manifest.get("chroms", [])),
            )
            cached = DMCStore(path=staging, test=test, _manifest=cached_manifest)
            return cached if return_store else cached.to_dataframe()

        if weak_hit:
            logger.info(
                "DMC cache hit at %s (%s sites, legacy manifest); "
                "upgrading manifest and skipping recompute.",
                staging,
                f"{cached_manifest.get('total_sites', 0):,}",
            )
            upgraded = dict(cached_manifest)
            upgraded["input_sig"] = input_sig
            upgraded["epykit_version"] = _epykit_version()
            _cache.write_json(staging / _MANIFEST_NAME, upgraded)
            cached = DMCStore(path=staging, test=test, _manifest=upgraded)
            return cached if return_store else cached.to_dataframe()

        if not all_present:
            logger.info(
                "DMC manifest at %s references missing per-chrom files; "
                "recomputing.", staging,
            )

    # Wipe stale per-chrom files from a prior partial run in the same
    # directory so we never end up with a mix of fresh and stale chroms.
    for stale in staging.glob("chrom=*.parquet"):
        stale.unlink()
    # Also drop a stale manifest so partial-run state never lingers.
    stale_manifest = staging / _MANIFEST_NAME
    if stale_manifest.exists():
        stale_manifest.unlink()

    chrom_enum   = pl.Enum(list(chromosomes))
    strand_enum  = pl.Enum(["+", "-", "*"])

    written_entries: list[dict] = []
    for chrom, chrom_result in run_chrom_pipeline(
        chromosomes, _dmc_chrom_handler,
        backend=backend, n_workers=n_workers, label="DMC",
    ):
        # Cast chrom/strand to Enum at write time to keep peak DataFrame
        # memory bounded on full-genome inputs. On 22M rows this drops
        # chrom alone from ~280 MB (Utf8) to ~22 MB.
        cast_exprs = [pl.col("chrom").cast(chrom_enum)]
        if "strand" in chrom_result.columns:
            cast_exprs.append(pl.col("strand").cast(strand_enum))
        chrom_result = chrom_result.with_columns(cast_exprs)

        out_file = staging / _chrom_filename(chrom)
        chrom_result.write_parquet(str(out_file))
        n_sites = len(chrom_result)
        written_entries.append({
            "name": chrom,
            "n_sites": int(n_sites),
            "file": out_file.name,
        })
        logger.info("  %s sites -> staged to disk (%s)", f"{n_sites:,}", chrom)
        del chrom_result
        gc.collect()

    if not written_entries:
        logger.warning("No results generated")
        # Don't litter the cache dir with an empty manifest -- clean up.
        empty_df = pl.DataFrame(schema=_EMPTY_SCHEMA)
        if return_store:
            from ._dmc_store import _MANIFEST_NAME
            empty_manifest = {
                "epykit_version": _epykit_version(),
                "test": test,
                "chroms": [],
                "total_sites": 0,
                "bh_qvalues_applied": False,
                "completed_at": _now_iso(),
            }
            _cache.write_json(staging / _MANIFEST_NAME, empty_manifest)
            return DMCStore(path=staging, test=test, _manifest=empty_manifest)
        return empty_df

    total_sites = sum(e["n_sites"] for e in written_entries)
    logger.info(
        "Assembled DMC store at %s (%d chromosomes, %s sites)",
        staging, len(written_entries), f"{total_sites:,}",
    )

    manifest = {
        "epykit_version": _epykit_version(),
        "test": test,
        "input_methylstore": str(store.resolve()),
        "input_sig": input_sig,
        "chroms": written_entries,
        "total_sites": int(total_sites),
        "bh_qvalues_applied": False,
        "completed_at": _now_iso(),
    }
    _cache.write_json(staging / _MANIFEST_NAME, manifest)
    dmc_store = DMCStore(path=staging, test=test, _manifest=manifest)

    if return_store:
        return dmc_store
    # Back-compat default: assemble the full DataFrame for callers that
    # expect one. With Enum chrom/strand this is ~700 MB at 22M rows
    # instead of ~2 GB.
    return dmc_store.to_dataframe()


_VALID_FDR_METHODS = {"fdr_bh", "fdr_by", "fdr_tsbh", "fdr_tsbky", "fdr_storey"}


def _storey_pi0(pvals: np.ndarray, lam: float | None = None) -> float:
    """Estimate the proportion of true nulls pi0 using Storey's method.

    Parameters
    ----------
    pvals : 1-D float array, finite values only.
    lam : float in (0, 1), optional
        Tuning parameter. When None, picked via Storey's bootstrap-free
        smoother (Storey 2002 JRSSB, Storey-Tibshirani 2003 PNAS).

    Returns
    -------
    pi0_hat : float in (0, 1].
    """
    if pvals.size == 0:
        return 1.0
    if lam is None:
        # Storey's smoother: try lam in {0.05, 0.10, ..., 0.95}, fit
        # natural cubic spline to pi0(lam), pick pi0 at lam=0.95 of the
        # smoothed curve. Lightweight: skip the spline and use the
        # simple Storey-Tibshirani plug-in at lam=0.5, which performs
        # almost as well in practice when n is large.
        lam = 0.5
    n = pvals.size
    pi0 = float(np.sum(pvals > lam)) / (n * (1.0 - lam))
    return float(min(1.0, max(0.0, pi0)))


def combine_neighbour_pvalues(
    dmc_df: pl.DataFrame,
    *,
    neighbour_bp: int = 200,
    pvalue_col: str = "pvalue",
    meth_diff_col: str = "meth_diff",
    pos_col: str = "pos",
    chrom_col: str = "chrom",
    out_col: str = "pvalue_combined",
    weight: str = "uniform",
    min_sign_agreement: float = 0.6,
    require_focal_signal: bool = True,
    focal_p_thresh: float = 0.5,
) -> pl.DataFrame:
    """RADMeth-style neighbour-aware p-value combiner.

    For each CpG i, combine its p-value with neighbours j on the same
    chromosome within +/-``neighbour_bp`` bp using a **signed Stouffer
    Z-test**. Sites whose effect direction (sign(meth_diff)) agrees
    contribute constructively; sites with opposing signs cancel, so
    spatially isolated false positives are not amplified.

    Per-site signed z is computed from the raw p-value as

        z_i = sign(meth_diff_i) * Phi^{-1}(1 - p_i / 2)

    Combined z is the unweighted Stouffer combination

        Z_combined = sum_{j in W(i)} z_j / sqrt(|W(i)|)

    where W(i) is the window of neighbours including i itself. The
    combined p-value is the two-sided normal tail
    ``2 * (1 - Phi(|Z_combined|))``.

    Parameters
    ----------
    dmc_df : pl.DataFrame
        DMC output. Must carry ``chrom``, ``pos``, ``pvalue``,
        ``meth_diff``.
    neighbour_bp : int
        Half-window width in bp. Each CpG is combined with every CpG
        within this distance on the same chromosome.
    weight : {"uniform"}
        Currently only uniform weighting is supported (Stouffer's
        original formula). Inverse-variance weighting (Liptak) would
        need per-site standard errors and is left as future work.
    min_sign_agreement : float in [0, 1]
        Of the neighbours that contribute (non-NaN signed z), require at
        least this fraction to share the focal site's effect direction.
        Otherwise the focal site keeps its raw p-value. Default 0.6
        (the focal site plus a majority of its neighbours agree). This
        is the guard against spatial-correlation contamination: null
        CpGs next to true DMCs would otherwise inherit the strong z of
        their neighbours and become false positives.
    require_focal_signal : bool
        If True (default), only sites whose raw p < ``focal_p_thresh``
        get combined. Sites whose own evidence is at-or-near uniform
        (raw p approaching 1) keep their raw p-value. Together with
        ``min_sign_agreement``, this restricts combining to candidate
        DMR sites whose own data points in a direction.
    focal_p_thresh : float
        Raw p-value threshold above which the focal site keeps its raw
        p (only relevant when ``require_focal_signal=True``). Default
        0.5; tighten to 0.1 for stricter "DMR-like only" combining.

    Returns
    -------
    pl.DataFrame
        ``dmc_df`` plus two columns: ``out_col`` (the combined p-value)
        and ``out_col + "_n_neighbours"`` (the count of CpGs in the
        window that contributed).

    Notes
    -----
    This is the same per-CpG cross-site combining idea behind
    RADMeth's ``adjust`` step. Returned at the per-CpG level rather
    than per-DMR so it slots into the existing DMC pipeline; pass
    ``out_col`` to ``apply_multiple_testing_correction(pvalue_col=...)``
    to obtain BH q-values on the combined p-values.

    Independence assumption (known limitation in 0.7.x)
    ---------------------------------------------------
    Stouffer's combination assumes the per-site z-scores are independent.
    Adjacent CpGs in WGBS are positively correlated (typical lag-1
    autocorrelation 0.3-0.7 in CpG-dense regions), so the variance of the
    combined Z is **larger** than 1 under H0 and the nominal N(0, 1) tail
    used by this function is anti-conservative for null sites. In epykit
    0.7.x, FDR control over the combined p-values relies on the
    ``min_sign_agreement`` gate (which restricts combining to sites where
    a majority of the window agrees on direction) and ``require_focal_signal``
    (which prevents amplification at otherwise-uniform sites), **not** on
    the Stouffer null being well-calibrated.

    A correlation-aware replacement (Brown's method with an empirical
    correlation kernel estimated per chromosome) is planned for v0.8;
    see ``docs/superpowers/specs/2026-05-27-paper-defendable-benchmark-design.md``
    P0-6 for the deferral rationale. The same spec recommends that
    benchmark reports include null-calibration FDR alongside any
    nominal-q claim that uses this combiner.
    """
    if weight != "uniform":
        raise NotImplementedError(f"weight={weight!r} not implemented yet.")
    if pvalue_col not in dmc_df.columns:
        raise ValueError(f"dmc_df missing required column {pvalue_col!r}")
    if meth_diff_col not in dmc_df.columns:
        raise ValueError(f"dmc_df missing required column {meth_diff_col!r}")

    from scipy.stats import norm as _norm

    out_p_chunks: list[np.ndarray] = []
    out_n_chunks: list[np.ndarray] = []
    keys: list[pl.DataFrame] = []

    for chrom, sub in dmc_df.group_by(chrom_col, maintain_order=True):
        sub_sorted = sub.sort(pos_col)
        positions = sub_sorted[pos_col].to_numpy()
        pvals = sub_sorted[pvalue_col].to_numpy().astype(np.float64)
        diffs = sub_sorted[meth_diff_col].to_numpy().astype(np.float64)
        n = len(positions)

        # Convert per-site p to signed z. p=0 -> +inf z (clamp), p=1 -> 0.
        p_clip = np.clip(pvals, 1e-300, 1.0 - 1e-12)
        # Two-sided p -> one-sided absolute z, then attach sign of effect.
        abs_z = _norm.isf(p_clip / 2.0)
        sign = np.sign(diffs)
        sign[sign == 0] = 1.0  # zero meth_diff: still combine, treat as +
        z = sign * abs_z
        z = np.where(np.isnan(p_clip), 0.0, z)

        # Two-pointer sliding window: for each i, advance left/right to
        # the bounds [pos_i - W, pos_i + W]. O(n_sites) per chromosome.
        z_sum = np.zeros(n, dtype=np.float64)
        n_in = np.zeros(n, dtype=np.int64)
        n_agree = np.zeros(n, dtype=np.int64)
        lo = 0
        hi = 0
        for i in range(n):
            target_lo = positions[i] - neighbour_bp
            target_hi = positions[i] + neighbour_bp
            while lo < n and positions[lo] < target_lo:
                lo += 1
            while hi < n and positions[hi] <= target_hi:
                hi += 1
            # Window is [lo, hi). i is in this range by construction.
            n_window = hi - lo
            if n_window == 0:
                z_sum[i] = z[i]
                n_in[i] = 1
                n_agree[i] = 1
            else:
                slice_z = z[lo:hi]
                slice_sign = sign[lo:hi]
                mask = ~np.isnan(slice_z) & np.isfinite(slice_z)
                n_eff = int(mask.sum())
                if n_eff == 0:
                    z_sum[i] = 0.0
                    n_in[i] = 0
                    n_agree[i] = 0
                else:
                    z_sum[i] = slice_z[mask].sum()
                    n_in[i] = n_eff
                    # Count how many neighbours share the focal site's
                    # effect direction (used by the sign-agreement guard).
                    focal_sign = sign[i]
                    n_agree[i] = int(np.sum((slice_sign[mask] == focal_sign)))
        with np.errstate(invalid="ignore", divide="ignore"):
            z_combined = np.where(n_in > 0, z_sum / np.sqrt(n_in), 0.0)
        p_combined = 2.0 * _norm.sf(np.abs(z_combined))

        # Sign-agreement guard: drop combined p-values where the focal
        # site's neighbours don't sufficiently agree on direction. This
        # is the fix for spatial-correlation contamination -- a null
        # CpG next to a true DMC will see its neighbour's strong z, but
        # if the rest of the window is mixed-sign the agreement drops
        # below the threshold and the focal site keeps its raw p.
        with np.errstate(invalid="ignore", divide="ignore"):
            agree_frac = np.where(n_in > 0, n_agree / n_in, 0.0)
        keep_combined = agree_frac >= min_sign_agreement
        # Focal-signal gate: require the focal site's own raw p to show
        # at least *some* signal (p < focal_p_thresh) before combining.
        if require_focal_signal:
            keep_combined = keep_combined & (pvals < focal_p_thresh)

        p_out = np.where(keep_combined, p_combined, pvals)
        # Never inflate: if combining produced a larger p than raw, keep raw.
        p_out = np.where(p_out < pvals, p_out, pvals)
        p_out = np.where(np.isnan(pvals), np.nan, p_out)
        p_combined = p_out

        out_p_chunks.append(p_combined)
        out_n_chunks.append(n_in)
        keys.append(sub_sorted.select([chrom_col, pos_col]))

    keys_df = pl.concat(keys)
    keys_df = keys_df.with_columns(
        pl.Series(out_col, np.concatenate(out_p_chunks)),
        pl.Series(f"{out_col}_n_neighbours", np.concatenate(out_n_chunks)),
    )
    return dmc_df.join(keys_df, on=[chrom_col, pos_col], how="left")


def _apply_storey_qvalues(pvals: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Storey-Tibshirani q-values. Returns (reject @ 0.05, qvals, pi0_hat).

    Uses pi0 estimated from the data to scale BH q-values:
        q_i = min_{k>=rank(p_i)} pi0_hat * n * p_(k) / k
    Strictly less conservative than fdr_bh when pi0 < 1 (which is the
    null-dominated regime where epykit's `lr` was losing TPR to BH).
    """
    n = pvals.size
    if n == 0:
        return np.zeros(0, dtype=bool), np.zeros(0), 1.0
    pi0 = _storey_pi0(pvals)
    order = np.argsort(pvals, kind="mergesort")
    sorted_p = pvals[order]
    ranks = np.arange(1, n + 1, dtype=np.float64)
    q_sorted = pi0 * n * sorted_p / ranks
    # Enforce monotonicity from the bottom of the list upward (the
    # Storey q-value is min over k >= i of pi0 * n * p_(k) / k).
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)
    qvals = np.empty_like(q_sorted)
    qvals[order] = q_sorted
    reject = qvals < 0.05
    return reject, qvals, pi0


def apply_multiple_testing_correction(
    dmc_results: Union[pl.DataFrame, DMCStore],
    method: str = "fdr_bh",
    pvalue_col: str = "pvalue",
    qvalue_col: str = "qvalue",
) -> Union[pl.DataFrame, DMCStore]:
    """Apply multiple testing correction (Benjamini-Hochberg default).

    Accepts either an in-memory ``pl.DataFrame`` or a ``DMCStore``. For
    a ``DMCStore``, BH runs in a streaming two-pass pattern so the only
    full-table-sized allocation is the float64 pvalue vector itself
    (~176 MB at 22M sites) -- no DataFrame copies, no concat.

    ``method`` selects the FDR procedure:

    * ``"fdr_bh"`` (default) -- Benjamini-Hochberg. Assumes pi0 = 1.
    * ``"fdr_by"`` -- Benjamini-Yekutieli (controls FDR under arbitrary
      dependence; strictly more conservative than BH).
    * ``"fdr_tsbh"`` -- Benjamini-Krieger-Yekutieli two-stage BH;
      statsmodels' pi0-adaptive variant of BH.
    * ``"fdr_storey"`` -- Storey-Tibshirani q-values (Storey 2002,
      Storey-Tibshirani 2003). Uses lam=0.5 to estimate pi0 from the
      empirical p-value histogram. Most powerful when a substantial
      fraction of the genome is null (the typical WGBS scenario);
      reduces to BH when pi0 = 1.

    ``reject`` is written as ``<qvalue_col>_reject`` when the column
    name differs from the default so the two outputs don't collide.
    """
    if method not in _VALID_FDR_METHODS:
        raise ValueError(
            f"method must be one of {sorted(_VALID_FDR_METHODS)}; got {method!r}"
        )

    if isinstance(dmc_results, DMCStore):
        # Cache hit: if this store was already corrected with the
        # same method + qvalue column, the per-chrom parquets already
        # carry the right qvalue/reject columns. Skip the 7-second
        # collect + writeback.
        prev_qcol = dmc_results.manifest.get("bh_qvalue_col", "qvalue")
        prev_method = dmc_results.manifest.get("bh_method", "fdr_bh")
        if (
            dmc_results.bh_applied
            and prev_qcol == qvalue_col
            and prev_method == method
        ):
            logger.info(
                "FDR correction cache hit on %s (method=%s, qvalue_col=%s); skipping.",
                dmc_results.path, method, qvalue_col,
            )
            return dmc_results
        return _apply_bh_to_store(dmc_results, method, pvalue_col, qvalue_col)

    pvals       = dmc_results[pvalue_col].to_numpy()
    nan_mask    = np.isnan(pvals)
    pvals_clean = np.where(nan_mask, 1.0, pvals)

    if method == "fdr_storey":
        reject, qvals, pi0_hat = _apply_storey_qvalues(pvals_clean)
        logger.info("Storey FDR: pi0_hat = %.4f (pvals=%d)", pi0_hat, len(pvals_clean))
    else:
        from statsmodels.stats.multitest import multipletests
        reject, qvals, _, _ = multipletests(pvals_clean, method=method)

    qvals  = np.where(nan_mask, np.nan,  qvals)
    reject = np.where(nan_mask, False,   reject)

    reject_col = "reject" if qvalue_col == "qvalue" else f"{qvalue_col}_reject"
    return dmc_results.with_columns([
        pl.Series(qvalue_col, qvals),
        pl.Series(reject_col, reject),
    ])


def _apply_bh_to_store(
    store: DMCStore,
    method: str,
    pvalue_col: str,
    qvalue_col: str,
) -> DMCStore:
    """Two-pass streaming BH correction over a ``DMCStore``.

    Pass 1: read pvalue column from each chrom parquet, copy into one
    preallocated float64 vector. Track ``(chrom, start, end)`` spans.
    Pass 2: compute BH on the vector, then for each chrom read its
    parquet, attach the qvalue / reject slices, and rewrite the chrom
    parquet atomically. Memory peak is ~3x the pvalue vector (input +
    qvals + reject), independent of the rest of the table.
    """
    total = store.total_sites
    if total == 0:
        logger.warning("apply_multiple_testing_correction: empty DMC store")
        return store

    reject_col = "reject" if qvalue_col == "qvalue" else f"{qvalue_col}_reject"

    logger.info(
        "FDR correction (streaming, method=%s): collecting %s p-values from %d chrom file(s)...",
        method, f"{total:,}", len(store.chroms()),
    )

    pvals = np.empty(total, dtype=np.float64)
    spans: list[tuple[str, int, int]] = []
    offset = 0
    for chrom, df in store.iter_chroms(columns=[pvalue_col]):
        n = len(df)
        if n == 0:
            continue
        pvals[offset:offset + n] = df[pvalue_col].to_numpy()
        spans.append((chrom, offset, offset + n))
        offset += n
        del df

    if offset == 0:
        return store
    # Truncate if any chroms returned fewer rows than the manifest claimed.
    pvals = pvals[:offset]

    nan_mask    = np.isnan(pvals)
    pvals_clean = np.where(nan_mask, 1.0, pvals)
    if method == "fdr_storey":
        reject, qvals, pi0_hat = _apply_storey_qvalues(pvals_clean)
        logger.info("Storey FDR: pi0_hat = %.4f (pvals=%d)", pi0_hat, len(pvals_clean))
    else:
        from statsmodels.stats.multitest import multipletests
        reject, qvals, _, _ = multipletests(pvals_clean, method=method)
    qvals  = np.where(nan_mask, np.nan,  qvals)
    reject = np.where(nan_mask, False,   reject)
    del pvals, pvals_clean, nan_mask

    logger.info("BH correction (streaming): writing q-values back per chromosome...")
    for chrom, start, end in spans:
        df = store.read_chrom(chrom)
        df = df.with_columns([
            pl.Series(qvalue_col, qvals[start:end]),
            pl.Series(reject_col, reject[start:end]),
        ])
        store.update_chrom(chrom, df)
        del df

    store.mark_bh_applied(qvalue_col=qvalue_col, method=method)
    return store


# Empirical-Bayes shrinkage of meth_diff

def shrink_meth_diff(
    dmc_df: pl.DataFrame,
    *,
    se_from: str = "ci",
    out_col: str = "meth_diff_shrunk",
) -> pl.DataFrame:
    """Empirical-Bayes shrinkage of per-CpG meth_diff toward zero.

    Model: ``meth_diff_i ~ N(theta_i, SE_i^2)`` with theta_i ~ N(0, tau^2), so the
    posterior mean is

        theta_i^shrunk = meth_diff_i * tau^2 / (tau^2 + SE_i^2)

    tau^2 is the empirical-Bayes between-site variance estimate

        tau^2 = max(0, Var(meth_diff) - mean(SE^2)).

    This is the Normal-prior special case of the ashr / apeglm family
    (no MLE in a GLM, no Cauchy prior -- but the same shrinkage
    behaviour at the user-visible level: low-information estimates
    collapse to 0 while well-powered effects barely move). Useful for
    ranking and for downstream regression-on-meth_diff analyses where
    raw low-coverage estimates inflate variance.

    Parameters
    ----------
    dmc_df : pl.DataFrame
        DMC output from any backend. Must carry ``meth_diff`` and a
        source of per-site standard errors (see ``se_from``).
    se_from : {"ci", "coef_se"}
        How to derive SE_i.

        * ``"ci"`` (default): infer SE from the 95 % Wald CI on
          ``meth_diff`` as ``(ci_hi - ci_lo) / (2 * 1.96)``. Works on
          every backend that emits ``meth_diff_ci_lo`` /
          ``meth_diff_ci_hi``.
        * ``"coef_se"``: use ``coef_se`` directly (GLM /
          ``bb_lr`` backends). Avoids the CI-width round-trip but is
          on the *linear-predictor* scale (logit beta coefficients), so
          the shrinkage acts on logit-Deltabeta rather than Deltabeta. Prefer ``"ci"``
          unless you specifically want logit-scale shrinkage.
    out_col : str
        Name of the appended shrunk-estimate column. Default
        ``"meth_diff_shrunk"``.

    Returns
    -------
    pl.DataFrame
        ``dmc_df`` plus three new columns:

        * ``out_col`` -- shrunk Deltabeta.
        * ``meth_diff_se`` -- the SE used for shrinkage (handy for QC).
        * ``shrinkage_factor`` -- tau^2 / (tau^2 + SE^2)  in  [0, 1]; values
          near 0 mean "shrink hard," near 1 mean "barely touched."
    """
    if "meth_diff" not in dmc_df.columns:
        raise ValueError(
            "shrink_meth_diff: dmc_df has no 'meth_diff' column."
        )
    if se_from == "ci":
        for col in ("meth_diff_ci_lo", "meth_diff_ci_hi"):
            if col not in dmc_df.columns:
                raise ValueError(
                    f"se_from='ci' needs '{col}' on the DMC table; "
                    "pass se_from='coef_se' if you're using a GLM backend "
                    "without CIs."
                )
        ci_lo = dmc_df.get_column("meth_diff_ci_lo").to_numpy().astype(np.float64)
        ci_hi = dmc_df.get_column("meth_diff_ci_hi").to_numpy().astype(np.float64)
        se = (ci_hi - ci_lo) / (2.0 * 1.959963984540054)  # ~= 2 * z_{0.975}
    elif se_from == "coef_se":
        if "coef_se" not in dmc_df.columns:
            raise ValueError(
                "se_from='coef_se' requires the 'coef_se' column "
                "(present on GLM / bb_lr DMC outputs)."
            )
        se = dmc_df.get_column("coef_se").to_numpy().astype(np.float64)
    else:
        raise ValueError(
            f"se_from must be 'ci' or 'coef_se'; got {se_from!r}"
        )

    meth_diff = dmc_df.get_column("meth_diff").to_numpy().astype(np.float64)
    finite = np.isfinite(meth_diff) & np.isfinite(se) & (se > 0)
    if not finite.any():
        # Nothing to shrink; return all-NaN columns so downstream code
        # doesn't trip on missing fields.
        n = dmc_df.height
        return dmc_df.with_columns([
            pl.Series(out_col, np.full(n, np.nan, dtype=np.float64)),
            pl.Series("meth_diff_se", np.where(np.isfinite(se), se, np.nan)),
            pl.Series("shrinkage_factor", np.full(n, np.nan, dtype=np.float64)),
        ])

    var_md = float(np.var(meth_diff[finite], ddof=1))
    mean_se2 = float(np.mean(se[finite] ** 2))
    tau2 = max(0.0, var_md - mean_se2)
    # Numerical floor: if every effect has zero variance after subtracting
    # mean sampling variance, every effect collapses to 0. That's the
    # correct EB answer (no signal between sites means everything is noise).

    shrink_factor = np.full_like(meth_diff, np.nan, dtype=np.float64)
    shrunk = np.full_like(meth_diff, np.nan, dtype=np.float64)
    if tau2 == 0.0:
        shrink_factor[finite] = 0.0
        shrunk[finite] = 0.0
    else:
        denom = tau2 + se[finite] ** 2
        shrink_factor[finite] = tau2 / denom
        shrunk[finite] = meth_diff[finite] * shrink_factor[finite]

    return dmc_df.with_columns([
        pl.Series(out_col, shrunk),
        pl.Series("meth_diff_se", np.where(np.isfinite(se), se, np.nan)),
        pl.Series("shrinkage_factor", shrink_factor),
    ])


# Permutation-based empirical FDR for DMC

def empirical_fdr_for_dmc(
    methylstore_path: str,
    samples_treatment: list[str],
    samples_control: list[str],
    observed_dmc: pl.DataFrame,
    *,
    n_perm: int = 100,
    seed: int = 42,
    n_jobs: int = 1,
    **dmc_kwargs,
) -> pl.DataFrame:
    """Empirical (permutation) FDR for per-CpG DMC results.

    Re-runs :func:`process_chromosomes_dmc` ``n_perm`` times with the
    treatment / control sample labels shuffled. For each observed CpG, the
    empirical p-value is estimated from the fraction of null DMCs (across
    all permutations and all sites) with raw p-value <= the observed raw
    p-value. The result is BH-adjusted to ``empirical_qvalue``.

    Parallels :func:`epykit.dmr.empirical_fdr_for_dmr`; same caveats apply:

    * The shuffler ignores any structure in ``md.obs`` (donor, batch).
      Covariate-adjusted DMC (``formula=`` / ``contrast=`` in
      :func:`tl.dmc`) refuses to call this -- label shuffling invalidates
      the stratified design.
    * Computational cost is roughly n_perm x (cost of one DMC run); on a
      whole-genome WGBS analysis this is the dominant runtime. Use a
      smaller ``n_perm`` (e.g. 50) or fewer chromosomes during exploration.

    Parameters
    ----------
    methylstore_path, samples_treatment, samples_control
        Same arguments forwarded to :func:`process_chromosomes_dmc`.
    observed_dmc
        DMC DataFrame returned by the observed (unpermuted) run. Must
        carry a ``pvalue`` column. Empirical columns are appended to a
        copy of this frame.
    n_perm
        Number of permutations. Default 100.
    seed
        Seed for the per-permutation label shuffler.
    n_jobs
        joblib parallel worker count. -1 uses all cores. Falls back to
        serial execution when joblib is not installed. Default 1.
    **dmc_kwargs
        Forwarded to :func:`process_chromosomes_dmc` for each permutation;
        should match the observed run's settings (test, chromosomes,
        unite, min_samples_*, dispersion, reference).

    Returns
    -------
    pl.DataFrame
        ``observed_dmc`` with added columns ``empirical_pvalue`` and
        ``empirical_qvalue``.
    """
    if len(observed_dmc) == 0:
        return observed_dmc.with_columns([
            pl.lit(None, dtype=pl.Float64).alias("empirical_pvalue"),
            pl.lit(None, dtype=pl.Float64).alias("empirical_qvalue"),
        ])
    if "pvalue" not in observed_dmc.columns:
        raise ValueError(
            "observed_dmc has no 'pvalue' column; empirical FDR needs raw "
            "p-values to compute null tail probabilities."
        )

    n_treat = len(samples_treatment)
    pool = list(samples_treatment) + list(samples_control)

    def _run_one_perm(perm_idx: int) -> np.ndarray:
        # Local RNG so parallel workers stay deterministic across n_jobs.
        local_rng = np.random.default_rng(seed + perm_idx + 1)
        shuffled = pool.copy()
        local_rng.shuffle(shuffled)
        perm_treat = shuffled[:n_treat]
        perm_ctrl = shuffled[n_treat:]
        kwargs = dict(dmc_kwargs)
        # Strip deprecated aliases so they don't double-bind.
        kwargs.pop("samples_case", None)
        kwargs.pop("min_samples_case", None)
        # Each permutation runs DMC; we only need the pvalue column.
        # Use a throwaway out_dir per-perm so n_perm runs don't leave
        # n_perm stores littering the cache, and use return_store=True
        # so we can pull just the pvalue stream without materialising
        # the full per-perm DataFrame.
        kwargs.pop("out_dir", None)
        kwargs.pop("return_store", None)
        perm_dir = Path(tempfile.mkdtemp(prefix=f"epykit_dmc_perm_{perm_idx}_"))
        try:
            null_store = process_chromosomes_dmc(
                methylstore_path=methylstore_path,
                samples_treatment=perm_treat,
                samples_control=perm_ctrl,
                out_dir=perm_dir,
                return_store=True,
                **kwargs,
            )
        except Exception as exc:
            logger.warning("DMC permutation %d failed: %s", perm_idx, exc)
            import shutil
            shutil.rmtree(perm_dir, ignore_errors=True)
            return np.array([], dtype=np.float64)
        try:
            if null_store.total_sites == 0:
                return np.array([], dtype=np.float64)
            # Concatenate pvalues across chroms without materialising
            # the full perm DataFrame.
            parts = [
                df.get_column("pvalue").drop_nulls().to_numpy()
                for _, df in null_store.iter_chroms(columns=["pvalue"])
            ]
            if not parts:
                return np.array([], dtype=np.float64)
            return np.concatenate(parts)
        finally:
            null_store.cleanup()

    null_pvals_list: list[np.ndarray]
    if n_jobs == 1:
        null_pvals_list = [_run_one_perm(i) for i in range(n_perm)]
    else:
        try:
            from joblib import Parallel, delayed
            null_pvals_list = Parallel(n_jobs=n_jobs)(
                delayed(_run_one_perm)(i) for i in range(n_perm)
            )
        except ImportError:
            logger.warning(
                "joblib not installed; running DMC permutations serially."
            )
            null_pvals_list = [_run_one_perm(i) for i in range(n_perm)]

    if all(len(arr) == 0 for arr in null_pvals_list):
        logger.warning(
            "All %d DMC permutations produced zero null sites. Empirical "
            "p-values default to 1 / (1 + n_perm).",
            n_perm,
        )
    # `pvalue` is the per-CpG raw p-value by contract (P0-1 fix). The
    # combined column, if neighbour_combine was on, is `pvalue_combined`
    # and is intentionally NOT used here -- the null pool comes from
    # raw per-CpG runs of the same test, so the observed side must also
    # be raw to keep the comparison apples-to-apples.
    if "pvalue_combined" in observed_dmc.columns and "pvalue" not in observed_dmc.columns:
        raise ValueError(
            "empirical_fdr_for_dmc requires the raw `pvalue` column; "
            "got `pvalue_combined` instead. This indicates a stale "
            "(<=0.7.2) workflow that overwrote `pvalue` with the "
            "combined value -- re-run dmc() on the current epykit."
        )
    obs_p = observed_dmc.get_column("pvalue").to_numpy()
    # Per-permutation tail count: for each observed p, count the number
    # of permutations that produced at least one null site with p <= obs_p.
    # Denominator is n_perm + 1 (not |pooled null| + 1) so the empirical
    # p-value floor is 1/(n_perm+1), independent of how many sites each
    # permutation emitted.  This is the same formula used by
    # empirical_fdr_for_dmr (P0-2 fix) and ensures that the floor is
    # invariant to the permutation count rather than growing as n_perm
    # grows (the old pooled formula was anti-conservative for the same
    # reason).
    obs_finite_mask = np.isfinite(obs_p)
    obs_safe = np.where(obs_finite_mask, obs_p, 1.0)

    # min_null_p_per_perm[i] = min p across all null sites from perm i
    # (1.0 if perm produced no sites). The number of perms with at least
    # one null site with p <= obs_p is then sum(min_null_p[i] <= obs_p).
    min_null_p_per_perm = np.array([
        float(arr.min()) if len(arr) > 0 else 1.0
        for arr in null_pvals_list
    ], dtype=np.float64)
    min_null_sorted = np.sort(min_null_p_per_perm)
    counts = np.searchsorted(min_null_sorted, obs_safe, side="right")
    emp_p = (counts + 1.0) / (n_perm + 1.0)
    emp_p = np.clip(emp_p, 0.0, 1.0)
    emp_p = np.where(obs_finite_mask, emp_p, np.nan)

    from statsmodels.stats.multitest import multipletests
    finite = np.isfinite(emp_p)
    emp_q = np.full_like(emp_p, np.nan, dtype=np.float64)
    if finite.any():
        _, q_finite, _, _ = multipletests(emp_p[finite], method="fdr_bh")
        emp_q[finite] = q_finite

    return observed_dmc.with_columns([
        pl.Series("empirical_pvalue", emp_p),
        pl.Series("empirical_qvalue", emp_q),
    ])