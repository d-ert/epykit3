"""Observed-data CAR/donor GLS through a small residual donor eigensystem.

Locus intercepts are eliminated with one normal-matrix Cholesky per CAR
length. The residual donor covariance has rank at most S-2, so profiling its
variance needs only O(S) work per evaluation rather than repeated KxK solves.
Missing cells never enter the design, likelihood, or residual degrees of freedom.
"""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True, nogil=True)
def _car_precision(p, length):
    n = len(p)
    diagonal = np.ones(n)
    off = np.zeros(max(n - 1, 0))
    logdet = 0.0
    if length > 0:
        for j in range(n - 1):
            gap = (p[j + 1] - p[j]) / length
            a, denominator = np.exp(-gap), -np.expm1(-2 * gap)
            diagonal[j] += a * a / denominator
            diagonal[j + 1] += a * a / denominator
            off[j] = -a / denominator
            logdet += np.log(denominator)
    return diagonal, off, logdet


@njit(cache=True, nogil=True)
def _apply_car(diagonal, off, values):
    result = diagonal * values
    for j in range(len(off)):
        result[j] += off[j] * values[j + 1]
        result[j + 1] += off[j] * values[j]
    return result


@njit(cache=True, nogil=True)
def _chol_solve(chol, rhs):
    n, columns = rhs.shape
    result = rhs.copy()
    for col in range(columns):
        for i in range(n):
            for j in range(i):
                result[i, col] -= chol[i, j] * result[j, col]
            result[i, col] /= chol[i, i]
        for i in range(n - 1, -1, -1):
            for j in range(i + 1, n):
                result[i, col] -= chol[j, i] * result[j, col]
            result[i, col] /= chol[i, i]
    return result


@njit(cache=True, nogil=True)
def _spectral_prepare(p, y, h, g, mask, length):
    samples, k = y.shape
    shared_locus = False
    for j in range(k):
        control_seen, case_seen = False, False
        for s in range(samples):
            if mask[s, j]:
                if g[s] == 0:
                    control_seen = True
                else:
                    case_seen = True
        if control_seen and case_seen:
            shared_locus = True
    if not shared_locus:
        raise ValueError("Observed regional treatment contrast is not identifiable")
    normal = np.zeros((k + 1, k + 1))
    rhs = np.zeros((k + 1, samples + 2))
    median_sqrt = np.sqrt(np.median(h * h))
    logdet = 0.0
    n_observed = 0
    for s in range(samples):
        indices = np.flatnonzero(mask[s])
        n_observed += len(indices)
        diagonal, off, det = _car_precision(p[indices], length)
        logdet += det
        hs, ys = h[indices], y[s, indices]
        qh, qy = _apply_car(diagonal, off, hs), _apply_car(diagonal, off, ys)
        information, response = np.dot(hs, qh), np.dot(hs, qy)
        for a in range(len(indices)):
            j = indices[a]
            normal[j, j] += hs[a] * hs[a] * diagonal[a]
            normal[j, k] += g[s] * hs[a] * qh[a]
            normal[k, j] = normal[j, k]
            rhs[j, 0] += hs[a] * qy[a]
            rhs[j, s + 1] += hs[a] * qh[a] / median_sqrt
            if a + 1 < len(indices):
                next_locus = indices[a + 1]
                value = hs[a] * hs[a + 1] * off[a]
                normal[j, next_locus] += value
                normal[next_locus, j] += value
        normal[k, k] += g[s] * g[s] * information
        rhs[k, 0] += g[s] * response
        rhs[k, s + 1] += g[s] * information / median_sqrt
    rhs[k, samples + 1] = 1.0
    df = n_observed - k - 1
    if df <= 0:
        raise ValueError("Observed regional design has no residual degrees of freedom")
    chol = np.linalg.cholesky(normal)
    solved = _chol_solve(chol, rhs)
    beta0, donor_coefficients = solved[:, 0], solved[:, 1:samples + 1]
    treatment_variance = solved[k, samples + 1]
    logdet += 2.0 * np.log(np.diag(chol)).sum()

    # Residual quadratics avoid cancellation in yy - b'G^-1b and F-D'G^-1D.
    donor_gram, donor_response = np.zeros((samples, samples)), np.zeros(samples)
    rss0 = 0.0
    for s in range(samples):
        indices = np.flatnonzero(mask[s])
        diagonal, off, _ = _car_precision(p[indices], length)
        residual = np.empty(len(indices))
        donor_residual = np.empty((len(indices), samples))
        for a in range(len(indices)):
            j = indices[a]
            residual[a] = y[s, j] - h[j] * (beta0[j] + g[s] * beta0[k])
            for t in range(samples):
                donor_residual[a, t] = h[j] * (
                    (1.0 if s == t else 0.0) / median_sqrt
                    - donor_coefficients[j, t] - g[s] * donor_coefficients[k, t]
                )
        qr = _apply_car(diagonal, off, residual)
        rss0 += np.dot(residual, qr)
        for t in range(samples):
            qv = _apply_car(diagonal, off, donor_residual[:, t])
            donor_response[t] += np.dot(donor_residual[:, t].copy(), qr)
            for u in range(samples):
                donor_gram[u, t] += np.dot(donor_residual[:, u].copy(), qv)
    eigenvalues, vectors = np.linalg.eigh((donor_gram + donor_gram.T) / 2.0)
    tolerance = 1e-10 * max(1.0, np.max(np.abs(eigenvalues)))
    for t in range(samples):
        if eigenvalues[t] < -tolerance:
            raise ValueError("Residual donor information is not positive semidefinite")
        if eigenvalues[t] < tolerance:
            eigenvalues[t] = 0.0
    z = vectors.T @ donor_response
    weights = donor_coefficients[k, :].copy() @ vectors
    energy = np.zeros(samples)
    for t in range(samples):
        if eigenvalues[t] > 0:
            energy[t] = z[t] * z[t] / eigenvalues[t]
        else:
            z[t] = 0.0
    perpendicular = rss0 - energy.sum()
    if perpendicular < -1e-8 * max(rss0, 1.0):
        raise ValueError("Residual spectral energy exceeds observed residual sum")
    perpendicular = max(perpendicular, 0.0)
    return eigenvalues, z, weights, energy, perpendicular, beta0[k], treatment_variance, logdet, df


@njit(cache=True, nogil=True)
def _spectral_evaluate(prepared, c):
    eigenvalues, z, weights, energy, perpendicular, beta0, base_variance, logdet, df = prepared
    ratio = c / (1.0 - c)
    rss = perpendicular
    effect, variance, determinant = beta0, base_variance, logdet
    for t in range(len(eigenvalues)):
        denominator = 1.0 + ratio * eigenvalues[t]
        rss += energy[t] / denominator
        effect -= weights[t] * ratio * z[t] / denominator
        # Confounded donor modes (lambda=0) still contribute to coefficient
        # variance, even though the restricted likelihood cannot see them.
        variance += weights[t] * weights[t] * ratio / denominator
        determinant += np.log1p(ratio * eigenvalues[t])
    if rss <= 0 or variance <= 0:
        return np.inf, effect, np.nan
    return df * np.log(rss / df) + determinant, effect, np.sqrt(rss / df * variance)


@njit(cache=True, nogil=True)
def _donor_golden(prepared, left, right):
    golden = .6180339887498949
    a, b = right-golden*(right-left), left+golden*(right-left)
    fa, fb = _spectral_evaluate(prepared, a)[0], _spectral_evaluate(prepared, b)[0]
    for _ in range(48):
        if right-left < 1e-8:
            break
        if fa < fb:
            right, b, fb = b, a, fa
            a = right-golden*(right-left)
            fa = _spectral_evaluate(prepared, a)[0]
        else:
            left, a, fa = a, b, fb
            b = left+golden*(right-left)
            fb = _spectral_evaluate(prepared, b)[0]
    return a if fa < fb else b


@njit(cache=True, nogil=True)
def _profile_masked(p, y, h, g, mask, length):
    prepared = _spectral_prepare(p, y, h, g, mask, length)
    fractions = np.linspace(0.0, .98, 9)
    objectives = np.empty(len(fractions))
    best, effect, se = _spectral_evaluate(prepared, 0.0)
    best_c = 0.0
    for j in range(len(fractions)):
        value, ee, ss = _spectral_evaluate(prepared, fractions[j])
        objectives[j] = value
        if value < best:
            best, best_c, effect, se = value, fractions[j], ee, ss
    for j in range(len(fractions)):
        left, right = max(0, j-1), min(len(fractions)-1, j+1)
        if objectives[j] <= objectives[left] and objectives[j] <= objectives[right]:
            c = _donor_golden(prepared, fractions[left], fractions[right])
            value, ee, ss = _spectral_evaluate(prepared, c)
            if value < best:
                best, best_c, effect, se = value, c, ee, ss
    return best, best_c, effect, se


@njit(cache=True, nogil=True)
def _length_golden(p, y, h, g, mask, left, right):
    golden = .6180339887498949
    a, b = right-golden*(right-left), left+golden*(right-left)
    fa = _profile_masked(p, y, h, g, mask, np.exp(a))[0]
    fb = _profile_masked(p, y, h, g, mask, np.exp(b))[0]
    for _ in range(36):
        if right-left < 1e-6:
            break
        if fa < fb:
            right, b, fb = b, a, fa
            a = right-golden*(right-left)
            fa = _profile_masked(p, y, h, g, mask, np.exp(a))[0]
        else:
            left, a, fa = a, b, fb
            b = left+golden*(right-left)
            fb = _profile_masked(p, y, h, g, mask, np.exp(b))[0]
    return a if fa < fb else b


@njit(cache=True, nogil=True)
def _fit_masked(p, y, h, g, mask):
    best, c, effect, se = _profile_masked(p, y, h, g, mask, 0.0)
    best_length = 0.0
    grid = np.linspace(np.log(.1), np.log(10000.0), 17)
    objectives = np.empty(len(grid))
    for j in range(len(grid)):
        length = np.exp(grid[j])
        value, cc, ee, ss = _profile_masked(p, y, h, g, mask, length)
        objectives[j] = value
        if value < best:
            best, best_length, c, effect, se = value, length, cc, ee, ss
    for j in range(1, len(grid)-1):
        if objectives[j] < objectives[j-1] and objectives[j] < objectives[j+1]:
            length = np.exp(_length_golden(p, y, h, g, mask, grid[j-1], grid[j+1]))
            value, cc, ee, ss = _profile_masked(p, y, h, g, mask, length)
            if value < best:
                best, best_length, c, effect, se = value, length, cc, ee, ss
    return best_length, c, effect, se, best
