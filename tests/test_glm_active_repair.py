"""Independent scalar IRLS stopping and heterogeneous active-row workloads."""

import numpy as np

from epykit import _glm


def scalar_coefficients(m, n, x, max_iter=25, tol=1e-6):
    good = n > 0
    prop = np.where(good, m/np.maximum(n, 1), .5)
    prop = np.clip(prop, 1e-6, 1-1e-6)
    z = np.log(prop/(1-prop))
    beta = np.linalg.lstsq(x[good]*np.sqrt(n[good, None]), z[good]*np.sqrt(n[good]), rcond=None)[0]
    for _ in range(max_iter):
        eta = np.clip(x @ beta, -30, 30)
        mu = np.clip(1/(1+np.exp(-eta)), 1e-6, 1-1e-6)
        var = mu*(1-mu)
        w = n*var
        z = eta + np.where(good, (m/np.maximum(n, 1)-mu)/var, 0)
        new = np.linalg.lstsq(x[good]*np.sqrt(w[good, None]), z[good]*np.sqrt(w[good]), rcond=None)[0]
        if np.max(np.abs(new-beta)) < tol:
            return new
        beta = new
    return np.full(x.shape[1], np.nan)


def test_easy_converged_rows_leave_active_batch(monkeypatch):
    x = np.column_stack([np.ones(6), [0., 0., 0., 1., 1., 1.]])
    n = np.full((4, 6), 100)
    m = np.array([[50]*6, [20, 30, 60, 40, 70, 80], [0, 0, 1, 97, 98, 100], [0]*6])
    rows = []
    original = _glm._solve_weighted_lsq

    def counted(x, z, w):
        rows.append(len(z))
        return original(x, z, w)

    monkeypatch.setattr(_glm, '_solve_weighted_lsq', counted)
    result = _glm.irls_binomial_batch(m, n, x, return_cov=True)
    assert any(0 < count < 4 for count in rows)
    for i in range(3):
        np.testing.assert_allclose(result[0][i], scalar_coefficients(m[i], n[i], x), atol=1e-8, rtol=1e-6)
    assert np.isnan(result[0][-1]).all()


def test_empty_and_all_invalid_batches_keep_shapes_and_nan_masks():
    x = np.column_stack([np.ones(6), [0., 0., 0., 1., 1., 1.]])
    for count in [0, 1, 7]:
        n = np.zeros((count, 6), int)
        result = _glm.irls_binomial_batch(n, n, x, return_cov=True)
        assert [a.shape for a in result] == [(count, 2), (count, 2), (count,), (count,), (count,), (count, 2, 2)]
        for a in [*result[:4], result[5]]:
            assert np.isnan(a).all()
        assert np.all(result[4] == 0)


def test_nonconverged_and_nonfinite_updates_do_not_emit_inference(monkeypatch):
    x = np.column_stack([np.ones(6), [0., 0., 0., 1., 1., 1.]])
    m = np.array([[1, 4, 30, 60, 90, 99]])
    n = np.full_like(m, 100)
    assert np.isnan(_glm.irls_binomial_batch(m, n, x, max_iter=0)[0]).all()
    original = _glm._solve_weighted_lsq
    calls = 0

    def fail_update(x, z, w):
        nonlocal calls
        calls += 1
        return original(x, z, w) if calls == 1 else np.full((len(z), 2), np.nan)

    monkeypatch.setattr(_glm, '_solve_weighted_lsq', fail_update)
    result = _glm.irls_binomial_batch(m, n, x, return_cov=True)
    assert np.isnan(result[0]).all()
    assert calls == 2
