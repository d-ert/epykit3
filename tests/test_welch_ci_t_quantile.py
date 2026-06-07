"""M13: welch_t's meth_diff CI must use the t(Satterthwaite df) critical
value that its p-value uses, not the hardcoded normal z=1.96. At n=3
(dof~4, t_.025=2.776) the z-based interval is ~30% too narrow vs its own
test, so a site can have p>0.05 yet a CI excluding 0 (CI/test disagreement)."""
import numpy as np
from scipy import stats as sp_stats

from epykit._glm import welch_meth_diff_ci


def test_welch_ci_uses_t_quantile_when_dof_given():
    # Means chosen away from the [-1,1] clamp so the half-width is exact.
    mean_case = np.array([0.30])
    mean_ctrl = np.array([0.20])
    vm_case = np.array([0.0025])   # var of the mean
    vm_ctrl = np.array([0.0025])
    dof = np.array([4.0])

    lo_t, hi_t = welch_meth_diff_ci(mean_case, vm_case, mean_ctrl, vm_ctrl, dof=dof)
    lo_z, hi_z = welch_meth_diff_ci(mean_case, vm_case, mean_ctrl, vm_ctrl)  # dof=None -> z

    half_t = (hi_t - lo_t) / 2.0
    half_z = (hi_z - lo_z) / 2.0
    # t(4) interval is wider than the normal interval by t_.025,4 / z_.025.
    expected_ratio = sp_stats.t.isf(0.025, 4) / sp_stats.norm.isf(0.025)
    np.testing.assert_allclose(half_t / half_z, expected_ratio, rtol=1e-9)
    assert expected_ratio > 1.4  # ~1.416 at dof=4


def test_welch_ci_t_half_width_matches_test_critical_value():
    mean_case = np.array([0.30])
    mean_ctrl = np.array([0.20])
    vm_case = np.array([0.0025])
    vm_ctrl = np.array([0.0025])
    dof = np.array([4.0])
    lo, hi = welch_meth_diff_ci(mean_case, vm_case, mean_ctrl, vm_ctrl, dof=dof)
    se = np.sqrt(vm_case + vm_ctrl)
    crit = sp_stats.t.isf(0.025, 4)
    np.testing.assert_allclose((hi - lo) / 2.0, crit * se, rtol=1e-9)


def test_welch_ci_default_is_legacy_normal():
    mean_case = np.array([0.30])
    mean_ctrl = np.array([0.20])
    vm = np.array([0.0025])
    lo, hi = welch_meth_diff_ci(mean_case, vm, mean_ctrl, vm)
    se = np.sqrt(vm + vm)
    z = sp_stats.norm.isf(0.025)
    np.testing.assert_allclose((hi - lo) / 2.0, z * se, rtol=1e-9)


def test_welch_ci_large_dof_approaches_normal():
    mean_case = np.array([0.30])
    mean_ctrl = np.array([0.20])
    vm = np.array([0.0025])
    lo_t, hi_t = welch_meth_diff_ci(mean_case, vm, mean_ctrl, vm, dof=np.array([10_000.0]))
    lo_z, hi_z = welch_meth_diff_ci(mean_case, vm, mean_ctrl, vm)
    np.testing.assert_allclose(hi_t - lo_t, hi_z - lo_z, rtol=1e-3)
