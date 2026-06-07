# tests/test_empirical_fdr_stratified.py
"""M1: stratified empirical FDR must permute within strata while preserving
each stratum's original treatment/control split. The pre-fix implementation
shuffled in dict order then split globally, so a paired 1T+1C design sent
ALL pairs' first elements to treatment every permutation."""
import numpy as np

from epykit.dmr import _stratified_permutation_assignment


def test_stratified_permutation_preserves_per_stratum_counts():
    # 4 paired strata, each 1 treatment + 1 control. n_treat = 4, n_ctrl = 4.
    strata_map = {
        "pair_A": ["A_T", "A_C"],
        "pair_B": ["B_T", "B_C"],
        "pair_C": ["C_T", "C_C"],
        "pair_D": ["D_T", "D_C"],
    }
    original_treatment = ["A_T", "B_T", "C_T", "D_T"]
    original_control = ["A_C", "B_C", "C_C", "D_C"]
    rng = np.random.default_rng(0)

    # Run many permutations; each must produce exactly one treatment sample
    # per stratum (k_treat = 1 per pair).
    for _ in range(200):
        perm_t, perm_c = _stratified_permutation_assignment(
            strata_map=strata_map,
            samples_treatment=original_treatment,
            samples_control=original_control,
            rng=rng,
        )
        assert len(perm_t) == 4
        assert len(perm_c) == 4
        for stratum, members in strata_map.items():
            n_in_treat = sum(1 for s in perm_t if s in members)
            assert n_in_treat == 1, (
                f"Stratum {stratum} must contribute exactly 1 sample to "
                f"treatment (its original count); got {n_in_treat}."
            )


def test_stratified_permutation_unequal_strata_preserves_counts():
    # Stratum sizes 2T+1C and 1T+2C.
    strata_map = {"S1": ["s1a", "s1b", "s1c"], "S2": ["s2a", "s2b", "s2c"]}
    original_treatment = ["s1a", "s1b", "s2a"]   # 2 from S1, 1 from S2
    original_control = ["s1c", "s2b", "s2c"]
    rng = np.random.default_rng(1)

    for _ in range(200):
        perm_t, _ = _stratified_permutation_assignment(
            strata_map=strata_map,
            samples_treatment=original_treatment,
            samples_control=original_control,
            rng=rng,
        )
        s1_in_t = sum(1 for s in perm_t if s in strata_map["S1"])
        s2_in_t = sum(1 for s in perm_t if s in strata_map["S2"])
        assert s1_in_t == 2 and s2_in_t == 1
