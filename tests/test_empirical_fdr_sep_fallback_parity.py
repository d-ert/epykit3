# tests/test_empirical_fdr_sep_fallback_parity.py
"""M3: empirical DMC FDR must run permutations with the SAME sep_fallback
and smoothing settings as the observed run; otherwise the Westfall-Young
statistic compares deflated observed p-values against an un-deflated null
pool, producing anti-conservative empirical p."""
import inspect

import epykit.dmc as ep_dmc


def test_empirical_fdr_for_dmc_signature_accepts_sep_and_smoothing_kwargs():
    """Pre-fix: empirical_fdr_for_dmc accepted no sep_fallback/sep_threshold/
    smoothing/smoothing_span_bp kwargs, so tl.dmc could not forward them
    even if it wanted to. The fix adds these as keyword-only args."""
    sig = inspect.signature(ep_dmc.empirical_fdr_for_dmc)
    names = set(sig.parameters)
    for required in (
        "sep_fallback",
        "sep_threshold",
        "smoothing",
        "smoothing_span_bp",
    ):
        assert required in names, (
            f"empirical_fdr_for_dmc must accept {required!r} as a kwarg "
            f"(M3 forwarding fix); current params: {sorted(names)}"
        )


def test_empirical_fdr_for_dmc_forwards_sep_and_smoothing(monkeypatch):
    """Patch the inner per-perm DMC engine call and verify each
    permutation receives sep_fallback / smoothing matching the
    empirical_fdr_for_dmc kwargs the user passed."""
    import polars as pl

    captured: list[dict] = []

    class _FakeDMCStore:
        """Minimal stand-in for DMCStore: zero null sites so the empirical
        worker treats each perm as having no contributions."""
        total_sites = 0

        def iter_chroms(self, columns=None):
            return iter(())

        def cleanup(self):
            pass

    def _fake_process_chromosomes_dmc(*args, **kwargs):
        captured.append(dict(kwargs))
        return _FakeDMCStore()

    # Patch the symbol as seen INSIDE epykit.dmc (where empirical_fdr_for_dmc
    # references it). The exact attribute path may be one of:
    #   epykit.dmc.process_chromosomes_dmc      (most likely)
    #   epykit.dmc._process_chromosomes_dmc     (if there's a wrapper)
    # We try both; one must exist.
    target_name = "process_chromosomes_dmc"
    assert hasattr(ep_dmc, target_name), (
        f"epykit.dmc must expose {target_name!r} for patch -- update test "
        f"target if the engine was renamed."
    )
    monkeypatch.setattr(ep_dmc, target_name, _fake_process_chromosomes_dmc)

    # Minimal observed DMC frame -- 1 'significant' site so empirical step runs.
    observed = pl.DataFrame({
        "chrom": ["chr1"],
        "pos": [100],
        "pvalue": [1e-3],
        "meth_diff": [0.4],
    })

    # We pass a methylstore path that the patched engine will ignore.
    ep_dmc.empirical_fdr_for_dmc(
        methylstore_path="/__unused__",
        samples_treatment=["s1", "s2", "s3"],
        samples_control=["c1", "c2", "c3"],
        observed_dmc=observed,
        n_perm=3,
        seed=0,
        n_jobs=1,
        test="lr",
        sep_fallback=True,
        sep_threshold=0.05,
        smoothing=True,
        smoothing_span_bp=500,
    )

    assert len(captured) >= 3, (
        f"Expected >=3 per-perm calls, got {len(captured)}. The patched "
        f"engine may have been bypassed."
    )
    for kwargs in captured:
        assert kwargs.get("sep_fallback") is True, (
            f"sep_fallback not forwarded: {kwargs}"
        )
        assert kwargs.get("sep_threshold") == 0.05, kwargs
        assert kwargs.get("smoothing") is True, kwargs
        assert kwargs.get("smoothing_span_bp") == 500, kwargs


def test_empirical_fdr_for_dmc_defaults_do_not_leak_into_perms(monkeypatch):
    """Negative-path companion to the forwarding test: when the user does
    NOT pass sep_fallback/smoothing kwargs, each per-perm call must see the
    documented defaults (False / False / 0.9 / 500). Guards against a future
    maintainer flipping a default and silently changing behaviour for every
    existing caller of empirical_fdr_for_dmc."""
    import polars as pl

    captured: list[dict] = []

    class _FakeDMCStore:
        total_sites = 0

        def iter_chroms(self, columns=None):
            return iter(())

        def cleanup(self):
            pass

    def _fake_process_chromosomes_dmc(*args, **kwargs):
        captured.append(dict(kwargs))
        return _FakeDMCStore()

    monkeypatch.setattr(
        ep_dmc, "process_chromosomes_dmc", _fake_process_chromosomes_dmc
    )

    observed = pl.DataFrame({
        "chrom": ["chr1"], "pos": [100],
        "pvalue": [1e-3], "meth_diff": [0.4],
    })

    # NB: deliberately omit sep_fallback / sep_threshold / smoothing /
    # smoothing_span_bp -- we want the helper's defaults to flow through.
    ep_dmc.empirical_fdr_for_dmc(
        methylstore_path="/__unused__",
        samples_treatment=["s1", "s2", "s3"],
        samples_control=["c1", "c2", "c3"],
        observed_dmc=observed,
        n_perm=3,
        seed=0,
        n_jobs=1,
        test="lr",
    )

    assert len(captured) >= 3
    for kwargs in captured:
        assert kwargs.get("sep_fallback") is False, kwargs
        assert kwargs.get("sep_threshold") == 0.9, kwargs
        assert kwargs.get("smoothing") is False, kwargs
        assert kwargs.get("smoothing_span_bp") == 500, kwargs
