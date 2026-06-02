"""P0-6: combine_neighbour_pvalues docstring must own the Stouffer
independence caveat. We are not fixing the math in 0.7.3 -- we are owning
the limitation so reviewers see it documented."""
from epykit.dmc import combine_neighbour_pvalues


def test_docstring_acknowledges_independence_violation():
    doc = combine_neighbour_pvalues.__doc__ or ""
    # The key phrases reviewers will look for:
    assert "independence" in doc.lower() or "correlated" in doc.lower(), (
        "Docstring must mention that adjacent CpGs are correlated and "
        "Stouffer assumes independence."
    )
    assert "min_sign_agreement" in doc or "sign agreement" in doc.lower(), (
        "Docstring must explain that the FDR safety net comes from the "
        "sign-agreement gate, not the Stouffer null."
    )
    assert "Brown" in doc or "v0.8" in doc or "future work" in doc.lower(), (
        "Docstring must point readers at the planned correlation-aware "
        "replacement (Brown's method) and note this is a 0.7.x limitation."
    )
