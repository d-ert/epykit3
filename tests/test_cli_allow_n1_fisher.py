"""D12: ``--allow-n1`` resolves --test to fisher at n=1.

``--allow-n1`` advertises a "Fisher fallback", but the CLI default ``--test lr``
(and explicit ``lr``) has no n=1 path -- so the advertised fallback never fired
and the n=1 run silently used lr. The fix: in ``_cli_n1_and_footgun_checks``,
when ``min(n_treat, n_ctrl) < 2`` and ``--allow-n1`` is set and the user is on
the lr/auto engine (no n=1 fallback), resolve ``args.test = "fisher"``. An
explicit non-lr engine choice (glm/welch_t) is respected, not overridden.
"""

from __future__ import annotations

import argparse
import warnings


def _checks_namespace(*, test, allow_n1, n_treat=1, n_ctrl=1):
    args = argparse.Namespace(
        test=test,
        allow_n1=allow_n1,
        unite=True,
        min_samples_treatment=2,
        min_samples_control=2,
    )
    args._samples = (
        [f"t{i}" for i in range(n_treat)],
        [f"c{i}" for i in range(n_ctrl)],
    )
    return args


def test_allow_n1_resolves_lr_to_fisher():
    """Bare ``--allow-n1`` (default --test lr) at n=1 must resolve to fisher."""
    from epykit.cli import _cli_n1_and_footgun_checks

    args = _checks_namespace(test="lr", allow_n1=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _cli_n1_and_footgun_checks(args, unit="sites")
    assert args.test == "fisher"


def test_allow_n1_explicit_welch_t_not_overridden():
    """An explicit ``--test welch_t`` at n=1 + --allow-n1 is left as the user
    chose it (fisher is NOT silently substituted)."""
    from epykit.cli import _cli_n1_and_footgun_checks

    args = _checks_namespace(test="welch_t", allow_n1=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _cli_n1_and_footgun_checks(args, unit="sites")
    assert args.test == "welch_t"


def test_allow_n1_explicit_glm_not_overridden():
    from epykit.cli import _cli_n1_and_footgun_checks

    args = _checks_namespace(test="glm", allow_n1=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _cli_n1_and_footgun_checks(args, unit="sites")
    assert args.test == "glm"


def test_n2_lr_not_changed():
    """At n>=2, lr is left untouched (no fisher substitution)."""
    from epykit.cli import _cli_n1_and_footgun_checks

    args = _checks_namespace(test="lr", allow_n1=False, n_treat=2, n_ctrl=2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _cli_n1_and_footgun_checks(args, unit="sites")
    assert args.test == "lr"


def test_allow_n1_resolution_warns():
    """The lr->fisher resolution must emit a user-facing warning."""
    from epykit.cli import _cli_n1_and_footgun_checks

    args = _checks_namespace(test="lr", allow_n1=True)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _cli_n1_and_footgun_checks(args, unit="sites")
    msgs = " ".join(str(w.message).lower() for w in caught)
    assert "fisher" in msgs
