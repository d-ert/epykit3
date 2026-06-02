"""DEPRECATED — renamed to ``epykit.dmr_segment`` in 0.7.5.

This shim re-exports the renamed function and emits a
``DeprecationWarning`` on import. Scheduled for removal in 1.1.

Rationale: the engine uses fixed state means and transition
probabilities (not Baum-Welch fitted), so calling it an HMM
was misleading.
"""

from __future__ import annotations

import logging
import warnings

from .dmr_segment import call_dmr_rule_segment

_msg = (
    "epykit.dmr_hmm is deprecated and will be removed in 1.1; "
    "use epykit.dmr_segment.call_dmr_rule_segment instead"
)
warnings.warn(_msg, DeprecationWarning, stacklevel=2)
logging.getLogger(__name__).warning(_msg)

# Legacy export name preserved.
call_dmr_hmm = call_dmr_rule_segment

__all__ = ["call_dmr_hmm", "call_dmr_rule_segment"]
