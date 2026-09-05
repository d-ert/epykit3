"""epykit: Parquet methylation pipeline (Bismark ingestion -> DMR annotation)

This package provides a complete WGBS analysis pipeline:

  - convert:   Bismark .cov -> partitioned Parquet methylstore
  - filter:    QC / coverage filtering / site intersection
  - dmc:       Differential methylation calling per CpG. Default test is
               ``lr`` (quasi-binomial likelihood-ratio with McCullagh-
               Nelder dispersion). Other backends: welch_t, fisher, glm.
  - dmr:       DMR calling (chain-merge default, with tile, sliding-window
               and HMM-segmentation alternatives) and Gaussian-kernel
               methylation smoothing
  - annotate:  Gene-feature and CpG-island context annotation
  - qc:        Bisulfite conversion rate, global methylation, coverage
               uniformity, plus opt-in clinical / cohort checks

DMC engine choice
-----------------
Bare ``lr`` is the recommended default and the engine the benchmark paper
characterises. ``welch_t``, ``fisher`` and ``glm`` are alternative
engines for specific situations (small samples, single replicates,
covariate-adjusted designs). The ``power_stack="lr+"`` tunable is an
*exploratory* research knob: it bundles four extensions (neighbour
combine, two-stage BH, separation fallback, EB dispersion shrinkage)
that were tuned on the Piao 2021 simulator and have **not** been shown
to improve on bare ``lr`` for real-world WGBS at realistic
overdispersion (real WGBS phi ~ 1.5-5; simulator phi ~ 0.4). On
GSE263850 ``lr+`` inflates DMC counts ~13x at the same q-threshold,
consistent with FPR drift. Use ``lr+`` only if you understand its
heuristics and are willing to validate them on your own data. See
``docs/architecture.md`` for the engine map.

Logging convention
------------------
Library modules (everything under ``epykit.*`` except ``epykit.cli``) emit
progress and diagnostics through the standard :mod:`logging` module via
``logger = logging.getLogger(__name__)`` -- they never call :func:`print`.
The CLI entry point (``epykit.cli``) reserves :func:`print` for the
final user-facing result lines on stdout; structured progress logs flow
through the same logging hierarchy and are controlled via ``-v`` / ``-q``.
This split lets host applications and notebooks consume epykit without
having their stdout polluted, while CLI users see the expected output.
"""

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _v

try:
    __version__ = _v("epykit")
except PackageNotFoundError:
    # editable install or running from source without install
    __version__ = "0.0.0+unknown"

from . import pl, pp, query, tl
from ._config import get_tmp_dir, set_tmp_dir
from ._dmc_store import DMCStore
from ._glm import build_design
from .anndata_io import to_anndata, to_mudata
from .annotate import (
    HOMER_FEATURES,
    annotate_cpg_islands,
    annotate_features,
)
from .clocks import age_clock, deconvolve
from .convert import convert_sample
from .dmr import (
    DMR_PRESETS,
    call_dmr_chain_merge,
    call_dmr_sliding_window,
    smooth_methylation_bsmooth,
    smooth_methylation_gaussian,
)
from .dvc import call_dvr_density, process_chromosomes_dvc

# Export / interop (lazy heavy deps inside)
from .export import dmcs_to_bed, dmrs_to_bed, export_tables, to_bedgraph, to_bigwig
from .impute import impute_knn_anndata, impute_knn_beta
from .io import (
    load,
    read_bismark,
    read_combined_strand_bed,
    read_methyldackel,
    read_nfcore_methylseq,
)
from .methyldata import MethylData
from .methylkit_io import to_methylkit_tabix
from .multiqc_export import report_multiqc
from .nfcore_qc import read_nfcore_methylseq_qc
from .qc import (
    bisulfite_conversion_rate,
    contamination_estimate,
    coverage_uniformity,
    global_methylation_report,
    sex_check,
)
from .qc import (
    power as power_calc,
)
from .qc import (
    sample_correlation as sample_correlation_qc,
)
from .report import generate_report

__all__ = [
    "DMR_PRESETS",
    "HOMER_FEATURES",
    # DMCStore is the streaming-store handle returned by tl.dmc(..., return_store=True).
    # The lower-level dmc engine functions (process_chromosomes_dmc, apply_multiple_testing_correction,
    # empirical_fdr_for_dmc, fisher_exact_vectorized, shrink_meth_diff) moved to
    # epykit.dmc submodule in 1.0. Import them via `from epykit.dmc import ...`.
    "DMCStore",
    "MethylData",
    "__version__",
    "age_clock",
    "annotate_cpg_islands",
    "annotate_features",
    "bisulfite_conversion_rate",
    "build_design",
    "call_dmr_chain_merge",
    "call_dmr_sliding_window",
    "call_dvr_density",
    "contamination_estimate",
    "convert_sample",
    "coverage_uniformity",
    "deconvolve",
    "dmcs_to_bed",
    "dmrs_to_bed",
    "export_tables",
    "generate_report",
    "get_tmp_dir",
    "global_methylation_report",
    "impute_knn_anndata",
    "impute_knn_beta",
    "load",
    "pl",
    "power_calc",
    "pp",
    "process_chromosomes_dvc",
    "query",
    "read_bismark",
    "read_combined_strand_bed",
    "read_methyldackel",
    "read_nfcore_methylseq",
    "read_nfcore_methylseq_qc",
    "report_multiqc",
    "sample_correlation_qc",
    "set_tmp_dir",
    "sex_check",
    "smooth_methylation_bsmooth",
    "smooth_methylation_gaussian",
    "tl",
    "to_anndata",
    "to_bedgraph",
    "to_bigwig",
    "to_methylkit_tabix",
    "to_mudata",
]

# --- Deprecation shim for demoted top-level names (1.0; removed in 1.2) ---
_DEMOTED_TO_DMC = frozenset(
    {
        "process_chromosomes_dmc",
        "apply_multiple_testing_correction",
        "empirical_fdr_for_dmc",
        "fisher_exact_vectorized",
        "shrink_meth_diff",
    }
)


def __getattr__(name: str):
    """Module-level __getattr__ provides backward compat for names removed
    from `__all__` at 1.0.

    Returns the function from the appropriate submodule and emits
    DeprecationWarning pointing users at the new import path. Shim is
    scheduled for removal in 1.2.
    """
    if name in _DEMOTED_TO_DMC:
        import warnings

        from . import dmc as _dmc_mod

        warnings.warn(
            f"epykit.{name} is no longer a top-level export; use "
            f"`from epykit.dmc import {name}` instead, or use the "
            f"recommended `epykit.tl.dmc` wrapper. This shim will be "
            f"removed in epykit 1.2.",
            DeprecationWarning,
            stacklevel=2,
        )
        return getattr(_dmc_mod, name)
    raise AttributeError(f"module 'epykit' has no attribute {name!r}")
