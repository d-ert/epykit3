"""epykit: Parquet methylation pipeline (Bismark ingestion -> DMR annotation)

This package provides a complete WGBS analysis pipeline:

  - convert:   Bismark .cov -> partitioned Parquet methylstore
  - filter:    QC / coverage filtering / site intersection
  - dmc:       Differential methylation calling per CpG. Default test is
               ``lr`` (quasi-binomial likelihood-ratio with McCullagh-
               Nelder dispersion). Other backends: score, glm, logit_t,
               welch_t, fisher.
  - dmr:       DMR calling (tile-based default, sliding-window legacy) and
               Gaussian-kernel methylation smoothing
  - annotate:  Gene-feature and CpG-island context annotation
  - qc:        Bisulfite conversion rate, global methylation, coverage
               uniformity, plus opt-in clinical / cohort checks

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

from importlib.metadata import version as _v, PackageNotFoundError

try:
    __version__ = _v("epykit")
except PackageNotFoundError:
    # editable install or running from source without install
    __version__ = "0.0.0+unknown"

from .methyldata import MethylData
from .io import (
    read_bismark,
    read_methyldackel,
    read_combined_strand_bed,
    read_nfcore_methylseq,
    load,
)
from . import pp, tl, pl, query
from ._config import set_tmp_dir, get_tmp_dir

from .convert import convert_sample
from ._dmc_store import DMCStore
from .dmc import (
    process_chromosomes_dmc,
    apply_multiple_testing_correction,
    empirical_fdr_for_dmc,
    fisher_exact_vectorized,
    shrink_meth_diff,
)
from .dmr import (
    DMR_PRESETS,
    call_dmr_chain_merge,
    call_dmr_sliding_window,
    smooth_methylation_gaussian,
    smooth_methylation_bsmooth,
)
from .annotate import (
    annotate_features,
    annotate_cpg_islands,
    HOMER_FEATURES,
)
from .qc import (
    bisulfite_conversion_rate,
    global_methylation_report,
    coverage_uniformity,
)
from ._glm import build_design

# Export / interop (lazy heavy deps inside)
from .export import to_bedgraph, to_bigwig, dmcs_to_bed, dmrs_to_bed
from .anndata_io import to_anndata, to_mudata
from .methylkit_io import to_methylkit_tabix
from .multiqc_export import report_multiqc
from .nfcore_qc import read_nfcore_methylseq_qc
from .report import generate_report
from .dvc import process_chromosomes_dvc, call_dvr_density  # noqa: F401
from .impute import impute_knn_beta, impute_knn_anndata
from .clocks import age_clock, deconvolve
from .qc import (
    sex_check,
    contamination_estimate,
    sample_correlation as sample_correlation_qc,
    power as power_calc,
)

__all__ = [
    # version
    "__version__",
    # data object
    "MethylData",
    # I/O
    "read_bismark", "read_methyldackel", "read_combined_strand_bed",
    "read_nfcore_methylseq", "load",
    # namespaces (scanpy-style)
    "pp", "tl", "pl", "query",
    # runtime config
    "set_tmp_dir", "get_tmp_dir",
    # ingestion
    "convert_sample",
    # DMC engines (advanced users; tl.dmc is the recommended entry)
    "process_chromosomes_dmc",
    "apply_multiple_testing_correction",
    "empirical_fdr_for_dmc",
    "fisher_exact_vectorized",
    "shrink_meth_diff",
    "DMCStore",
    # DMR engines
    "call_dmr_sliding_window",
    "call_dmr_chain_merge",
    "DMR_PRESETS",
    "smooth_methylation_gaussian",
    "smooth_methylation_bsmooth",
    # DVC / DVR engines
    "process_chromosomes_dvc",
    "call_dvr_density",
    # Imputation
    "impute_knn_beta",
    "impute_knn_anndata",
    # Clocks / deconvolution
    "age_clock",
    "deconvolve",
    # annotation
    "annotate_features",
    "HOMER_FEATURES",
    "annotate_cpg_islands",
    # QC
    "bisulfite_conversion_rate",
    "global_methylation_report",
    "coverage_uniformity",
    "sex_check",
    "contamination_estimate",
    "sample_correlation_qc",
    "power_calc",
    # GLM design
    "build_design",
    # Exports / interop
    "to_bedgraph",
    "to_bigwig",
    "dmcs_to_bed",
    "dmrs_to_bed",
    "to_anndata",
    "to_mudata",
    "to_methylkit_tabix",
    "report_multiqc",
    "read_nfcore_methylseq_qc",
    "generate_report",
]
