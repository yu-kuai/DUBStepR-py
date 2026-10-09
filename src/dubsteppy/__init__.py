"""dubsteppy — DUBStepR: correlation-based feature selection for single-cell RNA-seq data.

Python port of the R package DUBStepR (https://github.com/prabhakarlab/DUBStepR),
Ranjan et al., Nature Communications 2021 (doi:10.1038/s41467-021-26085-2).
"""

from .core import DUBStepRResult, dubstepr, dubstepr_anndata, log_normalize
from .datasets import pbmc_norm_small_data
from .filtering import get_filtered_data, load_gene_annotation
from .ggc import get_correlation_range, get_ggc
from .optimal import get_optimal_feature_set
from .stepwise import find_elbow, run_stepwise_reg

__version__ = "1.2.0"

__all__ = [
    "DUBStepRResult",
    "dubstepr",
    "dubstepr_anndata",
    "find_elbow",
    "get_correlation_range",
    "get_filtered_data",
    "get_ggc",
    "get_optimal_feature_set",
    "load_gene_annotation",
    "log_normalize",
    "pbmc_norm_small_data",
    "run_stepwise_reg",
]
