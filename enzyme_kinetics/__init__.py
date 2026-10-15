"""Enzymatic Kinetics Analyzer: Michaelis-Menten, dose-response and inhibition-model fitting
with t-based confidence intervals and replicate-aware batch processing.

Author: Oluwaseun O. Ajayi
License: MIT
"""

from .batch import batch_dose_response, batch_michaelis_menten, dose_response_warnings, summarize_replicates
from .fitting import (FitResult, compare_inhibition_models, fit_dose_response, fit_inhibition,
                      fit_michaelis_menten)
from .io import read_compound_table, read_replicate_table
from .models import competitive, four_parameter_logistic, michaelis_menten, mixed, uncompetitive

__version__ = "2.0.0"

__all__ = [
    "FitResult", "fit_michaelis_menten", "fit_dose_response", "fit_inhibition",
    "compare_inhibition_models", "batch_michaelis_menten", "batch_dose_response",
    "summarize_replicates", "dose_response_warnings", "read_replicate_table", "read_compound_table",
    "michaelis_menten", "four_parameter_logistic", "competitive", "uncompetitive", "mixed",
    "__version__",
]
