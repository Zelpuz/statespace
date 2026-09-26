"""ssfortran: linear Gaussian state space models (Durbin and Koopman 2012,
Part I), computed by a Fortran library.

The matrix-level API is `Representation`; built-in models and estimation
are added on top of it.
"""

from . import diagnostics
from ._lib import StateSpaceError, version
from .models import (
    ARIMA,
    Cycle,
    ContinuousLevel,
    ContinuousTrend,
    FitResults,
    Irregular,
    Level,
    MappedModel,
    MLEModel,
    Model,
    Regression,
    Seasonal,
    StructuralModel,
    Trend,
    fit_many,
)
from .representation import (
    DIFFUSE_MULTIVARIATE,
    DIFFUSE_UNIVARIATE,
    FILTER_CONVENTIONAL,
    FILTER_UNIVARIATE,
    INIT_APPROX_DIFFUSE,
    INIT_DIFFUSE,
    INIT_GENERAL,
    INIT_KNOWN,
    INIT_STATIONARY,
    FilterResults,
    Representation,
    SmootherResults,
)

__version__ = version()

__all__ = [
    "diagnostics", "Model", "StructuralModel", "MappedModel", "MLEModel", "FitResults", "fit_many",
    "Irregular", "Level", "Trend", "Seasonal", "Cycle", "Regression", "ARIMA",
    "ContinuousLevel", "ContinuousTrend",
    "Representation", "FilterResults", "SmootherResults", "StateSpaceError",
    "INIT_KNOWN", "INIT_APPROX_DIFFUSE", "INIT_STATIONARY", "INIT_DIFFUSE", "INIT_GENERAL",
    "FILTER_CONVENTIONAL", "FILTER_UNIVARIATE", "DIFFUSE_UNIVARIATE", "DIFFUSE_MULTIVARIATE",
]
