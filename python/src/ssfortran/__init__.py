"""Linear Gaussian state space models with a Fortran core.

``ssfortran`` covers Part I of Durbin and Koopman, *Time Series Analysis by
State Space Methods* (2nd ed., 2012), cited as DK:

* :class:`StructuralModel`, :class:`MappedModel`, :class:`MLEModel`:
  models with parameters, estimated by maximum likelihood;
* :class:`Representation`: the system at given parameter values, which
  each model builds, and the algorithms of DK ch. 4-7 (filtering,
  smoothing, likelihood, simulation, forecasting, diagnostics);
* :mod:`ssfortran.diagnostics`: residual tests.
"""

from . import diagnostics
from ._lib import StateSpaceError, version
from .models import (
    ARIMA,
    ContinuousLevel,
    ContinuousTrend,
    Cycle,
    FitResults,
    ForecastResults,
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
    AugmentedResults,
    FilterResults,
    Representation,
    SmootherResults,
)

__version__ = version()

__all__ = [
    "diagnostics",
    "Model",
    "StructuralModel",
    "MappedModel",
    "MLEModel",
    "FitResults",
    "ForecastResults",
    "fit_many",
    "Irregular",
    "Level",
    "Trend",
    "Seasonal",
    "Cycle",
    "Regression",
    "ARIMA",
    "ContinuousLevel",
    "ContinuousTrend",
    "Representation",
    "AugmentedResults",
    "FilterResults",
    "SmootherResults",
    "StateSpaceError",
    "INIT_KNOWN",
    "INIT_APPROX_DIFFUSE",
    "INIT_STATIONARY",
    "INIT_DIFFUSE",
    "INIT_GENERAL",
    "FILTER_CONVENTIONAL",
    "FILTER_UNIVARIATE",
    "DIFFUSE_UNIVARIATE",
    "DIFFUSE_MULTIVARIATE",
]
