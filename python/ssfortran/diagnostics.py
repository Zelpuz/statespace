"""Residual diagnostics (DK 2.12, 7.5) on a series, NaNs skipped."""

import ctypes

import numpy as np

from ._lib import call, ptr


def _x(x):
    return np.ascontiguousarray(np.asarray(x, dtype=np.float64).ravel())


def ljung_box(x, lags=10, model_df=0):
    """Ljung-Box statistics Q(k), k = 1..lags, and their p-values (with
    lags - model_df degrees of freedom)."""
    x = _x(x)
    stat, pvalue = np.empty(lags), np.empty(lags)
    call("ss_ljung_box", x.size, ptr(x), int(lags), int(model_df), ptr(stat), ptr(pvalue))
    return stat, pvalue


def jarque_bera(x):
    """Jarque-Bera normality test: (statistic, p-value, skewness, kurtosis)."""
    x = _x(x)
    out = [ctypes.c_double() for _ in range(4)]
    call("ss_jarque_bera", x.size, ptr(x), *(ctypes.byref(o) for o in out))
    return tuple(o.value for o in out)


def breakvar(x, h=None):
    """Test of equal variances in the first and last h observations
    (default a third of the sample): (statistic, two-sided p-value)."""
    x = _x(x)
    stat, pvalue = ctypes.c_double(), ctypes.c_double()
    call("ss_breakvar", x.size, ptr(x), 0 if h is None else int(h), ctypes.byref(stat),
         ctypes.byref(pvalue))
    return stat.value, pvalue.value
