"""Tests on residual series (DK §2.12, §7.5).

The tests apply to standardized one-step prediction errors
(:attr:`FilterResults.standardized_forecasts_error
<ssfortran.FilterResults.standardized_forecasts_error>`) after the diffuse
periods. NaNs are skipped.
"""

import ctypes

import numpy as np

from ._lib import call, ptr


def _x(x):
    return np.ascontiguousarray(np.asarray(x, dtype=np.float64).ravel())


def ljung_box(x, lags=10, model_df=0):
    """Test for serial correlation with the Ljung-Box statistic.

    :math:`Q(k) = n (n + 2) \\sum_{j=1}^k \\rho_j^2 / (n - j)`, referred to
    :math:`\\chi^2` with k - `model_df` degrees of freedom (DK §2.12.1).

    Parameters
    ----------
    x : array_like
        Residuals.
    lags : int, optional
        Largest lag k.
    model_df : int, optional
        Degrees of freedom to subtract, e.g. the number of ARMA
        coefficients.

    Returns
    -------
    stat : ndarray, shape (lags,)
        :math:`Q(1), \\dots, Q(k)`.
    pvalue : ndarray, shape (lags,)
        Their p-values.

    Notes
    -----
    Pairs with a NaN are skipped, and n is the number of non-NaN values.

    Examples
    --------
    >>> e = np.random.default_rng(0).standard_normal(200)
    >>> stat, pvalue = ss.diagnostics.ljung_box(e, lags=5)
    >>> stat.shape
    (5,)
    """
    x = _x(x)
    stat, pvalue = np.empty(lags), np.empty(lags)
    call("ss_ljung_box", x.size, ptr(x), int(lags), int(model_df), ptr(stat), ptr(pvalue))
    return stat, pvalue


def jarque_bera(x):
    """Test for normality with the Jarque-Bera statistic.

    :math:`JB = n (S^2 / 6 + (K - 3)^2 / 24)` with skewness S and kurtosis
    K, referred to :math:`\\chi^2_2` (DK §2.12.1).

    Parameters
    ----------
    x : array_like
        Residuals.

    Returns
    -------
    jb : float
        The statistic.
    pvalue : float
        Its p-value.
    skew : float
        Sample skewness.
    kurtosis : float
        Sample kurtosis (3 for the normal distribution).
    """
    x = _x(x)
    out = [ctypes.c_double() for _ in range(4)]
    call("ss_jarque_bera", x.size, ptr(x), *(ctypes.byref(o) for o in out))
    return tuple(o.value for o in out)


def breakvar(x, h=None):
    """Test for heteroskedasticity by comparing the two ends of the sample.

    :math:`H(h) = \\sum_{t=n-h+1}^n e_t^2 / \\sum_{t=1}^h e_t^2`, referred to
    the F(h, h) distribution (DK §2.12.1).

    Parameters
    ----------
    x : array_like
        Residuals.
    h : int, optional
        Number of observations at each end; default a third of the sample.

    Returns
    -------
    stat : float
        H(h).
    pvalue : float
        Two-sided p-value.
    """
    x = _x(x)
    stat, pvalue = ctypes.c_double(), ctypes.c_double()
    call("ss_breakvar", x.size, ptr(x), 0 if h is None else int(h), ctypes.byref(stat),
         ctypes.byref(pvalue))
    return stat.value, pvalue.value
