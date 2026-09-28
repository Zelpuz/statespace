"""Models with parameters and their estimation by maximum likelihood.

A model maps a parameter vector to the system matrices of a
:class:`~ssfortran.Representation`. Three kinds are provided:

* :class:`StructuralModel` assembles built-in components (DK ch. 3):
  irregular, level, trend, seasonal, cycle, regression, ARIMA and
  continuous-time components, for one series or several.
* :class:`MappedModel` is declared: each parameter sets given matrix
  entries, with a transform to the optimizer's unconstrained scale.
* :class:`MLEModel` is subclassed, as statsmodels' ``MLEModel``: Python code
  sets the matrices from the parameters.

The first two run entirely in Fortran and can be fitted in parallel with
:func:`fit_many`; the third calls Python at every likelihood evaluation.
Declared and Python models are set up on the model itself: fixed matrices
with ``mod[name] = value``, the initialization with the ``initialize_*``
methods.

Estimation (:meth:`Model.fit`) maximizes the log likelihood with L-BFGS-B
over the unconstrained parameters. The gradient uses the analytic score of
DK §7.3.3 for parameters in H, R, Q and a stationary initial variance, and
central differences for parameters in Z or T.
"""

import ctypes
import functools
import math
import warnings
import weakref
from dataclasses import dataclass, field
from statistics import NormalDist

import numpy as np

from ._lib import TRANSFORM_CB, UPDATE_CB, StateSpaceError, call, farray, ptr
from .representation import _ARR, FilterResults, Representation, SmootherResults

_COV = {None: 0, "none": 0, "diagonal": 1, "full": 2}
_SEASONAL = {"dummy": 1, "trig": 2, "trigonometric": 2, "harrison-stevens": 3, "hs": 3}
GRADIENT_AUTO, GRADIENT_NUMERICAL, GRADIENT_ANALYTIC = 0, 1, 2
_GRADIENT = {
    "auto": GRADIENT_AUTO,
    "numerical": GRADIENT_NUMERICAL,
    "analytic": GRADIENT_ANALYTIC,
}


def _cov(c):
    try:
        return _COV[c if c is None else c.lower()]
    except KeyError:
        raise ValueError(
            f"covariance must be one of None, 'diagonal', 'full', not {c!r}"
        )


# ---------------------------------------------------------------- components


@dataclass
class Irregular:
    """Observation disturbance :math:`\\varepsilon_t \\sim N(0, \\Sigma_\\varepsilon)`.

    Always acts on the observations.

    Parameters
    ----------
    cov : {"diagonal", "full", None}
        Form of :math:`\\Sigma_\\varepsilon` across series; ``"full"``
        estimates it through its Cholesky factor. ``None`` removes the
        irregular.
    weights : array_like, shape (n,), optional
        :math:`H_t = w_t \\Sigma_\\varepsilon`, for unequally spaced or
        aggregated observations (DK §3.8).
    """

    cov: str = "diagonal"
    weights: object = None


@dataclass
class Level:
    """Random walk level :math:`\\mu_{t+1} = \\mu_t + \\xi_t` (DK ch. 2, §3.2.1).

    Parameters
    ----------
    cov : {"diagonal", "full", None}
        Form of :math:`Var(\\xi_t)` across series (DK §3.3); ``None`` gives
        a fixed level.
    at_observations : bool
        Act on the observations rather than the signals when the model has
        loadings.
    """

    cov: str = "diagonal"
    at_observations: bool = False


@dataclass
class Trend:
    """Local linear trend (DK §3.2.1).

    :math:`\\mu_{t+1} = \\mu_t + \\nu_t + \\xi_t`,
    :math:`\\nu_{t+1} = \\nu_t + \\zeta_t`.

    Parameters
    ----------
    level_cov, slope_cov : {"diagonal", "full", None}
        Forms of :math:`Var(\\xi_t)` and :math:`Var(\\zeta_t)`.
        ``level_cov=None`` gives the smooth trend (integrated random walk);
        both ``None`` a deterministic linear trend.
    at_observations : bool
        Act on the observations rather than the signals.
    """

    level_cov: str = "diagonal"
    slope_cov: str = "diagonal"
    at_observations: bool = False


@dataclass
class Seasonal:
    """Seasonal component of a given period (DK §3.2.2).

    Parameters
    ----------
    period : int
        Number of periods per cycle, e.g. 12 for monthly data.
    form : {"dummy", "trig", "harrison-stevens"}
        Dummy (s - 1 states), trigonometric (s - 1 states, all harmonics
        with one variance) or Harrison-Stevens (s states).
    cov : {"diagonal", "full", None}
        Form of the disturbance variance; ``None`` for a fixed seasonal.
    at_observations : bool
        Act on the observations rather than the signals.
    """

    period: int
    form: str = "dummy"
    cov: str = "diagonal"
    at_observations: bool = False


@dataclass
class Cycle:
    """Stochastic cycle (DK §3.2.4).

    :math:`(c_{t+1}, c^*_{t+1})' = \\rho C(\\lambda) (c_t, c^*_t)' + w_t` with
    :math:`C(\\lambda)` the rotation by frequency :math:`\\lambda`.
    Parameters: the disturbance variance, :math:`\\lambda` and, if damped,
    :math:`\\rho`.

    Parameters
    ----------
    cov : {"diagonal", "full", None}
        Form of :math:`Var(w_t)`.
    damped : bool
        Estimate :math:`\\rho \\in (0, 1)`; the cycle is then stationary and
        starts from its unconditional distribution. Otherwise
        :math:`\\rho = 1` and the cycle is diffuse.
    period_bounds : tuple of float
        Bounds on the period :math:`2\\pi / \\lambda`, in observations;
        ``None`` as upper bound means the sample size.
    at_observations : bool
        Act on the observations rather than the signals.
    """

    cov: str = "diagonal"
    damped: bool = True
    period_bounds: tuple = (2.0, None)
    at_observations: bool = False


@dataclass
class Regression:
    """Regression effects :math:`x_t' \\beta_t` (DK §3.2.5, §3.6).

    The coefficients are diffuse states, fixed or random walks.

    Parameters
    ----------
    exog : array_like, shape (n,) or (n, k)
        Regressors; also interventions (step, pulse or slope dummies).
    random_walk : bool or sequence of bool
        Make all, or the flagged, coefficients random walks with estimated
        variances (DK eq. 3.15).
    series : int
        The series (0-based) the effects enter.
    at_observations : bool
        Act on the observations rather than the signals.
    """

    exog: object
    random_walk: object = False
    series: int = 0
    at_observations: bool = False


@dataclass
class ARIMA:
    """ARIMA(p, d, q)(P, D, Q)_s component (DK §3.4).

    The differences are states (diffuse), as in DK §3.4 and statsmodels'
    SARIMAX; the ARMA part is stationary. Parameters: the AR, seasonal AR,
    MA and seasonal MA coefficients, then the innovation variance.

    Parameters
    ----------
    order : tuple of int
        (p, d, q).
    seasonal_order : tuple of int
        (P, D, Q, s), with D at most 1.
    series : int
        The series (0-based).
    enforce_stationarity, enforce_invertibility : bool
        Keep the AR polynomials stationary and the MA polynomials
        invertible (Monahan 1984).
    at_observations : bool
        Act on the observations rather than the signals.
    """

    order: tuple = (1, 0, 0)
    seasonal_order: tuple = (0, 0, 0, 0)
    series: int = 0
    enforce_stationarity: bool = True
    enforce_invertibility: bool = True
    at_observations: bool = False


@dataclass
class ContinuousLevel:
    """Level in continuous time observed at arbitrary times (DK §3.8.1).

    Between observations the variance of the change is
    :math:`\\sigma^2 \\delta_t`, :math:`\\delta_t = \\tau_{t+1} - \\tau_t`.

    Parameters
    ----------
    times : array_like, shape (n,)
        Increasing observation times; repeated times are allowed.
    at_observations : bool
        Act on the observations rather than the signals.
    """

    times: object
    at_observations: bool = False


@dataclass
class ContinuousTrend:
    """Smooth trend in continuous time observed at arbitrary times (DK §3.8.2).

    With an :class:`Irregular`, the smoothed level is the cubic smoothing
    spline with smoothing parameter :math:`\\sigma^2_\\varepsilon / \\sigma^2`
    (DK §3.9.2).

    Parameters
    ----------
    times : array_like, shape (n,)
        Increasing observation times; repeated times are allowed.
    at_observations : bool
        Act on the observations rather than the signals.
    """

    times: object
    at_observations: bool = False


# ------------------------------------------------------------------ results


@dataclass
class ForecastResults:
    """Forecasts past the end of the sample (DK §4.11).

    Attributes
    ----------
    predicted_mean : ndarray, shape (steps,) or (steps, p)
        :math:`E(y_{n+j} | Y_n)`; one column per series.
    var_pred_mean : ndarray, shape (steps, p, p)
        :math:`Var(y_{n+j} | Y_n)`.
    """

    predicted_mean: np.ndarray
    var_pred_mean: np.ndarray

    @property
    def se_mean(self):
        """Standard errors of the forecasts, shaped as `predicted_mean`."""
        se = np.sqrt(np.diagonal(self.var_pred_mean, axis1=1, axis2=2))
        return se if se.shape[1] > 1 else se[:, 0]

    def conf_int(self, alpha=0.05):
        """Compute pointwise normal confidence intervals.

        Parameters
        ----------
        alpha : float, optional
            One minus the coverage.

        Returns
        -------
        lower, upper : ndarray
            Bounds, shaped as `predicted_mean`.
        """
        z = NormalDist().inv_cdf(1 - alpha / 2)
        return (
            self.predicted_mean - z * self.se_mean,
            self.predicted_mean + z * self.se_mean,
        )


@dataclass
class FitResults:
    """Maximum likelihood estimates and the fitted model (DK §7.3).

    Attributes
    ----------
    model : Model
        The model that was fitted.
    params : ndarray, shape (k,)
        Estimates, constrained.
    param_names : list of str
        Parameter names.
    bse : ndarray, shape (k,)
        Standard errors; NaN if the Hessian was not negative definite.
    cov_params : ndarray, shape (k, k)
        Covariance of the estimates: the inverse of the negative Hessian in
        the unconstrained parameters, mapped to the constrained ones by the
        delta method (DK §7.3.6).
    llf : float
        Maximized log likelihood.
    scale : float
        Estimated scale when it is concentrated out, else 1.
    aic, bic : float
        Information criteria, counting the diffuse states and a
        concentrated scale as parameters (DK §7.4, not divided by n).
    niter, nfev : int
        Optimizer iterations and log likelihood evaluations.
    converged : bool
        Whether L-BFGS-B met its convergence criterion.
    analytic_gradient : bool
        Whether the analytic score was used for at least one parameter.
    message : str
        The optimizer's final message.
    """

    model: "Model" = field(repr=False)
    params: np.ndarray
    param_names: list
    bse: np.ndarray
    cov_params: np.ndarray
    llf: float
    scale: float
    aic: float
    bic: float
    niter: int
    nfev: int
    converged: bool
    analytic_gradient: bool
    message: str

    def smooth(self):
        """Run the filter and smoother at the estimates.

        Returns
        -------
        SmootherResults
            Smoothed states and disturbances.
        """
        return self.model.smooth(self.params)

    def filter(self):
        """Run the filter at the estimates.

        Returns
        -------
        FilterResults
            The filter output.
        """
        return self.model.filter(self.params)

    def estimation_bias(self, ndraw=1000, antithetic=True, seed=None, variance=False):
        """Estimate the bias in the smoothed state from estimating the parameters.

        DK §7.3.7, eq. 7.20: treating :math:`\\hat\\psi` as the true value,
        draw :math:`\\psi^{(i)}` from :math:`N(\\hat\\psi, \\Omega)` on the
        unconstrained scale and average
        :math:`\\hat\\alpha(\\psi^{(i)}) - \\hat\\alpha(\\hat\\psi)`.

        Parameters
        ----------
        ndraw : int, optional
            Number of draws; even when `antithetic`.
        antithetic : bool, optional
            Pair each draw with its reflection about :math:`\\hat\\psi`.
        seed : int, optional
            Seed of the library's random numbers, for reproducible results.
        variance : bool, optional
            Also return the bias in the smoothed state variance.

        Returns
        -------
        bias_alpha : ndarray, shape (m, n)
            Bias in the smoothed state.
        bias_V : ndarray, shape (m, m, n)
            Bias in its variance; only with ``variance=True``.

        Raises
        ------
        ValueError
            If `cov_params` is not available.
        """
        if not np.all(np.isfinite(self.cov_params)):
            raise ValueError(
                "estimation_bias needs cov_params (fit with compute_cov=True)"
            )
        m, n = self.model.k_states, self.model.nobs
        ba = np.empty((m, n), order="F")
        bV = np.empty((m, m, n), order="F") if variance else None
        p = farray(self.params, (self.model.k_params,))
        c = farray(self.cov_params)
        failed = ctypes.c_int()
        self.model._call(
            "ss_estimation_bias",
            self.model._h,
            ptr(p),
            ptr(c),
            int(ndraw),
            int(bool(antithetic)),
            0 if seed is None else int(seed),
            ptr(ba),
            None if bV is None else ptr(bV),
            ctypes.byref(failed),
        )
        self.model.loglike(self.params)  # leave the model at the estimates
        return (ba, bV) if variance else ba

    # -- per-period output ------------------------------------------------

    @staticmethod
    def _series(a):
        """(p, n) -> (n,) for one series, else (n, p)."""
        a = np.asarray(a).T
        return a[:, 0] if a.shape[1] == 1 else a

    @property
    def fittedvalues(self):
        """One-step predictions of y, shape (n,) or (n, p)."""
        return self._series(self.filter().forecasts)

    @property
    def resid(self):
        """One-step prediction errors :math:`v_t`, shape (n,) or (n, p)."""
        return self._series(self.filter().forecasts_error)

    @property
    def standardized_residuals(self):
        """Standardized prediction errors (DK §7.5); NaN where missing or diffuse."""
        return self._series(self.filter().standardized_forecasts_error)

    def components(self, variance=False):
        """Compute the smoothed components at the estimates.

        Parameters
        ----------
        variance : bool, optional
            Also return their variances.

        Returns
        -------
        dict or tuple of dict
            See :meth:`Model.components`.
        """
        return self.model.components(self.params, variance=variance)

    def get_forecast(self, steps):
        """Forecast past the end of the sample (DK §4.11).

        Parameters
        ----------
        steps : int
            Number of periods ahead.

        Returns
        -------
        ForecastResults
            Forecasts and their variances.

        Raises
        ------
        StateSpaceError
            With code 4 for models with time-varying system matrices.
        """
        rep = self.model.representation(self.params)
        mean, cov = rep.forecast(int(steps))
        return ForecastResults(
            predicted_mean=mean.T if mean.shape[0] > 1 else mean[0],
            var_pred_mean=np.moveaxis(cov, 2, 0),
        )

    def forecast(self, steps):
        """Forecast past the end of the sample.

        Parameters
        ----------
        steps : int
            Number of periods ahead.

        Returns
        -------
        ndarray
            ``get_forecast(steps).predicted_mean``.
        """
        return self.get_forecast(steps).predicted_mean

    # -- summaries ---------------------------------------------------------

    def diagnostics(self, lags=None):
        """Test the standardized residuals (DK §2.12, §7.5).

        Uses the first series, after the diffuse periods.

        Parameters
        ----------
        lags : int, optional
            Lags of the Ljung-Box test; default min(10, n/5).

        Returns
        -------
        dict
            ``ljung_box``: (lags, Q, p-value); ``jarque_bera``: (statistic,
            p-value); ``skew``, ``kurtosis``; ``heteroskedasticity``: (H,
            p-value).
        """
        from . import diagnostics as dg

        f = self.filter()
        e = f.standardized_forecasts_error[0, f.diagnostic_start :]
        e = e[~np.isnan(e)]
        lags = lags or max(1, min(10, e.size // 5))
        q, qp = dg.ljung_box(e, lags)
        jb, jbp, skew, kurt = dg.jarque_bera(e)
        h, hp = dg.breakvar(e)
        return {
            "ljung_box": (lags, q[-1], qp[-1]),
            "jarque_bera": (jb, jbp),
            "skew": skew,
            "kurtosis": kurt,
            "heteroskedasticity": (h, hp),
        }

    def summary(self, alpha=0.05):
        """Summarize the estimates, the fit and the residual tests.

        Parameters
        ----------
        alpha : float, optional
            One minus the coverage of the confidence intervals.

        Returns
        -------
        str
            A table of estimates with standard errors, z statistics,
            p-values and confidence intervals; the log likelihood, AIC, BIC
            and number of diffuse periods; and the tests of
            :meth:`diagnostics`.
        """
        z = NormalDist().inv_cdf(1 - alpha / 2)
        f = self.filter()
        w = max(12, max(len(n) for n in self.param_names))
        lo_label, hi_label = f"[{alpha / 2:.3f}", f"{1 - alpha / 2:.3f}]"
        rule = "=" * (w + 72)
        lines = [
            rule,
            f"{'Model:':<16}{type(self.model).__name__:<24}"
            f"{'Log likelihood:':<18}{self.llf:>14.4f}",
            f"{'Observations:':<16}{self.model.nobs:<24}{'AIC:':<18}{self.aic:>14.4f}",
            f"{'Diffuse periods:':<16}{f.nobs_diffuse:<24}"
            f"{'BIC:':<18}{self.bic:>14.4f}",
        ]
        if self.model.concentrate_scale:
            lines.append(f"{'Scale:':<16}{self.scale:<24.6g}")
        lines += [
            rule,
            f"{'':{w}}  {'coef':>12} {'std err':>11} {'z':>8} {'P>|z|':>7} "
            f"{lo_label:>12} {hi_label:>12}",
            "-" * (w + 72),
        ]
        for n, p, s in zip(self.param_names, self.params, self.bse):
            if np.isfinite(s) and s > 0:
                zs = p / s
                pv = math.erfc(abs(zs) / math.sqrt(2))
                lines.append(
                    f"{n:{w}}  {p:12.6g} {s:11.4g} {zs:8.3f} {pv:7.3f} "
                    f"{p - z * s:12.6g} {p + z * s:12.6g}"
                )
            else:
                lines.append(
                    f"{n:{w}}  {p:12.6g} {'':11} {'':8} {'':7} {'':12} {'':12}"
                )
        d = self.diagnostics()
        lb, jb, het = d["ljung_box"], d["jarque_bera"], d["heteroskedasticity"]
        lines += [
            rule,
            f"Ljung-Box Q({lb[0]}): {lb[1]:.3f} (p = {lb[2]:.3f})   "
            f"Jarque-Bera: {jb[0]:.3f} (p = {jb[1]:.3f})   "
            f"H: {het[0]:.3f} (p = {het[1]:.3f})",
            f"Optimizer: {self.message}",
            rule,
        ]
        return "\n".join(lines)


# -------------------------------------------------------------------- models


class Model:
    """Base class of models with parameters, held by the Fortran library.

    Not instantiated directly; see :class:`StructuralModel`,
    :class:`MappedModel` and :class:`MLEModel`.

    Attributes
    ----------
    k_params : int
        Number of parameters.
    nobs, k_endog, k_states, k_posdef : int
        n, p, m and r of the representation.
    """

    def _attach(self, handle):
        self._h = handle
        self._finalizer = weakref.finalize(self, _free, "ss_model_free", handle)
        k, p, m, r, n = (ctypes.c_int() for _ in range(5))
        call("ss_model_info", handle, *(ctypes.byref(x) for x in (k, p, m, r, n)))
        self.k_params, self.k_endog, self.k_states = k.value, p.value, m.value
        self.k_posdef, self.nobs = r.value, n.value

    @property
    def ssm(self):
        """The model's own representation.

        Changes to it, such as ``filter_method`` or ``tol_steady``, apply to
        later likelihood evaluations. Its system matrices reflect the last
        parameters evaluated. :class:`MappedModel` and :class:`MLEModel` set
        its matrices and initialization through ``mod[name] = value`` and
        their ``initialize_*`` methods.
        """
        h = ctypes.c_void_p()
        self._call("ss_model_rep", self._h, ctypes.byref(h))
        return Representation._wrap(h, owner=self)

    @property
    def param_names(self):
        """Parameter names, e.g. ``"sigma2.level"``."""
        width = 32
        buf = ctypes.create_string_buffer(width * self.k_params)
        self._call("ss_model_param_names", self._h, width, buf)
        raw = buf.raw.decode()
        return [raw[i * width : (i + 1) * width].strip() for i in range(self.k_params)]

    @property
    def start_params(self):
        """Starting values for estimation, constrained."""
        out = np.empty(self.k_params)
        self._call("ss_model_start_params", self._h, ptr(out))
        return out

    @property
    def concentrate_scale(self):
        """Whether a scale is concentrated out of the likelihood.

        When true, H, Q and :math:`P_*` set by the parameters are relative
        to a scale :math:`\\sigma^2`, which is estimated in closed form
        (DK §2.10.2) and is not a parameter.
        """
        return getattr(self, "_concentrate", False)

    @concentrate_scale.setter
    def concentrate_scale(self, flag):
        self._call("ss_model_set_concentrate", self._h, int(bool(flag)))
        self._concentrate = bool(flag)

    def _params(self, params):
        a = farray(params, (self.k_params,))
        return a

    def _check_callbacks(self):
        """Re-raise an exception from a Python callback (MLEModel)."""

    def _call(self, name, *args):
        try:
            call(name, *args)
        except StateSpaceError:
            self._check_callbacks()
            raise
        self._check_callbacks()

    def transform_params(self, unconstrained):
        """Map unconstrained values to parameters.

        Parameters
        ----------
        unconstrained : array_like, shape (k,)
            Values on the optimizer's scale.

        Returns
        -------
        ndarray, shape (k,)
            Constrained parameters, e.g. variances :math:`\\exp(2x)`
            (DK §7.3.2).
        """
        x = self._params(unconstrained)
        out = np.empty(self.k_params)
        self._call("ss_model_transform", self._h, ptr(x), ptr(out))
        return out

    def untransform_params(self, constrained):
        """Map parameters to the optimizer's unconstrained scale.

        Parameters
        ----------
        constrained : array_like, shape (k,)
            Parameters.

        Returns
        -------
        ndarray, shape (k,)
            The inverse of :meth:`transform_params`.
        """
        p = self._params(constrained)
        out = np.empty(self.k_params)
        self._call("ss_model_untransform", self._h, ptr(p), ptr(out))
        return out

    def loglike(self, params):
        """Evaluate the log likelihood.

        Parameters
        ----------
        params : array_like, shape (k,)
            Constrained parameters.

        Returns
        -------
        float
            Log likelihood, concentrated when `concentrate_scale` is set; the
            scale is then in :attr:`scale`.
        """
        p = self._params(params)
        llf = ctypes.c_double()
        self._call("ss_model_loglike", self._h, ptr(p), ctypes.byref(llf))
        return llf.value

    @property
    def scale(self):
        """Concentrated scale at the last likelihood evaluation."""
        s = ctypes.c_double()
        self._call("ss_model_scale", self._h, ctypes.byref(s))
        return s.value

    def representation(self, params):
        """Build the representation at given parameters.

        Parameters
        ----------
        params : array_like, shape (k,)
            Constrained parameters.

        Returns
        -------
        Representation
            A copy, on the data's scale when the scale is concentrated out.
        """
        p = self._params(params)
        h = ctypes.c_void_p()
        self._call("ss_model_rep_at", self._h, ptr(p), ctypes.byref(h))
        return Representation._wrap(h, owner=None)

    def filter(self, params) -> FilterResults:
        """Run the Kalman filter at given parameters.

        Parameters
        ----------
        params : array_like, shape (k,)
            Constrained parameters.

        Returns
        -------
        FilterResults
            The filter output.
        """
        return self.representation(params).filter()

    def components(self, params, variance=False):
        """Compute the smoothed contribution of each component to the data.

        For a component with states :math:`b`, the contribution is
        :math:`Z_{t,b} \\hat\\alpha_{t,b}`; for the irregular it is
        :math:`\\hat\\varepsilon_t`. The contributions add up to the
        observations.

        Parameters
        ----------
        params : array_like, shape (k,)
            Constrained parameters.
        variance : bool, optional
            Also return the variances of the contributions.

        Returns
        -------
        components : dict
            ``{name: ndarray}``, shape (n,) or (n, p). Names are the
            component kinds (``"level"``, ``"seasonal"``, ...), numbered
            when repeated.
        variances : dict
            The same for the variances; only with ``variance=True``.

        Raises
        ------
        TypeError
            If the model is not a :class:`StructuralModel`.
        """
        maxc = 64
        first, size, at = (np.zeros(maxc, dtype=np.int32) for _ in range(3))
        nc = ctypes.c_int()
        self._call(
            "ss_model_components",
            self._h,
            maxc,
            ctypes.byref(nc),
            ptr(first),
            ptr(size),
            ptr(at),
        )
        if nc.value == 0:
            raise TypeError("components are defined for StructuralModel")
        rep = self.representation(params)
        sm = rep.smooth()
        Z = rep["design"]
        n = self.nobs
        if Z.ndim == 2:
            Z = np.repeat(Z[:, :, None], n, axis=2)
        names = _component_names(
            self.components_spec
            if hasattr(self, "components_spec")
            else [None] * nc.value
        )
        out, var = {}, {}
        for i in range(nc.value):
            b = slice(first[i], first[i] + size[i])
            if size[i] == 0:
                sig = sm.smoothed_measurement_disturbance
                v = np.diagonal(
                    sm.smoothed_measurement_disturbance_cov, axis1=0, axis2=1
                ).T
            else:
                Zb = Z[:, b, :]
                sig = np.einsum("pmt,mt->pt", Zb, sm.smoothed_state[b])
                v = np.einsum("pmt,mkt,pkt->pt", Zb, sm.smoothed_state_cov[b, b], Zb)
            sig, v = sig.T, v.T
            if sig.shape[1] == 1:
                sig, v = sig[:, 0], v[:, 0]
            out[names[i]] = sig
            var[names[i]] = v
        return (out, var) if variance else out

    def smooth(self, params) -> SmootherResults:
        """Run the filter and smoother at given parameters.

        Parameters
        ----------
        params : array_like, shape (k,)
            Constrained parameters.

        Returns
        -------
        SmootherResults
            Smoothed states and disturbances.
        """
        return self.representation(params).smooth()

    def fit(
        self,
        start_params=None,
        maxiter=500,
        m=10,
        factr=1e7,
        pgtol=1e-5,
        compute_cov=True,
        gradient="auto",
    ):
        """Estimate the parameters by maximum likelihood (DK §7.3).

        L-BFGS-B minimizes -loglike / n over the unconstrained parameters.
        The defaults are those of statsmodels.

        Parameters
        ----------
        start_params : array_like, optional
            Constrained starting values; default :attr:`start_params`.
        maxiter : int, optional
            Maximum number of iterations.
        m : int, optional
            Number of L-BFGS corrections.
        factr : float, optional
            Stop when the relative reduction of the objective is below
            ``factr`` times the machine epsilon; 10 gives a tight optimum.
        pgtol : float, optional
            Stop when the projected gradient is below `pgtol`.
        compute_cov : bool, optional
            Compute `cov_params` and `bse` from a numerical Hessian.
        gradient : {"auto", "analytic", "numerical"}, optional
            ``"auto"`` uses the analytic score where it applies (DK §7.3.3)
            and central differences elsewhere.

        Returns
        -------
        FitResults
            Estimates and fit statistics.

        Raises
        ------
        StateSpaceError
            If the likelihood cannot be evaluated at the start.

        Notes
        -----
        The score covers parameters in H, R, Q and :math:`P_*` through the
        smoothed disturbances (DK eq. 7.14, 7.16). Parameters that move Z or
        T get central differences, as DK recommend. Standard errors come
        from the Hessian in the unconstrained parameters and the delta
        method.

        Examples
        --------
        A local level with variances 1 (irregular) and 0.25 (level):

        >>> rng = np.random.default_rng(2)
        >>> y = np.cumsum(0.5 * rng.standard_normal(1000)) + rng.standard_normal(1000)
        >>> res = ss.StructuralModel(y, [ss.Irregular(), ss.Level()]).fit()
        >>> res.param_names
        ['sigma2.irregular', 'sigma2.level']
        >>> np.round(res.params, 2)
        array([1.04, 0.22])
        """
        h = ctypes.c_void_p()
        x0 = None if start_params is None else self._params(start_params)
        start = None if x0 is None else ptr(x0)
        code = self._lib_fit(start, maxiter, m, factr, pgtol, compute_cov, gradient, h)
        self._check_callbacks()
        return _fit_results(self, h, code)

    def _lib_fit(self, start, maxiter, m, factr, pgtol, compute_cov, gradient, h):
        return _call_code(
            "ss_model_fit",
            self._h,
            start,
            int(maxiter),
            int(m),
            float(factr),
            float(pgtol),
            int(bool(compute_cov)),
            _GRADIENT[gradient],
            ctypes.byref(h),
        )


class StructuralModel(Model):
    """A model assembled from structural components (DK ch. 3).

    The state vector stacks the components' states in the order given,
    and the system matrices are block diagonal (DK eq. 3.10). Each block is
    initialized as its component needs: diffuse for nonstationary
    components, from the unconditional distribution for stationary ones
    (DK §5.6).

    Parameters
    ----------
    endog : array_like, shape (n,) or (n, p)
        Observations, time first; NaN marks missing values.
    components : list
        :class:`Irregular`, :class:`Level`, :class:`Trend`,
        :class:`Seasonal`, :class:`Cycle`, :class:`Regression`,
        :class:`ARIMA`, :class:`ContinuousLevel` or :class:`ContinuousTrend`
        instances. With p > 1 each applies to every series (seemingly
        unrelated time series equations, DK §3.3), except
        :class:`Regression` and :class:`ARIMA`, which apply to one series.
    loading : array_like, shape (p, p_sig), optional
        Signal loadings :math:`\\Lambda`: components not marked
        ``at_observations`` describe p_sig signals, and
        :math:`Z = \\Lambda Z_{sig}` (DK §3.3.2, §3.7).
    loading_free : array_like of bool, shape (p, p_sig), optional
        Entries of `loading` to estimate; they become the last parameters.

    Notes
    -----
    Parameters are ordered component by component. Variances are mapped to
    the optimizer's scale by :math:`\\sigma^2 = \\exp(2\\psi)` (DK §7.3.2);
    full covariances through a Cholesky factor with log diagonal.

    Examples
    --------
    A basic structural model: level, trigonometric seasonal and irregular
    (DK §8.2):

    >>> rng = np.random.default_rng(2)
    >>> season = np.tile([1.0, -0.5, 0.3, -0.8], 30)
    >>> noise = 0.3 * rng.standard_normal(120)
    >>> y = np.cumsum(0.2 * rng.standard_normal(120)) + season + noise
    >>> comps = [ss.Irregular(), ss.Level(), ss.Seasonal(4, "trig")]
    >>> mod = ss.StructuralModel(y, comps)
    >>> mod.param_names
    ['sigma2.irregular', 'sigma2.level', 'sigma2.seasonal']
    >>> res = mod.fit()
    >>> sorted(res.components())
    ['irregular', 'level', 'seasonal']
    """

    def __init__(self, endog, components, loading=None, loading_free=None):
        y = np.asarray(endog, dtype=np.float64)
        if y.ndim == 1:
            y = y[:, None]
        n, p = y.shape
        yf = farray(y.T)
        b = ctypes.c_void_p()
        call("ss_struct_new", p, n, ptr(yf), ctypes.byref(b))
        try:
            keep = [yf]
            for c in components:
                self._add(b, c, n, keep)
            mh = ctypes.c_void_p()
            if loading is None:
                call("ss_struct_build", b, 0, None, None, ctypes.byref(mh))
            else:
                L = farray(loading)
                if L.ndim != 2 or L.shape[0] != p:
                    raise ValueError(f"loading must have shape ({p}, p_sig)")
                Lf = None
                if loading_free is not None:
                    Lf = np.asfortranarray(np.asarray(loading_free, dtype=np.int32))
                call(
                    "ss_struct_build",
                    b,
                    L.shape[1],
                    ptr(L),
                    None if Lf is None else ptr(Lf),
                    ctypes.byref(mh),
                )
        finally:
            call("ss_struct_free", b)
        self.components_spec = list(components)
        self._attach(mh)

    @staticmethod
    def _add(b, c, n, keep):
        at = int(bool(getattr(c, "at_observations", False)))
        if isinstance(c, Irregular):
            w = None
            if c.weights is not None:
                w = farray(c.weights, (n,))
                keep.append(w)
            call(
                "ss_struct_add_irregular", b, _cov(c.cov), None if w is None else ptr(w)
            )
        elif isinstance(c, Level):
            call("ss_struct_add_level", b, _cov(c.cov), at)
        elif isinstance(c, Trend):
            call("ss_struct_add_trend", b, _cov(c.level_cov), _cov(c.slope_cov), at)
        elif isinstance(c, Seasonal):
            call(
                "ss_struct_add_seasonal",
                b,
                int(c.period),
                _SEASONAL[c.form.lower()],
                _cov(c.cov),
                at,
            )
        elif isinstance(c, Cycle):
            lo, hi = c.period_bounds
            call(
                "ss_struct_add_cycle",
                b,
                _cov(c.cov),
                int(bool(c.damped)),
                float(lo),
                0.0 if hi is None else float(hi),
                at,
            )
        elif isinstance(c, Regression):
            x = np.asarray(c.exog, dtype=np.float64)
            if x.ndim == 1:
                x = x[:, None]
            if x.shape[0] != n:
                raise ValueError(f"exog must have {n} rows")
            x = np.asfortranarray(x)
            rw = np.broadcast_to(
                np.asarray(c.random_walk, dtype=np.int32), (x.shape[1],)
            )
            rw = np.ascontiguousarray(rw)
            keep += [x, rw]
            call(
                "ss_struct_add_regression",
                b,
                x.shape[1],
                ptr(x),
                ptr(rw),
                c.series + 1,
                at,
            )
        elif isinstance(c, ARIMA):
            p_, d, q = c.order
            P, D, Q, s = c.seasonal_order
            call(
                "ss_struct_add_arima",
                b,
                p_,
                d,
                q,
                P,
                D,
                Q,
                s,
                c.series + 1,
                int(c.enforce_stationarity),
                int(c.enforce_invertibility),
                at,
            )
        elif isinstance(c, (ContinuousLevel, ContinuousTrend)):
            t = farray(c.times, (n,))
            keep.append(t)
            call(
                "ss_struct_add_continuous",
                b,
                1 if isinstance(c, ContinuousLevel) else 2,
                ptr(t),
                at,
            )
        else:
            raise TypeError(f"unknown component {c!r}")


def _ssm_method(name):
    """A model method that calls the method of the same name of its
    representation, with that method's signature and documentation."""
    target = getattr(Representation, name)

    @functools.wraps(target)
    def method(self, *args, **kwargs):
        return getattr(self.ssm, name)(*args, **kwargs)

    return method


class _SystemAccess:
    """System matrices and initialization, set on the model itself."""

    def __setitem__(self, name, value):
        """Set a system matrix of the model's representation."""
        self.ssm[name] = value

    def __getitem__(self, name):
        """Copy of a system matrix of the model's representation."""
        return self.ssm[name]

    initialize_known = _ssm_method("initialize_known")
    initialize_diffuse = _ssm_method("initialize_diffuse")
    initialize_approximate_diffuse = _ssm_method("initialize_approximate_diffuse")
    initialize_stationary = _ssm_method("initialize_stationary")
    initialize = _ssm_method("initialize")
    initialize_block = _ssm_method("initialize_block")


class MappedModel(_SystemAccess, Model):
    """A model declared by the matrix entries each parameter sets.

    The declaration is held by the library, so no Python runs during
    estimation and :func:`fit_many` can fit the model in parallel.

    Parameters
    ----------
    endog : array_like, shape (n,) or (n, p), or Representation
        Observations, time first; NaN marks missing values. A
        :class:`Representation` is copied instead, with its matrices and
        initialization; pass `k_params` and the rest by keyword.
    k_states : int
        Number of states m; not used with a Representation.
    k_params : int
        Number of parameters.
    k_posdef : int, optional
        Number of state disturbances r; default `k_states`.
    param_names : list of str, optional
        Parameter names; default ``param1``, ``param2``, ...
    start_params : array_like, optional
        Constrained starting values; default 0.1 each.

    See Also
    --------
    StructuralModel : Models assembled from built-in components.
    MLEModel : Models whose matrices are set by Python code.

    Notes
    -----
    Set the fixed matrices with ``mod[name] = value`` and the
    initialization with the ``initialize_*`` methods, as on a
    :class:`Representation`; then declare the map with :meth:`map` (single
    entries), :meth:`cov` (covariance blocks) and :meth:`constrain`
    (transforms). Parameters in no transform group are unconstrained. Mapped
    entries are overwritten at every evaluation. A map to a single period of
    a time-varying matrix needs that matrix set first.

    Examples
    --------
    The local level model:

    >>> y = np.cumsum(np.random.default_rng(3).standard_normal(100))
    >>> mod = ss.MappedModel(y, k_states=1, k_params=2,
    ...                      param_names=["sigma2.irregular", "sigma2.level"],
    ...                      start_params=[0.5, 0.5])
    >>> mod["design"] = mod["transition"] = mod["selection"] = [[1.0]]
    >>> mod.initialize_diffuse()
    >>> _ = mod.map(0, "obs_cov", 0, 0).map(1, "state_cov", 0, 0)
    >>> _ = mod.constrain([0, 1], "positive")
    >>> res = mod.fit()
    >>> res.converged
    True
    """

    _GROUPS = {"positive": 1, "interval": 2, "stationary": 3, "invertible": 4}

    def __init__(
        self,
        endog,
        k_states=None,
        k_params=None,
        k_posdef=None,
        param_names=None,
        start_params=None,
    ):
        if isinstance(endog, Representation):
            rep = endog
            if k_states is not None:
                # The 0.1 signature, (representation, k_params, param_names,
                # start_params), called positionally.
                warnings.warn(
                    "MappedModel(representation, k_params, ...) with positional "
                    "arguments is deprecated; pass k_params, param_names and "
                    "start_params by keyword",
                    DeprecationWarning,
                    stacklevel=2,
                )
                k_states, k_params, k_posdef, param_names, start_params = (
                    None,
                    k_states,
                    None,
                    param_names if k_params is None else k_params,
                    start_params if k_posdef is None else k_posdef,
                )
        else:
            if k_states is None:
                raise TypeError("MappedModel needs k_states")
            rep = Representation(endog, k_states, k_posdef)
        if k_params is None:
            raise TypeError("MappedModel needs k_params")
        h = ctypes.c_void_p()
        call("ss_mapped_new", rep._h, int(k_params), ctypes.byref(h))
        self._attach(h)
        if param_names is not None:
            if len(param_names) != self.k_params:
                raise ValueError("need one name per parameter")
            width = 32
            raw = "".join(f"{n[:width]:<{width}}" for n in param_names).encode()
            call("ss_model_set_names", self._h, width, raw)
        if start_params is not None:
            self.start_params = start_params

    @Model.start_params.setter
    def start_params(self, values):
        """Set the starting values (constrained)."""
        v = self._params(values)
        call("ss_mapped_start", self._h, ptr(v))

    def map(self, param, matrix, i, j=0, t=None, coef=1.0):
        """Let a parameter set a matrix entry.

        ``matrix[i, j, t] = coef * params[param]``, all indices 0-based.

        Parameters
        ----------
        param : int
            Parameter index.
        matrix : str
            ``design``, ``obs_cov``, ``transition``, ``selection``,
            ``state_cov``, ``state_intercept`` or ``obs_intercept``.
        i, j : int
            Row and column; `j` is ignored for the intercepts.
        t : int, optional
            Period; ``None`` sets every period of a time-varying matrix.
        coef : float, optional
            Multiplier.

        Returns
        -------
        MappedModel
            This model, to chain calls.
        """
        if matrix not in _ARR or matrix == "endog":
            raise KeyError(matrix)
        call(
            "ss_mapped_entry",
            self._h,
            int(param) + 1,
            _ARR[matrix],
            int(i) + 1,
            int(j) + 1,
            0 if t is None else int(t) + 1,
            float(coef),
        )
        return self

    def cov(self, first_param, matrix, offset, dim):
        """Let parameters set a covariance block.

        Parameters ``first_param, ..., first_param + dim (dim + 1) / 2 - 1``
        are the lower triangle, by columns, of the block. They are estimated
        through a Cholesky factor with log diagonal, which keeps the block
        positive definite.

        Parameters
        ----------
        first_param : int
            First parameter index.
        matrix : {"obs_cov", "state_cov"}
            The matrix.
        offset : int
            First row and column of the block (0-based).
        dim : int
            Size of the block.

        Returns
        -------
        MappedModel
            This model, to chain calls.
        """
        call(
            "ss_mapped_block",
            self._h,
            _ARR[matrix],
            int(offset) + 1,
            int(dim),
            int(first_param) + 1,
        )
        return self

    def constrain(self, params, kind, bounds=(0.0, 0.0)):
        """Set the transform of consecutive parameters.

        Parameters
        ----------
        params : int or sequence of int
            Consecutive parameter indices.
        kind : {"positive", "interval", "stationary", "invertible"}
            ``positive``: :math:`\\exp(2x)`, for variances (DK §7.3.2).
            ``interval``: :math:`m + h x / \\sqrt{1 + x^2}` within `bounds`.
            ``stationary``: AR coefficients of a stationary polynomial
            (Monahan 1984). ``invertible``: MA coefficients of an invertible
            polynomial.
        bounds : tuple of float, optional
            (lower, upper) for ``interval``.

        Returns
        -------
        MappedModel
            This model, to chain calls.

        Raises
        ------
        ValueError
            If the parameters are not consecutive.
        """
        params = sorted(int(k) for k in np.atleast_1d(params))
        if params != list(range(params[0], params[-1] + 1)):
            raise ValueError("constrain needs consecutive parameters")
        lo, hi = bounds
        call(
            "ss_mapped_group",
            self._h,
            self._GROUPS[kind],
            params[0] + 1,
            params[-1] + 1,
            float(lo),
            float(hi),
        )
        return self


class MLEModel(_SystemAccess, Model):
    """A model whose system matrices are set by Python code.

    Subclass it as statsmodels' ``MLEModel``: set the fixed matrices
    (``self[name] = value``) and the initialization (``self.initialize_*``)
    in ``__init__``, and override

    * ``update(params)``, which sets the matrices that depend on the
      constrained parameters with ``self[name] = value``;
    * ``start_params``, a property;
    * optionally ``transform_params`` and ``untransform_params`` (the
      default is no transform) and ``param_names``.

    Parameters
    ----------
    endog : array_like, shape (n,) or (n, p)
        Observations, time first; NaN marks missing values.
    k_states : int
        Number of states m.
    k_params : int
        Number of parameters.
    k_posdef : int, optional
        Number of state disturbances r; default `k_states`.

    See Also
    --------
    MappedModel : Declared models, without Python in the loop.

    Notes
    -----
    The library calls ``update`` at every likelihood evaluation, and the
    analytic score calls it 2k more times per gradient for its matrix
    derivatives. Python therefore runs inside the optimization; declared
    and built-in models are faster, and only they work with
    :func:`fit_many`. An exception in ``update`` is raised again when the
    library returns, and the model stays usable.

    Examples
    --------
    >>> class LocalLevel(ss.MLEModel):
    ...     def __init__(self, y):
    ...         super().__init__(y, k_states=1, k_params=2)
    ...         self["design"] = self["transition"] = self["selection"] = [[1.0]]
    ...         self.initialize_diffuse()
    ...     @property
    ...     def start_params(self):
    ...         return np.array([1.0, 1.0])
    ...     def transform_params(self, x):
    ...         return np.exp(2 * np.asarray(x))
    ...     def untransform_params(self, p):
    ...         return 0.5 * np.log(p)
    ...     def update(self, params):
    ...         self["obs_cov"] = [[params[0]]]
    ...         self["state_cov"] = [[params[1]]]
    >>> y = np.cumsum(np.random.default_rng(4).standard_normal(100))
    >>> LocalLevel(y).fit().converged
    True
    """

    _python_callbacks = True

    def __init__(self, endog, k_states, k_params, k_posdef=None):
        template = Representation(endog, k_states, k_posdef)
        self._cb_error = None
        self._update_cb = UPDATE_CB(self._update_trampoline)
        cls = type(self)
        tr = cls.transform_params is not Model.transform_params
        ut = cls.untransform_params is not Model.untransform_params
        self._transform_cb = TRANSFORM_CB(self._transform_trampoline) if tr else None
        self._untransform_cb = (
            TRANSFORM_CB(self._untransform_trampoline) if ut else None
        )
        h = ctypes.c_void_p()
        call(
            "ss_callback_new",
            template._h,
            int(k_params),
            ctypes.cast(self._update_cb, ctypes.c_void_p),
            None if not tr else ctypes.cast(self._transform_cb, ctypes.c_void_p),
            None if not ut else ctypes.cast(self._untransform_cb, ctypes.c_void_p),
            ctypes.byref(h),
        )
        self._attach(h)
        self._ssm = Model.ssm.fget(self)

    @property
    def ssm(self):
        """The model's representation, which ``update`` fills."""
        return self._ssm

    def update(self, params):
        """Set the system matrices from the parameters; override this.

        Parameters
        ----------
        params : ndarray, shape (k,)
            Constrained parameters.
        """
        raise NotImplementedError("MLEModel subclasses must implement update(params)")

    @property
    def param_names(self):
        """Parameter names; override to name them."""
        return [f"param{i}" for i in range(self.k_params)]

    def fit(self, start_params=None, **kwargs):
        """Estimate the parameters by maximum likelihood.

        Parameters
        ----------
        start_params : array_like, optional
            Constrained starting values; default :attr:`start_params`.
        **kwargs
            As for :meth:`Model.fit`.

        Returns
        -------
        FitResults
            Estimates and fit statistics.
        """
        if start_params is None:
            start_params = self.start_params
        return super().fit(start_params=start_params, **kwargs)

    # -- callbacks ---------------------------------------------------------

    def _check_callbacks(self):
        err, self._cb_error = self._cb_error, None
        if err is not None:
            raise err

    def _run(self, fn):
        if self._cb_error is not None:
            return 1
        try:
            fn()
            return 0
        except BaseException as e:  # reported when control returns to Python
            self._cb_error = e
            return 1

    def _update_trampoline(self, k, params, rep):
        p = np.ctypeslib.as_array(params, (k,)).copy()
        return self._run(lambda: self.update(p))

    def _transform_trampoline(self, k, x, out):
        def f():
            np.ctypeslib.as_array(out, (k,))[:] = self.transform_params(
                np.ctypeslib.as_array(x, (k,)).copy()
            )

        return self._run(f)

    def _untransform_trampoline(self, k, x, out):
        def f():
            np.ctypeslib.as_array(out, (k,))[:] = self.untransform_params(
                np.ctypeslib.as_array(x, (k,)).copy()
            )

        return self._run(f)


_KIND = {
    "Irregular": "irregular",
    "Level": "level",
    "Trend": "trend",
    "Seasonal": "seasonal",
    "Cycle": "cycle",
    "Regression": "regression",
    "ARIMA": "arima",
    "ContinuousLevel": "level",
    "ContinuousTrend": "trend",
}


def _component_names(specs):
    names, seen = [], {}
    for c in specs:
        base = _KIND.get(type(c).__name__, "component")
        seen[base] = seen.get(base, 0) + 1
        names.append(base if seen[base] == 1 else f"{base}_{seen[base]}")
    return names


# ------------------------------------------------------------------- fitting


def fit_many(
    models, maxiter=500, m=10, factr=1e7, pgtol=1e-5, compute_cov=True, gradient="auto"
):
    """Fit independent models in parallel.

    The library fits the models on OpenMP threads, with the GIL released.

    Parameters
    ----------
    models : sequence of StructuralModel or MappedModel
        The models; each is left at its estimates.
    maxiter : int, optional
        Maximum number of iterations.
    m : int, optional
        Number of L-BFGS corrections.
    factr : float, optional
        Relative reduction tolerance, in units of the machine epsilon.
    pgtol : float, optional
        Projected gradient tolerance.
    compute_cov : bool, optional
        Compute standard errors.
    gradient : {"auto", "analytic", "numerical"}, optional
        As for :meth:`Model.fit`.

    Returns
    -------
    list of FitResults or None
        One result per model; None where the likelihood could not be
        evaluated.

    Raises
    ------
    TypeError
        If a model is an :class:`MLEModel`: Python callbacks cannot run on
        the library's threads.

    Notes
    -----
    The library must be built with OpenMP (the default with CMake). For
    large models set ``OPENBLAS_NUM_THREADS=1`` to avoid nested threading
    in BLAS.

    Examples
    --------
    >>> rng = np.random.default_rng(5)
    >>> ys = [np.cumsum(rng.standard_normal(80)) for _ in range(4)]
    >>> results = ss.fit_many([ss.StructuralModel(y, [ss.Irregular(), ss.Level()])
    ...                        for y in ys])
    >>> [r.converged for r in results]
    [True, True, True, True]
    """
    models = list(models)
    for mod in models:
        if not isinstance(mod, Model) or getattr(mod, "_python_callbacks", False):
            raise TypeError("fit_many needs models without Python callbacks")
    nm = len(models)
    hs = (ctypes.c_void_p * nm)(*[mod._h.value for mod in models])
    fhs = (ctypes.c_void_p * nm)()
    infos = (ctypes.c_int * nm)()
    call(
        "ss_fit_many",
        nm,
        hs,
        int(maxiter),
        int(m),
        float(factr),
        float(pgtol),
        int(bool(compute_cov)),
        _GRADIENT[gradient],
        fhs,
        infos,
    )
    out = []
    for mod, fh, code in zip(models, fhs, infos):
        out.append(None if not fh else _fit_results(mod, ctypes.c_void_p(fh), code))
    return out


def _call_code(name, *args):
    from ._lib import lib

    return getattr(lib, name)(*args)


def _fit_results(model, h, code):
    if not h:
        raise StateSpaceError(code, "fit")
    try:
        k = model.k_params
        d = [ctypes.c_double() for _ in range(4)]
        ints = [ctypes.c_int() for _ in range(5)]
        call("ss_fit_scalars", h, *(ctypes.byref(x) for x in d + ints))
        llf, scale, aic, bic = (x.value for x in d)
        niter, nfev, conv, ag, has_cov = (x.value for x in ints)
        params = np.empty(k)
        call("ss_fit_get", h, 1, ptr(params))
        bse = np.full(k, np.nan)
        cov = np.full((k, k), np.nan)
        if has_cov:
            call("ss_fit_get", h, 2, ptr(bse))
            cov = np.empty((k, k), order="F")
            call("ss_fit_get", h, 3, ptr(cov))
        buf = ctypes.create_string_buffer(64)
        call("ss_fit_message", h, len(buf), buf)
    finally:
        call("ss_fit_free", h)
    return FitResults(
        model=model,
        params=params,
        param_names=model.param_names,
        bse=bse,
        cov_params=cov,
        llf=llf,
        scale=scale,
        aic=aic,
        bic=bic,
        niter=niter,
        nfev=nfev,
        converged=bool(conv),
        analytic_gradient=bool(ag),
        message=buf.value.decode().strip(),
    )


def _free(routine, handle):
    if handle:
        call(routine, handle)
