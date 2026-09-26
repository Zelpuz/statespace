"""Models with parameters, and maximum likelihood estimation.

`StructuralModel` assembles the built-in components (Durbin and Koopman
2012, ch. 3): irregular, level, trend, seasonal, cycle, regression, ARIMA
and continuous-time components, for one series or several (SUTSE), with
optional signal loadings. Everything per likelihood evaluation runs in
Fortran; `fit_many` fits many models in parallel.
"""

import ctypes
import weakref
from dataclasses import dataclass, field

import numpy as np

from ._lib import TRANSFORM_CB, UPDATE_CB, StateSpaceError, call, farray, ptr
from .representation import _ARR, Representation, SmootherResults, FilterResults

_COV = {None: 0, "none": 0, "diagonal": 1, "full": 2}
_SEASONAL = {"dummy": 1, "trig": 2, "trigonometric": 2, "harrison-stevens": 3, "hs": 3}
GRADIENT_AUTO, GRADIENT_NUMERICAL, GRADIENT_ANALYTIC = 0, 1, 2
_GRADIENT = {"auto": GRADIENT_AUTO, "numerical": GRADIENT_NUMERICAL,
             "analytic": GRADIENT_ANALYTIC}


def _cov(c):
    try:
        return _COV[c if c is None else c.lower()]
    except KeyError:
        raise ValueError(f"covariance must be one of None, 'diagonal', 'full', not {c!r}")


# ---------------------------------------------------------------- components

@dataclass
class Irregular:
    """eps_t ~ N(0, Sigma_eps); `weights` (n) scale it per period (DK 3.8)."""
    cov: str = "diagonal"
    weights: object = None


@dataclass
class Level:
    """Random walk level (DK 2, 3.3)."""
    cov: str = "diagonal"
    at_observations: bool = False


@dataclass
class Trend:
    """Local linear trend (DK 3.2); level_cov=None gives the smooth trend."""
    level_cov: str = "diagonal"
    slope_cov: str = "diagonal"
    at_observations: bool = False


@dataclass
class Seasonal:
    """Seasonal of `period`: 'dummy', 'trig' or 'harrison-stevens' (DK 3.2.2);
    cov=None for a fixed seasonal."""
    period: int
    form: str = "dummy"
    cov: str = "diagonal"
    at_observations: bool = False


@dataclass
class Cycle:
    """Cycle (DK 3.2.4), damped (stationary) by default; the period is kept
    within `period_bounds` (in observations; upper None = the sample size)."""
    cov: str = "diagonal"
    damped: bool = True
    period_bounds: tuple = (2.0, None)
    at_observations: bool = False


@dataclass
class Regression:
    """Regression effects x_t' beta_t on `series` (0-based) (DK 3.2.5, 3.6);
    `random_walk` (bool or one per column) makes coefficients time-varying."""
    exog: object
    random_walk: object = False
    series: int = 0
    at_observations: bool = False


@dataclass
class ARIMA:
    """ARIMA(p, d, q)(P, D, Q)_s (DK 3.4), same state layout as SARIMAX."""
    order: tuple = (1, 0, 0)
    seasonal_order: tuple = (0, 0, 0, 0)
    series: int = 0
    enforce_stationarity: bool = True
    enforce_invertibility: bool = True
    at_observations: bool = False


@dataclass
class ContinuousLevel:
    """Continuous-time level observed at `times` (DK 3.8.1)."""
    times: object
    at_observations: bool = False


@dataclass
class ContinuousTrend:
    """Continuous-time smooth trend at `times` (DK 3.8.2); with an irregular
    this is the cubic smoothing spline (DK 3.9.2)."""
    times: object
    at_observations: bool = False


# ------------------------------------------------------------------ results

@dataclass
class FitResults:
    """Maximum likelihood estimates (DK 7.3)."""
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
        """Filter and smoother at the estimates."""
        return self.model.smooth(self.params)

    def filter(self):
        return self.model.filter(self.params)

    def summary(self):
        w = max(len(n) for n in self.param_names)
        lines = [f"log likelihood {self.llf:.4f}   AIC {self.aic:.4f}   BIC {self.bic:.4f}",
                 f"{'':{w}}  {'estimate':>14}  {'std err':>12}"]
        for n, p, s in zip(self.param_names, self.params, self.bse):
            lines.append(f"{n:{w}}  {p:14.6g}  {s:12.4g}")
        lines.append(self.message)
        return "\n".join(lines)


# -------------------------------------------------------------------- models

class Model:
    """A model with parameters, held by the Fortran library (base class)."""

    def _attach(self, handle):
        self._h = handle
        self._finalizer = weakref.finalize(self, _free, "ss_model_free", handle)
        k, p, m, r, n = (ctypes.c_int() for _ in range(5))
        call("ss_model_info", handle, *(ctypes.byref(x) for x in (k, p, m, r, n)))
        self.k_params, self.k_endog, self.k_states = k.value, p.value, m.value
        self.k_posdef, self.nobs = r.value, n.value

    @property
    def ssm(self):
        """The model's own representation (borrowed), e.g. to set
        filter_method or tol_steady before fitting."""
        h = ctypes.c_void_p()
        self._call("ss_model_rep", self._h, ctypes.byref(h))
        return Representation._wrap(h, owner=self)

    @property
    def param_names(self):
        width = 32
        buf = ctypes.create_string_buffer(width * self.k_params)
        self._call("ss_model_param_names", self._h, width, buf)
        raw = buf.raw.decode()
        return [raw[i * width:(i + 1) * width].strip() for i in range(self.k_params)]

    @property
    def start_params(self):
        out = np.empty(self.k_params)
        self._call("ss_model_start_params", self._h, ptr(out))
        return out

    @property
    def concentrate_scale(self):
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
        x = self._params(unconstrained)
        out = np.empty(self.k_params)
        self._call("ss_model_transform", self._h, ptr(x), ptr(out))
        return out

    def untransform_params(self, constrained):
        p = self._params(constrained)
        out = np.empty(self.k_params)
        self._call("ss_model_untransform", self._h, ptr(p), ptr(out))
        return out

    def loglike(self, params):
        """Log likelihood at constrained params (concentrated if
        concentrate_scale; the scale is then in `self.scale`)."""
        p = self._params(params)
        llf = ctypes.c_double()
        self._call("ss_model_loglike", self._h, ptr(p), ctypes.byref(llf))
        return llf.value

    @property
    def scale(self):
        s = ctypes.c_double()
        self._call("ss_model_scale", self._h, ctypes.byref(s))
        return s.value

    def representation(self, params):
        """A new Representation at params (on the data's scale)."""
        p = self._params(params)
        h = ctypes.c_void_p()
        self._call("ss_model_rep_at", self._h, ptr(p), ctypes.byref(h))
        return Representation._wrap(h, owner=None)

    def filter(self, params) -> FilterResults:
        return self.representation(params).filter()

    def smooth(self, params) -> SmootherResults:
        return self.representation(params).smooth()

    def fit(self, start_params=None, maxiter=500, m=10, factr=1e7, pgtol=1e-5,
            compute_cov=True, gradient="auto"):
        """Maximum likelihood with L-BFGS-B (statsmodels' defaults). The
        gradient is DK's analytic score where it applies (7.3.3) and central
        differences for parameters in Z or T."""
        h = ctypes.c_void_p()
        x0 = None if start_params is None else self._params(start_params)
        start = None if x0 is None else ptr(x0)
        code = self._lib_fit(start, maxiter, m, factr, pgtol, compute_cov, gradient, h)
        self._check_callbacks()
        return _fit_results(self, h, code)

    def _lib_fit(self, start, maxiter, m, factr, pgtol, compute_cov, gradient, h):
        return _call_code("ss_model_fit", self._h, start, int(maxiter), int(m), float(factr),
                          float(pgtol), int(bool(compute_cov)), _GRADIENT[gradient],
                          ctypes.byref(h))


class StructuralModel(Model):
    """Components assembled into one model (DK ch. 3).

    Parameters
    ----------
    endog : array_like, (n,) or (n, p)
        Observations, time first; NaN for missing values.
    components : list
        Irregular, Level, Trend, Seasonal, Cycle, Regression, ARIMA,
        ContinuousLevel, ContinuousTrend, in the order of the state vector.
    loading, loading_free : array_like (p, p_sig), optional
        Signal loadings Z = Lambda Z_sig (DK 3.3.2, 3.7) and which of their
        entries are estimated. Components not marked `at_observations`
        then describe the p_sig signals.
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
                call("ss_struct_build", b, L.shape[1], ptr(L), None if Lf is None else ptr(Lf),
                     ctypes.byref(mh))
        finally:
            call("ss_struct_free", b)
        self.components = list(components)
        self._attach(mh)

    @staticmethod
    def _add(b, c, n, keep):
        at = int(bool(getattr(c, "at_observations", False)))
        if isinstance(c, Irregular):
            w = None
            if c.weights is not None:
                w = farray(c.weights, (n,))
                keep.append(w)
            call("ss_struct_add_irregular", b, _cov(c.cov), None if w is None else ptr(w))
        elif isinstance(c, Level):
            call("ss_struct_add_level", b, _cov(c.cov), at)
        elif isinstance(c, Trend):
            call("ss_struct_add_trend", b, _cov(c.level_cov), _cov(c.slope_cov), at)
        elif isinstance(c, Seasonal):
            call("ss_struct_add_seasonal", b, int(c.period), _SEASONAL[c.form.lower()],
                 _cov(c.cov), at)
        elif isinstance(c, Cycle):
            lo, hi = c.period_bounds
            call("ss_struct_add_cycle", b, _cov(c.cov), int(bool(c.damped)), float(lo),
                 0.0 if hi is None else float(hi), at)
        elif isinstance(c, Regression):
            x = np.asarray(c.exog, dtype=np.float64)
            if x.ndim == 1:
                x = x[:, None]
            if x.shape[0] != n:
                raise ValueError(f"exog must have {n} rows")
            x = np.asfortranarray(x)
            rw = np.broadcast_to(np.asarray(c.random_walk, dtype=np.int32), (x.shape[1],))
            rw = np.ascontiguousarray(rw)
            keep += [x, rw]
            call("ss_struct_add_regression", b, x.shape[1], ptr(x), ptr(rw), c.series + 1, at)
        elif isinstance(c, ARIMA):
            p_, d, q = c.order
            P, D, Q, s = c.seasonal_order
            call("ss_struct_add_arima", b, p_, d, q, P, D, Q, s, c.series + 1,
                 int(c.enforce_stationarity), int(c.enforce_invertibility), at)
        elif isinstance(c, (ContinuousLevel, ContinuousTrend)):
            t = farray(c.times, (n,))
            keep.append(t)
            call("ss_struct_add_continuous", b, 1 if isinstance(c, ContinuousLevel) else 2,
                 ptr(t), at)
        else:
            raise TypeError(f"unknown component {c!r}")


class MappedModel(Model):
    """A model given by which matrix entries each parameter sets; no Python
    runs during estimation, and `fit_many` can fit it in parallel.

    Parameters
    ----------
    representation : Representation
        The fixed entries and the initialization (copied).
    k_params : int
    param_names : list of str, optional
    start_params : array_like, optional
        Constrained starting values (default 0.1).

    Declare the map with `map`, `cov` and `constrain`, e.g. for the local
    level model::

        mod = MappedModel(rep, 2, ["sigma2.irregular", "sigma2.level"])
        mod.map(0, "obs_cov", 0, 0)
        mod.map(1, "state_cov", 0, 0)
        mod.constrain([0, 1], "positive")
    """

    _GROUPS = {"positive": 1, "interval": 2, "stationary": 3, "invertible": 4}

    def __init__(self, representation, k_params, param_names=None, start_params=None):
        h = ctypes.c_void_p()
        call("ss_mapped_new", representation._h, int(k_params), ctypes.byref(h))
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
        v = self._params(values)
        call("ss_mapped_start", self._h, ptr(v))

    def map(self, param, matrix, i, j=0, t=None, coef=1.0):
        """Parameter `param` sets matrix[i, j, t] = coef * param (0-based;
        t=None for every period). matrix is one of design, obs_cov,
        transition, selection, state_cov, state_intercept, obs_intercept."""
        if matrix not in _ARR or matrix == "endog":
            raise KeyError(matrix)
        call("ss_mapped_entry", self._h, int(param) + 1, _ARR[matrix], int(i) + 1, int(j) + 1,
             0 if t is None else int(t) + 1, float(coef))
        return self

    def cov(self, first_param, matrix, offset, dim):
        """Parameters first_param.. are the lower triangle (by columns) of the
        dim x dim covariance block of obs_cov or state_cov at offset
        (0-based); estimated through its Cholesky factor."""
        call("ss_mapped_block", self._h, _ARR[matrix], int(offset) + 1, int(dim),
             int(first_param) + 1)
        return self

    def constrain(self, params, kind, bounds=(0.0, 0.0)):
        """Transform for consecutive parameters: 'positive' (variances,
        exp(2x)), 'interval' (within bounds), 'stationary' (AR
        coefficients), 'invertible' (MA coefficients)."""
        params = sorted(int(k) for k in np.atleast_1d(params))
        if params != list(range(params[0], params[-1] + 1)):
            raise ValueError("constrain needs consecutive parameters")
        lo, hi = bounds
        call("ss_mapped_group", self._h, self._GROUPS[kind], params[0] + 1, params[-1] + 1,
             float(lo), float(hi))
        return self


class MLEModel(Model):
    """A model written in Python, as with statsmodels' MLEModel: set the
    fixed matrices and the initialization in __init__ (``self["design"] =
    ...``, ``self.ssm.initialize_diffuse()``), and override

    * ``update(params)``: set the matrices that depend on the constrained
      params, with ``self[name] = value``;
    * ``start_params`` (property);
    * optionally ``transform_params``/``untransform_params`` (unconstrained
      <-> constrained) and ``param_names``.

    The library calls ``update`` for every likelihood evaluation (and the
    analytic score calls it 2k times per gradient for the matrix
    derivatives), so Python runs inside the optimization; for speed, prefer
    MappedModel or the built-in models. Cannot be used with fit_many.
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
        self._untransform_cb = TRANSFORM_CB(self._untransform_trampoline) if ut else None
        h = ctypes.c_void_p()
        call("ss_callback_new", template._h, int(k_params),
             ctypes.cast(self._update_cb, ctypes.c_void_p),
             None if not tr else ctypes.cast(self._transform_cb, ctypes.c_void_p),
             None if not ut else ctypes.cast(self._untransform_cb, ctypes.c_void_p),
             ctypes.byref(h))
        self._attach(h)
        self._ssm = Model.ssm.fget(self)

    @property
    def ssm(self):
        return self._ssm

    def __setitem__(self, name, value):
        self._ssm[name] = value

    def __getitem__(self, name):
        return self._ssm[name]

    def update(self, params):
        raise NotImplementedError("MLEModel subclasses must implement update(params)")

    @property
    def param_names(self):
        return [f"param{i}" for i in range(self.k_params)]

    def fit(self, start_params=None, **kwargs):
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
        except BaseException as e:   # reported when control returns to Python
            self._cb_error = e
            return 1

    def _update_trampoline(self, k, params, rep):
        p = np.ctypeslib.as_array(params, (k,)).copy()
        return self._run(lambda: self.update(p))

    def _transform_trampoline(self, k, x, out):
        def f():
            np.ctypeslib.as_array(out, (k,))[:] = self.transform_params(
                np.ctypeslib.as_array(x, (k,)).copy())
        return self._run(f)

    def _untransform_trampoline(self, k, x, out):
        def f():
            np.ctypeslib.as_array(out, (k,))[:] = self.untransform_params(
                np.ctypeslib.as_array(x, (k,)).copy())
        return self._run(f)


# ------------------------------------------------------------------- fitting

def fit_many(models, maxiter=500, m=10, factr=1e7, pgtol=1e-5, compute_cov=True,
             gradient="auto"):
    """Fit independent models in parallel (OpenMP threads in the library;
    the GIL is released). Returns a list of FitResults, with None where a fit
    failed. For large models set OPENBLAS_NUM_THREADS=1 to avoid nested
    threading in BLAS. Models defined by Python callbacks cannot be fitted
    this way; use Model.fit."""
    models = list(models)
    for mod in models:
        if not isinstance(mod, Model) or getattr(mod, "_python_callbacks", False):
            raise TypeError("fit_many needs models without Python callbacks")
    nm = len(models)
    hs = (ctypes.c_void_p * nm)(*[mod._h.value for mod in models])
    fhs = (ctypes.c_void_p * nm)()
    infos = (ctypes.c_int * nm)()
    call("ss_fit_many", nm, hs, int(maxiter), int(m), float(factr), float(pgtol),
         int(bool(compute_cov)), _GRADIENT[gradient], fhs, infos)
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
    return FitResults(model=model, params=params, param_names=model.param_names, bse=bse,
                      cov_params=cov, llf=llf, scale=scale, aic=aic, bic=bic, niter=niter,
                      nfev=nfev, converged=bool(conv), analytic_gradient=bool(ag),
                      message=buf.value.decode().strip())


def _free(routine, handle):
    if handle:
        call(routine, handle)
