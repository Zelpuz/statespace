"""Matrix-level API: a state space representation filled from numpy.

The model is

    y_t = d_t + Z_t alpha_t + eps_t,          eps_t ~ N(0, H_t)
    alpha_t+1 = c_t + T_t alpha_t + R_t eta_t,  eta_t ~ N(0, Q_t)

(Durbin and Koopman 2012, 3.1). The matrices use statsmodels' names and
shapes, time last: ``design`` Z (p, m[, n]), ``obs_cov`` H (p, p[, n]),
``transition`` T (m, m[, n]), ``selection`` R (m, r[, n]), ``state_cov``
Q (r, r[, n]), ``state_intercept`` c (m[, n]), ``obs_intercept`` d (p[, n]).
A matrix without the time dimension (or with a time dimension of 1) is
time-invariant.
"""

import ctypes
import weakref
from dataclasses import dataclass

import numpy as np

from ._lib import call, farray, ptr

# Initialization kinds (statespace_rep)
INIT_KNOWN = 1
INIT_APPROX_DIFFUSE = 2
INIT_STATIONARY = 3
INIT_DIFFUSE = 4
INIT_GENERAL = 5
# Filter and diffuse methods
FILTER_CONVENTIONAL = 0
FILTER_UNIVARIATE = 1
DIFFUSE_UNIVARIATE = 0
DIFFUSE_MULTIVARIATE = 1

_ARR = {"endog": 1, "design": 2, "obs_cov": 3, "transition": 4, "selection": 5,
        "state_cov": 6, "state_intercept": 7, "obs_intercept": 8}
_OPT_INT = {"filter_method": 1, "diffuse_method": 2, "loglikelihood_burn": 3,
            "marginal_likelihood": 4}
_OPT_REAL = {"tol_diffuse": 1, "tol_steady": 2}


@dataclass
class FilterResults:
    """Kalman filter output (DK 4.3, 5.2), in statsmodels' naming; time is the
    last axis. Periods are 0-based: predicted_state[:, t] is a_{t+1}."""
    llf: float
    llf_obs: np.ndarray
    nobs_diffuse: int
    k_diffuse: int
    t_steady: int
    predicted_state: np.ndarray
    predicted_state_cov: np.ndarray
    predicted_diffuse_state_cov: np.ndarray
    filtered_state: np.ndarray
    filtered_state_cov: np.ndarray
    forecasts: np.ndarray
    forecasts_error: np.ndarray
    forecasts_error_cov: np.ndarray
    forecasts_error_diffuse_cov: np.ndarray
    forecasts_error_cov_inv: np.ndarray
    kalman_gain: np.ndarray
    standardized_forecasts_error: np.ndarray
    diagnostic_start: int


@dataclass
class SmootherResults:
    """State and disturbance smoother output (DK 4.4-4.5, 5.3).
    scaled_smoothed_estimator[:, t] is r_t for t = 0..n."""
    filter: FilterResults
    smoothed_state: np.ndarray
    smoothed_state_cov: np.ndarray
    scaled_smoothed_estimator: np.ndarray
    scaled_smoothed_estimator_cov: np.ndarray
    smoothed_measurement_disturbance: np.ndarray
    smoothed_measurement_disturbance_cov: np.ndarray
    smoothed_state_disturbance: np.ndarray
    smoothed_state_disturbance_cov: np.ndarray


@dataclass
class AugmentedResults:
    """Augmented Kalman filter and smoother (DK 5.7): the delta = 0 filter,
    the GLS estimate of the diffuse elements delta and its variance, and
    the diffuse (7.2.2), fixed-but-unknown (7.2.4) and marginal (7.2.6)
    log likelihoods."""
    filter: FilterResults
    delta: np.ndarray
    delta_cov: np.ndarray
    llf: float
    llf_fixed: float
    llf_marginal: float
    smoothed_state: np.ndarray
    smoothed_state_cov: np.ndarray


class Representation:
    """A linear Gaussian state space model held by the Fortran library.

    Parameters
    ----------
    endog : array_like, shape (n,) or (n, p)
        Observations, time first as in statsmodels; NaN marks missing values.
    k_states : int
        Number of states m.
    k_posdef : int, optional
        Number of state disturbances r (default m).
    """

    def __init__(self, endog, k_states, k_posdef=None):
        y = np.asarray(endog, dtype=np.float64)
        if y.ndim == 1:
            y = y[:, None]
        nobs, k_endog = y.shape
        k_posdef = int(k_states) if k_posdef is None else int(k_posdef)
        h = ctypes.c_void_p()
        call("ss_rep_new", k_endog, int(k_states), k_posdef, nobs, ctypes.byref(h))
        self._attach(h, owner=None)
        self["endog"] = y.T

    @classmethod
    def _wrap(cls, handle, owner):
        """A Representation over an existing handle. With owner=None the
        handle is owned (freed with this object); otherwise it is borrowed
        from `owner`, which is kept alive."""
        self = object.__new__(cls)
        self._attach(handle, owner)
        return self

    def _attach(self, handle, owner):
        p, m, r, n = (ctypes.c_int() for _ in range(4))
        nts = np.zeros(8, dtype=np.int32)
        call("ss_rep_info", handle, ctypes.byref(p), ctypes.byref(m), ctypes.byref(r),
             ctypes.byref(n), ptr(nts))
        d = object.__setattr__
        d(self, "_h", handle)
        d(self, "_owner", owner)
        d(self, "k_endog", p.value)
        d(self, "k_states", m.value)
        d(self, "k_posdef", r.value)
        d(self, "nobs", n.value)
        if owner is None:
            d(self, "_finalizer", weakref.finalize(self, _free, "ss_rep_free", handle))

    # -- system matrices -----------------------------------------------

    def _shape(self, name):
        p, m, r = self.k_endog, self.k_states, self.k_posdef
        return {"endog": (p, self.nobs), "design": (p, m), "obs_cov": (p, p),
                "transition": (m, m), "selection": (m, r), "state_cov": (r, r),
                "state_intercept": (m,), "obs_intercept": (p,)}[name]

    def __setitem__(self, name, value):
        if name not in _ARR:
            raise KeyError(name)
        base = self._shape(name)
        a = np.asarray(value, dtype=np.float64)
        if name == "endog":
            if a.shape != base:
                raise ValueError(f"endog must have shape {base} after transposing")
            nt = self.nobs
            a = farray(a)
        elif a.shape == base or a.shape == base + (1,):
            nt = 1
            a = farray(a, base + (1,))
        elif a.shape == base + (self.nobs,):
            nt = self.nobs
            a = farray(a)
        else:
            raise ValueError(f"{name} must have shape {base} or {base + (self.nobs,)}, "
                             f"not {a.shape}")
        call("ss_rep_set", self._h, _ARR[name], nt, ptr(a))

    def __getitem__(self, name):
        """A copy of the array from the library; time-invariant matrices come
        without the time dimension."""
        if name not in _ARR:
            raise KeyError(name)
        p, m, r, n = (ctypes.c_int() for _ in range(4))
        nts = np.zeros(8, dtype=np.int32)
        call("ss_rep_info", self._h, ctypes.byref(p), ctypes.byref(m), ctypes.byref(r),
             ctypes.byref(n), ptr(nts))
        base = self._shape(name)
        nt = int(nts[_ARR[name] - 1])
        shape = base if name == "endog" else base + (nt,)
        out = np.empty(shape, dtype=np.float64, order="F")
        call("ss_rep_get", self._h, _ARR[name], ptr(out))
        return out[..., 0] if name != "endog" and nt == 1 else out

    @property
    def endog(self):
        """Observations, shape (p, n)."""
        return self["endog"]

    # -- options ---------------------------------------------------------

    def __setattr__(self, name, value):
        if name in _OPT_INT:
            call("ss_rep_set_int", self._h, _OPT_INT[name], int(value))
        elif name in _OPT_REAL:
            call("ss_rep_set_real", self._h, _OPT_REAL[name], float(value))
        object.__setattr__(self, name, value)

    # -- initialization (DK 5.1) ---------------------------------------

    def initialize_known(self, a1, P1):
        m = self.k_states
        a = farray(a1, (m,))
        P = farray(P1, (m, m))
        call("ss_rep_init_known", self._h, ptr(a), ptr(P))

    def initialize_diffuse(self):
        """Exact diffuse initialization (DK 5.2)."""
        call("ss_rep_init_diffuse", self._h)

    def initialize_approximate_diffuse(self, variance=1e6):
        call("ss_rep_init_approx_diffuse", self._h, float(variance))

    def initialize_stationary(self):
        call("ss_rep_init_stationary", self._h)

    def initialize(self, a1, Pstar, Pinf):
        """General initialization alpha_1 ~ N(a1, Pstar + kappa Pinf), kappa -> inf."""
        m = self.k_states
        a, Ps, Pi = farray(a1, (m,)), farray(Pstar, (m, m)), farray(Pinf, (m, m))
        call("ss_rep_init_general", self._h, ptr(a), ptr(Ps), ptr(Pi))

    def initialize_block(self, start, stop, kind, a1=None, P1=None, variance=1e6):
        """Initialize states start..stop-1 (0-based, like a slice) with `kind`
        (INIT_*); a1 and P1 are needed for INIT_KNOWN."""
        nb = stop - start
        if a1 is not None:
            a, P = farray(a1, (nb,)), farray(P1, (nb, nb))
            call("ss_rep_init_block", self._h, start + 1, stop, kind, ptr(a), ptr(P),
                 float(variance))
        else:
            call("ss_rep_init_block", self._h, start + 1, stop, kind, None, None,
                 float(variance))

    # -- computations ----------------------------------------------------

    def loglike(self):
        """Log likelihood (diffuse log likelihood with a diffuse initialization)."""
        llf = ctypes.c_double()
        call("ss_loglike", self._h, ctypes.byref(llf))
        return llf.value

    def loglike_concentrated(self):
        """Log likelihood with the scale concentrated out; returns (llf, scale)."""
        llf, scale = ctypes.c_double(), ctypes.c_double()
        call("ss_loglike_concentrated", self._h, ctypes.byref(llf), ctypes.byref(scale))
        return llf.value, scale.value

    def filter(self, method="conventional"):
        """Run the Kalman filter and return FilterResults; method="sqrt"
        uses the square root filter (DK 6.3)."""
        fh = ctypes.c_void_p()
        call("ss_sqrt_filter" if method == "sqrt" else "ss_filter", self._h, ctypes.byref(fh))
        try:
            return self._filter_results(fh)
        finally:
            _free("ss_filter_free", fh)

    def smooth(self, method="conventional"):
        """Run the filter and the state and disturbance smoother;
        method="sqrt" uses the square root filter and smoother (DK 6.3)."""
        fh, sh = ctypes.c_void_p(), ctypes.c_void_p()
        call("ss_sqrt_filter" if method == "sqrt" else "ss_filter", self._h, ctypes.byref(fh))
        try:
            fres = self._filter_results(fh)
            if method == "sqrt":
                call("ss_sqrt_smoother", self._h, fh, ctypes.byref(sh))
            else:
                call("ss_smooth", self._h, fh, ctypes.byref(sh))
            try:
                p, m, r, n = self.k_endog, self.k_states, self.k_posdef, self.nobs
                get = lambda code, shape: _get("ss_smoother_get", sh, code, shape)
                return SmootherResults(
                    filter=fres,
                    smoothed_state=get(1, (m, n)),
                    smoothed_state_cov=get(2, (m, m, n)),
                    scaled_smoothed_estimator=get(3, (m, n + 1)),
                    scaled_smoothed_estimator_cov=get(4, (m, m, n + 1)),
                    smoothed_measurement_disturbance=get(5, (p, n)),
                    smoothed_measurement_disturbance_cov=get(6, (p, p, n)),
                    smoothed_state_disturbance=get(7, (r, n)),
                    smoothed_state_disturbance_cov=get(8, (r, r, n)),
                )
            finally:
                _free("ss_smoother_free", sh)
        finally:
            _free("ss_filter_free", fh)

    def _filter_results(self, fh):
        p, m, n = self.k_endog, self.k_states, self.nobs
        llf = ctypes.c_double()
        nd, kd, ts = ctypes.c_int(), ctypes.c_int(), ctypes.c_int()
        call("ss_filter_scalars", fh, ctypes.byref(llf), ctypes.byref(nd), ctypes.byref(kd),
             ctypes.byref(ts))
        get = lambda code, shape: _get("ss_filter_get", fh, code, shape)
        return FilterResults(
            llf=llf.value, llf_obs=get(12, (n,)), nobs_diffuse=nd.value,
            k_diffuse=kd.value, t_steady=ts.value,
            predicted_state=get(1, (m, n + 1)),
            predicted_state_cov=get(2, (m, m, n + 1)),
            predicted_diffuse_state_cov=get(3, (m, m, n + 1)),
            filtered_state=get(4, (m, n)),
            filtered_state_cov=get(5, (m, m, n)),
            forecasts=get(6, (p, n)),
            forecasts_error=get(7, (p, n)),
            forecasts_error_cov=get(8, (p, p, n)),
            forecasts_error_diffuse_cov=get(9, (p, p, n)),
            forecasts_error_cov_inv=get(10, (p, p, n)),
            kalman_gain=get(11, (m, p, n)),
            standardized_forecasts_error=self._standardized(fh),
            diagnostic_start=self._diagnostic_start(fh),
        )

    def _standardized(self, fh):
        out = np.empty((self.k_endog, self.nobs), order="F")
        call("ss_standardized_residuals", fh, ptr(out))
        return out

    def _diagnostic_start(self, fh):
        t0 = ctypes.c_int()
        call("ss_diagnostic_start", self._h, fh, ctypes.byref(t0))
        return t0.value - 1   # 0-based

    # -- forecasting and simulation -----------------------------------------

    def forecast(self, steps):
        """Forecasts of y for `steps` periods past the sample (DK 4.11;
        time-invariant models): mean (p, steps) and covariance (p, p, steps)."""
        fh = ctypes.c_void_p()
        call("ss_filter", self._h, ctypes.byref(fh))
        try:
            mean = np.empty((self.k_endog, steps), order="F")
            cov = np.empty((self.k_endog, self.k_endog, steps), order="F")
            call("ss_forecast", self._h, fh, int(steps), ptr(mean), ptr(cov))
        finally:
            _free("ss_filter_free", fh)
        return mean, cov

    def _variates(self, rng, variates):
        p, m, r, n = self.k_endog, self.k_states, self.k_posdef, self.nobs
        if variates is not None:
            ui, ue, un = variates
            return farray(ui, (m,)), farray(ue, (p, n)), farray(un, (r, n))
        rng = np.random.default_rng(rng)
        return (rng.standard_normal(m), np.asfortranarray(rng.standard_normal((p, n))),
                np.asfortranarray(rng.standard_normal((r, n))))

    def simulate(self, rng=None, variates=None):
        """Simulate (y, alpha, eps, eta) from the model; rng is a numpy
        Generator or seed, or give the standard normal `variates` (u_init (m),
        u_eps (p, n), u_eta (r, n)). Diffuse states start at a1."""
        ui, ue, un = self._variates(rng, variates)
        p, m, r, n = self.k_endog, self.k_states, self.k_posdef, self.nobs
        y, alpha = np.empty((p, n), order="F"), np.empty((m, n), order="F")
        eps, eta = np.empty((p, n), order="F"), np.empty((r, n), order="F")
        call("ss_simulate", self._h, ptr(ui), ptr(ue), ptr(un), ptr(y), ptr(alpha), ptr(eps),
             ptr(eta))
        return y, alpha, eps, eta

    def simulation_smoother(self, rng=None, variates=None):
        """One draw of (alpha, eps, eta) given the data (DK 4.9.2, mean
        corrections); rng and variates as for `simulate`."""
        ui, ue, un = self._variates(rng, variates)
        p, m, r, n = self.k_endog, self.k_states, self.k_posdef, self.nobs
        state, eps = np.empty((m, n), order="F"), np.empty((p, n), order="F")
        eta = np.empty((r, n), order="F")
        call("ss_simulation_smoother", self._h, ptr(ui), ptr(ue), ptr(un), ptr(state),
             ptr(eps), ptr(eta))
        return state, eps, eta

    # -- DK ch. 4 extras ------------------------------------------------------

    def _with_filter(self, fn, smoother=False):
        fh, sh = ctypes.c_void_p(), ctypes.c_void_p()
        call("ss_filter", self._h, ctypes.byref(fh))
        try:
            if smoother:
                call("ss_smooth", self._h, fh, ctypes.byref(sh))
                try:
                    return fn(fh, sh)
                finally:
                    _free("ss_smoother_free", sh)
            return fn(fh)
        finally:
            _free("ss_filter_free", fh)

    def _mn(self):
        return (self.k_states, self.nobs), (self.k_states, self.k_states, self.nobs)

    def fast_smoother(self):
        """Smoothed state by the fast state smoother (DK 4.6.2), (m, n)."""
        a = np.empty(self._mn()[0], order="F")
        self._with_filter(lambda fh: call("ss_fast_smoother", self._h, fh, ptr(a)))
        return a

    def classical_smoother(self):
        """Classical (Rauch-Tung-Striebel) smoother (DK 4.6.1): (alphahat, V)."""
        a, V = np.empty(self._mn()[0], order="F"), np.empty(self._mn()[1], order="F")
        self._with_filter(lambda fh: call("ss_classical_smoother", self._h, fh, ptr(a), ptr(V)))
        return a, V

    def two_filter_smoother(self):
        """Two-filter formula for smoothing (DK 4.6.4): (alphahat, V)."""
        a, V = np.empty(self._mn()[0], order="F"), np.empty(self._mn()[1], order="F")
        self._with_filter(lambda fh: call("ss_two_filter_smoother", self._h, fh, ptr(a), ptr(V)))
        return a, V

    def whittle_smoother(self):
        """Smoothed state from the Whittle relation (DK 4.6.3), (m, n). The
        recursion runs backwards through T^-1 and can be numerically unstable
        (as DK note); prefer `smooth` in practice."""
        a = np.empty(self._mn()[0], order="F")
        self._with_filter(lambda fh: call("ss_whittle_smoother", self._h, fh, ptr(a)))
        return a

    def fixed_point_smoother(self, t):
        """Estimates of alpha_t (0-based t) given y up to k = t..n-1 (DK 4.4.6):
        (means (m, n-t), variances (m, m, n-t))."""
        k = self.nobs - t
        a = np.empty((self.k_states, k), order="F")
        V = np.empty((self.k_states, self.k_states, k), order="F")
        self._with_filter(lambda fh: call("ss_fixed_point_smoother", self._h, fh, int(t) + 1,
                                          ptr(a), ptr(V)))
        return a, V

    def fixed_lag_smoother(self, lag):
        """alpha_s given y up to s + lag (DK 4.4.6): (alphahat, V), NaN in the
        last `lag` periods."""
        a, V = np.empty(self._mn()[0], order="F"), np.empty(self._mn()[1], order="F")
        self._with_filter(lambda fh: call("ss_fixed_lag_smoother", self._h, fh, int(lag),
                                          ptr(a), ptr(V)))
        return a, V

    def update_smoothed(self, alphahat, V, nobs_old):
        """Update smoothed estimates for the first `nobs_old` periods (given
        those observations only) to all n observations (DK 4.4.5). alphahat
        (m, >= nobs_old) and V (m, m, >= nobs_old); returns full (m, n) and
        (m, m, n) arrays."""
        a = np.full(self._mn()[0], np.nan, order="F")
        W = np.full(self._mn()[1], np.nan, order="F")
        a[:, :nobs_old] = np.asarray(alphahat)[:, :nobs_old]
        W[:, :, :nobs_old] = np.asarray(V)[:, :, :nobs_old]
        self._with_filter(lambda fh: call("ss_update_smoothed", self._h, fh, int(nobs_old),
                                          ptr(a), ptr(W)))
        return a, W

    def smoothed_state_autocov(self):
        """Cov(alpha_t+1, alpha_t | Y_n), (m, m, n-1) (DK 4.7)."""
        out = np.empty((self.k_states, self.k_states, self.nobs - 1), order="F")
        self._with_filter(lambda fh, sh: call("ss_smoothed_state_autocov", self._h, fh, sh,
                                              ptr(out)), smoother=True)
        return out

    def smoothed_state_cov_between(self, t, j):
        """Cov(alpha_t, alpha_j | Y_n) (0-based periods; DK 4.7)."""
        out = np.empty((self.k_states, self.k_states), order="F")
        self._with_filter(lambda fh, sh: call("ss_smoothed_state_cov_between", self._h, fh, sh,
                                              int(t) + 1, int(j) + 1, ptr(out)), smoother=True)
        return out

    def filtered_state_weights(self):
        """Weights of past observations in a_t and a_t|t (DK 4.8.2): (Wa, Watt),
        each (m, p, n, n) with [:, :, t, j] the weight of y_j."""
        m, p, n = self.k_states, self.k_endog, self.nobs
        Wa, Watt = np.empty((m, p, n, n), order="F"), np.empty((m, p, n, n), order="F")
        self._with_filter(lambda fh: call("ss_filtered_state_weights", self._h, fh, ptr(Wa),
                                          ptr(Watt)))
        return Wa, Watt

    def smoothed_state_weights(self):
        """Weights of the smoothed state (DK 4.8.3): alphahat_t = sum_j
        W[:, :, t, j] (y_j - d_j) + sum_j C[:, :, t, j] c_j + A[:, :, t] a1."""
        m, p, n = self.k_states, self.k_endog, self.nobs
        W = np.empty((m, p, n, n), order="F")
        C = np.empty((m, m, n, n), order="F")
        A = np.empty((m, m, n), order="F")
        call("ss_smoothed_state_weights", self._h, ptr(W), ptr(C), ptr(A))
        return W, C, A

    def innovation_transition(self, t):
        """L_t (0-based t), the transition of the state prediction error (DK 4.3)."""
        out = np.empty((self.k_states, self.k_states), order="F")
        self._with_filter(lambda fh: call("ss_innovation_transition", self._h, fh, int(t) + 1,
                                          ptr(out)))
        return out

    def augmented(self):
        """Augmented Kalman filter and smoother (DK 5.7)."""
        ah = ctypes.c_void_p()
        call("ss_augmented_filter", self._h, ctypes.byref(ah))
        try:
            k = ctypes.c_int()
            llf, llf_fixed, llf_marg = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
            call("ss_augmented_info", ah, ctypes.byref(k), ctypes.byref(llf),
                 ctypes.byref(llf_fixed), ctypes.byref(llf_marg))
            delta = np.empty(k.value)
            delta_cov = np.empty((k.value, k.value), order="F")
            if k.value > 0:
                call("ss_augmented_get", ah, 1, ptr(delta))
                call("ss_augmented_get", ah, 2, ptr(delta_cov))
            fh = ctypes.c_void_p()
            call("ss_augmented_filter_result", ah, ctypes.byref(fh))
            fres = self._filter_results(fh)
            a, V = np.empty(self._mn()[0], order="F"), np.empty(self._mn()[1], order="F")
            call("ss_augmented_smoother", self._h, ah, ptr(a), ptr(V))
        finally:
            _free("ss_augmented_free", ah)
        return AugmentedResults(filter=fres, delta=delta, delta_cov=delta_cov, llf=llf.value,
                                llf_fixed=llf_fixed.value, llf_marginal=llf_marg.value,
                                smoothed_state=a, smoothed_state_cov=V)

    def em(self, maxiter=500, tol=1e-8, diagonal_obs_cov=False, diagonal_state_cov=False):
        """EM for obs_cov (H) and state_cov (Q) (DK 7.3.4), updating this
        representation. Returns (llf, niter, llf_path)."""
        llf, niter = ctypes.c_double(), ctypes.c_int()
        path = np.empty(int(maxiter))
        call("ss_em", self._h, int(maxiter), float(tol), int(bool(diagonal_obs_cov)),
             int(bool(diagonal_state_cov)), ctypes.byref(llf), ctypes.byref(niter), ptr(path))
        return llf.value, niter.value, path[:niter.value]

    def collapse(self):
        """The collapsed representation (DK 6.5) and the log likelihood
        adjustment: loglike() == collapsed.loglike() + adjust.sum()."""
        h = ctypes.c_void_p()
        adj = np.empty(self.nobs)
        call("ss_collapse", self._h, ctypes.byref(h), ptr(adj))
        return Representation._wrap(h, owner=None), adj

    def add_state_restrictions(self, restriction, value):
        """A representation with R*_t alpha_t = r*_t imposed (DK 6.6):
        restriction (q, m) or (q, m, n), value (q, n) with NaN where a
        restriction is inactive. Use filter_method=FILTER_UNIVARIATE."""
        Rm = np.asarray(restriction, dtype=np.float64)
        if Rm.ndim == 2:
            Rm = Rm[:, :, None]
        q, nt = Rm.shape[0], Rm.shape[2]
        Rm = np.asfortranarray(Rm)
        r = farray(value, (q, self.nobs))
        h = ctypes.c_void_p()
        call("ss_add_restrictions", self._h, q, nt, ptr(Rm), ptr(r), ctypes.byref(h))
        return Representation._wrap(h, owner=None)

    # -- more diagnostics ------------------------------------------------------

    def auxiliary_residuals(self, vector=False):
        """Standardized smoothed disturbances (DK 7.5): (eps (p, n), eta (r, n));
        vector=True standardizes each period's vector (Cholesky)."""
        e = np.empty((self.k_endog, self.nobs), order="F")
        u = np.empty((self.k_posdef, self.nobs), order="F")
        self._with_filter(lambda fh, sh: call("ss_auxiliary_residuals", self._h, sh,
                                              int(bool(vector)), ptr(e), ptr(u)), smoother=True)
        return e, u

    def de_jong_penzer(self):
        """de Jong and Penzer (1998) shock statistics (DK 7.5):
        (state (m, n), observation (p, n))."""
        r = np.empty((self.k_states, self.nobs), order="F")
        e = np.empty((self.k_endog, self.nobs), order="F")
        self._with_filter(lambda fh, sh: call("ss_de_jong_penzer", self._h, fh, sh, ptr(r), ptr(e)),
                          smoother=True)
        return r, e

    def least_squares_residuals(self, start, stop):
        """Least squares residuals (DK 6.2.4) with the regression coefficients
        in states start..stop-1 (0-based slice), (p, n)."""
        out = np.empty((self.k_endog, self.nobs), order="F")
        call("ss_least_squares_residuals", self._h, int(start) + 1, int(stop), ptr(out))
        return out

    def r2_diffuse(self, series=0):
        """R^2_D against a random walk with drift (Harvey 1989)."""
        r2 = ctypes.c_double()
        self._with_filter(lambda fh: call("ss_r2_diffuse", self._h, fh, int(series) + 1,
                                          ctypes.byref(r2)))
        return r2.value

    def djs_simulation_smoother(self, rng=None, variates=None):
        """de Jong-Shephard simulation smoother (DK 4.9.3): (alpha, eps, eta)."""
        ui, ue, un = self._variates(rng, variates)
        p, m, r, n = self.k_endog, self.k_states, self.k_posdef, self.nobs
        eps, eta = np.empty((p, n), order="F"), np.empty((r, n), order="F")
        alpha = np.empty((m, n), order="F")
        self._with_filter(lambda fh: call("ss_djs_simulation_smoother", self._h, fh, ptr(ue),
                                          ptr(un), ptr(ui), ptr(eps), ptr(eta), ptr(alpha)))
        return alpha, eps, eta

    def steady_state(self):
        """Steady-state (P, F) of a time-invariant model (DK 4.3.4)."""
        P = np.empty((self.k_states, self.k_states), order="F")
        F = np.empty((self.k_endog, self.k_endog), order="F")
        call("ss_steady_state", self._h, ptr(P), ptr(F))
        return P, F


def _get(routine, handle, code, shape):
    out = np.empty(shape, dtype=np.float64, order="F")
    call(routine, handle, code, ptr(out))
    return out


def _free(routine, handle):
    if handle:
        call(routine, handle)
