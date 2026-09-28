"""State space representations filled from numpy.

A :class:`Representation` holds the linear Gaussian state space model

.. math::

    y_t &= d_t + Z_t \\alpha_t + \\varepsilon_t, &
    \\varepsilon_t &\\sim N(0, H_t), \\\\
    \\alpha_{t+1} &= c_t + T_t \\alpha_t + R_t \\eta_t, &
    \\eta_t &\\sim N(0, Q_t), \\qquad t = 1, \\dots, n,

(DK eq. 3.1, with intercepts) together with the distribution of
:math:`\\alpha_1`, and runs the algorithms of DK Part I on it: filtering,
smoothing, likelihood evaluation, simulation and forecasting. Every model
builds one (:meth:`~ssfortran.Model.representation`); construct one
directly to run the algorithms at fixed matrices.

The system matrices take statsmodels' names and shapes, with time as the
last axis:

=====================  ========  ==================
name                   symbol    shape
=====================  ========  ==================
``design``             Z         (p, m) or (p, m, n)
``obs_cov``            H         (p, p) or (p, p, n)
``transition``         T         (m, m) or (m, m, n)
``selection``          R         (m, r) or (m, r, n)
``state_cov``          Q         (r, r) or (r, r, n)
``state_intercept``    c         (m,) or (m, n)
``obs_intercept``      d         (p,) or (p, n)
=====================  ========  ==================

Here p is the number of series, m the number of states, r the number of
state disturbances and n the number of periods. A matrix without the time
axis, or with a time axis of length 1, is time-invariant.

All arrays are copied into the Fortran library; reading an item returns a
copy.

References
----------
.. [1] Durbin, J. and Koopman, S. J. (2012). *Time Series Analysis by State
   Space Methods*, 2nd ed. Oxford University Press. Cited as DK.
"""

import ctypes
import weakref
from dataclasses import dataclass

import numpy as np

from ._lib import call, farray, ptr

#: Initialization of a block of states: known mean and variance.
INIT_KNOWN = 1
#: Initialization: approximate diffuse, variance ``kappa * I``.
INIT_APPROX_DIFFUSE = 2
#: Initialization: the unconditional distribution of a stationary block.
INIT_STATIONARY = 3
#: Initialization: exact diffuse (DK §5.2).
INIT_DIFFUSE = 4
#: Initialization: general ``N(a1, Pstar + kappa Pinf)``, kappa to infinity.
INIT_GENERAL = 5
#: Filter: the observation vector at once (DK §4.3).
FILTER_CONVENTIONAL = 0
#: Filter: one element of the observation vector at a time (DK §6.4).
FILTER_UNIVARIATE = 1
#: Diffuse periods: univariate exact initial filter (DK §5.2.5, the default).
DIFFUSE_UNIVARIATE = 0
#: Diffuse periods: multivariate exact initial filter where possible (DK §5.2).
DIFFUSE_MULTIVARIATE = 1

_ARR = {
    "endog": 1,
    "design": 2,
    "obs_cov": 3,
    "transition": 4,
    "selection": 5,
    "state_cov": 6,
    "state_intercept": 7,
    "obs_intercept": 8,
}
_OPT_INT = {
    "filter_method": 1,
    "diffuse_method": 2,
    "loglikelihood_burn": 3,
    "marginal_likelihood": 4,
}
_OPT_REAL = {"tol_diffuse": 1, "tol_steady": 2}


@dataclass
class FilterResults:
    """Output of the Kalman filter (DK §4.3, §5.2).

    Arrays have time as the last axis and 0-based periods: column ``t`` of
    ``predicted_state`` is :math:`a_{t+1}` in DK's 1-based notation, and
    its last column is the one-step prediction past the sample.

    Attributes
    ----------
    llf : float
        Log likelihood (DK eq. 7.2); the diffuse log likelihood (DK §7.2.2)
        with diffuse states.
    llf_obs : ndarray, shape (n,)
        Contribution of each period to `llf`.
    nobs_diffuse : int
        Number of periods with a diffuse state (DK's *d*).
    k_diffuse : int
        Rank of the diffuse part of the initial variance.
    t_steady : int
        First period (1-based) that used the steady-state shortcut, or 0.
    predicted_state : ndarray, shape (m, n+1)
        :math:`a_t = E(\\alpha_t | y_1, \\dots, y_{t-1})`.
    predicted_state_cov : ndarray, shape (m, m, n+1)
        :math:`P_t`, or :math:`P_{*,t}` in diffuse periods.
    predicted_diffuse_state_cov : ndarray, shape (m, m, n+1)
        :math:`P_{\\infty,t}`; zero after the diffuse periods.
    filtered_state : ndarray, shape (m, n)
        :math:`a_{t|t}` (DK eq. 4.24).
    filtered_state_cov : ndarray, shape (m, m, n)
        :math:`P_{t|t}`.
    forecasts : ndarray, shape (p, n)
        One-step predictions :math:`d_t + Z_t a_t`.
    forecasts_error : ndarray, shape (p, n)
        Innovations :math:`v_t`; NaN where y is missing.
    forecasts_error_cov : ndarray, shape (p, p, n)
        :math:`F_t`, in the original coordinates also in univariate periods.
    forecasts_error_diffuse_cov : ndarray, shape (p, p, n)
        :math:`F_{\\infty,t} = Z_t P_{\\infty,t} Z_t'`.
    forecasts_error_cov_inv : ndarray, shape (p, p, n)
        :math:`F_t^{-1}` over the observed elements, zero-padded; NaN in
        periods processed element by element.
    kalman_gain : ndarray, shape (m, p, n)
        :math:`K_t = T_t P_t Z_t' F_t^{-1}`; NaN in periods processed element
        by element.
    standardized_forecasts_error : ndarray, shape (p, n)
        :math:`L_t^{-1} v_t` with :math:`F_t = L_t L_t'` (DK §7.5); NaN
        where missing or diffuse.
    diagnostic_start : int
        First 0-based period after the burn-in and the diffuse periods.
    """

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
    """Output of the state and disturbance smoother (DK §4.4-4.5, §5.3).

    Attributes
    ----------
    filter : FilterResults
        The filter run the smoother used.
    smoothed_state : ndarray, shape (m, n)
        :math:`\\hat\\alpha_t = E(\\alpha_t | Y_n)` (DK eq. 4.39).
    smoothed_state_cov : ndarray, shape (m, m, n)
        :math:`V_t = Var(\\alpha_t | Y_n)` (DK eq. 4.43).
    scaled_smoothed_estimator : ndarray, shape (m, n+1)
        :math:`r_t` for t = 0, ..., n.
    scaled_smoothed_estimator_cov : ndarray, shape (m, m, n+1)
        :math:`N_t` for t = 0, ..., n.
    smoothed_measurement_disturbance : ndarray, shape (p, n)
        :math:`\\hat\\varepsilon_t = E(\\varepsilon_t | Y_n)`. For a missing
        element this is its conditional mean given the observed elements,
        not zero.
    smoothed_measurement_disturbance_cov : ndarray, shape (p, p, n)
        :math:`Var(\\varepsilon_t | Y_n)`.
    smoothed_state_disturbance : ndarray, shape (r, n)
        :math:`\\hat\\eta_t = Q_t R_t' r_t` (DK eq. 4.63).
    smoothed_state_disturbance_cov : ndarray, shape (r, r, n)
        :math:`Var(\\eta_t | Y_n)`.
    """

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
    """Output of the augmented Kalman filter and smoother (DK §5.7).

    The diffuse elements :math:`\\delta` of the initial state
    (DK eq. 5.2) are treated as unknown constants: the filter runs once with
    :math:`\\delta = 0` and carries their effect as extra columns.

    Attributes
    ----------
    filter : FilterResults
        The filter with :math:`\\delta = 0`.
    delta : ndarray, shape (k,)
        Generalized least squares estimate of :math:`\\delta`.
    delta_cov : ndarray, shape (k, k)
        Its variance.
    llf : float
        Diffuse log likelihood (DK §7.2.3); equal to that of the exact
        initial filter.
    llf_fixed : float
        Log likelihood with :math:`\\delta` fixed but unknown (DK §7.2.4).
    llf_marginal : float
        Marginal log likelihood (DK §7.2.6; Francke, Koopman and de Vos
        2010).
    smoothed_state : ndarray, shape (m, n)
        Smoothed state.
    smoothed_state_cov : ndarray, shape (m, m, n)
        Its variance, including the uncertainty about :math:`\\delta`.
    """

    filter: FilterResults
    delta: np.ndarray
    delta_cov: np.ndarray
    llf: float
    llf_fixed: float
    llf_marginal: float
    smoothed_state: np.ndarray
    smoothed_state_cov: np.ndarray


class Representation:
    """A linear Gaussian state space model and its algorithms.

    The system matrices start at zero, except ``selection``, which starts as
    the leading m x r identity; they are set by item assignment. See the
    module description for their names and shapes. The initial state
    must be set with one of the ``initialize_*`` methods before filtering.

    Parameters
    ----------
    endog : array_like, shape (n,) or (n, p)
        Observations, time first as in statsmodels. NaN marks a missing
        value; any subset of a period's elements may be missing.
    k_states : int
        Number of states m.
    k_posdef : int, optional
        Number of state disturbances r. The default is `k_states`.

    Attributes
    ----------
    nobs, k_endog, k_states, k_posdef : int
        n, p, m and r.
    filter_method : int
        :data:`FILTER_CONVENTIONAL` (default) or :data:`FILTER_UNIVARIATE`.
        Set-only.
    diffuse_method : int
        :data:`DIFFUSE_UNIVARIATE` (default) or :data:`DIFFUSE_MULTIVARIATE`.
        Set-only.
    loglikelihood_burn : int
        Number of leading periods left out of the log likelihood. Set-only.
    marginal_likelihood : bool
        Report the marginal instead of the diffuse log likelihood (DK
        §7.2.6). Set-only.
    tol_diffuse : float
        Threshold on :math:`F_\\infty` and :math:`\\|P_\\infty\\|_F^2` in the
        diffuse periods; default 1e-10, as in statsmodels. Set-only.
    tol_steady : float
        Relative change in :math:`P_t` below which the filter holds P, F and
        K fixed (DK §4.3.4); default 1e-15. Negative turns the shortcut off.
        Set-only.

    See Also
    --------
    StructuralModel, MappedModel, MLEModel : Models with parameters.

    Notes
    -----
    The filter treats diffuse periods element by element (DK §5.2.5,
    §6.4) and reports :math:`v_t` and :math:`F_t` in the original
    coordinates; the per-element quantities stay internal. With
    time-invariant system matrices the conventional filter stops updating
    :math:`P_t` once it converges and resumes after a missing observation.

    Examples
    --------
    The local level model (DK ch. 2) with exact diffuse initialization:

    >>> import numpy as np
    >>> import ssfortran as ss
    >>> rng = np.random.default_rng(0)
    >>> y = rng.standard_normal(100).cumsum() + rng.standard_normal(100)
    >>> rep = ss.Representation(y, k_states=1)
    >>> rep["design"] = rep["transition"] = rep["selection"] = [[1.0]]
    >>> rep["obs_cov"] = [[1.0]]
    >>> rep["state_cov"] = [[1.0]]
    >>> rep.initialize_diffuse()
    >>> round(rep.loglike(), 4)
    -189.1596
    >>> res = rep.smooth()
    >>> res.smoothed_state.shape
    (1, 100)
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
        call(
            "ss_rep_info",
            handle,
            ctypes.byref(p),
            ctypes.byref(m),
            ctypes.byref(r),
            ctypes.byref(n),
            ptr(nts),
        )
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
        return {
            "endog": (p, self.nobs),
            "design": (p, m),
            "obs_cov": (p, p),
            "transition": (m, m),
            "selection": (m, r),
            "state_cov": (r, r),
            "state_intercept": (m,),
            "obs_intercept": (p,),
        }[name]

    def __setitem__(self, name, value):
        """Set a system matrix, or ``endog`` with shape (p, n)."""
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
            raise ValueError(
                f"{name} must have shape {base} or {base + (self.nobs,)}, not {a.shape}"
            )
        call("ss_rep_set", self._h, _ARR[name], nt, ptr(a))

    def __getitem__(self, name):
        """Copy of a system matrix; time-invariant ones without the time axis."""
        if name not in _ARR:
            raise KeyError(name)
        p, m, r, n = (ctypes.c_int() for _ in range(4))
        nts = np.zeros(8, dtype=np.int32)
        call(
            "ss_rep_info",
            self._h,
            ctypes.byref(p),
            ctypes.byref(m),
            ctypes.byref(r),
            ctypes.byref(n),
            ptr(nts),
        )
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
        """Initialize with a known mean and variance.

        Parameters
        ----------
        a1 : array_like, shape (m,)
            Mean of the initial state.
        P1 : array_like, shape (m, m)
            Variance of the initial state.

        See Also
        --------
        initialize_diffuse : Exact diffuse initialization.
        initialize_stationary : The unconditional distribution.
        initialize_block : Different initializations for blocks of states.
        """
        m = self.k_states
        a = farray(a1, (m,))
        P = farray(P1, (m, m))
        call("ss_rep_init_known", self._h, ptr(a), ptr(P))

    def initialize_diffuse(self):
        """Initialize every state as exact diffuse (DK §5.2).

        The initial variance is :math:`\\kappa I` with
        :math:`\\kappa \\to \\infty`, handled exactly by the exact initial
        Kalman filter; the log likelihood is then the diffuse log likelihood
        (DK §7.2.2).

        See Also
        --------
        initialize_approximate_diffuse : A large finite variance instead.
        """
        call("ss_rep_init_diffuse", self._h)

    def initialize_approximate_diffuse(self, variance=1e6):
        """Initialize with mean zero and a large variance, ``variance * I``.

        Parameters
        ----------
        variance : float, optional
            Variance of each state; default 1e6, as in statsmodels.

        See Also
        --------
        initialize_diffuse : The exact treatment (DK §5.2).
        """
        call("ss_rep_init_approx_diffuse", self._h, float(variance))

    def initialize_stationary(self):
        """Initialize with the unconditional distribution of the state.

        The mean is :math:`(I - T)^{-1} c` and the variance solves
        :math:`P = T P T' + R Q R'` (DK §5.6.2), recomputed from the current
        matrices at every filter run. The model must be time-invariant in T,
        R, Q and c at the first period.

        Raises
        ------
        StateSpaceError
            From the filter, with code 6, if T has an eigenvalue on or
            outside the unit circle.
        """
        call("ss_rep_init_stationary", self._h)

    def initialize(self, a1, Pstar, Pinf):
        """Initialize with :math:`\\alpha_1 \\sim N(a_1, P_* + \\kappa P_\\infty)`.

        The general form of DK §5.1, with :math:`\\kappa \\to \\infty`
        handled exactly.

        Parameters
        ----------
        a1 : array_like, shape (m,)
            Mean.
        Pstar : array_like, shape (m, m)
            Variance of the proper part.
        Pinf : array_like, shape (m, m)
            Variance direction of the diffuse part; usually a selection
            :math:`A A'`.
        """
        m = self.k_states
        a, Ps, Pi = farray(a1, (m,)), farray(Pstar, (m, m)), farray(Pinf, (m, m))
        call("ss_rep_init_general", self._h, ptr(a), ptr(Ps), ptr(Pi))

    def initialize_block(self, start, stop, kind, a1=None, P1=None, variance=1e6):
        """Initialize a block of states, leaving the others as they are.

        Blocks combine, for example a diffuse trend with a stationary ARMA
        block (DK §5.6).

        Parameters
        ----------
        start, stop : int
            The block is states ``start:stop`` (0-based, like a slice).
        kind : int
            :data:`INIT_KNOWN`, :data:`INIT_APPROX_DIFFUSE`,
            :data:`INIT_STATIONARY` or :data:`INIT_DIFFUSE`.
        a1, P1 : array_like, optional
            Mean (stop - start,) and variance for :data:`INIT_KNOWN`.
        variance : float, optional
            Variance for :data:`INIT_APPROX_DIFFUSE`.

        Notes
        -----
        A stationary block must not depend on states outside it through T.

        Examples
        --------
        A diffuse level with a stationary AR(1) disturbance:

        >>> rep = ss.Representation(np.zeros(20), k_states=2)
        >>> rep["transition"] = [[1.0, 0.0], [0.0, 0.5]]
        >>> rep.initialize_block(0, 1, ss.INIT_DIFFUSE)
        >>> rep.initialize_block(1, 2, ss.INIT_STATIONARY)
        """
        nb = stop - start
        if a1 is not None:
            a, P = farray(a1, (nb,)), farray(P1, (nb, nb))
            call(
                "ss_rep_init_block",
                self._h,
                start + 1,
                stop,
                kind,
                ptr(a),
                ptr(P),
                float(variance),
            )
        else:
            call(
                "ss_rep_init_block",
                self._h,
                start + 1,
                stop,
                kind,
                None,
                None,
                float(variance),
            )

    # -- computations ----------------------------------------------------

    def loglike(self):
        """Evaluate the log likelihood without storing the filter output.

        Returns
        -------
        float
            The log likelihood (DK eq. 7.2), or the diffuse log likelihood
            (DK §7.2.2) with diffuse states, or the marginal log likelihood
            when `marginal_likelihood` is set.

        Raises
        ------
        StateSpaceError
            If :math:`F_t` is not positive definite (code 2), the
            initialization is missing (3) or not stationary (6).

        See Also
        --------
        loglike_concentrated : With the scale concentrated out.
        filter : The full filter output.
        """
        llf = ctypes.c_double()
        call("ss_loglike", self._h, ctypes.byref(llf))
        return llf.value

    def loglike_concentrated(self):
        """Evaluate the log likelihood with the scale concentrated out.

        H, Q and :math:`P_*` are taken relative to a scale :math:`\\sigma^2`,
        which is replaced by its maximum likelihood estimate (DK §2.10.2,
        §7.3).

        Returns
        -------
        llf : float
            Concentrated log likelihood.
        scale : float
            Estimate of :math:`\\sigma^2`.
        """
        llf, scale = ctypes.c_double(), ctypes.c_double()
        call("ss_loglike_concentrated", self._h, ctypes.byref(llf), ctypes.byref(scale))
        return llf.value, scale.value

    def filter(self, method="conventional"):
        """Run the Kalman filter.

        Parameters
        ----------
        method : {"conventional", "sqrt"}, optional
            ``"sqrt"`` propagates square roots of :math:`P_t` (DK §6.3);
            it gives the same output and does not support diffuse states.

        Returns
        -------
        FilterResults
            The filter output.

        See Also
        --------
        smooth : Filter and smoother.
        loglike : Log likelihood without storing the output.
        """
        fh = ctypes.c_void_p()
        call(
            "ss_sqrt_filter" if method == "sqrt" else "ss_filter",
            self._h,
            ctypes.byref(fh),
        )
        try:
            return self._filter_results(fh)
        finally:
            _free("ss_filter_free", fh)

    def smooth(self, method="conventional"):
        """Run the filter and the state and disturbance smoother.

        Parameters
        ----------
        method : {"conventional", "sqrt"}, optional
            ``"sqrt"`` uses the square root filter and smoother (DK §6.3).

        Returns
        -------
        SmootherResults
            Smoothed states and disturbances with their variances.

        See Also
        --------
        fast_smoother, classical_smoother, two_filter_smoother :
            Alternative smoothers (DK §4.6).

        Examples
        --------
        >>> rep = ss.Representation(np.arange(5.0), k_states=1)
        >>> rep["design"] = rep["transition"] = rep["selection"] = [[1.0]]
        >>> rep["obs_cov"] = rep["state_cov"] = [[1.0]]
        >>> rep.initialize_diffuse()
        >>> res = rep.smooth()
        >>> np.round(res.smoothed_state, 3)
        array([[0.6, 1.2, 2. , 2.8, 3.4]])
        """
        fh, sh = ctypes.c_void_p(), ctypes.c_void_p()
        call(
            "ss_sqrt_filter" if method == "sqrt" else "ss_filter",
            self._h,
            ctypes.byref(fh),
        )
        try:
            fres = self._filter_results(fh)
            if method == "sqrt":
                call("ss_sqrt_smoother", self._h, fh, ctypes.byref(sh))
            else:
                call("ss_smooth", self._h, fh, ctypes.byref(sh))
            try:
                p, m, r, n = self.k_endog, self.k_states, self.k_posdef, self.nobs

                def get(code, shape):
                    return _get("ss_smoother_get", sh, code, shape)

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
        call(
            "ss_filter_scalars",
            fh,
            ctypes.byref(llf),
            ctypes.byref(nd),
            ctypes.byref(kd),
            ctypes.byref(ts),
        )

        def get(code, shape):
            return _get("ss_filter_get", fh, code, shape)

        return FilterResults(
            llf=llf.value,
            llf_obs=get(12, (n,)),
            nobs_diffuse=nd.value,
            k_diffuse=kd.value,
            t_steady=ts.value,
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
        return t0.value - 1  # 0-based

    # -- forecasting and simulation -----------------------------------------

    def forecast(self, steps):
        """Forecast the observations past the end of the sample (DK §4.11).

        Parameters
        ----------
        steps : int
            Number of periods ahead.

        Returns
        -------
        mean : ndarray, shape (p, steps)
            :math:`E(y_{n+j} | Y_n)`.
        cov : ndarray, shape (p, p, steps)
            :math:`Var(y_{n+j} | Y_n)`.

        Raises
        ------
        StateSpaceError
            With code 4 if a system matrix varies over time. Append NaN
            observations, with the future matrices, to forecast such models.
        """
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
        return (
            rng.standard_normal(m),
            np.asfortranarray(rng.standard_normal((p, n))),
            np.asfortranarray(rng.standard_normal((r, n))),
        )

    def simulate(self, rng=None, variates=None):
        """Simulate observations and states from the model.

        Parameters
        ----------
        rng : numpy.random.Generator or int, optional
            Source of the standard normal variates.
        variates : tuple of array_like, optional
            The variates themselves, ``(u_init (m,), u_eps (p, n),
            u_eta (r, n))``; overrides `rng`.

        Returns
        -------
        y : ndarray, shape (p, n)
            Observations.
        alpha : ndarray, shape (m, n)
            States.
        eps : ndarray, shape (p, n)
            Observation disturbances.
        eta : ndarray, shape (r, n)
            State disturbances.

        Notes
        -----
        :math:`\\alpha_1 = a_1 + S(P_*) u_{init}`, so diffuse states start at
        their mean :math:`a_1`.
        """
        ui, ue, un = self._variates(rng, variates)
        p, m, r, n = self.k_endog, self.k_states, self.k_posdef, self.nobs
        y, alpha = np.empty((p, n), order="F"), np.empty((m, n), order="F")
        eps, eta = np.empty((p, n), order="F"), np.empty((r, n), order="F")
        call(
            "ss_simulate",
            self._h,
            ptr(ui),
            ptr(ue),
            ptr(un),
            ptr(y),
            ptr(alpha),
            ptr(eps),
            ptr(eta),
        )
        return y, alpha, eps, eta

    def simulation_smoother(self, rng=None, variates=None):
        """Draw states and disturbances from their distribution given the data.

        Uses the mean-correction method of Durbin and Koopman (2002)
        (DK §4.9.2): simulate from the model, smooth the simulated data and
        shift by the smoothed values of the actual data.

        Parameters
        ----------
        rng : numpy.random.Generator or int, optional
            Source of the standard normal variates.
        variates : tuple of array_like, optional
            ``(u_init (m,), u_eps (p, n), u_eta (r, n))``; overrides `rng`.

        Returns
        -------
        alpha : ndarray, shape (m, n)
            Draw of the states.
        eps : ndarray, shape (p, n)
            Draw of the observation disturbances.
        eta : ndarray, shape (r, n)
            Draw of the state disturbances.

        See Also
        --------
        djs_simulation_smoother : The de Jong-Shephard method (DK §4.9.3).
        """
        ui, ue, un = self._variates(rng, variates)
        p, m, r, n = self.k_endog, self.k_states, self.k_posdef, self.nobs
        state, eps = np.empty((m, n), order="F"), np.empty((p, n), order="F")
        eta = np.empty((r, n), order="F")
        call(
            "ss_simulation_smoother",
            self._h,
            ptr(ui),
            ptr(ue),
            ptr(un),
            ptr(state),
            ptr(eps),
            ptr(eta),
        )
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
        """Compute the smoothed state by the fast state smoother (DK §4.6.2).

        A backward pass for :math:`r_t` only, then the forward recursion
        :math:`\\hat\\alpha_{t+1} = c_t + T_t \\hat\\alpha_t + R_t Q_t R_t' r_t`;
        no variances.

        Returns
        -------
        ndarray, shape (m, n)
            Smoothed state, equal to that of `smooth`.
        """
        a = np.empty(self._mn()[0], order="F")
        self._with_filter(lambda fh: call("ss_fast_smoother", self._h, fh, ptr(a)))
        return a

    def classical_smoother(self):
        """Compute the smoothed state by the classical fixed-interval smoother.

        The Rauch-Tung-Striebel form (DK §4.6.1), from the filtered states.
        Not available with diffuse states.

        Returns
        -------
        alphahat : ndarray, shape (m, n)
            Smoothed state.
        V : ndarray, shape (m, m, n)
            Its variance.
        """
        a, V = np.empty(self._mn()[0], order="F"), np.empty(self._mn()[1], order="F")
        self._with_filter(
            lambda fh: call("ss_classical_smoother", self._h, fh, ptr(a), ptr(V))
        )
        return a, V

    def two_filter_smoother(self):
        """Compute the smoothed state by the two-filter formula (DK §4.6.4).

        A backward information filter combined with the forward filter. Not
        available with diffuse states.

        Returns
        -------
        alphahat : ndarray, shape (m, n)
            Smoothed state.
        V : ndarray, shape (m, m, n)
            Its variance.
        """
        a, V = np.empty(self._mn()[0], order="F"), np.empty(self._mn()[1], order="F")
        self._with_filter(
            lambda fh: call("ss_two_filter_smoother", self._h, fh, ptr(a), ptr(V))
        )
        return a, V

    def whittle_smoother(self):
        """Compute the smoothed state from the Whittle relation (DK §4.6.3).

        Returns
        -------
        ndarray, shape (m, n)
            Smoothed state.

        Notes
        -----
        The recursion runs backwards through :math:`T_t^{-1}` and loses
        accuracy with the length of the series, as DK note. It serves to
        illustrate the relation; use `smooth` in practice.
        """
        a = np.empty(self._mn()[0], order="F")
        self._with_filter(lambda fh: call("ss_whittle_smoother", self._h, fh, ptr(a)))
        return a

    def fixed_point_smoother(self, t):
        """Track the estimate of one state as observations arrive (DK §4.4.6).

        Parameters
        ----------
        t : int
            The period (0-based).

        Returns
        -------
        mean : ndarray, shape (m, n - t)
            Column k is :math:`E(\\alpha_t | y_1, \\dots, y_{t+k})`.
        var : ndarray, shape (m, m, n - t)
            The corresponding variances.
        """
        k = self.nobs - t
        a = np.empty((self.k_states, k), order="F")
        V = np.empty((self.k_states, self.k_states, k), order="F")
        self._with_filter(
            lambda fh: call(
                "ss_fixed_point_smoother", self._h, fh, int(t) + 1, ptr(a), ptr(V)
            )
        )
        return a, V

    def fixed_lag_smoother(self, lag):
        """Estimate each state from the observations up to `lag` periods later.

        DK §4.4.6.

        Parameters
        ----------
        lag : int
            Number of later observations used.

        Returns
        -------
        alphahat : ndarray, shape (m, n)
            :math:`E(\\alpha_s | y_1, \\dots, y_{s+lag})`; NaN in the last
            `lag` periods.
        V : ndarray, shape (m, m, n)
            Its variance.
        """
        a, V = np.empty(self._mn()[0], order="F"), np.empty(self._mn()[1], order="F")
        self._with_filter(
            lambda fh: call(
                "ss_fixed_lag_smoother", self._h, fh, int(lag), ptr(a), ptr(V)
            )
        )
        return a, V

    def update_smoothed(self, alphahat, V, nobs_old):
        """Update smoothed estimates for newly arrived observations (DK §4.4.5).

        Parameters
        ----------
        alphahat : array_like, shape (m, nobs_old)
            Smoothed state given the first `nobs_old` observations.
        V : array_like, shape (m, m, nobs_old)
            Its variance.
        nobs_old : int
            Number of observations the estimates were based on.

        Returns
        -------
        alphahat : ndarray, shape (m, n)
            Smoothed state given all n observations.
        V : ndarray, shape (m, m, n)
            Its variance.
        """
        a = np.full(self._mn()[0], np.nan, order="F")
        W = np.full(self._mn()[1], np.nan, order="F")
        a[:, :nobs_old] = np.asarray(alphahat)[:, :nobs_old]
        W[:, :, :nobs_old] = np.asarray(V)[:, :, :nobs_old]
        self._with_filter(
            lambda fh: call(
                "ss_update_smoothed", self._h, fh, int(nobs_old), ptr(a), ptr(W)
            )
        )
        return a, W

    def smoothed_state_autocov(self):
        """Compute the lag-one covariances of the smoothed state (DK §4.7).

        Returns
        -------
        ndarray, shape (m, m, n - 1)
            Slice t is :math:`Cov(\\alpha_{t+1}, \\alpha_t | Y_n)`; NaN where
            the diffuse periods are involved.
        """
        out = np.empty((self.k_states, self.k_states, self.nobs - 1), order="F")
        self._with_filter(
            lambda fh, sh: call("ss_smoothed_state_autocov", self._h, fh, sh, ptr(out)),
            smoother=True,
        )
        return out

    def smoothed_state_cov_between(self, t, j):
        """Compute the covariance of two states given the data (DK §4.7).

        Parameters
        ----------
        t, j : int
            Periods (0-based), past the diffuse periods.

        Returns
        -------
        ndarray, shape (m, m)
            :math:`Cov(\\alpha_t, \\alpha_j | Y_n)`.
        """
        out = np.empty((self.k_states, self.k_states), order="F")
        self._with_filter(
            lambda fh, sh: call(
                "ss_smoothed_state_cov_between",
                self._h,
                fh,
                sh,
                int(t) + 1,
                int(j) + 1,
                ptr(out),
            ),
            smoother=True,
        )
        return out

    def filtered_state_weights(self):
        """Compute the weights of past observations in the filtered state.

        DK §4.8.2; intercepts and the initial mean aside.

        Returns
        -------
        Wa : ndarray, shape (m, p, n, n)
            ``Wa[:, :, t, j]`` is the weight of :math:`y_j` in :math:`a_t`.
        Watt : ndarray, shape (m, p, n, n)
            The same for :math:`a_{t|t}`.
        """
        m, p, n = self.k_states, self.k_endog, self.nobs
        Wa, Watt = np.empty((m, p, n, n), order="F"), np.empty((m, p, n, n), order="F")
        self._with_filter(
            lambda fh: call(
                "ss_filtered_state_weights", self._h, fh, ptr(Wa), ptr(Watt)
            )
        )
        return Wa, Watt

    def smoothed_state_weights(self):
        """Compute the weights of the data in the smoothed state (DK §4.8.3).

        The smoothed state is linear in the data, the state intercepts and
        the initial mean:
        :math:`\\hat\\alpha_t = \\sum_j W_{tj} (y_j - d_j) + \\sum_j C_{tj} c_j
        + A_t a_1`.

        Returns
        -------
        W : ndarray, shape (m, p, n, n)
            ``W[:, :, t, j]`` is :math:`W_{tj}`.
        C : ndarray, shape (m, m, n, n)
            ``C[:, :, t, j]`` is :math:`C_{tj}`.
        A : ndarray, shape (m, m, n)
            ``A[:, :, t]`` is :math:`A_t`.

        Notes
        -----
        Costs O(n) smoother runs.
        """
        m, p, n = self.k_states, self.k_endog, self.nobs
        W = np.empty((m, p, n, n), order="F")
        C = np.empty((m, m, n, n), order="F")
        A = np.empty((m, m, n), order="F")
        call("ss_smoothed_state_weights", self._h, ptr(W), ptr(C), ptr(A))
        return W, C, A

    def innovation_transition(self, t):
        """Compute :math:`L_t`, the transition of the state prediction error.

        :math:`a_{t+1} - \\alpha_{t+1} = L_t (a_t - \\alpha_t) + \\dots`
        with :math:`L_t = T_t - K_t Z_t` (DK §4.3).

        Parameters
        ----------
        t : int
            Period (0-based), past the diffuse periods.

        Returns
        -------
        ndarray, shape (m, m)
            :math:`L_t`.
        """
        out = np.empty((self.k_states, self.k_states), order="F")
        self._with_filter(
            lambda fh: call(
                "ss_innovation_transition", self._h, fh, int(t) + 1, ptr(out)
            )
        )
        return out

    def augmented(self):
        """Run the augmented Kalman filter and smoother (DK §5.7).

        Returns
        -------
        AugmentedResults
            The filter, the estimate of the diffuse elements and three log
            likelihoods.

        See Also
        --------
        marginal_likelihood : Option for the marginal log likelihood.
        """
        ah = ctypes.c_void_p()
        call("ss_augmented_filter", self._h, ctypes.byref(ah))
        try:
            k = ctypes.c_int()
            llf, llf_fixed, llf_marg = (
                ctypes.c_double(),
                ctypes.c_double(),
                ctypes.c_double(),
            )
            call(
                "ss_augmented_info",
                ah,
                ctypes.byref(k),
                ctypes.byref(llf),
                ctypes.byref(llf_fixed),
                ctypes.byref(llf_marg),
            )
            delta = np.empty(k.value)
            delta_cov = np.empty((k.value, k.value), order="F")
            if k.value > 0:
                call("ss_augmented_get", ah, 1, ptr(delta))
                call("ss_augmented_get", ah, 2, ptr(delta_cov))
            fh = ctypes.c_void_p()
            call("ss_augmented_filter_result", ah, ctypes.byref(fh))
            fres = self._filter_results(fh)
            a, V = (
                np.empty(self._mn()[0], order="F"),
                np.empty(self._mn()[1], order="F"),
            )
            call("ss_augmented_smoother", self._h, ah, ptr(a), ptr(V))
        finally:
            _free("ss_augmented_free", ah)
        return AugmentedResults(
            filter=fres,
            delta=delta,
            delta_cov=delta_cov,
            llf=llf.value,
            llf_fixed=llf_fixed.value,
            llf_marginal=llf_marg.value,
            smoothed_state=a,
            smoothed_state_cov=V,
        )

    def em(
        self, maxiter=500, tol=1e-8, diagonal_obs_cov=False, diagonal_state_cov=False
    ):
        """Estimate H and Q by the EM algorithm (DK §7.3.4).

        Updates ``obs_cov`` and ``state_cov`` in place. The E-step uses the
        smoothed disturbances and their variances; for missing elements
        these are the conditional moments given the observed elements.

        Parameters
        ----------
        maxiter : int, optional
            Maximum number of iterations.
        tol : float, optional
            Stop when the log likelihood improves by less than `tol`.
        diagonal_obs_cov, diagonal_state_cov : bool, optional
            Restrict H or Q to be diagonal.

        Returns
        -------
        llf : float
            Log likelihood at the start of the last iteration.
        niter : int
            Number of iterations.
        llf_path : ndarray, shape (niter,)
            Log likelihood at each iteration; never decreasing.

        Notes
        -----
        H and Q must be time-invariant, the initialization must not depend
        on them (not stationary), and there must be no burn-in.
        """
        llf, niter = ctypes.c_double(), ctypes.c_int()
        path = np.empty(int(maxiter))
        call(
            "ss_em",
            self._h,
            int(maxiter),
            float(tol),
            int(bool(diagonal_obs_cov)),
            int(bool(diagonal_state_cov)),
            ctypes.byref(llf),
            ctypes.byref(niter),
            ptr(path),
        )
        return llf.value, niter.value, path[: niter.value]

    def collapse(self):
        """Collapse the observations to the state dimension (DK §6.5).

        Each :math:`y_t` is replaced by its generalized least squares
        projection on :math:`\\alpha_t`, which gives the same filtered and
        smoothed states with p reduced to m (Jungbacker and Koopman 2015).

        Returns
        -------
        collapsed : Representation
            The collapsed model.
        adjust : ndarray, shape (n,)
            Log likelihood adjustment: ``self.loglike() ==
            collapsed.loglike() + adjust.sum()``.

        Raises
        ------
        StateSpaceError
            With code 2 if a period's :math:`Z_t` (observed rows) lacks full
            column rank.
        """
        h = ctypes.c_void_p()
        adj = np.empty(self.nobs)
        call("ss_collapse", self._h, ctypes.byref(h), ptr(adj))
        return Representation._wrap(h, owner=None), adj

    def add_state_restrictions(self, restriction, value):
        """Impose linear restrictions on the states (DK §6.6).

        The restrictions :math:`R^*_t \\alpha_t = r^*_t` enter as exact
        observations.

        Parameters
        ----------
        restriction : array_like, shape (q, m) or (q, m, n)
            :math:`R^*_t`.
        value : array_like, shape (q, n)
            :math:`r^*_t`; NaN where a restriction does not apply.

        Returns
        -------
        Representation
            The model with the restrictions as extra observations. Filter it
            with ``filter_method = FILTER_UNIVARIATE``: once a restricted
            direction has no state noise, :math:`F_t` is singular.
        """
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
        """Compute the standardized smoothed disturbances (DK §7.5).

        Large values point to outliers (observation disturbances) and
        structural breaks (state disturbances).

        Parameters
        ----------
        vector : bool, optional
            Standardize each period's disturbance vector by the Cholesky
            factor of its variance instead of element by element.

        Returns
        -------
        eps : ndarray, shape (p, n)
            Auxiliary residuals of the observation equation.
        eta : ndarray, shape (r, n)
            Auxiliary residuals of the state equation. NaN where the
            variance is zero.
        """
        e = np.empty((self.k_endog, self.nobs), order="F")
        u = np.empty((self.k_posdef, self.nobs), order="F")
        self._with_filter(
            lambda fh, sh: call(
                "ss_auxiliary_residuals", self._h, sh, int(bool(vector)), ptr(e), ptr(u)
            ),
            smoother=True,
        )
        return e, u

    def de_jong_penzer(self):
        """Compute the shock statistics of de Jong and Penzer (1998) (DK §7.5).

        Returns
        -------
        state : ndarray, shape (m, n)
            :math:`r_{i,t} / \\sqrt{N_{ii,t}}`.
        observation : ndarray, shape (p, n)
            :math:`e_{i,t} / \\sqrt{D_{ii,t}}`. NaN where undefined.
        """
        r = np.empty((self.k_states, self.nobs), order="F")
        e = np.empty((self.k_endog, self.nobs), order="F")
        self._with_filter(
            lambda fh, sh: call("ss_de_jong_penzer", self._h, fh, sh, ptr(r), ptr(e)),
            smoother=True,
        )
        return r, e

    def least_squares_residuals(self, start, stop):
        """Compute least squares residuals for regression effects (DK §6.2.4).

        Parameters
        ----------
        start, stop : int
            The regression coefficients are states ``start:stop`` (0-based);
            they must be diffuse and constant.

        Returns
        -------
        ndarray, shape (p, n)
            Residuals with the coefficients at their full-sample estimate.
        """
        out = np.empty((self.k_endog, self.nobs), order="F")
        call("ss_least_squares_residuals", self._h, int(start) + 1, int(stop), ptr(out))
        return out

    def r2_diffuse(self, series=0):
        """Compute the coefficient of determination against a random walk.

        :math:`R^2_D = 1 - SSE / \\sum_t (\\Delta y_t - \\overline{\\Delta y})^2`
        with SSE the sum of squared innovations after the diffuse periods
        (Harvey 1989; DK §7.4).

        Parameters
        ----------
        series : int, optional
            The series (0-based).

        Returns
        -------
        float
            :math:`R^2_D`; negative when the model forecasts worse than a
            random walk with drift.
        """
        r2 = ctypes.c_double()
        self._with_filter(
            lambda fh: call(
                "ss_r2_diffuse", self._h, fh, int(series) + 1, ctypes.byref(r2)
            )
        )
        return r2.value

    def djs_simulation_smoother(self, rng=None, variates=None):
        """Draw states and disturbances by the de Jong-Shephard method.

        DK §4.9.3: the disturbances are drawn backwards from their
        distribution given the data and the later draws.

        Parameters
        ----------
        rng : numpy.random.Generator or int, optional
            Source of the standard normal variates.
        variates : tuple of array_like, optional
            ``(u_init (m,), u_eps (p, n), u_eta (r, n))``; overrides `rng`.

        Returns
        -------
        alpha : ndarray, shape (m, n)
            Draw of the states.
        eps : ndarray, shape (p, n)
            Draw of the observation disturbances.
        eta : ndarray, shape (r, n)
            Draw of the state disturbances.

        See Also
        --------
        simulation_smoother : The mean-correction method (DK §4.9.2).
        """
        ui, ue, un = self._variates(rng, variates)
        p, m, r, n = self.k_endog, self.k_states, self.k_posdef, self.nobs
        eps, eta = np.empty((p, n), order="F"), np.empty((r, n), order="F")
        alpha = np.empty((m, n), order="F")
        self._with_filter(
            lambda fh: call(
                "ss_djs_simulation_smoother",
                self._h,
                fh,
                ptr(ue),
                ptr(un),
                ptr(ui),
                ptr(eps),
                ptr(eta),
                ptr(alpha),
            )
        )
        return alpha, eps, eta

    def steady_state(self):
        """Compute the steady state of the Kalman filter (DK §2.11, §4.3.4).

        The limit :math:`\\bar P` of the recursion for :math:`P_t`, by the
        structure-preserving doubling algorithm (Chu, Fan, Lin and Wang 2004),
        or by the recursion itself when H is singular.

        Returns
        -------
        P : ndarray, shape (m, m)
            :math:`\\bar P`.
        F : ndarray, shape (p, p)
            :math:`\\bar F = Z \\bar P Z' + H`, the prediction error variance
            (DK §7.4).

        Raises
        ------
        StateSpaceError
            With code 4 if a system matrix varies over time, or 7 if the
            recursion does not converge.

        Examples
        --------
        The local level model with signal-to-noise ratio q has
        :math:`\\bar P = x \\sigma^2_\\varepsilon`,
        :math:`x = (q + \\sqrt{q^2 + 4q}) / 2` (DK §2.11):

        >>> rep = ss.Representation(np.zeros(10), k_states=1)
        >>> rep["design"] = rep["transition"] = rep["selection"] = [[1.0]]
        >>> rep["obs_cov"], rep["state_cov"] = [[1.0]], [[2.0]]
        >>> P, F = rep.steady_state()
        >>> q = 2.0
        >>> float(P[0, 0]), float((q + np.sqrt(q**2 + 4 * q)) / 2)
        (2.732050807568877, 2.732050807568877)
        """
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
