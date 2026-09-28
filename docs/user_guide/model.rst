The state space model
=====================

The library works with the linear Gaussian state space model (DK §3.1)

.. math::

   y_t &= d_t + Z_t \alpha_t + \varepsilon_t, &
   \varepsilon_t &\sim N(0, H_t), \\
   \alpha_{t+1} &= c_t + T_t \alpha_t + R_t \eta_t, &
   \eta_t &\sim N(0, Q_t), \qquad t = 1, \dots, n,

with :math:`y_t` the p observations at period t, :math:`\alpha_t` the m
states and :math:`\eta_t` the r state disturbances. The disturbances are
independent across periods and of each other and of :math:`\alpha_1`.

The representation
------------------

A :class:`~ssfortran.Representation` (Fortran ``ssm_rep_t``) holds the data,
the system matrices and the distribution of :math:`\alpha_1`: the model at
given parameter values. Models build and own one; a
:class:`~ssfortran.MappedModel` or :class:`~ssfortran.MLEModel` sets its
matrices with ``mod[name] = value``. Working with a Representation directly
serves to run the filters, smoothers and other algorithms at fixed values
(:doc:`filtering`). The names follow statsmodels:

=================== ======== ========== ===================
Python              Fortran  symbol     shape
=================== ======== ========== ===================
``design``          ``Z``    Z          (p, m) or (p, m, n)
``obs_cov``         ``H``    H          (p, p) or (p, p, n)
``transition``      ``T``    T          (m, m) or (m, m, n)
``selection``       ``R``    R          (m, r) or (m, r, n)
``state_cov``       ``Q``    Q          (r, r) or (r, r, n)
``state_intercept`` ``c``    c          (m,) or (m, n)
``obs_intercept``   ``d``    d          (p,) or (p, n)
=================== ======== ========== ===================

Time is the last axis. A matrix without it is time-invariant; one with a
time axis of length n varies over time, and each can be chosen separately.
In Fortran the time dimension is always present, with length 1 or n.

>>> rng = np.random.default_rng(0)
>>> rep = ss.Representation(rng.standard_normal((50, 2)), k_states=3, k_posdef=2)
>>> rep.nobs, rep.k_endog, rep.k_states, rep.k_posdef
(50, 2, 3, 2)
>>> rep["transition"] = np.diag([0.9, 0.5, 1.0])
>>> rep["design"] = rng.standard_normal((2, 3, 50))      # time-varying Z
>>> rep["design"].shape, rep["transition"].shape
((2, 3, 50), (3, 3))

Observations are passed time first, (n, p), as in statsmodels, and stored
as (p, n).

Missing observations
--------------------

NaN marks a missing value, and any subset of a period's elements may be
missing (DK §4.10). The filter uses the observed elements only; a period
with all elements missing is a prediction step. Forecasting is filtering
with missing observations at the end.

For a missing element, the smoothed observation disturbance is its
conditional mean given the observed elements of the same period, which is
nonzero when :math:`H_t` has correlations (see :doc:`../design/missing_data`).

Dimensions and cost
-------------------

A filter step costs :math:`O(m^3 + p^3 + m^2 p)`. When p is large relative
to m, collapsing the observations (DK §6.5,
:meth:`~ssfortran.Representation.collapse`) reduces p to m without changing
the states or, up to a known term, the likelihood. The univariate treatment
(DK §6.4) avoids inverting :math:`F_t` altogether.
