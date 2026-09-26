``statespace_structural``
=========================

The structural components of DK §3.2-3.3, §3.8 and §3.9.

Each component applies to all p series of the model at once, as a
seemingly unrelated time series equations (SUTSE) model (DK §3.3): each
series has its own states, and the disturbances are correlated across
series through a p × p variance matrix of type

``COV_DIAGONAL`` (1)
   one variance per series;
``COV_FULL`` (2)
   a full matrix, estimated through its Cholesky factor with log diagonal;
``COV_NONE`` (0)
   no disturbance: a fixed component.

Within a component the states are ordered series by series, and within a
series in DK's order. Variances use :math:`\sigma^2 = \exp(2\psi)`
(DK §7.3.2); parameter names follow ``sigma2.<component>``, with the series
number appended when p > 1 and ``cov.<component>.i.j`` for full matrices.

Components
----------

.. code-block:: fortran

   type, extends(component_t) :: irregular_t
     integer :: cov = COV_DIAGONAL
     real(dp), allocatable :: weights(:)       ! (n), optional
   end type

The irregular :math:`\varepsilon_t \sim N(0, \Sigma_\varepsilon)`; adds to
H and has no states. With ``weights``, :math:`H_t = w_t
\Sigma_\varepsilon`, for unequally spaced or aggregated data (DK §3.8).
Always at the observations.

.. code-block:: fortran

   type, extends(component_t) :: level_t
     integer :: cov = COV_DIAGONAL
   end type

The random walk level :math:`\mu_{t+1} = \mu_t + \xi_t` (DK ch. 2).

.. code-block:: fortran

   type, extends(component_t) :: trend_t
     integer :: cov_level = COV_DIAGONAL, cov_slope = COV_DIAGONAL
   end type

The local linear trend :math:`\mu_{t+1} = \mu_t + \nu_t + \xi_t`,
:math:`\nu_{t+1} = \nu_t + \zeta_t` (DK §3.2.1). ``cov_level = COV_NONE``
gives the smooth trend (integrated random walk); both ``COV_NONE`` a
deterministic linear trend.

.. code-block:: fortran

   type, extends(component_t) :: seasonal_t
     integer :: period = 4
     integer :: form = SEASONAL_DUMMY
     integer :: cov = COV_DIAGONAL
   end type

The seasonal of period s (DK §3.2.2): ``SEASONAL_DUMMY`` (1; s - 1 states,
DK eq. 3.3), ``SEASONAL_TRIG`` (2; s - 1 states, all harmonics with one
variance, DK eqs. 3.7-3.8) or ``SEASONAL_HS`` (3; Harrison-Stevens, s
states). ``COV_NONE`` gives a fixed seasonal.

.. code-block:: fortran

   type, extends(component_t) :: cycle_t
     integer :: cov = COV_DIAGONAL
     logical :: damped = .true.
     real(dp) :: period_min = 2.0_dp
     real(dp) :: period_max = 0.0_dp          ! 0: the number of observations
   end type

The cycle (DK §3.2.4, eq. 3.13): :math:`(c_{t+1}, c^*_{t+1})' = \rho
C(\lambda) (c_t, c^*_t)' + w_t` with :math:`C(\lambda)` the rotation by
:math:`\lambda` and :math:`Var(w_t) = \Sigma \otimes I_2`. Parameters: the
covariance, then :math:`\lambda`, bounded by the period bounds, then
:math:`\rho \in (0, 1)` if damped. A damped cycle is stationary and starts
from its unconditional distribution (DK §3.2.4); statsmodels starts it
diffuse (see :doc:`../../design/cycle_initialization`).

.. code-block:: fortran

   type, extends(component_t) :: regression_t
     real(dp), allocatable :: x(:, :)              ! (n, k_x)
     logical, allocatable :: random_walk(:)        ! (k_x), default all .false.
     integer :: series = 1
   end type

Regression effects :math:`x_t'\beta_t` for one series (DK §3.2.5, §3.6).
The coefficients are diffuse states, fixed or random walks (DK eq. 3.15)
with estimated variances named ``sigma2.beta.j``. With loadings, set
``at_observations`` for effects on an observed series, such as the
intercepts of the common levels model (DK §3.3.2).

.. code-block:: fortran

   type, extends(component_t) :: continuous_level_t
     real(dp), allocatable :: times(:)             ! (n)
   end type
   type, extends(component_t) :: continuous_trend_t
     real(dp), allocatable :: times(:)             ! (n)
   end type

Continuous-time components observed at times :math:`t_1 \le \dots \le
t_n` (DK §3.8). The level has :math:`Var(\eta_i) = \sigma^2 \delta_i`,
:math:`\delta_i = t_{i+1} - t_i`. The smooth trend :math:`d\nu = \sigma
dw`, :math:`d\mu = \nu\,dt` gives (DK eqs. 3.41-3.42)

.. math::

   T_i = \begin{bmatrix} 1 & \delta_i \\ 0 & 1 \end{bmatrix}, \quad
   Q_i = \sigma^2 \delta_i \begin{bmatrix} \delta_i^2/3 & \delta_i/2 \\
                                          \delta_i/2 & 1 \end{bmatrix};

with diffuse initial states and an irregular of variance
:math:`\sigma_\varepsilon^2`, its smoothed level is the cubic smoothing
spline with :math:`\lambda = \sigma_\varepsilon^2 / \sigma^2` (DK §3.9.2;
Wahba 1978). Both are univariate; repeated times are allowed.

Interventions
-------------

.. code-block:: fortran

   pure function step_intervention(n, tau) result(w)    ! 0 before tau, 1 from tau
   pure function pulse_intervention(n, tau) result(w)   ! 1 at tau, 0 elsewhere
   pure function slope_intervention(n, tau) result(w)   ! 0 before tau, 1 + t - tau from tau

Regressors for intervention effects (DK §3.2.5), to pass to
``regression_t``.

Covariance transforms
---------------------

.. code-block:: fortran

   pure function cov_constrain(cov, p, x) result(c)
   pure function cov_unconstrain(cov, p, c) result(x)

Map between the unconstrained values and the covariance parameters: for
``COV_DIAGONAL`` the variances, :math:`\exp(2x)`; for ``COV_FULL`` the lower
triangle of :math:`\Sigma = L L'` by columns, with L lower triangular and
diagonal :math:`\exp(x)`. Also used by :doc:`statespace_mapped`.
