``statespace_rep``
==================

The representation of a linear Gaussian state space model (DK eq. 3.1):

.. math::

   y_t &= d_t + Z_t \alpha_t + \varepsilon_t, &
   \varepsilon_t &\sim N(0, H_t), \\
   \alpha_{t+1} &= c_t + T_t \alpha_t + R_t \eta_t, &
   \eta_t &\sim N(0, Q_t),

with :math:`y_t` of length p, :math:`\alpha_t` of length m,
:math:`\eta_t` of length r, and t = 1, ..., n, together with the
distribution of :math:`\alpha_1`.

``ssm_rep_t`` is plain data: its components are public, and the
algorithms of the other modules take it as an argument. Missing
observations are NaN in ``y`` and may be any subset of a period's elements
(DK §4.10).

Example
-------

The local level model for the Nile data, with exact diffuse
initialization (from ``example/nile_local_level.f90``):

.. code-block:: fortran

   type(ssm_rep_t) :: rep
   type(filter_result_t) :: fres
   integer :: info

   rep = ssm_rep(y, m=1, r=1)        ! y(1, n)
   rep%Z = 1.0_dp
   rep%T = 1.0_dp
   rep%H = 15099.0_dp
   rep%Q = 1469.1_dp
   call rep%initialize_diffuse()
   call kalman_filter(rep, fres, info)

``ssm_rep_t``
-------------

.. code-block:: fortran

   type :: ssm_rep_t
     integer :: k_endog, k_states, k_posdef, nobs         ! p, m, r, n
     real(dp), allocatable :: y(:, :)                      ! (p, n)
     real(dp), allocatable :: Z(:, :, :), H(:, :, :)       ! (p, m, 1|n), (p, p, 1|n)
     real(dp), allocatable :: T(:, :, :), R(:, :, :)       ! (m, m, 1|n), (m, r, 1|n)
     real(dp), allocatable :: Q(:, :, :)                   ! (r, r, 1|n)
     real(dp), allocatable :: c(:, :), d(:, :)             ! (m, 1|n), (p, 1|n)
     integer :: filter_method = FILTER_CONVENTIONAL
     integer :: diffuse_method = DIFFUSE_UNIVARIATE
     real(dp) :: tol_diffuse = 1.0e-10_dp
     real(dp) :: tol_steady = 1.0e-15_dp
     integer :: loglikelihood_burn = 0
     logical :: marginal_likelihood = .false.
     ! initialization: blk_first, blk_last, blk_kind, a1, Pstar1, Pinf1
   end type

A time dimension of 1 makes a matrix time-invariant; assigning an array
with time dimension n (``rep%Z = Z_tv``) makes it time-varying.

``filter_method``
   ``FILTER_CONVENTIONAL`` processes :math:`y_t` as a vector (DK §4.3);
   ``FILTER_UNIVARIATE`` one element at a time (DK §6.4) in every period.
``diffuse_method``
   ``DIFFUSE_UNIVARIATE`` (default) filters the diffuse periods element by
   element (DK §5.2.5); ``DIFFUSE_MULTIVARIATE`` uses the multivariate
   exact initial filter (DK §5.2) where :math:`F_\infty` is zero or
   nonsingular and the univariate step where it is singular but nonzero.
``tol_diffuse``
   Threshold on :math:`F_\infty` and :math:`\|P_\infty\|_F^2` in the
   diffuse periods, as in statsmodels.
``tol_steady``
   In a time-invariant model the conventional filter holds P, F and K fixed
   once :math:`\|P_{t+1} - P_t\|_F \le` ``tol_steady`` :math:`\|P_{t+1}\|_F`
   (DK §4.3.4), for as long as :math:`y_t` is fully observed. A negative
   value turns this off. See :doc:`../../design/steady_state`.
``loglikelihood_burn``
   Number of leading periods left out of the log likelihood.
``marginal_likelihood``
   Report the marginal log likelihood of Francke, Koopman and de Vos (2010)
   instead of the diffuse log likelihood (DK §7.2.6). The two differ by a
   term in Z, T and the diffuse directions only.

Constants
---------

Initialization kinds, for ``initialize_block``: ``INIT_NONE`` (0),
``INIT_KNOWN`` (1), ``INIT_APPROX_DIFFUSE`` (2), ``INIT_STATIONARY`` (3),
``INIT_DIFFUSE`` (4), ``INIT_GENERAL`` (5).

Filter methods: ``FILTER_CONVENTIONAL`` (0), ``FILTER_UNIVARIATE`` (1).
Diffuse methods: ``DIFFUSE_UNIVARIATE`` (0), ``DIFFUSE_MULTIVARIATE`` (1).

``ssm_rep``
-----------

.. code-block:: fortran

   function ssm_rep(y, m, r) result(rep)
     real(dp), intent(in) :: y(:, :)
     integer, intent(in) :: m, r
     type(ssm_rep_t) :: rep

Create a time-invariant representation for data ``y(p, n)`` with m states
and r state disturbances. All system matrices start at zero except R, which
starts as the leading m × r identity. The initialization must be set before
filtering.

``tidx``
--------

.. code-block:: fortran

   pure integer function tidx(nt, t)

The slice of period t in an array whose time dimension has length ``nt``:
``min(t, nt)``.

Initialization
--------------

The initial state is :math:`\alpha_1 = a + A\delta + R_0 \eta_0` with
:math:`\delta` diffuse (DK eq. 5.2), that is
:math:`\alpha_1 \sim N(a, P_* + \kappa P_\infty)` with
:math:`\kappa \to \infty`. It is set for the whole state vector, or block
by block, for example a diffuse trend beside a stationary ARMA block
(DK §5.6).

.. code-block:: fortran

   subroutine initialize_known(self, a1, P1)
   subroutine initialize_approximate_diffuse(self, kappa, a1)   ! both optional
   subroutine initialize_stationary(self)
   subroutine initialize_diffuse(self, a1)                      ! a1 optional
   subroutine initialize_general(self, a1, Pstar, Pinf)
   subroutine initialize_block(self, first, last, kind, a1, P1, kappa)

``initialize_known``
   :math:`\alpha_1 \sim N(a_1, P_1)`.
``initialize_approximate_diffuse``
   :math:`\alpha_1 \sim N(a_1, \kappa I)`, :math:`\kappa` = 1e6 and
   :math:`a_1 = 0` unless given.
``initialize_stationary``
   The unconditional distribution: :math:`a_1 = (I - T)^{-1} c` and
   :math:`P_1 = T P_1 T' + R Q R'`, from the matrices at t = 1 and
   recomputed at each filter run, so it follows parameter updates.
``initialize_diffuse``
   Exact diffuse: :math:`P_\infty = I`, :math:`P_* = 0`. The diffuse
   log likelihood (DK §7.2.2) accounts for the diffuse states, so
   ``loglikelihood_burn`` should be 0.
``initialize_general``
   DK §5.1 in full, with given :math:`a_1, P_*, P_\infty`.
``initialize_block``
   States ``first..last`` get ``kind``. Blocks must not overlap and must
   together cover the state vector. A stationary block must not depend on
   states outside it through T.

``initial_state``
-----------------

.. code-block:: fortran

   subroutine initial_state(self, a1, Pstar, Pinf, info)
     real(dp), intent(out) :: a1(:), Pstar(:, :), Pinf(:, :)
     integer, intent(out) :: info

The mean and the two variance parts of :math:`\alpha_1` implied by the
initialization. ``info`` is ``SS_ERR_INIT`` if a state is not initialized
and ``SS_ERR_NOT_STATIONARY`` if a stationary block has a unit root.

``k_diffuse``
-------------

.. code-block:: fortran

   integer function k_diffuse(self)

Number of diffuse elements of :math:`\alpha_1`, the rank of
:math:`P_\infty`; counted as parameters in AIC and BIC (DK §7.4).

``scale_by``
------------

.. code-block:: fortran

   subroutine scale_by(self, s)
     real(dp), intent(in) :: s

Multiply H, Q and :math:`P_*` by s: puts a model whose scale was
concentrated out (``loglike_concentrated``) on the data's scale.
:math:`P_\infty` is not affected.

``validate``
------------

.. code-block:: fortran

   subroutine validate(self, info)
     integer, intent(out) :: info

Check that the dimensions are consistent, that each time dimension is 1 or
n, and that the model is initialized. The filters call it first.
