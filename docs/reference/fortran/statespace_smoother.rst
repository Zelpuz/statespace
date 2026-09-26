``statespace_smoother``
=======================

The state and disturbance smoother (DK §4.4, §4.5), with the exact initial
smoother for diffuse periods (DK §5.3).

The backward recursion computes :math:`r_{t-1}` and :math:`N_{t-1}`
(DK eqs. 4.39, 4.43), then

.. math::

   \hat\alpha_t = a_t + P_t r_{t-1}, \qquad V_t = P_t - P_t N_{t-1} P_t,

and the smoothed disturbances :math:`\hat\varepsilon_t = H_t u_t` and
:math:`\hat\eta_t = Q_t R_t' r_t` (DK eq. 4.63) with their variances. In
diffuse periods it adds the recursions for :math:`r^{(1)}, N^{(1)},
N^{(2)}`, element by element or in the multivariate form, matching how
the filter processed the period.

``smoother_result_t``
---------------------

.. code-block:: fortran

   type :: smoother_result_t
     integer :: k_endog, k_states, k_posdef, nobs
     real(dp), allocatable :: alphahat(:, :)   ! (m, n)       E(alpha_t | Y_n)
     real(dp), allocatable :: V(:, :, :)       ! (m, m, n)    Var(alpha_t | Y_n)
     real(dp), allocatable :: r(:, :)          ! (m, 0:n)
     real(dp), allocatable :: N(:, :, :)       ! (m, m, 0:n)
     real(dp), allocatable :: epshat(:, :)     ! (p, n)       E(eps_t | Y_n)
     real(dp), allocatable :: epsvar(:, :, :)  ! (p, p, n)
     real(dp), allocatable :: etahat(:, :)     ! (r, n)       E(eta_t | Y_n)
     real(dp), allocatable :: etavar(:, :, :)  ! (r, r, n)
   end type

``r`` and ``N`` use DK's indexing: 0, ..., n with :math:`r_n = 0` and
:math:`N_n = 0`; in diffuse periods they hold :math:`r^{(0)}, N^{(0)}`.

For a missing element of :math:`y_t`, ``epshat`` and ``epsvar`` are the
conditional moments of that element of :math:`\varepsilon_t` given the
observed data, nonzero when :math:`H_t` has correlations; statsmodels
reports 0 and :math:`H_t`. EM needs the conditional moments (see
:doc:`../../design/missing_data`). The moments are in the original
coordinates of :math:`y_t` in every period.

``state_smoother``
------------------

.. code-block:: fortran

   subroutine state_smoother(rep, fres, sres, info)
     type(ssm_rep_t), intent(in) :: rep
     type(filter_result_t), intent(in) :: fres
     type(smoother_result_t), intent(out) :: sres
     integer, intent(out) :: info

Smooth a filter run of ``rep`` (from ``kalman_filter``). ``info`` is
``SS_ERR_DIM`` if ``fres`` does not match ``rep``.

``diffuse_mv_gains``
--------------------

.. code-block:: fortran

   subroutine diffuse_mv_gains(rep, fres, t, Zo, vo, F1, F2, L0, L1)
     type(ssm_rep_t), intent(in) :: rep
     type(filter_result_t), intent(in) :: fres
     integer, intent(in) :: t
     real(dp), allocatable, intent(out) :: Zo(:, :), vo(:), F1(:, :), F2(:, :), L0(:, :), L1(:, :)

The quantities of the multivariate exact initial smoother at a period with
:math:`F_\infty` nonsingular over the observed elements o:
:math:`F^{(1)} = F_\infty^{-1}`, :math:`F^{(2)} = -F^{(1)} F_* F^{(1)}`,
:math:`L^{(0)} = T - K^{(0)} Z_o` and :math:`L^{(1)} = -K^{(1)} Z_o`
(DK §5.2-5.3). Used by the smoothing extras (:doc:`statespace_smoothing`).
