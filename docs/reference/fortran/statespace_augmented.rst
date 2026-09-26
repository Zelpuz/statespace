``statespace_augmented``
========================

The augmented Kalman filter and smoother (DK §5.7): an alternative to the
exact initial filter for diffuse states, and the basis of regression
estimation (DK §6.2).

With :math:`\alpha_1 = a + A\delta + R_0\eta_0` and :math:`\delta` unknown,
the filter runs once with :math:`\delta = 0`, giving :math:`v^*_t, F_t, K_t,
P_t`, and carries the effect of :math:`\delta` as extra columns:

.. math::

   A_1 = A, \quad V^A_t = -Z_t A_t, \quad A_{t+1} = T_t A_t + K_t V^A_t, \\
   s = \sum_t V^{A\prime}_t F_t^{-1} v^*_t, \quad
   S = \sum_t V^{A\prime}_t F_t^{-1} V^A_t,

so that :math:`v_t(\delta) = v^*_t + V^A_t \delta`. The generalized least
squares estimate is :math:`\hat\delta = -S^{-1} s` with variance
:math:`S^{-1}`, and the diffuse log likelihood (DK §7.2.2-7.2.3) is

.. math::

   \log L_d = \log L^* + \tfrac12 s' S^{-1} s - \tfrac12 \log|S|.

The smoother adds the uncertainty about :math:`\delta` to :math:`V_t` by the
law of total variance. The filter needs :math:`F_t` nonsingular over the
observed elements; the exact initial filter does not.

``augmented_result_t``
----------------------

.. code-block:: fortran

   type :: augmented_result_t
     integer :: k_diffuse                       ! k
     type(filter_result_t) :: filter            ! the delta = 0 filter
     real(dp), allocatable :: A(:, :, :)        ! (m, k, n+1)
     real(dp), allocatable :: VA(:, :, :)       ! (p, k, n)
     real(dp), allocatable :: s(:), S_mat(:, :)
     real(dp), allocatable :: delta(:), delta_cov(:, :)
     real(dp) :: llf, llf_fixed, llf_marginal
   end type

``llf``
   Diffuse log likelihood; equal to that of the exact initial filter.
``llf_fixed``
   Log likelihood with :math:`\delta` fixed but unknown, concentrated at
   :math:`\hat\delta` (DK §7.2.4): ``llf`` + :math:`\tfrac12 \log|S|`.
``llf_marginal``
   Marginal log likelihood (DK §7.2.6); see ``marginal_correction`` in
   :doc:`statespace_filter`.

``augmented_filter``
--------------------

.. code-block:: fortran

   subroutine augmented_filter(rep, res, info)
     type(ssm_rep_t), intent(in) :: rep
     type(augmented_result_t), intent(out) :: res
     integer, intent(out) :: info

``info`` is ``SS_ERR_NOT_PD`` if S is singular: the data do not identify
:math:`\delta`.

``augmented_smoother``
----------------------

.. code-block:: fortran

   subroutine augmented_smoother(rep, res, alphahat, V, info)
     real(dp), intent(out) :: alphahat(:, :), V(:, :, :)   ! (m, n), (m, m, n)

.. math::

   \hat\alpha_t = a^*_t + A_t \hat\delta + P_t (r^*_{t-1} + R^A_{t-1} \hat\delta), \quad
   V_t = P_t - P_t N_{t-1} P_t + B_t S^{-1} B_t', \quad
   B_t = A_t + P_t R^A_{t-1}.
