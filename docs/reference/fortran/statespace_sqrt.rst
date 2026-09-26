``statespace_sqrt``
===================

The square root filter and smoother (DK §6.3).

The filter propagates a square root :math:`\tilde P_t` of :math:`P_t`, so
:math:`P_t` stays positive semi-definite by construction. Each step
triangularizes the pre-array

.. math::

   \begin{bmatrix} Z_o \tilde P & \tilde H_{oo} & 0 \\
                   T \tilde P & 0 & R \tilde Q \end{bmatrix}
   \;\to\;
   \begin{bmatrix} \tilde F & 0 & 0 \\ \tilde K & \tilde P_{t+1} & 0 \end{bmatrix}

by an orthogonal transformation (an LQ factorization by Householder QR,
where DK use Givens rotations). Then :math:`F = \tilde F \tilde F'`,
:math:`K = \tilde K \tilde F^{-1}` and
:math:`P_{t+1} = \tilde P_{t+1} \tilde P_{t+1}'`. The smoother does the same
for :math:`N_{t-1}`.

The output uses the types of ``kalman_filter`` and ``state_smoother``, so
the rest of the library works on it. Diffuse states are not supported; use
the exact initial filter or an approximate diffuse initialization.

``sqrt_kalman_filter``
----------------------

.. code-block:: fortran

   subroutine sqrt_kalman_filter(rep, res, info, Pchol)
     type(ssm_rep_t), intent(in) :: rep
     type(filter_result_t), intent(out) :: res
     integer, intent(out) :: info
     real(dp), allocatable, intent(out), optional :: Pchol(:, :, :)   ! (m, m, n+1)

``Pchol`` receives the square roots :math:`\tilde P_t`.

``sqrt_state_smoother``
-----------------------

.. code-block:: fortran

   subroutine sqrt_state_smoother(rep, fres, sres, info)
     type(filter_result_t), intent(in) :: fres
     type(smoother_result_t), intent(out) :: sres

:math:`r_t` and :math:`\hat\alpha_t` as in ``state_smoother``, with
:math:`N_{t-1}` propagated through its square root and
:math:`V_t = P_t - (P_t \tilde N_{t-1})(P_t \tilde N_{t-1})'`.
