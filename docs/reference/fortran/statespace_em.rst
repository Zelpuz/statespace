``statespace_em``
=================

The EM algorithm for the variance matrices H and Q (DK §7.3.4).

With the disturbances as complete data, the M-step for time-invariant H and
Q is closed form:

.. math::

   H_{new} = \frac1n \sum_{t=1}^n \big(\hat\varepsilon_t \hat\varepsilon_t'
             + Var(\varepsilon_t | Y)\big), \qquad
   Q_{new} = \frac1{n-1} \sum_{t=1}^{n-1} \big(\hat\eta_t \hat\eta_t'
             + Var(\eta_t | Y)\big),

since :math:`\eta_n` does not affect the data. Missing elements of
:math:`y_t` enter through their conditional moments. No iteration decreases
the log likelihood, and a fixed point is a stationary point of it. EM needs
an initialization that does not depend on H or Q (no stationary blocks) and
no burn-in. It converges slowly near the optimum; a few EM steps are a good
start for ``fit``.

``em_step``
-----------

.. code-block:: fortran

   subroutine em_step(rep, llf, info, diagonal_H, diagonal_Q)
     type(ssm_rep_t), intent(inout) :: rep
     real(dp), intent(out) :: llf
     integer, intent(out) :: info
     logical, intent(in), optional :: diagonal_H, diagonal_Q

One iteration: the E-step at the current H and Q, the M-step replacing them.
``llf`` is the log likelihood before the update. ``diagonal_H`` or
``diagonal_Q`` keep only the diagonal, the restricted maximum for diagonal
matrices.

``em_variances``
----------------

.. code-block:: fortran

   subroutine em_variances(rep, maxiter, tol, llf, niter, info, diagonal_H, diagonal_Q, llf_path)
     integer, intent(in) :: maxiter
     real(dp), intent(in) :: tol
     real(dp), intent(out) :: llf
     integer, intent(out) :: niter, info
     real(dp), allocatable, intent(out), optional :: llf_path(:)

Iterate ``em_step`` until the log likelihood improves by less than ``tol``
or for ``maxiter`` iterations. ``llf`` is the log likelihood at the start of
the last iteration; ``llf_path`` holds it for every iteration.
