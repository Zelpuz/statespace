``statespace_restrict``
=======================

Filtering and smoothing under linear restrictions on the state (DK §6.6).

The restrictions :math:`R^*_t \alpha_t = r^*_t` enter the observation
equation as exact observations:

.. math::

   y^+_t = \begin{bmatrix} y_t \\ r^*_t \end{bmatrix}, \quad
   Z^+_t = \begin{bmatrix} Z_t \\ R^*_t \end{bmatrix}, \quad
   H^+_t = \begin{bmatrix} H_t & 0 \\ 0 & 0 \end{bmatrix}.

A NaN in :math:`r^*_t` leaves that restriction inactive in period t. The
filtered and smoothed states satisfy the active restrictions exactly. The
log likelihood of the result includes the restriction rows, so estimate
parameters with the unrestricted model. Filter the result with
``FILTER_UNIVARIATE``: once a restricted direction has no state noise,
:math:`F_t` is singular.

``add_state_restrictions``
--------------------------

.. code-block:: fortran

   subroutine add_state_restrictions(rep, Rmat, rval, rrep, info)
     type(ssm_rep_t), intent(in) :: rep
     real(dp), intent(in) :: Rmat(:, :, :)     ! (q, m, 1|n)  R*_t
     real(dp), intent(in) :: rval(:, :)        ! (q, n)       r*_t
     type(ssm_rep_t), intent(out) :: rrep
     integer, intent(out) :: info
