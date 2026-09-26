``statespace_dense``
====================

The matrix formulation of the model (DK §4.13): the whole sample as one
Gaussian vector. With :math:`\alpha` the stacked states and :math:`Y_n`
the stacked observed elements,

.. math::

   E(Y_n) = \mu = d + Z a^*, \quad Var(Y_n) = \Omega = Z V^* Z' + H, \\
   \log L = -\tfrac12 \big(N \log 2\pi + \log|\Omega|
            + (Y_n - \mu)' \Omega^{-1} (Y_n - \mu)\big), \\
   E(\alpha | Y_n) = a^* + V^* Z' \Omega^{-1} (Y_n - \mu), \quad
   Var(\alpha | Y_n) = V^* - V^* Z' \Omega^{-1} Z V^*.

The cost is :math:`O((nm)^3)`. It is a reference for small problems and a
check on the recursive algorithms, which the tests use.

``dense_loglike_smooth``
------------------------

.. code-block:: fortran

   subroutine dense_loglike_smooth(rep, llf, alphahat, V, info)
     type(ssm_rep_t), intent(in) :: rep
     real(dp), intent(out) :: llf
     real(dp), intent(out) :: alphahat(:, :)   ! (m, n)
     real(dp), intent(out) :: V(:, :, :)       ! (m, m, n)
     integer, intent(out) :: info

Known or approximate diffuse initialization only; ``info`` is
``SS_ERR_UNSUPPORTED`` with diffuse states.
