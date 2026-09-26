``statespace_collapse``
=======================

Collapsing large observation vectors (DK §6.5; Jungbacker and Koopman
2015).

When p is large relative to m, each :math:`y_t` is replaced by its
generalized least squares projection on the state:

.. math::

   \bar y_t = (Z' H^{-1} Z)^{-1} Z' H^{-1} (y_t - d_t) = \alpha_t + \bar\varepsilon_t,
   \qquad Var(\bar\varepsilon_t) = (Z' H^{-1} Z)^{-1},

over the observed elements. The residual
:math:`e_t = y_t - d_t - Z \bar y_t` is independent of :math:`\bar y_t`
and of the states, so the collapsed model (observation :math:`\bar y_t`,
Z = I) gives the same filtered and smoothed states, and

.. math::

   \log L = \log L^* + \sum_t \Big[ -\tfrac{p_t - m}{2} \log 2\pi
     - \tfrac12 \log \tfrac{|H_t|}{|\bar H_t|} - \tfrac12 e_t' H_t^{-1} e_t \Big].

``collapse_observations``
-------------------------

.. code-block:: fortran

   subroutine collapse_observations(rep, crep, llf_adjust, info)
     type(ssm_rep_t), intent(in) :: rep
     type(ssm_rep_t), intent(out) :: crep
     real(dp), intent(out) :: llf_adjust(:)   ! (n)
     integer, intent(out) :: info

Build the collapsed representation and the per-period adjustment, so that
``loglike(rep) = loglike(crep) + sum(llf_adjust)``. Every period with
observations needs the observed block of :math:`H_t` nonsingular and the
observed rows of :math:`Z_t` of full column rank m; otherwise ``info`` is
``SS_ERR_NOT_PD``.
