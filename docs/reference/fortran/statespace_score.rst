``statespace_score``
====================

The analytic score (DK §7.3.3) from one smoother pass.

By the Fisher identity, :math:`\partial \log L / \partial\psi =
E[\partial \log p(Y, \alpha; \psi) / \partial\psi \mid Y]` (DK eqs.
7.13-7.15). For parameters in :math:`H_t`, :math:`R_t Q_t R_t'` and the
:math:`P_*` part of the initialization this gives (DK eq. 7.16; Koopman and
Shephard 1992)

.. math::

   \frac{\partial \log L}{\partial \psi_i}
   = \tfrac12 \sum_t \mathrm{tr}(\dot H_t G^\varepsilon_t)
   + \tfrac12 \sum_{t<n} \mathrm{tr}\big(\dot W_t (r_t r_t' - N_t)\big)
   + \tfrac12 \mathrm{tr}(\dot P_* G^1),

with :math:`W_t = R_t Q_t R_t'`,
:math:`G^\varepsilon_t = H_t^{-1}(\hat\varepsilon_t \hat\varepsilon_t' +
Var(\varepsilon_t | Y) - H_t) H_t^{-1}` and
:math:`G^1 = P_*^+ (V_1 + (\hat\alpha_1 - a_1)(\hat\alpha_1 - a_1)' - P_*)
P_*^+`.

The derivatives :math:`\dot H, \dot W, \dot P_*` are central differences of
what the model's ``update`` sets, so no extra code is needed. The identity
holds for the diffuse log likelihood (DK §7.3.5) and for a concentrated
scale, but not with a burn-in. With missing data it needs the conditional
moments of the missing elements of :math:`\varepsilon_t`, which the
smoother provides.

``analytic_score``
------------------

.. code-block:: fortran

   subroutine analytic_score(model, params, score, llf, info, analytic)
     class(ssm_model_t), intent(inout) :: model
     real(dp), intent(in) :: params(:)
     real(dp), intent(out) :: score(:)
     real(dp), intent(out) :: llf
     integer, intent(out) :: info
     logical, intent(out), optional :: analytic(:)

The score at constrained ``params`` and the log likelihood from the same
filter run. A parameter that moves Z, T, c, d, :math:`a_1` or
:math:`P_\infty`, or a singular :math:`H_t`, is outside the formula: without
``analytic`` the call returns ``SS_ERR_UNSUPPORTED``; with it, such
parameters are marked false (score 0) and the others are computed, for the
hybrid gradient of ``fit``. An indefinite H or Q returns
``SS_ERR_NOT_PD``.
