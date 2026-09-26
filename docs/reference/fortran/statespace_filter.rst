``statespace_filter``
=====================

The Kalman filter (DK §4.3), the exact initial Kalman filter for diffuse
states (DK §5.2), the log likelihood (DK §7.2) and the steady state
(DK §4.3.4).

Each period is processed in one of three ways, recorded in
``res%method(t)``:

``METHOD_CONVENTIONAL`` (0)
   The observation vector at once (DK eq. 4.24). In a diffuse period this
   is the case :math:`F_\infty = 0`.
``METHOD_UNIVARIATE`` (1)
   One observed element at a time (DK §6.4). When the observed block of
   :math:`H_t` is not diagonal it is first factorized as :math:`L D L'`
   and the observation equation transformed (DK §6.4.3). Diffuse periods
   use this form of the exact initial filter by default.
``METHOD_DIFFUSE_MV`` (2)
   The multivariate exact initial update (DK §5.2), with
   ``diffuse_method = DIFFUSE_MULTIVARIATE``.

Missing observations are handled as in DK §4.10: only the observed elements
of :math:`y_t` enter the update. In conventional periods :math:`F_t^{-1}` is
stored zero-padded over the missing elements, which keeps the full-size
smoother formulas exact.

The output is in the original coordinates of :math:`y_t` in every period;
the per-element quantities of univariate periods are kept separately for the
smoother (see :doc:`../../design/output_coordinates`).

``filter_result_t``
-------------------

.. code-block:: fortran

   type :: filter_result_t
     integer :: k_endog, k_states, nobs
     integer :: nobs_diffuse          ! DK's d
     integer :: t_steady              ! first steady-state period, 0 if none
     integer :: k_diffuse             ! rank of P_inf,1
     integer, allocatable :: method(:)             ! (n)
     real(dp), allocatable :: a(:, :), P(:, :, :)  ! (m, n+1), (m, m, n+1)
     real(dp), allocatable :: Pinf(:, :, :)        ! (m, m, n+1)
     real(dp), allocatable :: att(:, :), Ptt(:, :, :)
     real(dp), allocatable :: yhat(:, :), v(:, :)  ! (p, n)
     real(dp), allocatable :: F(:, :, :), Finf(:, :, :), Finv(:, :, :)
     real(dp), allocatable :: K(:, :, :)           ! (m, p, n)
     real(dp), allocatable :: llf_obs(:)           ! (n)
     real(dp) :: llf
     ! uv_n, uv_idx, uv_Z, uv_sig2, uv_v, uv_Fstar, uv_Finf, uv_Mstar,
     ! uv_Minf: per-element quantities of univariate periods
   end type

``a``, ``P``
   Predicted state :math:`a_t` and variance :math:`P_t` (the :math:`P_*`
   part in diffuse periods) for t = 1, ..., n+1.
``Pinf``
   :math:`P_{\infty,t}`; zero after the diffuse periods.
``att``, ``Ptt``
   Filtered :math:`a_{t|t}` and :math:`P_{t|t}` (DK eq. 4.24).
``yhat``, ``v``
   One-step prediction :math:`d_t + Z_t a_t` and innovation :math:`v_t`
   (NaN where missing).
``F``, ``Finf``
   :math:`F_t = Z_t P_t Z_t' + H_t` and :math:`Z_t P_{\infty,t} Z_t'`, over
   all elements, also where some are missing.
``Finv``, ``K``
   :math:`F_t^{-1}` over the observed elements (zero-padded) and
   :math:`K_t = T_t P_t Z_t' F_t^{-1}`; NaN in periods without a
   conventional update.
``llf_obs``, ``llf``
   Log likelihood contributions and their sum after the burn-in.

``kalman_filter``
-----------------

.. code-block:: fortran

   subroutine kalman_filter(rep, res, info)
     type(ssm_rep_t), intent(in) :: rep
     type(filter_result_t), intent(out) :: res
     integer, intent(out) :: info

Run the filter and store what the smoother needs. ``info`` is
``SS_ERR_NOT_PD`` if some :math:`F_t` is not positive definite, or the
status of ``validate`` or ``initial_state``.

``loglike``
-----------

.. code-block:: fortran

   function loglike(rep, info) result(llf)
     type(ssm_rep_t), intent(in) :: rep
     integer, intent(out) :: info
     real(dp) :: llf

The log likelihood (DK eq. 7.2), or the diffuse log likelihood with diffuse
states (DK §7.2.2), without storing the filter history. Estimation calls
this one.

``loglike_concentrated``
------------------------

.. code-block:: fortran

   function loglike_concentrated(rep, scale, info) result(llf)
     type(ssm_rep_t), intent(in) :: rep
     real(dp), intent(out) :: scale
     integer, intent(out) :: info
     real(dp) :: llf

The log likelihood with a scale :math:`\sigma^2` concentrated out (DK
§2.10.2): H, Q and :math:`P_*` are taken relative to :math:`\sigma^2`. With
S the sum of the :math:`v_t' F_t^{-1} v_t` terms at :math:`\sigma^2 = 1`
over the :math:`n_*` elements outside the diffuse updates,

.. math::

   \hat\sigma^2 = S / n_*, \qquad
   \log L_c = \log L_1 + S/2 - n_* (\log \hat\sigma^2 + 1) / 2.

This equals the ordinary log likelihood of the model scaled by
:math:`\hat\sigma^2` (see ``scale_by``). It is computed from the part of
the log likelihood without the quadratic terms, which avoids cancellation.

``marginal_correction``
-----------------------

.. code-block:: fortran

   real(dp) function marginal_correction(rep, info) result(corr)

The marginal minus the diffuse log likelihood (Francke, Koopman and de Vos
2010, eqs. 16, 21; DK §7.2.6):

.. math::

   \log L_M = \log L_d + \tfrac{k}{2} \log 2\pi + \tfrac12 \log|S_*|, \qquad
   S_* = \sum_t V_t^{*\prime} V_t^*, \quad V_t^* = Z_{t,o} A_t^*, \quad
   A_{t+1}^* = T_t A_t^*,

with :math:`A_1^* = A` the diffuse directions of :math:`\alpha_1`. It
depends on Z, T and A only, so it changes estimates only when parameters
enter Z or T. Zero without diffuse states. Used when
``rep%marginal_likelihood`` is set.

``steady_state``
----------------

.. code-block:: fortran

   subroutine steady_state(rep, P, F, info, tol, maxiter, niter)
     type(ssm_rep_t), intent(in) :: rep
     real(dp), intent(out) :: P(:, :), F(:, :)
     integer, intent(out) :: info
     real(dp), intent(in), optional :: tol        ! default 1e-12
     integer, intent(in), optional :: maxiter     ! default 100000
     integer, intent(out), optional :: niter

The limit :math:`\bar P` of the recursion for :math:`P_t` (DK §2.11,
§4.3.4) started from :math:`P_1 = 0`, and :math:`\bar F = Z \bar P Z' + H`.
With H positive definite it uses the structure-preserving doubling
algorithm (Chu, Fan, Lin and Wang 2004): step k equals step
:math:`2^k` of the recursion, so the slow convergence of fixed regression
coefficients or nearly fixed seasonals takes some 40 steps. With H singular
it iterates the recursion itself. ``info`` is ``SS_ERR_UNSUPPORTED`` for
time-varying matrices and ``SS_ERR_NOT_CONVERGED`` if the change in P does
not fall below ``tol`` relative to P.
