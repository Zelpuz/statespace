``statespace_simsmooth``
========================

Simulation from the model and the simulation smoothers (DK §4.9, §5.5).

The random input is standard normal variates, so draws are reproducible:
``u_init`` (m), ``u_eps`` (p, n) and ``u_eta`` (r, n). They are scaled by
lower triangular square roots of :math:`P_*`, :math:`H_t` and :math:`Q_t`
(the Cholesky factor when the matrix is positive definite), as in
statsmodels, so the same variates give the same draws.

``simsmooth_result_t``
----------------------

.. code-block:: fortran

   type :: simsmooth_result_t
     real(dp), allocatable :: state(:, :)  ! (m, n) draw of alpha_t | Y_n
     real(dp), allocatable :: eps(:, :)    ! (p, n) draw of eps_t | Y_n
     real(dp), allocatable :: eta(:, :)    ! (r, n) draw of eta_t | Y_n
   end type

``simulate``
------------

.. code-block:: fortran

   subroutine simulate(rep, u_init, u_eps, u_eta, y, alpha, eps, eta, info)
     real(dp), intent(in) :: u_init(:), u_eps(:, :), u_eta(:, :)
     real(dp), intent(out) :: y(:, :), alpha(:, :), eps(:, :), eta(:, :)

Simulate from the model:

.. math::

   \alpha_1 = a_1 + S(P_*) u_{init}, \quad
   y_t = d_t + Z_t \alpha_t + S(H_t) u_{\varepsilon,t}, \quad
   \alpha_{t+1} = c_t + T_t \alpha_t + R_t S(Q_t) u_{\eta,t},

with S(A) a lower triangular square root of A. Missing values in
``rep%y`` are ignored. Diffuse states start at their mean.

``simulation_smoother``
-----------------------

.. code-block:: fortran

   subroutine simulation_smoother(rep, u_init, u_eps, u_eta, sim, info)
     type(simsmooth_result_t), intent(out) :: sim

Draw :math:`(\alpha, \varepsilon, \eta)` given the data by mean corrections
(Durbin and Koopman 2002; DK §4.9.2):

1. simulate :math:`(y^+, \alpha^+, \varepsilon^+, \eta^+)` from the model;
2. smooth :math:`y - y^+` in the model with zero intercepts and
   :math:`a_1 = 0`, which gives :math:`\hat\alpha(y) - \hat\alpha(y^+)`
   because smoothing is affine;
3. add the result to the simulated values.

Diffuse directions of :math:`\alpha_1` get no draw: the exact diffuse
smoother reproduces any shift in them, so they cancel in step 3 (DK §5.5).
For a missing element the drawn :math:`\varepsilon` is conditional on the
observed elements (statsmodels draws it unconditionally).

``djs_measurement_disturbances``
--------------------------------

.. code-block:: fortran

   subroutine djs_measurement_disturbances(rep, fres, u_eps, eps, info)
     type(filter_result_t), intent(in) :: fres
     real(dp), intent(in) :: u_eps(:, :)
     real(dp), intent(out) :: eps(:, :)

The de Jong-Shephard simulation smoother for :math:`\varepsilon` (DK §4.9.3,
eqs. 4.83-4.88): going back from t = n, draw :math:`\varepsilon_t` from its
distribution given the data and the later draws. Not defined with diffuse
states.

``djs_state_disturbances``
--------------------------

.. code-block:: fortran

   subroutine djs_state_disturbances(rep, fres, u_eta, u_init, eta, alpha, info)
     real(dp), intent(in) :: u_eta(:, :), u_init(:)
     real(dp), intent(out) :: eta(:, :), alpha(:, :)

The same for :math:`\eta` (DK eqs. 4.89-4.91), then the states forward
(DK eq. 4.92). DK start the states from the mean
:math:`a_1 + P_1 \tilde r_0` of :math:`\alpha_1` given the data and the
:math:`\eta` draws; this routine draws :math:`\alpha_1` from that
distribution, with variance :math:`P_1 - P_1 \tilde N_0 P_1`, using
``u_init``. Without that draw the state draws are too concentrated. Not
defined with diffuse states.

``psd_sqrt``
------------

.. code-block:: fortran

   function psd_sqrt(A) result(S)

Lower triangular S with :math:`S S' = A` for symmetric positive
semi-definite A: :math:`S = L D^{1/2}` from :math:`A = L D L'`; the Cholesky
factor when A is positive definite.

``draw_standard_normal``
------------------------

.. code-block:: fortran

   subroutine draw_standard_normal(x)
     real(dp), intent(out) :: x(:)

Fill x with independent N(0, 1) draws (Box-Muller on ``random_number``).
