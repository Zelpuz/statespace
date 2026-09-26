``statespace_mle``
==================

Maximum likelihood estimation (DK §7.3) with L-BFGS-B (Zhu, Byrd, Lu and
Nocedal 1997; Morales and Nocedal 2011).

The optimizer minimizes -loglike / n over the unconstrained parameters, as
statsmodels does, with statsmodels' default tolerances. The gradient is the
analytic score (DK §7.3.3; :doc:`statespace_score`) for the parameters it
covers and central differences for the others (parameters in Z or T, as
DK recommend). Standard errors come from a numerical Hessian in the
unconstrained parameters, mapped to the constrained ones by the delta
method (DK §7.3.6). Points where the likelihood cannot be evaluated get a
large objective value, so the line search backs off.

See :doc:`../../design/estimation` for the choices behind the gradient and
the standard errors.

Example
-------

.. code-block:: fortran

   type(fit_result_t) :: res
   type(fit_options_t) :: opts

   opts%factr = 10.0_dp        ! tight convergence
   opts%pgtol = 1.0e-9_dp
   call fit(model, res, options=opts, info=info)
   print *, res%params, res%bse, res%llf

``fit_options_t``
-----------------

.. code-block:: fortran

   type :: fit_options_t
     integer :: maxiter = 500
     integer :: m = 10                  ! L-BFGS corrections
     real(dp) :: factr = 1.0e7_dp       ! relative reduction tolerance / machine epsilon
     real(dp) :: pgtol = 1.0e-5_dp      ! projected gradient tolerance
     logical :: compute_cov = .true.
     integer :: iprint = -1             ! L-BFGS-B output; negative is silent
     integer :: gradient = GRADIENT_AUTO
   end type

``gradient`` is ``GRADIENT_AUTO`` (0; analytic where it applies),
``GRADIENT_NUMERICAL`` (1) or ``GRADIENT_ANALYTIC`` (2; fail with
``SS_ERR_UNSUPPORTED`` if no parameter is covered).

``fit_result_t``
----------------

.. code-block:: fortran

   type :: fit_result_t
     real(dp), allocatable :: params(:)          ! constrained estimates
     real(dp), allocatable :: cov_params(:, :)
     real(dp), allocatable :: bse(:)
     real(dp) :: llf, scale, aic, bic
     integer :: niter, nfev
     logical :: converged, analytic_gradient
     character(len=60) :: message
   end type

``aic`` and ``bic`` count the diffuse states and a concentrated scale as
parameters and use n minus the burn-in, as statsmodels does; DK §7.4 also
divide by n. ``analytic_gradient`` is true if the score covered at least one
parameter.

``fit``
-------

.. code-block:: fortran

   subroutine fit(model, res, start_params, options, info)
     class(ssm_model_t), intent(inout) :: model
     type(fit_result_t), intent(out) :: res
     real(dp), intent(in), optional :: start_params(:)
     type(fit_options_t), intent(in), optional :: options
     integer, intent(out) :: info

Estimate the parameters; the model is left at the estimates. If the
estimates exist but the Hessian is not negative definite, ``info`` is
``SS_ERR_NOT_PD`` and ``cov_params`` and ``bse`` are not allocated.

``fit_many``
------------

.. code-block:: fortran

   subroutine fit_many(models, res, info, options)
     class(ssm_model_t), intent(inout) :: models(:)
     type(fit_result_t), intent(out) :: res(:)
     integer, intent(out) :: info(:)

Fit independent models in parallel with OpenMP when the library is built
with ``-fopenmp``, in turn otherwise. Each model is touched by one thread
only; L-BFGS-B keeps its state in the caller's arrays, so it is reentrant.
For large models set ``OPENBLAS_NUM_THREADS=1`` to avoid nested threading.

``numerical_hessian``
---------------------

.. code-block:: fortran

   function numerical_hessian(model, params, info, transformed) result(hess)

Central-difference Hessian of the log likelihood at ``params``, constrained
unless ``transformed = .false.``.

``estimation_bias``
-------------------

.. code-block:: fortran

   subroutine estimation_bias(model, fres, ndraw, bias_alpha, info, bias_V, antithetic, failed)
     type(fit_result_t), intent(in) :: fres
     integer, intent(in) :: ndraw
     real(dp), intent(out) :: bias_alpha(:, :)                ! (m, n)
     real(dp), intent(out), optional :: bias_V(:, :, :)       ! (m, m, n)
     logical, intent(in), optional :: antithetic              ! default .true.
     integer, intent(out), optional :: failed

The bias of the smoothed state from estimating the parameters (DK §7.3.7,
eq. 7.20). Treating :math:`\hat\psi` as the true value, draw
:math:`\psi^{(i)} \sim N(\hat\psi, \Omega)` on the unconstrained scale, with
:math:`\Omega` from ``fres%cov_params`` by the delta method, and estimate

.. math::

   B = \frac1N \sum_i \hat\alpha(\psi^{(i)}) - \hat\alpha(\hat\psi).

With ``antithetic``, each draw is paired with its reflection
:math:`2\hat\psi - \psi^{(i)}`, and ``ndraw`` must be even. Draws where the
model cannot be evaluated are skipped with their pair and counted in
``failed``. Uses ``random_number``; seed it for reproducible results. The
model is left at :math:`\hat\psi`.
