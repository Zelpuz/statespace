``statespace_model``
====================

The abstract model with parameters, the Fortran counterpart of statsmodels'
``MLEModel``, and the standard parameter transforms.

A model extends ``ssm_model_t``. Its constructor sets up ``rep`` (data,
fixed matrices, initialization) and ``k_params``; it implements

* ``update(params)``: write the constrained parameters into the system
  matrices;
* ``start_params()``: constrained starting values;

and may override ``transform_params`` (unconstrained to constrained),
``untransform_params`` and ``param_names``. The optimizer works on
unconstrained values; ``update`` always receives constrained ones.

Example
-------

The local level model (``example/nile_mle.f90``):

.. literalinclude:: ../../../example/nile_mle.f90
   :language: fortran
   :lines: 6-78

``ssm_model_t``
---------------

.. code-block:: fortran

   type, abstract :: ssm_model_t
     type(ssm_rep_t) :: rep
     integer :: k_params = 0
     logical :: concentrate_scale = .false.
     real(dp) :: scale = 1.0_dp
   contains
     procedure(update_iface), deferred :: update
     procedure(start_params_iface), deferred :: start_params
     procedure :: transform_params, untransform_params, param_names
     procedure :: loglike, rep_at, filter, smooth
   end type

``concentrate_scale``
   H, Q and :math:`P_*` set by ``update`` are relative to a scale
   :math:`\sigma^2`, which is concentrated out of the likelihood
   (DK §2.10.2) and is not a parameter. ``scale`` holds its estimate at the
   last evaluation.

Deferred procedures:

.. code-block:: fortran

   subroutine update(self, params)
     class(ssm_model_t), intent(inout) :: self
     real(dp), intent(in) :: params(:)

   function start_params(self) result(params)
     class(ssm_model_t), intent(in) :: self
     real(dp), allocatable :: params(:)

Overridable procedures, with defaults:

.. code-block:: fortran

   function transform_params(self, unconstrained) result(constrained)     ! identity
   function untransform_params(self, constrained) result(unconstrained)   ! identity
   function param_names(self) result(names)            ! "param1", ...; character(32)

Evaluation at given parameters (constrained unless
``transformed = .false.``):

.. code-block:: fortran

   function loglike(self, params, info, transformed) result(llf)
   subroutine rep_at(self, params, rep, info, transformed)
   subroutine filter(self, params, fres, info, transformed)
   subroutine smooth(self, params, fres, sres, info, transformed)

``loglike`` calls ``update`` and returns the log likelihood, concentrated
when ``concentrate_scale`` is set. ``rep_at`` returns the representation at
``params`` on the data's scale (H, Q and :math:`P_*` multiplied by the
estimated scale); ``filter`` and ``smooth`` run on it.

Transforms
----------

.. code-block:: fortran

   elemental real(dp) function constrain_positive(x)      ! exp(2 x)
   elemental real(dp) function unconstrain_positive(x)    ! log(x) / 2
   pure function constrain_stationary(unconstrained) result(phi)
   pure function unconstrain_stationary(phi) result(unconstrained)

``constrain_positive``
   :math:`\sigma^2 = \exp(2\psi)` for variances (DK §7.3.2).
``constrain_stationary``
   Coefficients :math:`\phi` of a stationary AR polynomial
   :math:`1 - \phi_1 L - \dots - \phi_p L^p` (Monahan 1984): the partial
   autocorrelations are :math:`x / \sqrt{1 + x^2}`, DK's bounded transform
   with a = 1. The negated result gives invertible MA coefficients.
   Matches statsmodels' ``constrain_stationary_univariate`` to 1e-14.
