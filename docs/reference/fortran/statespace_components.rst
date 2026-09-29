``statespace_components``
=========================

Model components and their assembly into a state space model (DK ch. 3).

A component owns a block of the state vector: its columns of Z, its
diagonal blocks of T, R and Q, and optionally a contribution to H (the
irregular). ``structural_model`` stacks components as DK do (DK eq. 3.10):

.. math::

   Z_t = (Z^{[1]}_t, Z^{[2]}_t, \dots), \quad
   T_t = \mathrm{diag}(T^{[1]}_t, T^{[2]}_t, \dots),

and likewise R and Q, and initializes each block as the component asks:
diffuse for nonstationary components, stationary for ARMA blocks and damped
cycles (DK §5.6). The result is an ordinary ``ssm_model_t``.

**Signal loadings** (DK §3.3.2, §3.3.3, §3.7): with ``loading``
(p × p_sig), the components describe p_sig signals and the observations
load on them, :math:`Z_t = \Lambda Z^{sig}_t`; ``loading_free`` marks
entries of :math:`\Lambda` that are parameters, appended after the
components' parameters. Components whose ``observation_level()`` is true
(the irregular, intercepts) act on the observations directly.

The built-in components are in :doc:`statespace_structural` and
:doc:`statespace_arima`.

Example
-------

From ``example/dk_8_2_seatbelt.f90``: a basic structural model, then the
same with regression effects.

.. code-block:: fortran

   type(component_holder_t) :: comps(4)
   type(structural_model_t) :: model

   allocate (irregular_t :: comps(1)%c)
   allocate (level_t :: comps(2)%c)
   comps(3)%c = seasonal_t(period=12, form=SEASONAL_TRIG)
   model = structural_model(y, comps(1:3), info)
   call fit(model, res, info=info)

   comps(4)%c = regression_t(x=x)             ! petrol price and the seat belt law
   model = structural_model(y, comps, info)

``component_t``
---------------

.. code-block:: fortran

   type, abstract :: component_t
     integer :: p, n        ! series and observations (set by setup)
     integer :: m, r, k     ! states, disturbances, parameters
     logical :: tv_Z, tv_H, tv_T, tv_R, tv_Q    ! time-varying matrices
     logical :: at_observations = .false.
   contains
     procedure(setup_iface), deferred :: setup
     procedure(fill_iface), deferred :: fill
     procedure :: transform_params, untransform_params, start_params
     procedure :: param_names, init_blocks, observation_level
   end type

To write a component, implement

.. code-block:: fortran

   subroutine setup(self, y)                      ! set m, r, k and the tv_ flags
     class(component_t), intent(inout) :: self
     real(dp), intent(in) :: y(:, :)              ! the data (or the signals' stand-ins)

   subroutine fill(self, params, Z, H, T, R, Q)   ! write the blocks from constrained params
     class(component_t), intent(in) :: self
     real(dp), intent(in) :: params(:)
     real(dp), intent(inout) :: Z(:, :, :), H(:, :, :), T(:, :, :), R(:, :, :), Q(:, :, :)

``fill`` receives the component's blocks with time dimension 1 or n as its
flags say; ``H`` is the full observation variance, to which only
observation-level components add. The defaults of the other procedures:
no transform, start values 0.1, names ``param1, ...``, the whole block
diffuse (``init_blocks`` returns rows ``[first, last, kind]`` relative to
the block), acting at the signals.

``component_holder_t``
----------------------

.. code-block:: fortran

   type :: component_holder_t
     class(component_t), allocatable :: c
   end type

An array element holding any component.

``structural_model_t``
----------------------

.. code-block:: fortran

   type, extends(ssm_model_t) :: structural_model_t
     type(component_holder_t), allocatable :: comps(:)
     integer, allocatable :: s0(:), e0(:), k0(:)    ! offsets of states, disturbances, parameters
     real(dp), allocatable :: loading(:, :)         ! (p, p_sig)
     logical, allocatable :: loading_free(:, :)
     integer :: k_comp                              ! parameters of the components
     real(dp), allocatable :: Zsig(:, :, :)         ! (p_sig, m, 1|n)
   end type

Component i owns states ``s0(i)+1 .. s0(i)+comps(i)%c%m``. Parameters are
ordered component by component, then the free loadings.

``structural_model``
--------------------

.. code-block:: fortran

   function structural_model(y, comps, info, loading, loading_free) result(model)
     real(dp), intent(in) :: y(:, :)                        ! (p, n)
     type(component_holder_t), intent(in) :: comps(:)
     integer, intent(out) :: info
     real(dp), intent(in), optional :: loading(:, :)        ! (p, p_sig)
     logical, intent(in), optional :: loading_free(:, :)
     type(structural_model_t) :: model

Assemble the components; the model starts at its start parameters.

Helpers
-------

.. code-block:: fortran

   elemental real(dp) function constrain_interval(psi, lo, hi) result(chi)
   elemental real(dp) function unconstrain_interval(chi, lo, hi) result(psi)
   real(dp) function diff_variance(y) result(v)

``constrain_interval``
   DK's bounded transform (§7.3.2) shifted to an interval:
   :math:`\chi = m + h \psi / \sqrt{1 + \psi^2}` maps the real line onto
   (lo, hi), with m and h the midpoint and half-width.
``diff_variance``
   Sample variance of the first differences of the first series, NaNs
   skipped; the scale of the default start values.
