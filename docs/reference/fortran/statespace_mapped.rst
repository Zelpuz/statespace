``statespace_mapped``
=====================

A model declared as a map from parameters to entries of the system
matrices, for models written without Fortran code. The Python package's
:class:`~ssfortran.MappedModel` builds on it.

The representation holds every fixed entry and the initialization. The map
says which entries the parameters set:

* **entries**: parameter k sets :math:`X(i, j, t) = c\,\psi_k`, for X one
  of Z, H, T, R, Q, c, d (codes 2-8, as ``SS_ARR_*`` in the C interface)
  and t = 0 for every time slice;
* **covariance blocks**: parameters :math:`k, \dots, k + q - 1`,
  :math:`q = d(d+1)/2`, are the lower triangle, by columns, of a d × d
  block of H or Q starting at row and column ``offset``.

**Transform groups** map consecutive parameters to the optimizer's scale:

=====================  =====  ==================================================
kind                   value  transform
=====================  =====  ==================================================
``MAP_NONE``           0      none
``MAP_POSITIVE``       1      :math:`\exp(2x)` (DK §7.3.2)
``MAP_INTERVAL``       2      :math:`m + h x / \sqrt{1 + x^2}` onto (lo, hi)
``MAP_STATIONARY``     3      stationary AR coefficients (Monahan 1984)
``MAP_INVERTIBLE``     4      invertible MA coefficients
``MAP_COV``            5      Cholesky factor with log diagonal (covariance blocks)
=====================  =====  ==================================================

Parameters in no group are unconstrained. Since ``update`` only overwrites
the mapped entries, the fixed entries are set once.

``mapped_model_t``
------------------

.. code-block:: fortran

   type, extends(ssm_model_t) :: mapped_model_t
     type(map_entry_t), allocatable :: entries(:)
     type(map_block_t), allocatable :: blocks(:)
     type(map_group_t), allocatable :: groups(:)
     real(dp), allocatable :: start(:)
     character(len=32), allocatable :: names(:)
   contains
     procedure :: add_entry, add_block, add_group
   end type

.. code-block:: fortran

   subroutine add_entry(self, param, array, i, j, t, coef)
   subroutine add_block(self, array, offset, dim, first)     ! also adds a MAP_COV group
   subroutine add_group(self, kind, first, last, lo, hi)

``map_entry_t``, ``map_block_t`` and ``map_group_t`` hold these
declarations.

``mapped_model``
----------------

.. code-block:: fortran

   function mapped_model(rep, k) result(model)
     type(ssm_rep_t), intent(in) :: rep
     integer, intent(in) :: k
     type(mapped_model_t) :: model

A mapped model with k parameters over a copy of ``rep``; start values 0.1
and names ``param1, ...`` until set.

Example
-------

The local level model:

.. code-block:: fortran

   type(mapped_model_t) :: model
   rep = ssm_rep(y, m=1, r=1)
   rep%Z = 1.0_dp; rep%T = 1.0_dp
   call rep%initialize_diffuse()
   model = mapped_model(rep, 2)
   call model%add_entry(1, 3, 1, 1, 0, 1.0_dp)          ! H(1, 1) = params(1)
   call model%add_entry(2, 6, 1, 1, 0, 1.0_dp)          ! Q(1, 1) = params(2)
   call model%add_group(MAP_POSITIVE, 1, 2, 0.0_dp, 0.0_dp)
