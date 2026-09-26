``statespace_callback``
=======================

A model whose ``update``, and optionally its parameter transforms, are C
function pointers, for models defined in another language. The Python
package's :class:`~ssfortran.MLEModel` builds on it.

``update_cb(k, params, rep)`` fills the model's representation, passed as a
handle usable with the ``ss_rep_*`` routines of the C interface, from the k
constrained parameters, and returns 0 on success.
``transform_cb(k, x, out)`` and ``untransform_cb(k, x, out)`` map between
unconstrained and constrained parameters; without them the parameters are
unconstrained.

When ``update_cb`` fails, H is saved and set to NaN, so the likelihood
cannot be evaluated and the optimizer backs off; H is restored at the next
call. The caller reports the failure (the Python package re-raises the
exception). Callback models must not be fitted on several threads.

``callback_model_t``
--------------------

.. code-block:: fortran

   type, extends(ssm_model_t) :: callback_model_t
     type(c_funptr) :: update_cb, transform_cb, untransform_cb
     type(c_ptr) :: rep_handle          ! handle to this model's rep, passed to update_cb
     real(dp), allocatable :: start(:)
     integer :: failures = 0            ! failed callbacks so far
   end type

C prototypes of the callbacks:

.. code-block:: c

   int update_cb(int k, const double *params, void *rep);
   int transform_cb(int k, const double *x, double *out);
