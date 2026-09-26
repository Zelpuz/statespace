Interfaces
==========

Representation and model
------------------------

**Decision.** The representation (``ssm_rep_t``) is plain data with public
components, separate from the model with parameters (``ssm_model_t``), as
statsmodels separates ``Representation`` and ``KalmanFilter`` from
``MLEModel``.

**Why.** The algorithms need only the matrices; the separation lets other
languages fill the matrices directly, and lets one model produce many
representations (at different parameters) without copying its logic.

The C interface
---------------

**Decision.** The library exports ``bind(C)`` routines with opaque handles,
copies inputs in, copies outputs into caller-allocated buffers, returns a
status code from every routine, and uses 1-based indices.

**Why.** Opaque handles keep Fortran's derived types out of the interface.
Copying in and out makes ownership plain: the caller owns its buffers, the
library owns its objects, and nothing is shared across the boundary. Status
codes, not ``stop``, leave error handling to the caller.

**Alternatives.** Exposing pointers into the library's arrays would save
copies, but ties the caller to the library's memory layout and lifetime.
The copies cost little next to a filter run.

**Where.** ``statespace_capi``, ``statespace_capi_models``,
``statespace_capi_extras``; :doc:`../reference/c_api`.

Python binding and build
------------------------

**Decision.** The Python package calls the C interface through ctypes. CMake
builds the shared library and scikit-build-core packages it into a wheel;
fpm remains the build for Fortran users and for the tests.

**Why.** ctypes is part of the standard library, needs no compiled
extension, and releases the GIL during calls, which ``fit_many`` relies on.
The wheel is then independent of the Python version. fpm does not build
shared libraries; CMake does, and fetches L-BFGS-B at the same pinned
commit as fpm.

**Alternatives.** cffi adds a dependency whose wheels can lag new Python
versions; Cython or f2py add a compiled extension.

**Scope.** The package returns numpy arrays. pandas indexes and plotting
were added and then removed as outside the library's scope.

Three kinds of user model
-------------------------

**Decision.** Python users define models in three ways: built-in components
(``StructuralModel``), a declared map from parameters to matrix entries
(``MappedModel``), or Python code that sets the matrices
(``MLEModel``).

**Why.** Speed depends on whether Python runs inside the optimization. The
first two run entirely in Fortran and can be fitted in parallel. The third
is the most general and familiar to statsmodels users, but calls Python at
every likelihood evaluation, and 2k more times per gradient because the
score differentiates ``update`` numerically.

**Callback failures.** An exception in ``update`` must not leave the model
broken. The library then sets H to NaN, so the evaluation fails and the
optimizer backs off; it saves H first and restores it at the next call.
Python stores the exception and raises it when the library returns.
Callback models cannot run on ``fit_many``'s threads; both the Python layer
and the library reject them.

**Where.** ``statespace_mapped``, ``statespace_callback``;
``python/ssfortran/models.py``.
