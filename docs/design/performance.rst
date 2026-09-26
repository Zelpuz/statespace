Performance
===========

The filter in statsmodels is compiled Cython calling BLAS. The time around
it is spent in Python: building the matrices at each evaluation, finite
differences for the gradient, storing the full filter history when only
the likelihood is needed, and fitting series one at a time. The library
moves all of this into Fortran and adds four changes to the core.

Small matrices
--------------

**Decision.** ``gemm``, ``gemv`` and ``chol_inv`` use inline loops when the
operands are small (below about 10 × 10 × 10 for products and order 16 for
the inverse), and BLAS or LAPACK otherwise.

**Why.** State space models are small; m and p are often below 10. At that
size the cost of a call into OpenBLAS, which checks its threading and
dispatches on the operation, exceeds the arithmetic. The inline loops
halved the time of the local level and ARMA filters.

A copy on every call
--------------------

**Decision.** The array arguments of the numerical routines are declared
``contiguous``.

**Why.** gfortran 15 copies an assumed-shape argument when it passes it on
to a ``contiguous`` dummy, even when the array is contiguous at run time; a
debug build reported millions of such copies during the tests. Declaring
the dummies ``contiguous`` all the way down removed them and halved the
time per filter step again. The interfaces users implement (``update``,
``fill`` and the parameter transforms) are left unchanged, so that existing
models still compile.

The steady state and the likelihood
-----------------------------------

The log likelihood is computed without storing the filter history, and in
time-invariant models the filter switches to its steady state
(:doc:`steady_state`). For estimation, which evaluates only the likelihood
and the score, these two changes matter most.

Parallel fits
-------------

**Decision.** ``fit_many`` fits independent models on OpenMP threads, with
the GIL released in Python.

**Why.** Fitting many short series is a common use and parallelizes without
coordination. L-BFGS-B keeps its state in the caller's arrays and the
library has no global state, so the routines are reentrant; each model is
touched by one thread only. BLAS should then run single-threaded
(``OPENBLAS_NUM_THREADS=1``).

Numbers
-------

From ``example/bench.f90`` and ``bench/bench_statsmodels.py``, time for the
log likelihood at n = 100,000 in milliseconds, before and after these
changes, against statsmodels' Cython filter:

.. list-table::
   :header-rows: 1

   * - model
     - before
     - after
     - statsmodels
   * - local level (n = 1,000,000)
     - 728
     - 49
     - 415
   * - level and trigonometric seasonal
     - 141
     - 9.4
     - 251
   * - 8 series, 3 factors
     - 180
     - 8.5
     - 118
   * - ARMA(2, 1)
     - 79
     - 5.2
     - 44

Problems found
--------------

* An address sanitizer run found a heap overflow in
  ``collapse_observations``: gfortran at ``-O2`` sized an inlined
  ``matmul`` wrongly in an assignment whose shape changed between periods.
  The array is now allocated explicitly. The tests run under debug, release,
  OpenMP and address-sanitizer builds.
* gfortran 15 at ``-O3`` fails with an internal compiler error on
  ``names = mb%model%param_names()`` for a polymorphic ``model``; a helper
  routine avoids it (``model_names`` in ``statespace_capi_models``).
