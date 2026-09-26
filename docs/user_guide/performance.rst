Performance
===========

Timings in milliseconds, minimum over repeats, on one core of a 12-core
x86-64 machine with OpenBLAS; from ``bench/bench_python.py``
(ssfortran and statsmodels both called from Python):

.. list-table::
   :header-rows: 1

   * - case
     - ssfortran
     - statsmodels
   * - log likelihood, local level, n = 100,000
     - 5.0
     - 41
   * - smoother, local level, n = 100,000
     - 79
     - 284
   * - fit, level and trigonometric seasonal, n = 1000
     - 114
     - 528
   * - 1000 local level fits, n = 100, ``fit_many``
     - 201
     - 7185
   * - 1000 local level fits, n = 100, one by one
     - 1150
     - 7185
   * - fit, local level as a Python ``MLEModel``
     - 11
     - 23

``example/bench.f90`` and ``bench/bench_statsmodels.py`` compare the
Fortran core with statsmodels' Cython filter directly.

Where the time goes
-------------------

* **The likelihood.** Estimation evaluates the log likelihood without
  storing the filter output. In time-invariant models the steady-state
  shortcut (DK §4.3.4) then reduces a step to a few matrix-vector products.
* **The gradient.** The analytic score costs one filter and smoother pass,
  against 2k likelihood evaluations for central differences.
* **Small matrices.** State space models are small: m and p are often
  below 10. At that size the per-call cost of BLAS dominates, so the
  library uses inline loops below a size threshold.
* **Many series.** :func:`~ssfortran.fit_many` fits independent models on
  OpenMP threads without the GIL. Set ``OPENBLAS_NUM_THREADS=1`` so BLAS
  does not start threads of its own.
* **Python models.** An :class:`~ssfortran.MLEModel` calls Python at every
  evaluation; a :class:`~ssfortran.MappedModel` does not.

See :doc:`../design/performance` for the details.
