Installation
============

Requirements
------------

* A Fortran 2018 compiler; gfortran 13 or later is tested.
* LAPACK and BLAS, for example OpenBLAS.
* For the Python package: Python 3.10 or later and numpy. The build uses
  CMake and scikit-build-core, which ``pip`` installs.

Python
------

From the repository root::

   pip install .

This compiles ``libstatespace`` with CMake, with OpenMP when the compiler
supports it, and puts it inside the ``ssfortran`` package. The library is
loaded with ctypes, so the wheel does not depend on the Python version.

To develop without installing, build the library and point Python at the
source tree::

   cmake -S . -B build/cmake -G Ninja
   cmake --build build/cmake
   PYTHONPATH=python pytest python/tests

``SSFORTRAN_LIB`` overrides the location of the library.

Fortran
-------

The Fortran library builds with `fpm <https://fpm.fortran-lang.org>`_::

   fpm build --profile release
   fpm test --profile release
   fpm run --profile release --example nile_mle

``fit_many`` runs in parallel when OpenMP is enabled::

   fpm build --profile release --flag -fopenmp --link-flag -fopenmp

To use the library from another fpm project, add it as a dependency; the
modules are then available through ``use statespace``.

Data for the examples
---------------------

The Nile data are in ``data/nile.csv``. The chapter 8 examples need data
that are downloaded rather than stored in the repository::

   python data/fetch_dk_data.py

See ``data/README.md`` for their sources and terms.

Documentation
-------------

With Sphinx, numpydoc, pydata-sphinx-theme and myst-parser installed::

   make -C docs html        # into build/docs/html
   make -C docs doctest     # run the examples in the documentation

The examples print estimates rounded to the digits that optimizers and BLAS
libraries agree on; on another platform a last digit may still differ.
