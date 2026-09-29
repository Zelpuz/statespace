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

From PyPI, with pip or uv::

   pip install ssfortran
   uv add ssfortran             # in a uv project; or: uv pip install ssfortran

Wheels for Linux (x86_64 and aarch64, glibc 2.28 or later) include the
compiled library, gfortran's runtime and OpenBLAS, and need no compiler. On
other platforms pip and uv build from the source distribution, which needs
the compiler and libraries listed above.

From a clone of the repository::

   pip install .

This compiles ``libstatespace`` with CMake, with OpenMP when the compiler
supports it, and puts it inside the ``ssfortran`` package. The library is
loaded with ctypes, so the wheel does not depend on the Python version.

To develop without installing, build the library and point Python at the
source tree::

   cmake -S . -B build/cmake -G Ninja
   cmake --build build/cmake
   pytest                  # python/src and build/cmake, set in pyproject.toml

Code is formatted to 88 columns: ``ruff format`` and ``ruff check`` for the
Python (settings in ``pyproject.toml``), and `Fortitude
<https://fortitude.readthedocs.io>`_'s ``fortitude check`` for the Fortran
(settings in ``fpm.toml``), whose long lines are wrapped by hand.

``SSFORTRAN_LIB`` overrides the location of the library.

Linux wheels are built with `cibuildwheel <https://cibuildwheel.pypa.io>`_
(configuration in ``pyproject.toml``, workflow in
``.github/workflows/wheels.yml``). Each wheel carries its own gfortran runtime
and OpenBLAS, so it needs no compiler. To build one locally, with Docker::

   pipx run cibuildwheel --platform linux

The version is written once, in ``src/statespace_capi.f90``;
``python/tests/test_version.py`` checks the copies in ``fpm.toml``,
``CMakeLists.txt`` and ``docs/conf.py``.

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

The documentation is published at https://zelpuz.github.io/statespace/
with each release. To build it locally, with the pinned tools::

   pip install --group docs   # or: uv sync --group docs
   make -C docs html          # into build/docs/html
   make -C docs doctest       # run the examples in the documentation

The Makefile uses the tools in ``.venv``; ``VENV=<bin directory>`` points
it elsewhere.

The examples print estimates rounded to the digits that optimizers and BLAS
libraries agree on; on another platform a last digit may still differ.
