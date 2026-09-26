Fortran API
===========

The library is a set of Fortran 2018 modules. ``use statespace`` imports
the whole public interface; the modules can also be used one at a time.

Conventions
-----------

* Reals are ``real(dp)``, 64-bit (:doc:`statespace_kinds`).
* Arrays are column-major with time as the last dimension. A system matrix
  has a time dimension of 1 (time-invariant) or n (time-varying);
  ``tidx(nt, t)`` gives the slice for period t.
* Periods and indices are 1-based.
* Routines that can fail have an ``info`` argument set to a status code
  (``SS_OK`` or ``SS_ERR_*``, :doc:`statespace_kinds`). They do not stop the
  program.
* Missing observations are NaN in ``y``.
* References such as "DK §4.3" are to Durbin and Koopman (2012).

Modules
-------

Model and algorithms, in dependency order:

.. toctree::
   :maxdepth: 1

   statespace_kinds
   statespace_rep
   statespace_filter
   statespace_smoother
   statespace_forecast
   statespace_simsmooth
   statespace_smoothing
   statespace_augmented
   statespace_sqrt
   statespace_collapse
   statespace_restrict
   statespace_dense

Models, estimation and diagnostics:

.. toctree::
   :maxdepth: 1

   statespace_model
   statespace_mle
   statespace_score
   statespace_em
   statespace_components
   statespace_structural
   statespace_arima
   statespace_mapped
   statespace_callback
   statespace_diagnostics

Numerical support:

.. toctree::
   :maxdepth: 1

   statespace_linalg
   statespace_special

The C interface (``statespace_capi*``) is described in :doc:`../c_api`.
