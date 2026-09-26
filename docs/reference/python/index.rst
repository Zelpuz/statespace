Python API
==========

.. currentmodule:: ssfortran

The package ``ssfortran`` calls the Fortran library through its C interface
(:doc:`../c_api`). Arrays are numpy arrays with time as the last axis, as in
statsmodels; periods passed as arguments are 0-based.

Representation
--------------

.. autosummary::
   :toctree: generated

   Representation
   FilterResults
   SmootherResults
   AugmentedResults

Models and estimation
---------------------

.. autosummary::
   :toctree: generated

   StructuralModel
   MappedModel
   MLEModel
   Model
   FitResults
   ForecastResults
   fit_many

Components
----------

.. autosummary::
   :toctree: generated

   Irregular
   Level
   Trend
   Seasonal
   Cycle
   Regression
   ARIMA
   ContinuousLevel
   ContinuousTrend

Diagnostics
-----------

.. autosummary::
   :toctree: generated

   diagnostics.ljung_box
   diagnostics.jarque_bera
   diagnostics.breakvar

Errors and constants
--------------------

.. autosummary::
   :toctree: generated

   StateSpaceError
   version

The initialization kinds :data:`INIT_KNOWN`, :data:`INIT_APPROX_DIFFUSE`,
:data:`INIT_STATIONARY`, :data:`INIT_DIFFUSE` and :data:`INIT_GENERAL` are
passed to :meth:`Representation.initialize_block`;
:data:`FILTER_CONVENTIONAL`, :data:`FILTER_UNIVARIATE`,
:data:`DIFFUSE_UNIVARIATE` and :data:`DIFFUSE_MULTIVARIATE` are values of
the ``filter_method`` and ``diffuse_method`` options.

.. data:: INIT_KNOWN
.. data:: INIT_APPROX_DIFFUSE
.. data:: INIT_STATIONARY
.. data:: INIT_DIFFUSE
.. data:: INIT_GENERAL
.. data:: FILTER_CONVENTIONAL
.. data:: FILTER_UNIVARIATE
.. data:: DIFFUSE_UNIVARIATE
.. data:: DIFFUSE_MULTIVARIATE
