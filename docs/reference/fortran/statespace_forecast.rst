``statespace_forecast``
=======================

Forecasts past the end of the sample (DK §4.11).

Forecasting is filtering with missing observations. ``forecast`` runs the
prediction recursion from the end of a filter run, which must be past the
diffuse periods, and needs time-invariant system matrices. For
time-varying models, append NaN observations to ``y`` with the matrices
for the forecast horizon, and read ``yhat`` and ``F`` from the filter.

``forecast_result_t``
---------------------

.. code-block:: fortran

   type :: forecast_result_t
     integer :: horizon                          ! h
     real(dp), allocatable :: mean(:, :)         ! (p, h)     E(y_n+j | Y_n)
     real(dp), allocatable :: cov(:, :, :)       ! (p, p, h)  Var(y_n+j | Y_n)
     real(dp), allocatable :: state(:, :)        ! (m, h)     E(alpha_n+j | Y_n)
     real(dp), allocatable :: state_cov(:, :, :) ! (m, m, h)
   end type

``forecast``
------------

.. code-block:: fortran

   subroutine forecast(rep, fres, h, fc, info)
     type(ssm_rep_t), intent(in) :: rep
     type(filter_result_t), intent(in) :: fres
     integer, intent(in) :: h
     type(forecast_result_t), intent(out) :: fc
     integer, intent(out) :: info

Forecast j = 1, ..., h steps ahead. ``info`` is ``SS_ERR_UNSUPPORTED`` if a
system matrix or intercept varies over time.
