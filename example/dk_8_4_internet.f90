!> DK 8.4: ARMA models for the first differences of the number of users
!> logged on to an Internet server each minute (99 observations), fitted in
!> state space form (DK 3.4), with the AIC of DK 7.4,
!>
!>     AIC = n^-1 [-2 log L + 2 w],   w = p + q + 1,
!>
!> for p, q = 0..5 (Table 8.1); then again with 14 observations treated as
!> missing (Table 8.2); then forecasts from the ARMA(1, 1) model with 50%
!> intervals (Fig. 8.9).
!>
!> The log likelihoods agree with statsmodels' SARIMAX to 4 decimals, with
!> and without the missing observations (e.g. ARMA(1, 1): -254.1497 and
!> -225.7704). Higher-order models can have several local optima.
!> DK's tables use another scale: they are close to n^-1 [-log L + 2 (p + q)]
!> for the complete series, but no simple formula fits the table with missing
!> observations, and DK report optimizer failures for several models. With
!> the formula above, AIC chooses among the same few models: ARMA(3, 0),
!> ARMA(1, 1) and some higher-order models are within 0.03 of each other.
!>
!> Run from the repo root, after `.venv/bin/python data/fetch_dk_data.py`:
!>     fpm run --example dk_8_4_internet
program dk_8_4_internet
  use statespace
  use csv_io, only: read_csv, column
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none

  integer, parameter :: missing(14) = [6, 16, 26, 36, 46, 56, 66, 72, 73, 74, 75, 76, &
                                       86, 96]
  character(len=32), allocatable :: names(:)
  real(dp), allocatable :: data(:, :), users(:), dy(:, :)

  call read_csv("data/wwwusage.csv", names, data)
  users = column(names, data, "value")
  allocate (dy(1, size(users) - 1))
  dy(1, :) = users(2:) - users(:size(users) - 1)

  print '(a)', "Table 8.1: AIC for ARMA(p, q), differenced series"
  call aic_table(dy)
  dy(1, missing) = ieee_value(1.0_dp, ieee_quiet_nan)
  print '(/, a)', "Table 8.2: the same with 14 observations missing"
  call aic_table(dy)
  print '(/, a)', &
      "ARMA(1, 1) forecasts with 50% intervals, series with missing observations"
  call forecasts(dy, 10)

contains

  function arma(y, p, q) result(mod)
    real(dp), intent(in) :: y(:, :)
    integer, intent(in) :: p, q
    type(structural_model_t) :: mod
    type(component_holder_t) :: comps(1)
    integer :: info

    comps(1)%c = arima_t(ar=p, ma=q)
    mod = structural_model(y, comps, info)
    if (info /= SS_OK) error stop "model"
  end function arma

  subroutine aic_table(y)
    real(dp), intent(in) :: y(:, :)
    type(structural_model_t) :: mod
    type(fit_result_t) :: res
    type(fit_options_t) :: opts
    character(len=10) :: cell(0:5)
    integer :: p, q, info

    opts%compute_cov = .false.
    opts%factr = 10.0_dp
    opts%pgtol = 1.0e-8_dp
    print '(4x, 6(i10))', [(q, q=0, 5)]
    do p = 0, 5
      cell = ""
      do q = 0, 5
        if (p == 0 .and. q == 0) cycle
        mod = arma(y, p, q)
        call fit(mod, res, options=opts, info=info)
        if (info /= SS_OK) then
          cell(q) = "failed"
        else
          write (cell(q), '(f10.3)') (-2 * res%llf + 2 * (p + q + 1)) / size(y, 2)
        end if
      end do
      print '(i4, 6a10)', p, cell
    end do
  end subroutine aic_table

  subroutine forecasts(y, h)
    real(dp), intent(in) :: y(:, :)
    integer, intent(in) :: h
    type(structural_model_t) :: mod
    type(fit_result_t) :: res
    type(filter_result_t) :: fres
    type(forecast_result_t) :: fc
    real(dp), parameter :: z50 = 0.6744897501960817_dp    ! 75% normal quantile
    real(dp) :: sd
    integer :: j, info

    mod = arma(y, 1, 1)
    call fit(mod, res, info=info)
    if (info /= SS_OK) error stop "fit"
    call mod%filter(res%params, fres, info)
    call forecast(mod%rep, fres, h, fc, info)
    if (info /= SS_OK) error stop "forecast"
    print '(4x, a, 3a10)', "t", "forecast", "lower", "upper"
    do j = 1, h
      sd = sqrt(fc%cov(1, 1, j))
      print '(i5, 3f10.3)', size(y, 2) + j, fc%mean(1, j), fc%mean(1, j) - z50 * sd, &
          fc%mean(1, j) + z50 * sd
    end do
  end subroutine forecasts
end program dk_8_4_internet
