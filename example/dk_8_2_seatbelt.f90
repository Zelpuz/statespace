!> DK 8.2: basic structural model for the log of monthly UK car drivers
!> killed or seriously injured, January 1969 - December 1984, then with the
!> seat belt law (a level shift from February 1983) and the log real petrol
!> price as regression effects.
!>
!> DK's published values: variances 0.00341598 (irregular), 0.000935852
!> (level), 5.01096e-7 (seasonal); petrol -0.29140 (rmse 0.09832), law
!> -0.23773 (0.04632). These are reproduced. DK also print a log likelihood
!> of 435.295 and a prediction error variance of 0.00586717, which use
!> other conventions: the steady-state F here is 0.0062583 (DK 2.11, 4.3.4),
!> and 0.00586717 is that times (n - d) / n = 180 / 192. The log likelihood
!> here is DK's diffuse log L_d (7.2.2), the same as statsmodels'.
!>
!> Run from the repo root, after `.venv/bin/python data/fetch_dk_data.py`:
!>     fpm run --example dk_8_2_seatbelt
program dk_8_2_seatbelt
  use statespace
  use csv_io, only: read_csv, column
  implicit none

  integer, parameter :: n = 192, t_law = 170      ! February 1983
  character(len=32), allocatable :: names(:)
  real(dp), allocatable :: data(:, :), y(:, :), x(:, :)
  type(component_holder_t) :: comps(4)
  type(structural_model_t) :: mod
  type(fit_result_t) :: res
  type(fit_options_t) :: opts
  integer :: info

  call read_csv("data/seatbelts.csv", names, data)
  allocate (y(1, n))
  y(1, :) = log(column(names, data, "drivers"))
  opts%factr = 10.0_dp
  opts%pgtol = 1.0e-9_dp

  ! Basic structural model: level, trigonometric seasonal, irregular
  allocate (irregular_t :: comps(1)%c)
  allocate (level_t :: comps(2)%c)
  comps(3)%c = seasonal_t(period=12, form=SEASONAL_TRIG)
  mod = structural_model(y, comps(1:3), info)
  if (info /= SS_OK) error stop "model"
  call fit(mod, res, options=opts, info=info)
  if (info /= SS_OK) error stop "fit"
  print '(a)', "Basic structural model (DK 8.2)"
  call report(mod, res)

  ! With the seat belt law and the log petrol price
  allocate (x(n, 2))
  x(:, 1) = log(column(names, data, "PetrolPrice"))
  x(:, 2) = step_intervention(n, t_law)
  comps(4)%c = regression_t(x=x)
  mod = structural_model(y, comps, info)
  if (info /= SS_OK) error stop "model"
  call fit(mod, res, options=opts, info=info)
  if (info /= SS_OK) error stop "fit"
  print '(/, a)', "With petrol price and the seat belt law"
  call report(mod, res)
  call coefficients(mod, res)

contains

  subroutine report(mod, res)
    type(structural_model_t), intent(inout) :: mod
    type(fit_result_t), intent(in) :: res
    type(filter_result_t) :: fres
    character(len=32), allocatable :: pn(:)
    real(dp) :: F(1, 1)
    integer :: i, info

    pn = mod%param_names()
    do i = 1, mod%k_params
      print '(2x, a20, es14.6, a, f9.6, a)', trim(pn(i)), res%params(i), &
          "   (q-ratio ", res%params(i) / res%params(1), ")"
    end do
    call mod%filter(res%params, fres, info)
    print '(2x, a, f12.3, a, i0, a)', "log likelihood ", res%llf, "   (d = ", &
        fres%nobs_diffuse, ")"
    call prediction_error_variance(mod%rep, F, info)
    if (info == SS_OK) then
      print '(2x, a, f12.8)', "prediction error variance ", F(1, 1)
    else
      print '(2x, a)', &
          "prediction error variance: none (time-varying Z, no steady state)"
    end if
  end subroutine report

  !> Regression coefficients: the smoothed (= final filtered) states, with
  !> their root mean square errors.
  subroutine coefficients(mod, res)
    type(structural_model_t), intent(inout) :: mod
    type(fit_result_t), intent(in) :: res
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    character(len=8), parameter :: label(2) = ["petrol  ", "law 83.2"]
    integer :: i, j, info

    call mod%smooth(res%params, fres, sres, info)
    print '(2x, a8, 3a12)', "", "coef", "rmse", "t-value"
    do i = 1, 2
      j = mod%rep%k_states - 2 + i
      print '(2x, a8, 3f12.5)', label(i), sres%alphahat(j, n), sqrt(sres%V(j, j, n)), &
        sres%alphahat(j, n) / sqrt(sres%V(j, j, n))
    end do
  end subroutine coefficients
end program dk_8_2_seatbelt
