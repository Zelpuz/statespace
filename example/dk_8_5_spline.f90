!> DK 8.5: cubic smoothing spline for the motorcycle acceleration data (133
!> observations, irregularly spaced in time, some at the same time point),
!> as the continuous-time smooth trend model (DK 3.8.2, 3.9.2, eq. 8.1)
!> plus an irregular. The smoothness parameter lambda = sigma2_zeta /
!> sigma2_eps is estimated by maximum likelihood.
!>
!> DK report psi = log lambda = -3.59 (standard error 0.22), lambda =
!> 0.0275 with a 95% interval of 0.018 to 0.043, and AIC 9.43 (DK 7.4 with
!> diffuse initialization: n^-1 [-2 log L_d + 2 (q + w)], q = 2 diffuse
!> states, w = 2 parameters).
!>
!> Our AIC is 9.42, matching DK's 9.43, but our estimate is psi = -2.36
!> (se 0.44), lambda = 0.094. The profile log likelihood over psi has a
!> single maximum there; at DK's psi = -3.59 it is about 4.3 lower (-626.9),
!> which would give AIC 9.49. So DK's AIC is that of this optimum, and their
!> printed lambda seems to come from a different scaling that we couldn't
!> identify. The spline model itself is checked against SciPy's
!> smoothing spline in the tests.
!>
!> Run from the repo root, after `.venv/bin/python data/fetch_dk_data.py`:
!>     fpm run --example dk_8_5_spline
program dk_8_5_spline
  use statespace
  use csv_io, only: read_csv, column
  implicit none

  character(len=32), allocatable :: names(:)
  real(dp), allocatable :: data(:, :), times(:), y(:, :), g(:)
  type(component_holder_t) :: comps(2)
  type(structural_model_t) :: mod
  type(fit_result_t) :: res
  type(fit_options_t) :: opts
  type(filter_result_t) :: fres
  type(smoother_result_t) :: sres
  real(dp) :: psi, se, aic
  integer :: info, n, t

  call read_csv("data/mcycle.csv", names, data)
  times = column(names, data, "times")
  n = size(times)
  allocate (y(1, n))
  y(1, :) = column(names, data, "accel")

  allocate (irregular_t :: comps(1)%c)
  comps(2)%c = continuous_trend_t(times=times)
  mod = structural_model(y, comps, info)
  if (info /= SS_OK) error stop "model"
  opts%factr = 10.0_dp
  opts%pgtol = 1.0e-9_dp
  call fit(mod, res, options=opts, info=info)
  if (info /= SS_OK) error stop "fit"

  ! psi = log sigma2_zeta - log sigma2_eps, its variance by the delta method
  psi = log(res%params(2) / res%params(1))
  g = [-1.0_dp / res%params(1), 1.0_dp / res%params(2)]
  se = sqrt(dot_product(g, matmul(res%cov_params, g)))
  aic = (-2 * res%llf + 2 * (mod%rep%k_diffuse() + mod%k_params)) / n
  print '(a, 2es12.4)', "sigma2_eps, sigma2_zeta: ", res%params
  print '(a, f8.3, a, f6.3, a)', "psi = log lambda: ", psi, "  (se ", se, ")"
  print '(a, f8.4, a, f6.3, a, f6.3)', "lambda: ", exp(psi), "  95% interval ", &
    exp(psi - 1.96_dp * se), " to ", exp(psi + 1.96_dp * se)
  print '(a, f8.3)', "AIC: ", aic

  ! The spline (the smoothed level) with 95% intervals, every 10th point
  call mod%smooth(res%params, fres, sres, info)
  print '(/, a8, 4a10)', "time", "accel", "spline", "lower", "upper"
  do t = 1, n, 10
    print '(f8.1, 4f10.2)', times(t), y(1, t), sres%alphahat(1, t), &
      sres%alphahat(1, t) - 1.96_dp * sqrt(sres%V(1, 1, t)), &
      sres%alphahat(1, t) + 1.96_dp * sqrt(sres%V(1, 1, t))
  end do
end program dk_8_5_spline
