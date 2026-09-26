!> DK 3.10.3: modelling a survey series y_t (subject to survey error) and a
!> related, accurately measured series x_t together, as in Harvey and Chung
!> (2000) for UK unemployment: the bivariate local linear trend model (3.49)
!>
!>     (y_t, x_t)' = mu_t + eps_t,           eps_t ~ N(0, Sigma_eps),
!>     mu_t+1 = mu_t + nu_t + xi_t,          xi_t ~ N(0, Sigma_xi),
!>     nu_t+1 = nu_t + zeta_t,               zeta_t ~ N(0, Sigma_zeta),
!>
!> with full 2 x 2 covariance matrices. x_n is available a month before
!> y_n, which is handled as a missing observation (DK 4.10). The data are
!> simulated from the model; the example fits it by maximum likelihood and
!> compares the estimates of the level and slope of y with those from a
!> univariate local linear trend model for y alone.
!>
!> Run from the repo root:  fpm run --example dk_3_10_3_two_sources
program dk_3_10_3_two_sources
  use statespace
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none

  integer, parameter :: n = 120
  ! True parameters: lower triangles of Sigma_eps, Sigma_xi, Sigma_zeta.
  ! The survey is noisy, x accurate, and the trends closely related.
  real(dp), parameter :: truth(9) = [0.5_dp, 0.0_dp, 0.02_dp, &
                                     0.10_dp, 0.085_dp, 0.10_dp, &
                                     0.004_dp, 0.0035_dp, 0.004_dp]
  type(component_holder_t) :: comps(2)
  type(structural_model_t) :: biv, uni
  type(fit_result_t) :: res_b, res_u
  type(filter_result_t) :: fres
  type(smoother_result_t) :: s_b, s_u
  real(dp) :: y(2, n), alpha(4, n), eps(2), eta(4)
  real(dp), allocatable :: Hs(:, :), Qs(:, :)
  integer :: t, info, seed_size, i
  integer, allocatable :: seed(:)

  call random_seed(size=seed_size)
  seed = [(4242 + 11 * i, i=1, seed_size)]
  call random_seed(put=seed)

  ! The bivariate model, used first to simulate the data at the true values
  comps(1)%c = irregular_t(cov=COV_FULL)
  comps(2)%c = trend_t(cov_level=COV_FULL, cov_slope=COV_FULL)
  y = 0.0_dp
  biv = structural_model(y, comps, info)
  if (info /= SS_OK) error stop "model"
  call biv%update(truth)
  Hs = psd_sqrt(biv%rep%H(:, :, 1))
  Qs = psd_sqrt(biv%rep%Q(:, :, 1))
  alpha(:, 1) = [100.0_dp, 0.2_dp, 60.0_dp, 0.15_dp]   ! (mu_y, nu_y, mu_x, nu_x)
  do t = 1, n
    call draw_standard_normal(eps)
    call draw_standard_normal(eta)
    y(:, t) = matmul(biv%rep%Z(:, :, 1), alpha(:, t)) + matmul(Hs, eps)
    if (t < n) alpha(:, t + 1) = matmul(biv%rep%T(:, :, 1), alpha(:, t)) + &
                                 matmul(biv%rep%R(:, :, 1), matmul(Qs, eta))
  end do
  y(1, n) = ieee_value(1.0_dp, ieee_quiet_nan)   ! the survey is a month behind

  ! Bivariate fit
  biv = structural_model(y, comps, info)
  call fit(biv, res_b, info=info)
  if (info /= SS_OK) error stop "fit (bivariate)"
  call biv%smooth(res_b%params, fres, s_b, info)

  ! Univariate local linear trend for the survey alone
  uni = structural_model(y(1:1, :), comps, info)
  call fit(uni, res_u, info=info)
  if (info /= SS_OK) error stop "fit (univariate)"
  call uni%smooth(res_u%params, fres, s_u, info)

  print '(a)', "Bivariate estimates (true values in brackets):"
  call show("Sigma_eps  ", res_b%params(1:3), truth(1:3))
  call show("Sigma_xi   ", res_b%params(4:6), truth(4:6))
  call show("Sigma_zeta ", res_b%params(7:9), truth(7:9))
  print '(/, a)', "Smoothed level and slope of y, RMSE against the truth:"
  print '(2x, a, 2f10.4)', "bivariate (level, slope): ", rmse(s_b%alphahat(1, :), alpha(1, :)), &
    rmse(s_b%alphahat(2, :), alpha(2, :))
  print '(2x, a, 2f10.4)', "univariate:                ", rmse(s_u%alphahat(1, :), alpha(1, :)), &
    rmse(s_u%alphahat(2, :), alpha(2, :))
  print '(/, a, i0, a)', "Nowcast of y's level at t = ", n, " (survey not yet available):"
  print '(2x, a, f9.3)', "truth:      ", alpha(1, n)
  print '(2x, a, f9.3, a, f7.3)', "bivariate:  ", s_b%alphahat(1, n), " +- ", sqrt(s_b%V(1, 1, n))
  print '(2x, a, f9.3, a, f7.3)', "univariate: ", s_u%alphahat(1, n), " +- ", sqrt(s_u%V(1, 1, n))

contains

  subroutine show(label, est, tru)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: est(3), tru(3)

    print '(2x, a, 3(f8.4, " [", f6.4, "]"))', label, (est(i), tru(i), i=1, 3)
  end subroutine show

  real(dp) function rmse(a, b)
    real(dp), intent(in) :: a(:), b(:)

    rmse = sqrt(sum((a - b)**2) / size(a))
  end function rmse
end program dk_3_10_3_two_sources
