!> DK 3.10.2: benchmarking monthly survey data to error-free annual totals,
!> with DK's state space form (3.48):
!>
!>     y_t = mu_t + gamma_t + delta_t w_t + eps_t + sigma^s_t xi^s_t,
!>
!> an integrated random walk trend, the dummy seasonal (3.3), a calendar
!> effect delta_t w_t whose coefficient changes once a year (in January),
!> the irregular, and an AR(1) survey error with known standard deviation
!> sigma^s_t. The true values are y*_t = y_t - sigma^s_t xi^s_t, and the
!> benchmarks x_i = sum of y*_t over year i. The series is arranged as
!>
!>     y_1, ..., y_12, x_1, y_13, ..., y_24, x_2, ...,
!>
!> with the state
!>
!>     alpha_t = (mu_t..mu_t-11, gamma_t..gamma_t-11, delta_t,
!>                eps_t..eps_t-11, xi^s_t)',
!>
!> so that both kinds of observation are exact linear functions of the
!> state (H = 0). From month 12 to the benchmark point the transition is
!> the identity; from the benchmark point to January it is the monthly one
!> plus the yearly change of delta. The data are simulated from the model
!> itself, with known parameters.
!>
!> The example checks that the smoothed true values add up to the
!> benchmarks exactly, and compares their error with and without the
!> benchmarks.
!>
!> Run from the repo root:  fpm run --example dk_3_10_2_benchmarking
program dk_3_10_2_benchmarking
  use statespace
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none

  integer, parameter :: years = 10, n = 13 * years
  integer, parameter :: im = 0, ig = 12, id = 24, ie = 25, ix = 38, m = 38, r = 5
  real(dp), parameter :: s2_trend = 0.01_dp, s2_seas = 0.05_dp, s2_delta = 0.01_dp, &
                         s2_eps = 0.25_dp, phi = 0.7_dp
  type(ssm_rep_t) :: rep
  real(dp) :: alpha(m, n + 1), eta(r), w(n), sig(n), ystar(n), yhat(n, 2), bench_err
  real(dp) :: Pinf(m, m), Pstar(m, m)
  real(dp), allocatable :: y(:, :)
  type(filter_result_t) :: fres
  type(smoother_result_t) :: sres
  integer :: t, i, k, info, seed_size
  integer, allocatable :: seed(:)
  logical :: is_bench(n)

  call random_seed(size=seed_size)
  seed = [(12345 + 7 * i, i=1, seed_size)]
  call random_seed(put=seed)

  is_bench = [(mod(t, 13) == 0, t=1, n)]
  do t = 1, n
    call random_number(w(t))
    w(t) = 2 * w(t) - 1                       ! calendar variable
    sig(t) = 1.0_dp + 0.5_dp * sin(0.3_dp * t) ! known survey error sd
  end do

  allocate (y(1, n))
  y = 0.0_dp
  rep = ssm_rep(y, m=m, r=r)
  deallocate (rep%Z, rep%T, rep%R)
  allocate (rep%Z(1, m, n), rep%T(m, m, n), rep%R(m, r, n), source=0.0_dp)
  rep%H = 0.0_dp
  rep%Q(:, :, 1) = 0.0_dp
  rep%Q(1, 1, 1) = s2_trend; rep%Q(2, 2, 1) = s2_seas; rep%Q(3, 3, 1) = s2_delta
  rep%Q(4, 4, 1) = s2_eps; rep%Q(5, 5, 1) = 1 - phi**2

  do t = 1, n
    if (is_bench(t)) then
      ! x_i = sum over the year of mu + gamma + delta w + eps
      rep%Z(1, im + 1:im + 12, t) = 1.0_dp
      rep%Z(1, ig + 1:ig + 12, t) = 1.0_dp
      rep%Z(1, id + 1, t) = sum(w(t - 12:t - 1))
      rep%Z(1, ie + 1:ie + 12, t) = 1.0_dp
      call monthly(t, delta_changes=.true.)
    else
      rep%Z(1, im + 1, t) = 1.0_dp
      rep%Z(1, ig + 1, t) = 1.0_dp
      rep%Z(1, id + 1, t) = w(t)
      rep%Z(1, ie + 1, t) = 1.0_dp
      rep%Z(1, ix, t) = sig(t)
      if (t < n) then
        if (is_bench(t + 1)) then
          rep%T(:, :, t) = eye(m)                ! month 12 -> benchmark
        else
          call monthly(t, delta_changes=.false.)
        end if
      end if
    end if
  end do

  ! Simulate: alpha_1 with a trend, a seasonal pattern and random errors
  alpha(:, 1) = 0.0_dp
  do k = 0, 11
    alpha(im + 1 + k, 1) = 100.0_dp - 0.5_dp * k
    alpha(ig + 1 + k, 1) = 3.0_dp * cos(2 * acos(-1.0_dp) * k / 12)
  end do
  alpha(id + 1, 1) = 1.0_dp
  call draw_standard_normal(alpha(ie + 1:ie + 12, 1))
  alpha(ie + 1:ie + 12, 1) = sqrt(s2_eps) * alpha(ie + 1:ie + 12, 1)
  call draw_standard_normal(alpha(ix:ix, 1))
  do t = 1, n
    call draw_standard_normal(eta)
    eta = sqrt([s2_trend, s2_seas, s2_delta, s2_eps, 1 - phi**2]) * eta
    y(1, t) = dot_product(rep%Z(1, :, t), alpha(:, t))
    ystar(t) = y(1, t) - merge(0.0_dp, sig(t) * alpha(ix, t), is_bench(t))
    alpha(:, t + 1) = matmul(rep%T(:, :, t), alpha(:, t)) + matmul(rep%R(:, :, t), eta)
  end do
  rep%y = y

  ! Trend, seasonal and calendar states diffuse; eps lags and xi^s known
  ! (mean 0, variances sigma2_eps and 1)
  Pinf = 0.0_dp; Pstar = 0.0_dp
  do k = 1, id + 1
    Pinf(k, k) = 1.0_dp
  end do
  do k = ie + 1, ie + 12
    Pstar(k, k) = s2_eps
  end do
  Pstar(ix, ix) = 1.0_dp
  call rep%initialize_general(spread(0.0_dp, 1, m), Pstar, Pinf)

  ! With and without the benchmarks
  do k = 1, 2
    if (k == 2) where (is_bench) rep%y(1, :) = ieee_value(1.0_dp, ieee_quiet_nan)
    call kalman_filter(rep, fres, info)
    if (info /= SS_OK) error stop "filter"
    call state_smoother(rep, fres, sres, info)
    if (info /= SS_OK) error stop "smoother"
    do t = 1, n
      ! smoothed y*_t = Z_t alpha_hat_t without the survey error
      yhat(t, k) = dot_product(rep%Z(1, 1:ix - 1, t), sres%alphahat(1:ix - 1, t))
    end do
    if (k == 1) then
      bench_err = 0.0_dp
      do i = 1, years
        t = 13 * i
        bench_err = max(bench_err, abs(sum(yhat(t - 12:t - 1, 1)) - y(1, t)))
      end do
    end if
  end do

  print '(a, es10.2)', "largest |sum of smoothed y* over a year - benchmark|: ", &
      bench_err
  print '(a, f8.4)', "RMSE of smoothed y*, with benchmarks:    ", rmse(yhat(:, 1))
  print '(a, f8.4)', "RMSE of smoothed y*, without benchmarks: ", rmse(yhat(:, 2))
  print '(a, f8.4)', "RMSE of the raw survey values:           ", rmse(y(1, :))

contains

  !> Monthly transition from t to t+1, with the yearly change of delta
  !> (noise zeta) when leaving a benchmark point.
  subroutine monthly(t, delta_changes)
    integer, intent(in) :: t
    logical, intent(in) :: delta_changes
    integer :: j

    rep%T(im + 1, im + 1, t) = 2.0_dp            ! integrated random walk
    rep%T(im + 1, im + 2, t) = -1.0_dp
    rep%T(ig + 1, ig + 1:ig + 11, t) = -1.0_dp   ! dummy seasonal (3.3)
    do j = 1, 11
      rep%T(im + 1 + j, im + j, t) = 1.0_dp      ! lags
      rep%T(ig + 1 + j, ig + j, t) = 1.0_dp
      rep%T(ie + 1 + j, ie + j, t) = 1.0_dp
    end do
    rep%T(id + 1, id + 1, t) = 1.0_dp
    rep%T(ix, ix, t) = phi
    rep%R(im + 1, 1, t) = 1.0_dp
    rep%R(ig + 1, 2, t) = 1.0_dp
    if (delta_changes) rep%R(id + 1, 3, t) = 1.0_dp
    rep%R(ie + 1, 4, t) = 1.0_dp                 ! eps_t+1, a new draw
    rep%R(ix, 5, t) = 1.0_dp
  end subroutine monthly

  !> RMSE against the true y* over the survey months.
  real(dp) function rmse(v)
    real(dp), intent(in) :: v(:)

    rmse = sqrt(sum((v - ystar)**2, mask=.not. is_bench) / count(.not. is_bench))
  end function rmse
end program dk_3_10_2_benchmarking
