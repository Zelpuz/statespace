!> Timing benchmark, the Fortran side of bench/bench_statsmodels.py. Each
!> case is timed as the minimum over repeats, in milliseconds:
!>
!>   ll    local level, known init (a1 = 0, P1 = 1e6), H = 1, Q = 0.1
!>   bsm   level + trigonometric seasonal (12) + irregular, exact diffuse
!>   mv    p = 8, m = 3: Nelson-Siegel loadings, H = 0.01 I, T diagonal
!>   arma  ARMA(2, 1), stationary init
!>
!> for loglike, kalman_filter, and kalman_filter + state_smoother; then
!> 1000 local level fits of 100 observations each.
!>
!> Run from the repo root:  fpm run --profile release --example bench
!> (add --flag -fopenmp --link-flag -fopenmp for a parallel fit_many)
program bench
  use statespace
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none

  integer, parameter :: i8 = selected_int_kind(18)
  integer :: k
  integer, parameter :: sizes_small(2) = [10000, 100000], &
                        sizes_ll(3) = [10000, 100000, 1000000]

  call seed(1)
  print '(a8, a10, 3a14)', "case", "n", "loglike", "filter", "filter+smooth"
  do k = 1, size(sizes_ll)
    call run_ll(sizes_ll(k))
  end do
  do k = 1, size(sizes_small)
    call run_bsm(sizes_small(k))
  end do
  do k = 1, size(sizes_small)
    call run_mv(sizes_small(k))
  end do
  do k = 1, size(sizes_small)
    call run_arma(sizes_small(k))
  end do
  call run_many(1000, 100)

contains

  subroutine seed(s)
    integer, intent(in) :: s
    integer :: n, i
    integer, allocatable :: v(:)

    call random_seed(size=n)
    v = [(s + 31 * i, i=1, n)]
    call random_seed(put=v)
  end subroutine seed

  !> Random walk plus noise.
  function rw_noise(n) result(y)
    integer, intent(in) :: n
    real(dp) :: y(1, n)
    real(dp) :: e(n), u(n)
    integer :: t

    call draw_standard_normal(e)
    call draw_standard_normal(u)
    y(1, 1) = e(1)
    do t = 2, n
      y(1, t) = y(1, t - 1) + sqrt(0.1_dp) * e(t)
    end do
    y(1, :) = y(1, :) + u
  end function rw_noise

  real(dp) function now()
    integer(i8) :: c, r

    call system_clock(c, r)
    now = real(c, dp) / real(r, dp)
  end function now

  !> Time the three operations on rep, repeating each `reps` times.
  subroutine time_rep(name, rep, reps)
    character(len=*), intent(in) :: name
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: reps
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp) :: best(3), t0, llf
    integer :: i, info

    best = huge(1.0_dp)
    do i = 1, reps
      t0 = now()
      llf = loglike(rep, info)
      best(1) = min(best(1), now() - t0)
      t0 = now()
      call kalman_filter(rep, fres, info)
      best(2) = min(best(2), now() - t0)
      t0 = now()
      call kalman_filter(rep, fres, info)
      call state_smoother(rep, fres, sres, info)
      best(3) = min(best(3), now() - t0)
    end do
    if (info /= SS_OK) error stop "bench: "//name
    print '(a8, i10, 3f14.3)', name, rep%nobs, 1000 * best
  end subroutine time_rep

  integer function nreps(n)
    integer, intent(in) :: n

    nreps = max(3, min(20, 2000000 / n))
  end function nreps

  subroutine run_ll(n)
    integer, intent(in) :: n
    type(ssm_rep_t) :: rep

    rep = ssm_rep(rw_noise(n), m=1, r=1)
    rep%Z = 1.0_dp; rep%T = 1.0_dp; rep%R = 1.0_dp
    rep%H = 1.0_dp; rep%Q = 0.1_dp
    call rep%initialize_known([0.0_dp], reshape([1.0e6_dp], [1, 1]))
    call time_rep("ll", rep, nreps(n))
  end subroutine run_ll

  subroutine run_bsm(n)
    integer, intent(in) :: n
    type(component_holder_t) :: comps(3)
    type(structural_model_t) :: model
    integer :: info

    allocate (irregular_t :: comps(1)%c)
    allocate (level_t :: comps(2)%c)
    comps(3)%c = seasonal_t(period=12, form=SEASONAL_TRIG)
    model = structural_model(rw_noise(n), comps, info)
    call model%update([1.0_dp, 0.1_dp, 0.01_dp])
    call time_rep("bsm", model%rep, nreps(n))
  end subroutine run_bsm

  subroutine run_mv(n)
    integer, intent(in) :: n
    real(dp), parameter :: tau(8) = [3, 6, 12, 24, 36, 60, 84, 120], lambda = 0.06_dp
    type(ssm_rep_t) :: rep
    real(dp), allocatable :: y(:, :), e(:)
    real(dp) :: zi(8)

    allocate (y(8, n), e(8 * n))
    call draw_standard_normal(e)
    y = reshape(e, [8, n])
    rep = ssm_rep(y, m=3, r=3)
    zi = exp(-lambda * tau)
    rep%Z(:, 1, 1) = 1.0_dp
    rep%Z(:, 2, 1) = (1 - zi) / (lambda * tau)
    rep%Z(:, 3, 1) = rep%Z(:, 2, 1) - zi
    rep%H(:, :, 1) = 0.01_dp * eye(8)
    rep%T(:, :, 1) = 0.0_dp
    rep%T(1, 1, 1) = 0.99_dp; rep%T(2, 2, 1) = 0.95_dp; rep%T(3, 3, 1) = 0.9_dp
    rep%R(:, :, 1) = eye(3)
    rep%Q(:, :, 1) = 0.1_dp * eye(3)
    call rep%initialize_known(spread(0.0_dp, 1, 3), eye(3))
    call time_rep("mv", rep, nreps(n))
    rep%tol_steady = -1.0_dp
    call time_rep("mv-full", rep, nreps(n))
    rep%tol_steady = 1.0e-15_dp
    rep%filter_method = FILTER_UNIVARIATE
    call time_rep("mv-uv", rep, nreps(n))
    rep%filter_method = FILTER_CONVENTIONAL
    call random_number(e)
    where (reshape(e, [8, n]) < 0.1_dp) rep%y = ieee_value(1.0_dp, ieee_quiet_nan)
    call time_rep("mv-miss", rep, nreps(n))
  end subroutine run_mv

  subroutine run_arma(n)
    integer, intent(in) :: n
    type(component_holder_t) :: comps(1)
    type(structural_model_t) :: model
    real(dp), allocatable :: e(:)
    real(dp) :: y(1, n)
    integer :: info

    allocate (e(n))
    call draw_standard_normal(e)
    y(1, :) = e
    comps(1)%c = arima_t(ar=2, ma=1)
    model = structural_model(y, comps, info)
    call model%update([0.5_dp, 0.2_dp, 0.3_dp, 1.0_dp])
    call time_rep("arma", model%rep, nreps(n))
  end subroutine run_arma

  !> nseries local level fits (exact diffuse), n observations each.
  subroutine run_many(nseries, n)
    integer, intent(in) :: nseries, n
    type(component_holder_t) :: comps(2)
    type(structural_model_t) :: model
    type(fit_result_t) :: res
    real(dp) :: t0, y(1, n, nseries)
    integer :: i, info, nfail

    do i = 1, nseries
      y(:, :, i) = rw_noise(n)
    end do
    allocate (irregular_t :: comps(1)%c)
    allocate (level_t :: comps(2)%c)
    nfail = 0
    t0 = now()
    do i = 1, nseries
      model = structural_model(y(:, :, i), comps, info)
      call fit(model, res, info=info)
      if (info /= SS_OK) nfail = nfail + 1
    end do
    print '(/, a, i0, a, i0, a, f10.1, a, i0, a)', "fit ", nseries, &
        " local level series of n = ", n, ": ", 1000 * (now() - t0), " ms (", nfail, &
        " without standard errors)"
    call run_fit_many(y, comps)
  end subroutine run_many

  !> The same with fit_many (parallel when built with -fopenmp).
  subroutine run_fit_many(y, comps)
    real(dp), intent(in) :: y(:, :, :)
    type(component_holder_t), intent(in) :: comps(:)
    type(structural_model_t), allocatable :: models(:)
    type(fit_result_t), allocatable :: res(:)
    integer, allocatable :: info(:)
    real(dp) :: t0
    integer :: i, stat

    allocate (models(size(y, 3)), res(size(y, 3)), info(size(y, 3)))
    t0 = now()
    do i = 1, size(y, 3)
      models(i) = structural_model(y(:, :, i), comps, stat)
    end do
    call fit_many(models, res, info)
    print '(a, f10.1, a)', "fit_many, same series: ", 1000 * (now() - t0), " ms"
  end subroutine run_fit_many
end program bench
