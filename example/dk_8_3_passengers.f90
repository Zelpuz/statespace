!> DK 8.3: bivariate structural model for the log of monthly UK front and
!> rear seat passengers killed or seriously injured. Each series has a
!> local level, a fixed trigonometric seasonal and an irregular, with the
!> level and irregular disturbances correlated across the series (SUTSE,
!> DK 3.3).
!>
!>   1. Before 1983 (168 months), full-rank level covariance.
!>   2. The same with a rank-one level covariance: one common level with
!>      a free loading, plus an intercept for the rear seat series (DK 3.3.2).
!>   3. All 192 months with a level intervention for the seat belt law
!>      (February 1983) in both series.
!>   4. The intervention in the front seat series only; the rear seat
!>      series acts as a control (Harvey 1996).
!>   5. As 4 with a rank-one level covariance.
!>
!> DK's text also mentions km travelled and the petrol price as regressors,
!> but their printed estimates are those of the model without them: with
!> them the estimates differ, and without them all five models match DK to
!> every printed digit. DK's printed exponents are also one too small (the
!> irregular variances are x 1e-3, the level variances x 1e-4), and model
!> 1's level covariance is 2.933, not 2.993 (their correlation, 0.893,
!> implies 2.933); likewise model 3's irregular covariance is 4.593, not
!> 4.493 (correlation 0.660). DK's values, on those scales and corrected:
!>   1. Sigma_eps [5.006 4.569; 9.143], Sigma_eta [4.834 2.933; 2.234]
!>   2. Sigma_eps [5.062 4.791; 10.02], Sigma_eta [4.802 2.792; 1.623],
!>      log likelihood 2.689 below model 1
!>   3. Sigma_eps [5.135 4.593; 9.419], Sigma_eta [4.896 3.025; 2.317],
!>      law: front -0.32799 (rmse 0.05699), rear 0.03376 (0.05025)
!>   4. Sigma_eps [5.147 4.588; 9.380], Sigma_eta [4.754 2.926; 2.282],
!>      law: front -0.35630 (0.03655)
!>   5. Sigma_eps [5.206 4.789; 10.24], Sigma_eta [4.970 2.860; 1.646],
!>      law: front -0.41557 (0.02621)
!>
!> Run from the repo root, after `.venv/bin/python data/fetch_dk_data.py`:
!>     fpm run --example dk_8_3_passengers
program dk_8_3_passengers
  use statespace
  use csv_io, only: read_csv, column
  implicit none

  integer, parameter :: n_all = 192, n_pre = 168, t_law = 170
  character(len=32), allocatable :: names(:)
  real(dp), allocatable :: data(:, :), y(:, :)
  real(dp) :: llf1, llf2

  call read_csv("data/seatbelts.csv", names, data)
  allocate (y(2, n_all))
  y(1, :) = log(column(names, data, "front"))
  y(2, :) = log(column(names, data, "rear"))

  print '(a)', "1. Before 1983, full-rank level covariance"
  llf1 = estimate(n_pre, rank_one=.false., law=[.false., .false.])
  print '(/, a)', "2. Before 1983, rank-one level covariance"
  llf2 = estimate(n_pre, rank_one=.true., law=[.false., .false.])
  print '(2x, a, f10.3)', "log likelihood below model 1: ", llf1 - llf2
  print '(/, a)', "3. All data, law in both series"
  llf1 = estimate(n_all, rank_one=.false., law=[.true., .true.])
  print '(/, a)', "4. All data, law in the front seat series"
  llf1 = estimate(n_all, rank_one=.false., law=[.true., .false.])
  print '(/, a)', "5. As 4, rank-one level covariance"
  llf1 = estimate(n_all, rank_one=.true., law=[.true., .false.])

contains

  !> Fit the model to the first n months, print the covariances and the law
  !> coefficients, and return the log likelihood.
  real(dp) function estimate(n, rank_one, law) result(llf)
    integer, intent(in) :: n
    logical, intent(in) :: rank_one, law(2)
    type(component_holder_t), allocatable :: comps(:)
    type(structural_model_t) :: mod
    type(fit_result_t) :: res
    type(fit_options_t) :: opts
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp) :: Seps(2, 2), Seta(2, 2), lam
    integer :: info, i, j, nc, ireg(2)
    character(len=5), parameter :: label(2) = ["front", "rear "]

    ! Irregular, level, seasonal, then a regression component per series
    ! that needs one: the law, and in the rank-one model the rear seat
    ! intercept.
    nc = 3 + count(law .or. [.false., rank_one])
    allocate (comps(nc))
    comps(1)%c = irregular_t(cov=COV_FULL)
    if (rank_one) then
      ! One common level theta_t, on which the rear seat series loads with a
      ! free loading (DK 3.3.2)
      comps(2)%c = level_t()
      comps(3)%c = seasonal_t(period=12, form=SEASONAL_TRIG, cov=COV_NONE, at_observations=.true.)
    else
      comps(2)%c = level_t(cov=COV_FULL)
      comps(3)%c = seasonal_t(period=12, form=SEASONAL_TRIG, cov=COV_NONE)
    end if
    ireg = 0
    j = 3
    do i = 1, 2
      if (.not. (law(i) .or. (rank_one .and. i == 2))) cycle
      j = j + 1
      ireg(i) = j
      comps(j)%c = regression_t(x=regressors(n, law(i), rank_one .and. i == 2), series=i, &
                                at_observations=rank_one)
    end do
    if (rank_one) then
      mod = structural_model(y(:, 1:n), comps, info, loading=reshape([1.0_dp, 1.0_dp], [2, 1]), &
                             loading_free=reshape([.false., .true.], [2, 1]))
    else
      mod = structural_model(y(:, 1:n), comps, info)
    end if
    if (info /= SS_OK) error stop "model"
    opts%factr = 10.0_dp
    opts%pgtol = 1.0e-8_dp
    opts%compute_cov = .false.
    call fit(mod, res, options=opts, info=info)
    if (info /= SS_OK) error stop "fit"
    llf = res%llf

    Seps = full_cov(res%params(1:3))
    if (rank_one) then
      lam = res%params(mod%k_params)
      Seta = res%params(4) * reshape([1.0_dp, lam, lam, lam**2], [2, 2])
    else
      Seta = full_cov(res%params(4:6))
    end if
    print '(2x, a, 3f8.3, a, f6.3)', "Sigma_eps x 1e3: ", 1.0e3_dp * [Seps(1, 1), Seps(2, 1), Seps(2, 2)], &
      "   rho ", Seps(2, 1) / sqrt(Seps(1, 1) * Seps(2, 2))
    print '(2x, a, 3f8.3, a, f6.3)', "Sigma_eta x 1e4: ", 1.0e4_dp * [Seta(1, 1), Seta(2, 1), Seta(2, 2)], &
      "   rho ", Seta(2, 1) / sqrt(Seta(1, 1) * Seta(2, 2))
    print '(2x, a, f12.3)', "log likelihood ", llf

    if (.not. any(law)) return
    call mod%smooth(res%params, fres, sres, info)
    print '(2x, a8, 3a12)', "law", "coef", "rmse", "t-value"
    do i = 1, 2
      if (.not. law(i)) cycle
      j = mod%s0(ireg(i)) + mod%comps(ireg(i))%c%m     ! the law is the last regressor
      print '(2x, a8, 3f12.5)', label(i), sres%alphahat(j, n), sqrt(sres%V(j, j, n)), &
        sres%alphahat(j, n) / sqrt(sres%V(j, j, n))
    end do
  end function estimate

  !> The level intervention and/or an intercept.
  function regressors(n, law, intercept) result(x)
    integer, intent(in) :: n
    logical, intent(in) :: law, intercept
    real(dp), allocatable :: x(:, :)

    allocate (x(n, 0))
    if (intercept) x = reshape([spread(1.0_dp, 1, n)], [n, 1])
    if (law) x = reshape([pack(x, .true.), step_intervention(n, t_law)], [n, size(x, 2) + 1])
  end function regressors

  !> A 2 x 2 covariance from its lower triangle (s11, s21, s22).
  pure function full_cov(c) result(S)
    real(dp), intent(in) :: c(3)
    real(dp) :: S(2, 2)

    S = reshape([c(1), c(2), c(2), c(3)], [2, 2])
  end function full_cov
end program dk_8_3_passengers
