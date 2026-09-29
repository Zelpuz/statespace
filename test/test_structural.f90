!> Structural components (DK 3.2-3.3) against statsmodels' UnobservedComponents.
module test_structural
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use statespace
  use fixture_io, only: fixture_t, load_fixture
  use statespace_linalg, only: solve
  implicit none
  private

  public :: collect_structural

  real(dp), parameter :: rtol = 1.0e-9_dp

contains

  subroutine collect_structural(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ new_unittest("basic_structural_model", test_bsm), &
                 new_unittest("trigonometric_seasonal", test_trig), &
                 new_unittest("smooth_trend", test_smooth_trend), &
                 new_unittest("damped_cycle", test_cycle), &
                 new_unittest("regression_and_intervention", test_regression), &
                 new_unittest("random_walk_regression", test_rw_regression), &
                 new_unittest("sutse_independent_series", test_sutse_diagonal), &
                 new_unittest("full_covariance_transform", test_full_cov), &
                 new_unittest("fixed_seasonals_marginal", test_fixed_seasonals), &
                 new_unittest("common_levels", test_common_levels), &
                 new_unittest("latent_risk", test_latent_risk), &
                 new_unittest("dynamic_factor", test_dynamic_factor), &
                 new_unittest("continuous_level_unit_spacing", test_continuous_unit), &
                 new_unittest("continuous_level_interpolation", &
                              test_continuous_interp), &
                 new_unittest("continuous_spline_vs_scipy", test_cubic_spline), &
                 new_unittest("discrete_spline_penalized_ls", test_discrete_spline) ]
  end subroutine collect_structural

  subroutine check_close(error, actual, expected, label)
    type(error_type), allocatable, intent(out) :: error
    real(dp), intent(in) :: actual(:), expected(:)
    character(len=*), intent(in) :: label
    real(dp) :: err
    character(len=32) :: buf

    err = maxval(abs(actual - expected)) / max(1.0_dp, maxval(abs(expected)))
    write (buf, '(es10.3)') err
    call check(error, err <= rtol, label//": relative error "//trim(buf))
  end subroutine check_close

  !> Smoothed signal of component i: its block of Z times its block of alphahat.
  function signal(model, sres, i) result(s)
    type(structural_model_t), intent(in) :: model
    type(smoother_result_t), intent(in) :: sres
    integer, intent(in) :: i
    real(dp), allocatable :: s(:)
    integer :: t, a, m, iz

    a = model%s0(i)
    m = model%comps(i)%c%m
    allocate (s(model%rep%nobs))
    do t = 1, model%rep%nobs
      iz = min(t, size(model%rep%Z, 3))
      s(t) = dot_product(model%rep%Z(1, a + 1:a + m, iz), sres%alphahat(a + 1:a + m, t))
    end do
  end function signal

  !> Assemble, smooth at the fixture's parameters, and compare llf_obs.
  subroutine run(error, model, comps, fx, params, fres, sres)
    type(error_type), allocatable, intent(out) :: error
    type(structural_model_t), intent(out) :: model
    type(component_holder_t), intent(in) :: comps(:)
    type(fixture_t), intent(in) :: fx
    real(dp), intent(in) :: params(:)
    type(filter_result_t), intent(out) :: fres
    type(smoother_result_t), intent(out) :: sres
    integer :: info

    model = structural_model(fx%get2('y'), comps, info)
    call check(error, info, SS_OK, "structural_model info")
    if (allocated(error)) return
    call check(error, model%k_params, size(params), "number of parameters")
    if (allocated(error)) return
    call model%smooth(params, fres, sres, info)
    call check(error, info, SS_OK, "smooth info")
    if (allocated(error)) return
    call check_close(error, fres%llf_obs, fx%get1('llf_obs'), "llf_obs")
  end subroutine run

  subroutine test_bsm(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: comps(3)
    type(structural_model_t) :: model
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres

    fx = load_fixture("test/fixtures/uc_bsm_dummy.txt")
    allocate (irregular_t :: comps(1)%c)
    allocate (trend_t :: comps(2)%c)
    comps(3)%c = seasonal_t(period=4, form=SEASONAL_DUMMY)
    call run(error, model, comps, fx, fx%get1('params'), fres, sres)
    if (allocated(error)) return
    call check_close(error, sres%alphahat(model%s0(2) + 1, :), fx%get1('level'), &
                     "level")
    if (allocated(error)) return
    call check_close(error, sres%alphahat(model%s0(2) + 2, :), fx%get1('trend'), &
                     "slope")
    if (allocated(error)) return
    call check_close(error, signal(model, sres, 3), fx%get1('seasonal'), "seasonal")
  end subroutine test_bsm

  subroutine test_trig(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: comps(3)
    type(structural_model_t) :: model
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres

    fx = load_fixture("test/fixtures/uc_llevel_trig12.txt")
    allocate (irregular_t :: comps(1)%c)
    allocate (level_t :: comps(2)%c)
    comps(3)%c = seasonal_t(period=12, form=SEASONAL_TRIG)
    call run(error, model, comps, fx, fx%get1('params'), fres, sres)
    if (allocated(error)) return
    call check(error, model%comps(3)%c%m, 11, "s - 1 trigonometric states (DK 3.2.3)")
    if (allocated(error)) return
    call check_close(error, signal(model, sres, 3), fx%get1('seasonal'), "seasonal")
  end subroutine test_trig

  subroutine test_smooth_trend(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: comps(3)
    type(structural_model_t) :: model
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres

    fx = load_fixture("test/fixtures/uc_smooth_trend_dummy.txt")
    allocate (irregular_t :: comps(1)%c)
    comps(2)%c = trend_t(cov_level=COV_NONE)
    comps(3)%c = seasonal_t(period=4)
    call run(error, model, comps, fx, fx%get1('params'), fres, sres)
    if (allocated(error)) return
    call check_close(error, sres%alphahat(model%s0(2) + 1, :), fx%get1('level'), &
                     "level")
    if (allocated(error)) return
    call check_close(error, signal(model, sres, 3), fx%get1('seasonal'), "seasonal")
  end subroutine test_smooth_trend

  subroutine test_cycle(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: comps(3)
    type(structural_model_t) :: model
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres

    fx = load_fixture("test/fixtures/uc_llevel_cycle.txt")
    allocate (irregular_t :: comps(1)%c)
    allocate (level_t :: comps(2)%c)
    allocate (cycle_t :: comps(3)%c)
    call run(error, model, comps, fx, fx%get1('params'), fres, sres)
    if (allocated(error)) return
    call check_close(error, signal(model, sres, 3), fx%get1('cycle'), "cycle")
  end subroutine test_cycle

  subroutine test_regression(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: comps(3)
    type(structural_model_t) :: model
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: x(:, :)

    fx = load_fixture("test/fixtures/uc_llevel_regression.txt")
    x = fx%get2('x')
    allocate (irregular_t :: comps(1)%c)
    allocate (level_t :: comps(2)%c)
    ! The intervention is rebuilt with step_intervention (DK 3.2.5).
    comps(3)%c = regression_t(x=reshape([x(:, 1), step_intervention(size(x, 1), 50)], &
                                        shape(x)))
    call run(error, model, comps, fx, fx%get1('params'), fres, sres)
    if (allocated(error)) return
    call check_close(error, sres%alphahat(model%s0(2) + 1, :), fx%get1('level'), &
                     "level")
    if (allocated(error)) return
    call check_close(error, &
                     sres%alphahat(model%s0(3) + 1:model%s0(3) + 2, model%rep%nobs), &
                     fx%get1('beta'), "coefficients")
  end subroutine test_regression

  subroutine test_rw_regression(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: comps(3)
    type(structural_model_t) :: model
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres

    fx = load_fixture("test/fixtures/uc_llevel_rw_regression.txt")
    allocate (irregular_t :: comps(1)%c)
    allocate (level_t :: comps(2)%c)
    comps(3)%c = regression_t(x=fx%get2('x'), random_walk=[.true.])
    call run(error, model, comps, fx, [0.25_dp, 0.09_dp, 0.01_dp], fres, sres)
    if (allocated(error)) return
    call check_close(error, pack(sres%alphahat, .true.), fx%get1('alphahat'), &
                     "alphahat")
  end subroutine test_rw_regression

  !> SUTSE (DK 3.3) with diagonal covariances: independent series, so the
  !> log likelihood is the sum of the univariate ones.
  subroutine test_sutse_diagonal(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fa, fb
    type(component_holder_t) :: comps(3)
    type(structural_model_t) :: mv, ua, ub
    real(dp), allocatable :: ya(:, :), yb(:, :), y2(:, :)
    real(dp) :: la, lb, l2
    integer :: info

    fa = load_fixture("test/fixtures/uc_bsm_dummy.txt")
    fb = load_fixture("test/fixtures/uc_smooth_trend_dummy.txt")
    ya = fa%get2('y')
    yb = ya * 0.5_dp + 1.0_dp
    yb(1, 30) = ieee_value_nan()
    allocate (y2(2, size(ya, 2)))
    y2(1, :) = ya(1, :)
    y2(2, :) = yb(1, :)
    allocate (irregular_t :: comps(1)%c)
    allocate (trend_t :: comps(2)%c)
    comps(3)%c = seasonal_t(period=4)
    mv = structural_model(y2, comps, info)
    ua = structural_model(ya, comps, info)
    ub = structural_model(yb, comps, info)
    ! Parameters: per component, series 1 then 2.
    l2 = mv%loglike([0.3_dp, 0.2_dp, 0.1_dp, 0.05_dp, 0.003_dp, 0.002_dp, 0.01_dp, &
                     0.02_dp], info)
    la = ua%loglike([0.3_dp, 0.1_dp, 0.003_dp, 0.01_dp], info)
    lb = ub%loglike([0.2_dp, 0.05_dp, 0.002_dp, 0.02_dp], info)
    call check_close(error, [l2], [la + lb], "sum of univariate llf")

  contains

    real(dp) function ieee_value_nan()
      ieee_value_nan = 0.0_dp
      ieee_value_nan = ieee_value_nan / ieee_value_nan
    end function ieee_value_nan
  end subroutine test_sutse_diagonal

  !> Full covariance parameters (lower triangle of Sigma) round-trip through
  !> the Cholesky transform, and the assembled Q is that Sigma.
  subroutine test_full_cov(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: comps(2)
    type(structural_model_t) :: model
    real(dp), allocatable :: y2(:, :), p(:)
    integer :: info

    fx = load_fixture("test/fixtures/mv_invariant.txt")
    y2 = fx%get2('y')
    comps(1)%c = irregular_t(cov=COV_FULL)
    comps(2)%c = level_t(cov=COV_FULL)
    model = structural_model(y2, comps, info)
    p = [1.0_dp, 0.3_dp, 2.0_dp, 0.5_dp, -0.2_dp, 0.4_dp]
    call check_close(error, model%transform_params(model%untransform_params(p)), p, &
                     "roundtrip")
    if (allocated(error)) return
    call model%update(p)
    call check_close(error, [model%rep%H(:, :, 1)], [1.0_dp, 0.3_dp, 0.3_dp, 2.0_dp], &
                     "H")
    if (allocated(error)) return
    call check_close(error, [model%rep%Q(:, :, 1)], &
                     [0.5_dp, -0.2_dp, -0.2_dp, 0.4_dp], "Q")
    if (allocated(error)) return
    call check(error, info, SS_OK, "info")
  end subroutine test_full_cov

  !> Fixed (deterministic) dummy and trigonometric seasonals span the same
  !> patterns, so their marginal likelihoods (Francke, Koopman and de Vos
  !> 2010) are equal; their diffuse likelihoods differ by a constant.
  subroutine test_fixed_seasonals(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: cd(3), ct(3)
    type(structural_model_t) :: md, mt
    real(dp) :: ld, lt, dd, dt
    integer :: info

    fx = load_fixture("test/fixtures/uc_bsm_dummy.txt")
    allocate (irregular_t :: cd(1)%c, ct(1)%c)
    allocate (trend_t :: cd(2)%c, ct(2)%c)
    cd(3)%c = seasonal_t(period=4, form=SEASONAL_DUMMY, cov=COV_NONE)
    ct(3)%c = seasonal_t(period=4, form=SEASONAL_TRIG, cov=COV_NONE)
    md = structural_model(fx%get2('y'), cd, info)
    mt = structural_model(fx%get2('y'), ct, info)
    dd = md%loglike([0.25_dp, 0.09_dp, 0.0025_dp], info)
    dt = mt%loglike([0.25_dp, 0.09_dp, 0.0025_dp], info)
    md%rep%marginal_likelihood = .true.
    mt%rep%marginal_likelihood = .true.
    ld = md%loglike([0.25_dp, 0.09_dp, 0.0025_dp], info)
    lt = mt%loglike([0.25_dp, 0.09_dp, 0.0025_dp], info)
    call check_close(error, [ld], [lt], "marginal llf")
    if (allocated(error)) return
    call check(error, abs(dd - dt) > 1.0e-6_dp, &
               "diffuse llf differ between representations")
  end subroutine test_fixed_seasonals
  !> DK 3.3.2: y_t = a + A mu*_t + eps_t with one common level, A = (1, a2)'
  !> (a2 a free loading) and a = (0, c2)' (a diffuse intercept on series 2).
  !> Must equal the hand-built representation with state (mu*, c2).
  subroutine test_common_levels(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: comps(3)
    type(structural_model_t) :: model
    type(ssm_rep_t) :: rep
    real(dp), allocatable :: y(:, :)
    real(dp) :: a2, llf
    integer :: info

    fx = load_fixture("test/fixtures/mv_invariant.txt")
    y = fx%get2('y')
    a2 = 0.7_dp
    allocate (irregular_t :: comps(1)%c)
    allocate (level_t :: comps(2)%c)
    comps(3)%c = regression_t(x=reshape(spread(1.0_dp, 1, size(y, 2)), &
                                        [size(y, 2), 1]), series=2, &
                              at_observations=.true.)
    model = structural_model(y, comps, info, &
                             loading=reshape([1.0_dp, 0.0_dp], [2, 1]), &
                             loading_free=reshape([.false., .true.], [2, 1]))
    call check(error, info, SS_OK, "model info")
    if (allocated(error)) return
    call check(error, model%k_params, 4, "irregular (2), level (1), loading (1)")
    if (allocated(error)) return
    llf = model%loglike([1.0_dp, 2.0_dp, 0.25_dp, a2], info)

    rep = ssm_rep(y, 2, 1)
    rep%Z(:, :, 1) = reshape([1.0_dp, a2, 0.0_dp, 1.0_dp], [2, 2])
    rep%T(:, :, 1) = eye(2)
    rep%R(:, 1, 1) = [1.0_dp, 0.0_dp]
    rep%Q = 0.25_dp
    rep%H(:, :, 1) = reshape([1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp], [2, 2])
    call rep%initialize_diffuse()
    call check_close(error, [llf], [loglike(rep, info)], "llf")
  end subroutine test_common_levels

  !> DK 3.3.3: three log signals (exposure, risk, severity), here local
  !> levels, loaded by Z = [1 0 0; 1 1 0; 1 1 1] S_t.
  subroutine test_latent_risk(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: comps(2)
    type(structural_model_t) :: model
    type(ssm_rep_t) :: rep
    real(dp) :: L(3, 3), llf
    integer :: info

    fx = load_fixture("test/fixtures/mv_missing.txt")
    L = reshape([1, 1, 1, 0, 1, 1, 0, 0, 1], [3, 3])
    allocate (irregular_t :: comps(1)%c)
    allocate (level_t :: comps(2)%c)
    model = structural_model(fx%get2('y'), comps, info, loading=L)
    llf = model%loglike([1.0_dp, 1.5_dp, 0.5_dp, 0.1_dp, 0.2_dp, 0.3_dp], info)
    call check_close(error, pack(model%rep%Z(:, :, 1), .true.), pack(L, .true.), &
                     "Z = L")
    if (allocated(error)) return
    rep = ssm_rep(fx%get2('y'), 3, 3)
    rep%Z(:, :, 1) = L
    rep%T(:, :, 1) = eye(3)
    rep%Q(:, :, 1) = reshape([0.1_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.2_dp, 0.0_dp, 0.0_dp, &
                              0.0_dp, 0.3_dp], [3, 3])
    rep%H(:, :, 1) = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.5_dp, 0.0_dp, 0.0_dp, &
                              0.0_dp, 0.5_dp], [3, 3])
    call rep%initialize_diffuse()
    call check_close(error, [llf], [loglike(rep, info)], "llf")
  end subroutine test_latent_risk

  !> DK 3.7: y_t = Lambda f_t + eps_t with one AR(1) factor; the first
  !> loading is fixed at 1 for identification. Must equal the hand-built
  !> model, and its collapsed form (DK 6.5, the p >> m case).
  subroutine test_dynamic_factor(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: comps(2)
    type(structural_model_t) :: model
    type(ssm_rep_t) :: rep, crep
    real(dp), allocatable :: adj(:)
    real(dp) :: llf
    integer :: info

    fx = load_fixture("test/fixtures/mv_missing.txt")
    allocate (irregular_t :: comps(1)%c)
    comps(2)%c = arima_t(ar=1)
    model = structural_model(fx%get2('y'), comps, info, &
                           loading=reshape([1.0_dp, 0.0_dp, 0.0_dp], [3, 1]), &
                           loading_free=reshape([.false., .true., .true.], [3, 1]))
    call check(error, model%k_params, 7, "irregular (3), AR (2), loadings (2)")
    if (allocated(error)) return
    llf = model%loglike([1.0_dp, 1.5_dp, 0.5_dp, 0.6_dp, 0.8_dp, -0.4_dp, 1.3_dp], info)

    rep = ssm_rep(fx%get2('y'), 1, 1)
    rep%Z(:, 1, 1) = [1.0_dp, -0.4_dp, 1.3_dp]
    rep%T = 0.6_dp
    rep%Q = 0.8_dp
    rep%H(:, :, 1) = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.5_dp, 0.0_dp, 0.0_dp, &
                              0.0_dp, 0.5_dp], [3, 3])
    call rep%initialize_stationary()
    call check_close(error, [llf], [loglike(rep, info)], "llf")
    if (allocated(error)) return
    allocate (adj(rep%nobs))
    call collapse_observations(model%rep, crep, adj, info)
    call check_close(error, [loglike(crep, info) + sum(adj)], [llf], "collapsed")
  end subroutine test_dynamic_factor
  !> DK 3.8.1 with unit spacing is the discrete local level.
  subroutine test_continuous_unit(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: cc(2), cd(2)
    type(structural_model_t) :: mc, md
    integer :: info, t

    fx = load_fixture("test/fixtures/nile_llevel_exact.txt")
    allocate (irregular_t :: cc(1)%c, cd(1)%c)
    cc(2)%c = continuous_level_t(times=[(real(t, dp), t=1, 100)])
    allocate (level_t :: cd(2)%c)
    mc = structural_model(fx%get2('y'), cc, info)
    md = structural_model(fx%get2('y'), cd, info)
    call check_close(error, [mc%loglike([15099.0_dp, 1469.1_dp], info)], &
                     [md%loglike([15099.0_dp, 1469.1_dp], info)], "llf")
  end subroutine test_continuous_unit

  !> DK (3.36): estimating alpha between observation times by inserting a
  !> missing observation leaves the likelihood and the estimates at the
  !> observation times unchanged.
  subroutine test_continuous_interp(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: ca(2), cb(2)
    type(structural_model_t) :: ma, mb
    type(filter_result_t) :: fa, fb
    type(smoother_result_t) :: sa, sb
    real(dp), allocatable :: y(:, :), yb(:, :), ta(:), tb(:)
    integer :: info, t

    fx = load_fixture("test/fixtures/nile_llevel_exact.txt")
    y = fx%get2('y')
    y = y(:, 1:30)
    ta = [(real(t, dp) + 0.3_dp * sin(real(t, dp)), t=1, 30)]
    ! Insert a time halfway between observations 10 and 11.
    tb = [ta(1:10), 0.5_dp * (ta(10) + ta(11)), ta(11:)]
    allocate (yb(1, 31))
    yb(1, :) = [y(1, 1:10), ieee_nan(), y(1, 11:)]
    allocate (irregular_t :: ca(1)%c, cb(1)%c)
    ca(2)%c = continuous_level_t(times=ta)
    cb(2)%c = continuous_level_t(times=tb)
    ma = structural_model(y, ca, info)
    mb = structural_model(yb, cb, info)
    call ma%smooth([15099.0_dp, 1469.1_dp], fa, sa, info)
    call mb%smooth([15099.0_dp, 1469.1_dp], fb, sb, info)
    call check_close(error, [fb%llf], [fa%llf], "llf")
    if (allocated(error)) return
    call check_close(error, [sb%alphahat(1, 1:10), sb%alphahat(1, 12:)], &
                     sa%alphahat(1, :), "alphahat")

  contains

    real(dp) function ieee_nan()
      real(dp) :: z
      z = 0.0_dp
      ieee_nan = z / z
    end function ieee_nan
  end subroutine test_continuous_interp

  !> DK 3.9.2: continuous smooth trend + irregular, lambda = s_eps^2 / s^2,
  !> gives the cubic smoothing spline (Wahba 1978); against scipy.
  subroutine test_cubic_spline(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: comps(2)
    type(structural_model_t) :: model
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp) :: lam
    integer :: info

    fx = load_fixture("test/fixtures/smoothing_spline.txt")
    lam = sum(fx%get1('lam'))
    allocate (irregular_t :: comps(1)%c)
    comps(2)%c = continuous_trend_t(times=fx%get1('x'))
    model = structural_model(fx%get2('y'), comps, info)
    call model%smooth([1.0_dp, 1.0_dp / lam], fres, sres, info)
    call check(error, info, SS_OK, "smooth info")
    if (allocated(error)) return
    call check_close(error, sres%alphahat(1, :), fx%get1('fitted'), "spline")
  end subroutine test_cubic_spline

  !> DK 3.9.1: y = alpha + eps, Delta^2 alpha = zeta, Var(zeta) = sigma^2 /
  !> lambda: alphahat minimizes sum (y - alpha)^2 + lambda sum (Delta^2 alpha)^2,
  !> i.e. solves (I + lambda D'D) alpha = y (diffuse initial level and slope).
  subroutine test_discrete_spline(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: comps(2)
    type(structural_model_t) :: model
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: y(:, :), A(:, :), D(:, :), b(:, :)
    real(dp) :: lam
    integer :: info, n, t

    fx = load_fixture("test/fixtures/nile_llevel_exact.txt")
    y = fx%get2('y')
    y = y(:, 1:40)
    n = 40
    lam = 25.0_dp
    allocate (irregular_t :: comps(1)%c)
    comps(2)%c = trend_t(cov_level=COV_NONE)
    model = structural_model(y, comps, info)
    call model%smooth([1.0_dp, 1.0_dp / lam], fres, sres, info)
    allocate (D(n - 2, n), source=0.0_dp)
    do t = 1, n - 2
      D(t, t:t + 2) = [1.0_dp, -2.0_dp, 1.0_dp]
    end do
    A = eye(n) + lam * matmul(transpose(D), D)
    b = transpose(y)
    call solve(A, b, info)
    call check_close(error, sres%alphahat(1, :), b(:, 1), "penalized least squares")
  end subroutine test_discrete_spline
end module test_structural
