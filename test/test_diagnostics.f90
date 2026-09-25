!> Goodness of fit and diagnostics (DK 7.4-7.5) against statsmodels.
module test_diagnostics
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use statespace
  use statespace_special, only: chi2_sf, f_cdf
  use statespace_linalg, only: chol_inv, solve
  use statespace_diagnostics
  use fixture_io, only: fixture_t, load_fixture, rep_from_fixture
  implicit none
  private

  public :: collect_diagnostics

contains

  subroutine collect_diagnostics(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("special_functions", test_special), &
                new_unittest("nile_residual_tests", test_nile), &
                new_unittest("multivariate_standardized_residuals", test_mv), &
                new_unittest("auxiliary_residuals", test_auxiliary), &
                new_unittest("r2_diffuse", test_r2), &
                new_unittest("de_jong_penzer", test_djp), &
                new_unittest("vector_auxiliary_residuals", test_vector_aux), &
                new_unittest("likelihood_fixed_unknown_initial", test_llf_fixed), &
                new_unittest("marginal_ar1_constant", test_marginal_ar1), &
                new_unittest("marginal_common_trend_invariance", test_marginal_common_trend), &
                new_unittest("marginal_augmented_matches_exact", test_marginal_augmented), &
                new_unittest("recursive_and_least_squares_residuals", test_ls_residuals) &
                ]
  end subroutine collect_diagnostics

  subroutine check_rel(error, actual, expected, tol, label)
    type(error_type), allocatable, intent(out) :: error
    real(dp), intent(in) :: actual(:), expected(:), tol
    character(len=*), intent(in) :: label
    real(dp) :: err
    character(len=32) :: buf

    err = maxval(abs(actual - expected) / max(1.0_dp, abs(expected)))
    write (buf, '(es10.3)') err
    call check(error, err <= tol, label//": relative error "//trim(buf))
  end subroutine check_rel

  subroutine test_special(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    real(dp), allocatable :: x(:)

    fx = load_fixture("test/fixtures/diagnostics.txt")
    x = fx%get1('chi2_sf_in')
    call check_rel(error, [chi2_sf(x(1), x(2)), chi2_sf(x(3), x(4))], fx%get1('chi2_sf_out'), &
                   1.0e-12_dp, "chi2 survival")
    if (allocated(error)) return
    x = fx%get1('f_cdf_in')
    call check_rel(error, [f_cdf(x(1), x(2), x(3))], fx%get1('f_cdf_out'), 1.0e-12_dp, "F cdf")
  end subroutine test_special

  !> Nile local level at DK's parameters, exact diffuse: statsmodels' tests on
  !> the standardized one-step errors after the diffuse period.
  subroutine test_nile(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx, dx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    real(dp), allocatable :: e(:, :), expected(:, :)
    real(dp) :: lb(10), lbp(10), jb, jbp, skew, kurt, hstat, hp
    integer :: info, d

    fx = load_fixture("test/fixtures/nile_llevel_exact.txt")
    dx = load_fixture("test/fixtures/diagnostics.txt")
    rep = rep_from_fixture(fx)
    call kalman_filter(rep, fres, info)
    allocate (e(1, rep%nobs))
    call standardized_residuals(fres, e)
    d = diagnostic_start(rep, fres)
    call check(error, d, 2, "diagnostics start after the diffuse period")
    if (allocated(error)) return
    call check(error, ieee_is_nan(e(1, 1)), "NaN in the diffuse period")
    if (allocated(error)) return
    expected = dx%get2('std_err')
    call check_rel(error, e(1, d:), expected(1, d:), 1.0e-9_dp, "standardized residuals")
    if (allocated(error)) return

    call ljung_box(e(1, d:), lb, lbp)
    call check_rel(error, lb, dx%get1('lb_stat'), 1.0e-9_dp, "Ljung-Box statistics")
    if (allocated(error)) return
    call check_rel(error, lbp, dx%get1('lb_pvalue'), 1.0e-9_dp, "Ljung-Box p-values")
    if (allocated(error)) return
    call jarque_bera(e(1, d:), jb, jbp, skew, kurt)
    call check_rel(error, [jb, jbp, skew, kurt], dx%get1('jb'), 1.0e-9_dp, "Jarque-Bera")
    if (allocated(error)) return
    call breakvar_test(e(1, d:), hstat, hp)
    call check_rel(error, [hstat, hp], dx%get1('het'), 1.0e-9_dp, "heteroskedasticity H(h)")
  end subroutine test_nile

  !> p = 3, correlated H, partly missing: e_t = L^-1 v_t over the observed
  !> block (statsmodels writes 0 for missing elements; we write NaN).
  subroutine test_mv(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx, dx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    real(dp), allocatable :: e(:, :), expected(:, :)
    logical, allocatable :: obs(:, :)
    integer :: info

    fx = load_fixture("test/fixtures/mv_missing.txt")
    dx = load_fixture("test/fixtures/diagnostics.txt")
    rep = rep_from_fixture(fx)
    call kalman_filter(rep, fres, info)
    allocate (e(rep%k_endog, rep%nobs))
    call standardized_residuals(fres, e)
    obs = .not. ieee_is_nan(rep%y)
    call check(error, all(ieee_is_nan(e) .eqv. .not. obs), "NaN exactly at missing elements")
    if (allocated(error)) return
    expected = dx%get2('mv_std_err')
    call check_rel(error, pack(e, obs), pack(expected, obs), 1.0e-9_dp, "standardized residuals")
  end subroutine test_mv

  !> Auxiliary residuals equal DK's u_t / sqrt(D_t) and
  !> Q R' r_t / sqrt(Q R' N_t R Q) (univariate, conventional periods).
  subroutine test_auxiliary(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: eps_std(:, :), eta_std(:, :), u(:), Dt(:), a(:), b(:)
    integer :: info, t, n

    fx = load_fixture("test/fixtures/nile_llevel_known.txt")
    rep = rep_from_fixture(fx)
    n = rep%nobs
    call kalman_filter(rep, fres, info)
    call state_smoother(rep, fres, sres, info)
    allocate (eps_std(1, n), eta_std(1, n), u(n), Dt(n), a(n - 1), b(n - 1))
    call auxiliary_residuals(rep, sres, eps_std, eta_std)
    do t = 1, n
      u(t) = fres%Finv(1, 1, t) * fres%v(1, t) - fres%K(1, 1, t) * sres%r(1, t)
      Dt(t) = fres%Finv(1, 1, t) + fres%K(1, 1, t)**2 * sres%N(1, 1, t)
    end do
    call check_rel(error, eps_std(1, :), u / sqrt(Dt), 1.0e-9_dp, "measurement")
    if (allocated(error)) return
    ! eta_n has r_n = N_n = 0: no information, NaN.
    a = eta_std(1, :n - 1)
    b = sres%r(1, 1:n - 1) / sqrt(sres%N(1, 1, 1:n - 1))
    call check_rel(error, a, b, 1.0e-9_dp, "state")
    if (allocated(error)) return
    call check(error, ieee_is_nan(eta_std(1, n)), "last state residual undefined")
  end subroutine test_auxiliary

  !> R2_D = 1 - SSE / sum (Delta y - mean Delta y)^2 for the Nile local level.
  subroutine test_r2(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    real(dp), allocatable :: dy(:)
    real(dp) :: sse, sst
    integer :: info, n

    fx = load_fixture("test/fixtures/nile_llevel_exact.txt")
    rep = rep_from_fixture(fx)
    n = rep%nobs
    call kalman_filter(rep, fres, info)
    sse = sum(fres%v(1, 2:)**2)
    dy = rep%y(1, 2:) - rep%y(1, :n - 1)
    sst = sum((dy - sum(dy) / (n - 1))**2)
    call check_rel(error, [r2_diffuse(rep, fres)], [1.0_dp - sse / sst], 1.0e-12_dp, "R2_D")
  end subroutine test_r2
  !> For the local level model the de Jong-Penzer statistics are the
  !> auxiliary residuals: u/sqrt(D) and r/sqrt(N) (DK 7.5).
  subroutine test_djp(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: eps_std(:, :), eta_std(:, :), r_stat(:, :), e_stat(:, :)
    integer :: info, n

    fx = load_fixture("test/fixtures/nile_llevel_known.txt")
    rep = rep_from_fixture(fx)
    n = rep%nobs
    call kalman_filter(rep, fres, info)
    call state_smoother(rep, fres, sres, info)
    allocate (eps_std(1, n), eta_std(1, n), r_stat(1, n), e_stat(1, n))
    call auxiliary_residuals(rep, sres, eps_std, eta_std)
    call de_jong_penzer(rep, fres, sres, r_stat, e_stat)
    call check_rel(error, e_stat(1, :), eps_std(1, :), 1.0e-9_dp, "e / sqrt(D)")
    if (allocated(error)) return
    call check_rel(error, r_stat(1, :n - 1), eta_std(1, :n - 1), 1.0e-9_dp, "r / sqrt(N)")
  end subroutine test_djp

  !> Vector-standardized auxiliary residuals: |B epshat|^2 = epshat' Var^-1 epshat,
  !> and for p = 1 they equal the element-wise ones.
  subroutine test_vector_aux(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: ev(:, :), hv(:, :), es(:, :), hs(:, :), Vinv(:, :)
    real(dp) :: logdet
    integer :: info, t

    fx = load_fixture("test/fixtures/mv_invariant.txt")
    rep = rep_from_fixture(fx)
    call kalman_filter(rep, fres, info)
    call state_smoother(rep, fres, sres, info)
    allocate (ev(rep%k_endog, rep%nobs), hv(rep%k_posdef, rep%nobs))
    call auxiliary_residuals_vector(rep, sres, ev, hv)
    do t = 1, rep%nobs, 7
      Vinv = rep%H(:, :, 1) - sres%epsvar(:, :, t)
      call chol_inv(Vinv, logdet, info)
      call check_rel(error, [sum(ev(:, t)**2)], &
                     [dot_product(sres%epshat(:, t), matmul(Vinv, sres%epshat(:, t)))], 1.0e-9_dp, &
                     "quadratic form")
      if (allocated(error)) return
    end do

    fx = load_fixture("test/fixtures/nile_llevel_known.txt")
    rep = rep_from_fixture(fx)
    call kalman_filter(rep, fres, info)
    call state_smoother(rep, fres, sres, info)
    deallocate (ev, hv)
    allocate (ev(1, rep%nobs), hv(1, rep%nobs), es(1, rep%nobs), hs(1, rep%nobs))
    call auxiliary_residuals_vector(rep, sres, ev, hv)
    call auxiliary_residuals(rep, sres, es, hs)
    call check_rel(error, ev(1, :), es(1, :), 1.0e-12_dp, "p = 1")
  end subroutine test_vector_aux

  !> DK 7.2.4: with delta fixed but unknown, the concentrated log likelihood
  !> equals the ordinary one with delta = delta_hat (known initialization
  !> a + A delta_hat, P_star).
  subroutine test_llf_fixed(error)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), parameter :: paths(2) = [character(len=44) :: &
                                               "test/fixtures/nile_lltrend_exact.txt", &
                                               "test/fixtures/uc_level_ar1_mixed.txt"]
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep, fixed
    type(augmented_result_t) :: ares
    real(dp), allocatable :: a1(:), Pstar(:, :), Pinf(:, :)
    integer :: info, k, m

    do k = 1, size(paths)
      fx = load_fixture(trim(paths(k)))
      rep = rep_from_fixture(fx)
      m = rep%k_states
      call augmented_filter(rep, ares, info)
      call check(error, info, SS_OK, "augmented info")
      if (allocated(error)) return
      allocate (a1(m), Pstar(m, m), Pinf(m, m))
      call rep%initial_state(a1, Pstar, Pinf, info)
      fixed = rep
      call fixed%initialize_general(a1 + matmul(ares%A(:, :, 1), ares%delta), Pstar, 0.0_dp * Pinf)
      call check_rel(error, [ares%llf_fixed], [loglike(fixed, info)], 1.0e-9_dp, trim(paths(k)))
      if (allocated(error)) return
      deallocate (a1, Pstar, Pinf)
    end do
  end subroutine test_llf_fixed
  !> FKV (2010) 4.1: y_t = mu + u_t, u_t+1 = rho u_t + eta_t, mu diffuse. The
  !> marginal likelihood is the density of an orthonormal transform of y
  !> orthogonal to 1, i.e. of Delta y with |D D'| = n, so it equals the
  !> likelihood of FKV's state space form for Delta y plus log(n) / 2.
  subroutine test_marginal_ar1(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: level, diff
    real(dp), allocatable :: y(:, :)
    real(dp) :: rho, s2, lm, ld
    integer :: info, k, n

    fx = load_fixture("test/fixtures/nile_llevel_known.txt")
    y = fx%get2('y')
    n = size(y, 2)
    s2 = 15000.0_dp
    do k = 1, 3
      rho = 0.3_dp * k
      level = ssm_rep(y, 2, 1)
      level%Z(1, :, 1) = [1.0_dp, 1.0_dp]
      level%T(:, :, 1) = reshape([1.0_dp, 0.0_dp, 0.0_dp, rho], [2, 2])
      level%R(:, 1, 1) = [0.0_dp, 1.0_dp]
      level%Q = s2
      call level%initialize_block(1, 1, INIT_DIFFUSE)
      call level%initialize_block(2, 2, INIT_STATIONARY)
      level%marginal_likelihood = .true.
      lm = loglike(level, info)
      call check(error, info, SS_OK, "level model")
      if (allocated(error)) return

      ! alpha_t = (Delta u_t, eta_t-1)', Delta y_t = Delta u_t (FKV 4.1)
      diff = ssm_rep(y(:, 2:) - y(:, :n - 1), 2, 1)
      diff%Z(1, :, 1) = [1.0_dp, 0.0_dp]
      diff%T(:, :, 1) = reshape([rho, 0.0_dp, -1.0_dp, 0.0_dp], [2, 2])
      diff%R(:, 1, 1) = [1.0_dp, 1.0_dp]
      diff%Q = s2
      call diff%initialize_known([0.0_dp, 0.0_dp], &
                                 s2 * reshape([2.0_dp / (1.0_dp + rho), 1.0_dp, 1.0_dp, 1.0_dp], [2, 2]))
      ld = loglike(diff, info)
      call check_rel(error, [lm], [ld + 0.5_dp * log(real(n, dp))], 1.0e-10_dp, "marginal = Delta y")
      if (allocated(error)) return
    end do
  end subroutine test_marginal_ar1

  !> FKV (2010) 4.2: y_t = c + K mu_t + eps_t with a common random walk mu_t,
  !> K = (kappa, 1)', c = (0, c_2)'. Representation A has state (mu, c_2), so
  !> Z depends on kappa; representation B has state c + K mu with Z = I,
  !> R = K. The marginal likelihoods differ by a constant across kappa; the
  !> diffuse ones do not.
  subroutine test_marginal_common_trend(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: ra, rb
    real(dp) :: kappa, dm(2), dd(2)
    integer :: info, k
    character(len=64) :: buf

    fx = load_fixture("test/fixtures/mv_invariant.txt")
    do k = 1, 2
      kappa = 0.5_dp * k
      ra = ssm_rep(fx%get2('y'), 2, 1)
      ra%Z(:, :, 1) = reshape([kappa, 1.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      ra%T(:, :, 1) = eye(2)
      ra%R(:, 1, 1) = [1.0_dp, 0.0_dp]
      ra%Q = 0.25_dp
      ra%H(:, :, 1) = eye(2)
      call ra%initialize_diffuse()
      rb = ssm_rep(fx%get2('y'), 2, 1)
      rb%Z(:, :, 1) = eye(2)
      rb%T(:, :, 1) = eye(2)
      rb%R(:, 1, 1) = [kappa, 1.0_dp]
      rb%Q = 0.25_dp
      rb%H(:, :, 1) = eye(2)
      call rb%initialize_diffuse()
      dd(k) = loglike(ra, info) - loglike(rb, info)
      ra%marginal_likelihood = .true.
      rb%marginal_likelihood = .true.
      dm(k) = loglike(ra, info) - loglike(rb, info)
    end do
    write (buf, '(4es12.4)') dm, dd
    call check(error, abs(dm(1) - dm(2)) < 1.0e-9_dp .and. abs(dd(1) - dd(2)) > 1.0e-3_dp, &
               "marginal (A - B) and diffuse (A - B) at two kappas: "//trim(buf))
  end subroutine test_marginal_common_trend

  subroutine test_marginal_augmented(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(augmented_result_t) :: ares
    integer :: info

    fx = load_fixture("test/fixtures/mv_diffuse.txt")
    rep = rep_from_fixture(fx)
    call augmented_filter(rep, ares, info)
    rep%marginal_likelihood = .true.
    call check_rel(error, [ares%llf_marginal], [loglike(rep, info)], 1.0e-9_dp, "llf_marginal")
  end subroutine test_marginal_augmented
  !> DK 6.2.4 for a pure regression y_t = x_t' beta + eps_t: the recursive
  !> residuals are y_t - x_t' beta_hat_t-1 (OLS on y_1..y_t-1) and the least
  !> squares residuals y_t - x_t' beta_hat (OLS on all data).
  subroutine test_ls_residuals(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: n = 30, k = 2
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    real(dp) :: y(1, n), x(n, k), vplus(1, n), beta(k, 1), XtX(k, k), rec(n)
    integer :: info, t

    do t = 1, n
      x(t, :) = [1.0_dp, sin(0.7_dp * t)]
      y(1, t) = 2.0_dp + 0.5_dp * x(t, 2) + 0.3_dp * cos(2.1_dp * t)
    end do
    rep = ssm_rep(y, k, 1)
    rep%Z = reshape(transpose(x), [1, k, n])
    rep%T(:, :, 1) = eye(k)
    rep%R = 0.0_dp
    rep%H = 1.0_dp
    call rep%initialize_diffuse()
    call kalman_filter(rep, fres, info)
    do t = k + 1, n
      XtX = matmul(transpose(x(1:t - 1, :)), x(1:t - 1, :))
      beta = matmul(transpose(x(1:t - 1, :)), reshape(y(1, 1:t - 1), [t - 1, 1]))
      call solve(XtX, beta, info)
      rec(t) = y(1, t) - dot_product(x(t, :), beta(:, 1))
    end do
    call check_rel(error, fres%v(1, k + 1:), rec(k + 1:), 1.0e-9_dp, "recursive residuals")
    if (allocated(error)) return

    call least_squares_residuals(rep, 1, k, vplus, info)
    call check(error, info, SS_OK, "info")
    if (allocated(error)) return
    XtX = matmul(transpose(x), x)
    beta = matmul(transpose(x), transpose(y))
    call solve(XtX, beta, info)
    call check_rel(error, vplus(1, :), y(1, :) - matmul(x, beta(:, 1)), 1.0e-9_dp, &
                   "least squares residuals")
  end subroutine test_ls_residuals
end module test_diagnostics
