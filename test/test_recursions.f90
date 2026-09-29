!> DK chapter 4 recursions: matrix form (4.13), updating and fixed-point /
!> fixed-lag smoothing (4.4.5-4.4.6), filtering weights (4.8.2), the Whittle
!> relation (4.6.3), de Jong-Shephard simulation (4.9.3) and the steady
!> state (2.11, 4.3.4).
module test_recursions
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use statespace
  use fixture_io, only: fixture_t, load_fixture, rep_from_fixture
  implicit none
  private

  public :: collect_recursions

contains

  subroutine collect_recursions(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ new_unittest("dense_matches_recursions", test_dense), &
                 new_unittest("dense_rejects_diffuse", test_dense_diffuse), &
                 new_unittest("update_smoothed", test_update), &
                 new_unittest("fixed_point_smoother", test_fixed_point), &
                 new_unittest("fixed_lag_smoother", test_fixed_lag), &
                 new_unittest("filtering_weights", test_filter_weights), &
                 new_unittest("smoothing_weights_table_4_6", &
                              test_smoothing_weights_dk), &
                 new_unittest("whittle_recursion", test_whittle), &
                 new_unittest("de_jong_shephard", test_djs), &
                 new_unittest("steady_state_local_level", test_steady_nile), &
                 new_unittest("steady_state_long_filter", test_steady_filter), &
                 new_unittest("steady_state_arma_no_noise", test_steady_arma), &
                 new_unittest("steady_state_deterministic_level", test_steady_fixed) ]
  end subroutine collect_recursions

  subroutine check_rel(error, actual, expected, tol, label)
    type(error_type), allocatable, intent(out) :: error
    real(dp), intent(in) :: actual(:), expected(:), tol
    character(len=*), intent(in) :: label
    real(dp) :: err
    character(len=32) :: buf

    err = maxval(abs(actual - expected)) / max(1.0_dp, maxval(abs(expected)))
    write (buf, '(es10.3)') err
    call check(error, err <= tol, label//": relative error "//trim(buf))
  end subroutine check_rel

  !> The first k observations of a time-invariant fixture model.
  function truncated(rep, k) result(r)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: k
    type(ssm_rep_t) :: r

    r = rep
    r%y = rep%y(:, 1:k)
    r%nobs = k
  end function truncated

  subroutine smooth(rep, fres, sres)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(out) :: fres
    type(smoother_result_t), intent(out) :: sres
    integer :: info

    call kalman_filter(rep, fres, info)
    call state_smoother(rep, fres, sres, info)
  end subroutine smooth

  subroutine test_dense(error)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), parameter :: paths(4) = [character(len=40) :: &
                                               "test/fixtures/nile_llevel_known.txt", &
                                               "test/fixtures/mv_missing.txt", &
                                               "test/fixtures/mv_timevarying.txt", &
                                               "test/fixtures/mv_stationary.txt"]
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: a(:, :), V(:, :, :)
    real(dp) :: llf
    integer :: info, k

    do k = 1, size(paths)
      fx = load_fixture(trim(paths(k)))
      rep = rep_from_fixture(fx)
      call smooth(rep, fres, sres)
      allocate (a(rep%k_states, rep%nobs), V(rep%k_states, rep%k_states, rep%nobs))
      call dense_loglike_smooth(rep, llf, a, V, info)
      call check(error, info, SS_OK, "dense info")
      if (allocated(error)) return
      call check_rel(error, [llf], [fres%llf], 1.0e-9_dp, trim(paths(k))//" llf")
      if (allocated(error)) return
      call check_rel(error, pack(a, .true.), pack(sres%alphahat, .true.), 1.0e-8_dp, &
                     trim(paths(k))//" alphahat")
      if (allocated(error)) return
      call check_rel(error, pack(V, .true.), pack(sres%V, .true.), 1.0e-8_dp, &
                     trim(paths(k))//" V")
      if (allocated(error)) return
      deallocate (a, V)
    end do
  end subroutine test_dense

  subroutine test_dense_diffuse(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    real(dp) :: llf, a(1, 100), V(1, 1, 100)
    integer :: info

    fx = load_fixture("test/fixtures/nile_llevel_exact.txt")
    rep = rep_from_fixture(fx)
    call dense_loglike_smooth(rep, llf, a, V, info)
    call check(error, info, SS_ERR_UNSUPPORTED, "exact diffuse")
  end subroutine test_dense_diffuse

  !> Smooth the first 25 observations, update to all 40 (DK 4.4.5), and
  !> compare with smoothing all 40; for both filter methods.
  subroutine test_update(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep, part
    type(filter_result_t) :: fres, fpart
    type(smoother_result_t) :: sres, spart
    real(dp), allocatable :: a(:, :), V(:, :, :)
    integer :: info, k, n0

    fx = load_fixture("test/fixtures/mv_missing.txt")
    n0 = 25
    do k = 1, 2
      rep = rep_from_fixture(fx)
      if (k == 2) rep%filter_method = FILTER_UNIVARIATE
      part = truncated(rep, n0)
      call smooth(rep, fres, sres)
      call smooth(part, fpart, spart)
      allocate (a(rep%k_states, rep%nobs), V(rep%k_states, rep%k_states, rep%nobs), &
                source=0.0_dp)
      a(:, 1:n0) = spart%alphahat
      V(:, :, 1:n0) = spart%V
      call update_smoothed(rep, fres, n0, a, V, info)
      call check(error, info, SS_OK, "update info")
      if (allocated(error)) return
      call check_rel(error, pack(a, .true.), pack(sres%alphahat, .true.), 1.0e-9_dp, &
                     "alphahat")
      if (allocated(error)) return
      call check_rel(error, pack(V, .true.), pack(sres%V, .true.), 1.0e-9_dp, "V")
      if (allocated(error)) return
      deallocate (a, V)
    end do
  end subroutine test_update

  !> alphahat_10|k for k = 10..40 equals the smoother on y_1..y_k at t = 10.
  subroutine test_fixed_point(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: t0 = 10
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres, fk
    type(smoother_result_t) :: sk
    real(dp), allocatable :: pa(:, :), pV(:, :, :)
    integer :: info, k

    fx = load_fixture("test/fixtures/mv_missing.txt")
    rep = rep_from_fixture(fx)
    call kalman_filter(rep, fres, info)
    allocate (pa(rep%k_states, rep%nobs - t0 + 1), &
              pV(rep%k_states, rep%k_states, rep%nobs - t0 + 1))
    call fixed_point_smoother(rep, fres, t0, pa, pV, info)
    call check(error, info, SS_OK, "info")
    if (allocated(error)) return
    do k = t0, rep%nobs, 7
      call smooth(truncated(rep, k), fk, sk)
      call check_rel(error, pa(:, k - t0 + 1), sk%alphahat(:, t0), 1.0e-9_dp, &
                     "alphahat")
      if (allocated(error)) return
      call check_rel(error, pack(pV(:, :, k - t0 + 1), .true.), &
                     pack(sk%V(:, :, t0), .true.), 1.0e-9_dp, "V")
      if (allocated(error)) return
    end do
  end subroutine test_fixed_point

  !> alphahat_s|s+3 equals the smoother on y_1..y_s+3 at s.
  subroutine test_fixed_lag(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: j = 3
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres, fk
    type(smoother_result_t) :: sk
    real(dp), allocatable :: a(:, :), V(:, :, :)
    integer :: info, s

    fx = load_fixture("test/fixtures/mv_missing.txt")
    rep = rep_from_fixture(fx)
    call kalman_filter(rep, fres, info)
    allocate (a(rep%k_states, rep%nobs), V(rep%k_states, rep%k_states, rep%nobs))
    call fixed_lag_smoother(rep, fres, j, a, V, info)
    call check(error, info, SS_OK, "info")
    if (allocated(error)) return
    do s = 1, rep%nobs - j, 5
      call smooth(truncated(rep, s + j), fk, sk)
      call check_rel(error, a(:, s), sk%alphahat(:, s), 1.0e-9_dp, "alphahat")
      if (allocated(error)) return
      call check_rel(error, pack(V(:, :, s), .true.), pack(sk%V(:, :, s), .true.), &
                     1.0e-9_dp, "V")
      if (allocated(error)) return
    end do
    call check(error, all(ieee_is_nan(a(:, rep%nobs - j + 1:))), &
               "undefined columns are NaN")
  end subroutine test_fixed_lag

  !> DK 4.8.2 weights against filtering unit inputs (the filter is linear in
  !> y when c = d = a_1 = 0).
  subroutine test_filter_weights(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep, base, unit
    type(filter_result_t) :: fres, fu
    real(dp), allocatable :: Wa(:, :, :, :), Watt(:, :, :, :)
    integer :: info, i, j, m, p, n

    fx = load_fixture("test/fixtures/mv_missing.txt")
    rep = rep_from_fixture(fx)
    m = rep%k_states; p = rep%k_endog; n = rep%nobs
    call kalman_filter(rep, fres, info)
    allocate (Wa(m, p, n, n), Watt(m, p, n, n))
    call filtered_state_weights(rep, fres, Wa, Watt, info)
    call check(error, info, SS_OK, "info")
    if (allocated(error)) return
    base = rep
    base%y = merge(rep%y, 0.0_dp, ieee_is_nan(rep%y))
    base%c = 0.0_dp
    base%d = 0.0_dp
    call base%initialize_known(0.0_dp * fres%a(:, 1), fres%P(:, :, 1))
    do j = 1, n, 3
      do i = 1, p
        if (ieee_is_nan(rep%y(i, j))) cycle
        unit = base
        unit%y(i, j) = 1.0_dp
        call kalman_filter(unit, fu, info)
        call check_rel(error, pack(Wa(:, i, :, j), .true.), &
                       pack(fu%a(:, 1:n), .true.), 1.0e-10_dp, "a_t weights")
        if (allocated(error)) return
        call check_rel(error, pack(Watt(:, i, :, j), .true.), pack(fu%att, .true.), &
                       1.0e-10_dp, "a_t|t weights")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_filter_weights

  !> Our smoothing weights (by linearity) against DK Table 4.6:
  !> weight of y_j in alphahat_t = (I - P_t N_t-1) L_t-1 ... L_j+1 K_j, j < t.
  subroutine test_smoothing_weights_dk(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: W(:, :, :, :), C(:, :, :, :), A(:, :, :), prod(:, :), &
                             L(:, :)
    real(dp), allocatable :: K(:, :), Finv(:, :)
    integer :: info, t, j, s, m, p, n

    fx = load_fixture("test/fixtures/mv_missing.txt")
    rep = rep_from_fixture(fx)
    m = rep%k_states; p = rep%k_endog; n = rep%nobs
    call smooth(rep, fres, sres)
    allocate (W(m, p, n, n), C(m, m, n, n), A(m, m, n), K(m, p), Finv(p, p), L(m, m))
    call smoothed_state_weights(rep, W, C, A, info)
    do t = 5, n, 11
      do j = 1, t - 1, 2
        call conventional_gain(rep, fres, j, K, Finv, info)
        prod = K
        do s = j + 1, t - 1
          call innovation_transition(rep, fres, s, L, info)
          prod = matmul(L, prod)
        end do
        prod = matmul(eye(m) - matmul(fres%P(:, :, t), sres%N(:, :, t - 1)), prod)
        call check_rel(error, pack(W(:, :, t, j), .true.), pack(prod, .true.), &
                       1.0e-9_dp, "Table 4.6 weights")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_smoothing_weights_dk

  !> The Whittle recursion on a short Nile series (it is unstable over long
  !> ones, as DK note), for known and exact diffuse initialization.
  subroutine test_whittle(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: a(:, :)
    integer :: info, k

    fx = load_fixture("test/fixtures/nile_lltrend_exact.txt")
    do k = 1, 2
      rep = truncated(rep_from_fixture(fx), 25)
      if (k == 2) call rep%initialize_known([1000.0_dp, 0.0_dp], 1.0e4_dp * eye(2))
      call smooth(rep, fres, sres)
      allocate (a(rep%k_states, rep%nobs))
      call whittle_smoother(rep, fres, a, info)
      call check(error, info, SS_OK, "info")
      if (allocated(error)) return
      call check_rel(error, pack(a, .true.), pack(sres%alphahat, .true.), 1.0e-6_dp, &
                     "alphahat")
      if (allocated(error)) return
      deallocate (a)
    end do
  end subroutine test_whittle

  !> de Jong-Shephard draws: mean and variance of eps, eta and the state
  !> against the smoothed moments (5 standard errors), including missing
  !> elements with correlated H.
  subroutine test_djs(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: ndraw = 3000
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: u(:, :), u0(:), eps(:, :), eta(:, :), alpha(:, :)
    real(dp), allocatable :: se(:, :), se2(:, :), sh(:, :), sh2(:, :), sa(:, :), &
                             sa2(:, :)
    integer, allocatable :: seed(:)
    integer :: info, k, nseed, p, m, r, n
    real(dp) :: z
    character(len=32) :: buf

    call random_seed(size=nseed)
    seed = [(4242 + 3 * k, k=1, nseed)]
    call random_seed(put=seed)
    fx = load_fixture("test/fixtures/mv_missing.txt")
    rep = rep_from_fixture(fx)
    p = rep%k_endog; m = rep%k_states; r = rep%k_posdef; n = rep%nobs
    call smooth(rep, fres, sres)
    allocate (eps(p, n), eta(r, n), alpha(m, n))
    allocate (se(p, n), se2(p, n), sh(r, n), sh2(r, n), sa(m, n), sa2(m, n), &
              source=0.0_dp)
    do k = 1, ndraw
      allocate (u(p, n))
      do info = 1, n
        call draw_standard_normal(u(:, info))
      end do
      call djs_measurement_disturbances(rep, fres, u, eps, info)
      deallocate (u)
      allocate (u(r, n), u0(m))
      do info = 1, n
        call draw_standard_normal(u(:, info))
      end do
      call draw_standard_normal(u0)
      call djs_state_disturbances(rep, fres, u, u0, eta, alpha, info)
      deallocate (u0)
      deallocate (u)
      se = se + eps; se2 = se2 + eps**2
      sh = sh + eta; sh2 = sh2 + eta**2
      sa = sa + alpha; sa2 = sa2 + alpha**2
    end do
    se = se / ndraw; se2 = se2 / ndraw - se**2
    sh = sh / ndraw; sh2 = sh2 / ndraw - sh**2
    sa = sa / ndraw; sa2 = sa2 / ndraw - sa**2
    z = max(zmax(se, se2, sres%epshat, diag3(sres%epsvar)), &
            zmax(sh, sh2, sres%etahat, diag3(sres%etavar)), &
            zmax(sa, sa2, sres%alphahat, diag3(sres%V)))
    write (buf, '(f8.2)') z
    call check(error, z < 5.0_dp, "largest standardized error "//trim(buf))

  contains

    function diag3(A) result(D)
      real(dp), intent(in) :: A(:, :, :)
      real(dp) :: D(size(A, 1), size(A, 3))
      integer :: i

      do i = 1, size(A, 1)
        D(i, :) = A(i, i, :)
      end do
    end function diag3

    !> Largest |z| for the sample means and variances of the draws.
    real(dp) function zmax(mean, var, true_mean, true_var)
      real(dp), intent(in) :: mean(:, :), var(:, :), true_mean(:, :), true_var(:, :)

      zmax = max(maxval(abs(mean - true_mean) / sqrt(true_var / ndraw), &
                        mask=true_var > 1.0e-12_dp), maxval(abs(var - true_var) &
          / (sqrt(2.0_dp / ndraw) * true_var), mask=true_var > 1.0e-12_dp))
    end function zmax
  end subroutine test_djs

  !> DK 2.11: the local level model's steady state is P = x sigma2_eps with
  !> x = (q + sqrt(q^2 + 4 q)) / 2, q = sigma2_eta / sigma2_eps.
  subroutine test_steady_nile(error)
    type(error_type), allocatable, intent(out) :: error
    type(ssm_rep_t) :: r
    real(dp) :: y(1, 5), P(1, 1), F(1, 1), q, x
    integer :: info

    y = 0.0_dp
    r = ssm_rep(y, m=1, r=1)
    r%Z = 1.0_dp; r%T = 1.0_dp; r%R = 1.0_dp
    r%H = 15099.0_dp; r%Q = 1469.1_dp
    call steady_state(r, P, F, info)
    call check(error, info, SS_OK, "info")
    if (allocated(error)) return
    q = 1469.1_dp / 15099.0_dp
    x = (q + sqrt(q**2 + 4 * q)) / 2
    call check_rel(error, [P(1, 1), F(1, 1)], 15099.0_dp * [x, 1 + x], 1.0e-12_dp, &
                   "P, F")
  end subroutine test_steady_nile

  !> Correlated H, p = 2: the filter's P_t and F_t after many periods.
  subroutine test_steady_filter(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: r
    type(filter_result_t) :: fres
    real(dp), allocatable :: P(:, :), F(:, :)
    integer :: info, n

    fx = load_fixture("test/fixtures/mv_invariant.txt")
    r = rep_from_fixture(fx)
    n = 2000
    r%y = spread(spread(0.0_dp, 1, r%k_endog), 2, n)
    r%nobs = n
    call kalman_filter(r, fres, info)
    call check(error, info, SS_OK, "filter info")
    if (allocated(error)) return
    allocate (P(r%k_states, r%k_states), F(r%k_endog, r%k_endog))
    call steady_state(r, P, F, info)
    call check(error, info, SS_OK, "info")
    if (allocated(error)) return
    call check_rel(error, pack(P, .true.), pack(fres%P(:, :, n + 1), .true.), &
                   1.0e-10_dp, "P")
    if (allocated(error)) return
    call check_rel(error, pack(F, .true.), pack(fres%F(:, :, n), .true.), 1.0e-10_dp, &
                   "F")
  end subroutine test_steady_filter

  !> Singular H (the plain recursion): an invertible ARMA(1, 1) without
  !> measurement error has steady-state F = sigma2, the innovation variance.
  subroutine test_steady_arma(error)
    type(error_type), allocatable, intent(out) :: error
    type(component_holder_t) :: comps(1)
    type(structural_model_t) :: model
    real(dp) :: y(1, 5), P(2, 2), F(1, 1)
    integer :: info

    y = 0.0_dp
    comps(1)%c = arima_t(ar=1, ma=1)
    model = structural_model(y, comps, info)
    call model%update([0.6_dp, 0.3_dp, 2.0_dp])
    call steady_state(model%rep, P, F, info)
    call check(error, info, SS_OK, "info")
    if (allocated(error)) return
    call check_rel(error, [F(1, 1)], [2.0_dp], 1.0e-10_dp, "F = sigma2")
  end subroutine test_steady_arma

  !> A fixed level is learned like 1/t, so P -> 0 and F -> sigma2_eps; the
  !> doubling takes about 40 steps where the recursion would take ~1e12.
  subroutine test_steady_fixed(error)
    type(error_type), allocatable, intent(out) :: error
    type(ssm_rep_t) :: r
    real(dp) :: y(1, 5), P(1, 1), F(1, 1)
    integer :: info, niter

    y = 0.0_dp
    r = ssm_rep(y, m=1, r=1)
    r%Z = 1.0_dp; r%T = 1.0_dp; r%R = 1.0_dp
    r%H = 2.0_dp; r%Q = 0.0_dp
    call steady_state(r, P, F, info, niter=niter)
    call check(error, info, SS_OK, "info")
    if (allocated(error)) return
    call check(error, niter <= 60, "doubling steps")
    if (allocated(error)) return
    call check_rel(error, [P(1, 1), F(1, 1)], [0.0_dp, 2.0_dp], 1.0e-11_dp, "P, F")
  end subroutine test_steady_fixed
end module test_recursions
