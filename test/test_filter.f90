!> Kalman filter and state smoother against statsmodels fixtures.
module test_filter
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use statespace
  use fixture_io, only: fixture_t, load_fixture, rep_from_fixture
  implicit none
  private

  public :: collect_filter

  real(dp), parameter :: rtol = 1.0e-9_dp

contains

  subroutine collect_filter(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ new_unittest("nile_local_level_known", test_nile_known), &
                 new_unittest("multivariate_invariant", test_mv_invariant), &
                 new_unittest("multivariate_time_varying", test_mv_timevarying), &
                 new_unittest("nile_missing", test_nile_missing), &
                 new_unittest("multivariate_missing", test_mv_missing), &
                 new_unittest("multivariate_stationary", test_mv_stationary), &
                 new_unittest("forecast", test_forecast), &
                 new_unittest("exact_diffuse_nile_level", test_diffuse_nile_level), &
                 new_unittest("exact_diffuse_nile_trend", test_diffuse_nile_trend), &
                 new_unittest("exact_diffuse_multivariate", test_diffuse_mv), &
                 new_unittest("exact_diffuse_trend_seasonal_missing", &
                              test_diffuse_uc_seasonal), &
                 new_unittest("mixed_init_level_ar1", test_mixed_level_ar1), &
                 new_unittest("exact_diffuse_mv_time_varying", test_diffuse_mv_tv), &
                 new_unittest("exact_diffuse_vs_approximate", &
                              test_diffuse_vs_approximate), &
                 new_unittest("univariate_diffuse_nile_trend", test_uv_diffuse_trend), &
                 new_unittest("univariate_diffuse_mv", test_uv_diffuse_mv), &
                 new_unittest("univariate_diffuse_seasonal", &
                              test_uv_diffuse_seasonal), &
                 new_unittest("general_init_mv_diffuse", test_general_mv_diffuse), &
                 new_unittest("general_init_mixed", test_general_mixed), &
                 new_unittest("diffuse_mv_nile_level", test_dmv_nile_level), &
                 new_unittest("diffuse_mv_nile_trend", test_dmv_nile_trend), &
                 new_unittest("diffuse_mv_mv", test_dmv_mv), &
                 new_unittest("diffuse_mv_mv_time_varying", test_dmv_mv_tv), &
                 new_unittest("diffuse_mv_seasonal_missing", test_dmv_seasonal), &
                 new_unittest("diffuse_mv_mixed", test_dmv_mixed), &
                 new_unittest("diffuse_mv_singular_fallback", test_dmv_fallback), &
                 new_unittest("diffuse_mv_vs_univariate_time_varying", &
                              test_dmv_vs_uv_tv), &
                 new_unittest("diffuse_mv_branches_used", test_dmv_branches), &
                 new_unittest("sqrt_nile_known", test_sqrt_nile_known), &
                 new_unittest("sqrt_nile_missing", test_sqrt_nile_missing), &
                 new_unittest("sqrt_mv_invariant", test_sqrt_mv_invariant), &
                 new_unittest("sqrt_mv_time_varying", test_sqrt_mv_tv), &
                 new_unittest("sqrt_mv_missing", test_sqrt_mv_missing), &
                 new_unittest("sqrt_mv_stationary", test_sqrt_mv_stationary), &
                 new_unittest("sqrt_rejects_diffuse", test_sqrt_diffuse), &
                 new_unittest("collapse_known", test_collapse_known), &
                 new_unittest("collapse_diffuse", test_collapse_diffuse), &
                 new_unittest("linear_restrictions", test_restrictions), &
                 new_unittest("augmented_nile_level", test_aug_nile_level), &
                 new_unittest("augmented_nile_trend", test_aug_nile_trend), &
                 new_unittest("augmented_mv", test_aug_mv), &
                 new_unittest("augmented_mv_time_varying", test_aug_mv_tv), &
                 new_unittest("augmented_seasonal_missing", test_aug_seasonal), &
                 new_unittest("augmented_mixed", test_aug_mixed), &
                 new_unittest("univariate_nile_known", test_uv_nile_known), &
                 new_unittest("univariate_nile_missing", test_uv_nile_missing), &
                 new_unittest("univariate_mv_invariant", test_uv_mv_invariant), &
                 new_unittest("univariate_mv_time_varying", test_uv_mv_timevarying), &
                 new_unittest("univariate_mv_missing", test_uv_mv_missing), &
                 new_unittest("univariate_mv_stationary", test_uv_mv_stationary), &
                 new_unittest("stationary_init_rejects_unit_root", &
                              test_nonstationary), new_unittest("loglikelihood_burn", &
        test_burn), new_unittest("validate_errors", test_validate), &
                 new_unittest("steady_state_matches_full", test_steady_mv), &
                 new_unittest("steady_state_missing", test_steady_missing), &
                 new_unittest("steady_state_after_diffuse", test_steady_diffuse), &
                 new_unittest("steady_state_time_varying_c_d", test_steady_cd) ]
  end subroutine collect_filter

  !> Largest absolute difference relative to the scale of `expected`, over the
  !> elements where `mask` is true.
  pure real(dp) function relerr(actual, expected, mask)
    real(dp), intent(in) :: actual(:), expected(:)
    logical, intent(in) :: mask(:)

    relerr = maxval(abs(actual - expected), mask=mask) &
             / max(1.0_dp, maxval(abs(expected), mask=mask))
  end function relerr

  !> NaNs must match exactly; other elements must agree to `rtol`. With
  !> `mask`, only the elements where it is true are compared.
  subroutine check_close(error, actual, expected, label, mask)
    type(error_type), allocatable, intent(out) :: error
    real(dp), intent(in) :: actual(:), expected(:)
    character(len=*), intent(in) :: label
    logical, intent(in), optional :: mask(:)
    logical, allocatable :: sel(:), use(:)
    character(len=32) :: buf

    call check(error, size(actual) == size(expected), label//": size mismatch")
    if (allocated(error)) return
    allocate (sel(size(actual)), source=.true.)
    if (present(mask)) sel = mask
    call check(error, all((ieee_is_nan(actual) .eqv. ieee_is_nan(expected)) &
                          .or. .not. sel), label//": NaN pattern differs")
    if (allocated(error)) return
    use = sel .and. .not. ieee_is_nan(expected)
    if (.not. any(use)) return
    write (buf, '(es10.3)') relerr(actual, expected, use)
    call check(error, relerr(actual, expected, use) <= rtol, &
               label//": relative error "//trim(buf))
  end subroutine check_close

  !> Compare the full filter and smoother output with a fixture.
  !>
  !> In diffuse periods statsmodels reports v, F, Finf, yhat and the eps
  !> variances element by element in LDL-transformed coordinates, where ours
  !> are in original coordinates. For p = 1 the two coincide. For p > 1 the
  !> original-coordinate arrays are compared after the diffuse period only,
  !> and our per-element `uv_` arrays against statsmodels' within it.
  subroutine check_fixture(error, path, univariate, general_init, skip_smoother_until, &
                           diffuse_mv)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), intent(in) :: path
    logical, intent(in), optional :: univariate
    !> Use the multivariate exact initial filter (DIFFUSE_MULTIVARIATE).
    logical, intent(in), optional :: diffuse_mv
    !> Skip smoother comparisons for t <= this (where the fixture is known to
    !> be wrong).
    integer, intent(in), optional :: skip_smoother_until
    !> Re-specify the initialization with `initialize_general`, using the
    !> (a1, P_star, P_inf) implied by the fixture's initialization.
    logical, intent(in), optional :: general_init
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    integer :: info, p, m, n, nd, i, j, t
    real(dp) :: llf
    logical, allocatable :: obs(:, :), obs2(:, :, :), late(:), late2(:, :), &
                            late3(:, :, :)
    logical, allocatable :: conv3(:, :, :)
    real(dp), allocatable :: v_fx(:, :), F_fx(:, :, :), Finf_fx(:, :, :)
    real(dp), allocatable :: a1(:), Pstar1(:, :), Pinf1(:, :), blocks(:, :)
    integer :: k_diffuse
    logical, allocatable :: sm2m(:, :), sm3m(:, :, :), sm2p(:, :), sm3p(:, :, :)
    logical, allocatable :: sm2r(:, :), sm3r(:, :, :)
    logical :: keep

    fx = load_fixture(path)
    rep = rep_from_fixture(fx)
    if (present(univariate)) then
      if (univariate) rep%filter_method = FILTER_UNIVARIATE
    end if
    if (present(diffuse_mv)) then
      if (diffuse_mv) rep%diffuse_method = DIFFUSE_MULTIVARIATE
    end if
    if (present(general_init)) then
      if (general_init) then
        allocate (a1(rep%k_states), Pstar1(rep%k_states, rep%k_states), &
                  Pinf1(rep%k_states, rep%k_states))
        call rep%initial_state(a1, Pstar1, Pinf1, info)
        call check(error, info, SS_OK, "initial_state info")
        if (allocated(error)) return
        call rep%initialize_general(a1, Pstar1, Pinf1)
      end if
    end if
    call kalman_filter(rep, fres, info)
    call check(error, info, SS_OK, "kalman_filter info")
    if (allocated(error)) return
    call state_smoother(rep, fres, sres, info)
    call check(error, info, SS_OK, "state_smoother info")
    if (allocated(error)) return

    p = rep%k_endog; m = rep%k_states; n = rep%nobs
    nd = 0
    if (fx%has('nobs_diffuse')) nd = nint(sum(fx%get1('nobs_diffuse')))
    call check(error, fres%nobs_diffuse, nd, "nobs_diffuse")
    if (allocated(error)) return
    k_diffuse = 0
    if (fx%has('init_blocks')) then
      blocks = fx%get2('init_blocks')
      do i = 1, size(blocks, 1)
        if (nint(blocks(i, 3)) == INIT_DIFFUSE) then
          k_diffuse = k_diffuse + nint(blocks(i, 2) - blocks(i, 1)) + 1
        end if
      end do
    end if
    call check(error, fres%k_diffuse, k_diffuse, "k_diffuse")
    if (allocated(error)) return

    ! Time masks: original-coordinate arrays, and conventional-gain arrays.
    late = [(t > nd .or. p == 1, t=1, n)]
    allocate (late2(p, n), late3(p, p, n), conv3(m, p, n))
    do t = 1, n
      late2(:, t) = late(t)
      late3(:, :, t) = late(t)
      conv3(:, :, t) = fres%method(t) == METHOD_CONVENTIONAL .and. t > nd
    end do

    call check_close(error, pack(fres%yhat, .true.), fx%get1('yhat'), "yhat", &
                     pack(late2, .true.))
    if (allocated(error)) return
    call check_close(error, pack(fres%a, .true.), fx%get1('a'), "a")
    if (allocated(error)) return
    call check_close(error, pack(fres%P, .true.), fx%get1('P'), "P")
    if (allocated(error)) return
    call check_close(error, pack(fres%v, .true.), fx%get1('v'), "v", &
                     pack(late2, .true.))
    if (allocated(error)) return
    call check_close(error, pack(fres%F, .true.), fx%get1('F'), "F", &
                     pack(late3, .true.))
    if (allocated(error)) return
    call check_close(error, pack(fres%att, .true.), fx%get1('att'), "att")
    if (allocated(error)) return
    call check_close(error, pack(fres%Ptt, .true.), fx%get1('Ptt'), "Ptt")
    if (allocated(error)) return
    call check_close(error, pack(fres%K, .true.), fx%get1('K'), "K", &
                     pack(conv3, .true.))
    if (allocated(error)) return
    call check_close(error, fres%llf_obs, fx%get1('llf_obs'), "llf_obs")
    if (allocated(error)) return
    call check_close(error, [fres%llf], fx%get1('llf'), "llf")
    if (allocated(error)) return

    obs = .not. ieee_is_nan(rep%y)
    allocate (obs2(p, p, n))
    do i = 1, p
      do j = 1, p
        obs2(i, j, :) = obs(i, :) .and. obs(j, :)
      end do
    end do

    ! statsmodels reports Finf = 0 for missing y; ours is still Z P_inf Z'.
    if (fx%has('Pinf')) then
      call check_close(error, pack(fres%Pinf, .true.), fx%get1('Pinf'), "Pinf")
      if (allocated(error)) return
      call check_close(error, pack(fres%Finf, .true.), fx%get1('Finf'), "Finf", &
                       pack(late3 .and. obs2, .true.))
      if (allocated(error)) return
    end if

    ! Per-element univariate quantities in diffuse periods (p > 1).
    if (p > 1 .and. nd > 0) then
      v_fx = fx%get2('v')
      F_fx = fx%get3('F')
      Finf_fx = fx%get3('Finf')
      do t = 1, nd
        do i = 1, fres%uv_n(t)
          j = fres%uv_idx(i, t)
          call check_close(error, [fres%uv_v(i, t), fres%uv_Fstar(i, t), &
                                   fres%uv_Finf(i, t)], [v_fx(j, t), F_fx(j, j, t), &
              Finf_fx(j, j, t)], "diffuse elements")
          if (allocated(error)) return
        end do
      end do
    end if

    ! Smoother time masks
    allocate (sm2m(m, n), sm3m(m, m, n), sm2p(p, n), sm2r(rep%k_posdef, n), &
              sm3r(rep%k_posdef, rep%k_posdef, n), sm3p(p, p, n))
    do t = 1, n
      keep = .true.
      if (present(skip_smoother_until)) keep = t > skip_smoother_until
      sm2m(:, t) = keep; sm3m(:, :, t) = keep; sm2p(:, t) = keep; sm3p(:, :, t) = keep
      sm2r(:, t) = keep; sm3r(:, :, t) = keep
    end do

    call check_close(error, pack(sres%alphahat, .true.), fx%get1('alphahat'), &
                     "alphahat", pack(sm2m, .true.))
    if (allocated(error)) return
    call check_close(error, pack(sres%V, .true.), fx%get1('V'), "V", pack(sm3m, .true.))
    if (allocated(error)) return
    ! For missing y, statsmodels reports eps moments of 0 and H rather than the
    ! conditional moments, so compare observed elements only.
    call check_close(error, pack(sres%epshat, .true.), fx%get1('epshat'), "epshat", &
                     pack(obs .and. sm2p, .true.))
    if (allocated(error)) return
    call check_close(error, pack(sres%epsvar, .true.), fx%get1('epsvar'), "epsvar", &
                     pack(obs2 .and. late3 .and. sm3p, .true.))
    if (allocated(error)) return
    call check_close(error, pack(sres%etahat, .true.), fx%get1('etahat'), "etahat", &
                     pack(sm2r, .true.))
    if (allocated(error)) return
    call check_close(error, pack(sres%etavar, .true.), fx%get1('etavar'), "etavar", &
                     pack(sm3r, .true.))
    if (allocated(error)) return

    ! The likelihood-only path must agree with the full filter.
    llf = loglike(rep, info)
    call check(error, info, SS_OK, "loglike info")
    if (allocated(error)) return
    call check_close(error, [llf], [fres%llf], "loglike vs kalman_filter")
  end subroutine check_fixture

  subroutine test_nile_known(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/nile_llevel_known.txt")
  end subroutine test_nile_known

  subroutine test_mv_invariant(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/mv_invariant.txt")
  end subroutine test_mv_invariant

  subroutine test_mv_timevarying(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/mv_timevarying.txt")
  end subroutine test_mv_timevarying

  subroutine test_diffuse_nile_level(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/nile_llevel_exact.txt")
  end subroutine test_diffuse_nile_level

  subroutine test_diffuse_nile_trend(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/nile_lltrend_exact.txt")
  end subroutine test_diffuse_nile_trend

  subroutine test_diffuse_mv(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/mv_diffuse.txt")
  end subroutine test_diffuse_mv

  subroutine test_diffuse_uc_seasonal(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/uc_trend_seasonal_exact.txt")
  end subroutine test_diffuse_uc_seasonal

  subroutine test_mixed_level_ar1(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/uc_level_ar1_mixed.txt")
  end subroutine test_mixed_level_ar1

  ! statsmodels' univariate diffuse smoother uses the wrong transition matrix
  ! when T varies within the diffuse period (its source has a TODO there), so
  ! its smoothed output for t <= nobs_diffuse is not a reference. Those periods
  ! are checked against the large-kappa limit in test_diffuse_vs_approximate.
  subroutine test_diffuse_mv_tv(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/mv_diffuse_timevarying.txt", &
                       skip_smoother_until=2)
  end subroutine test_diffuse_mv_tv

  !> Exact diffuse smoothing is the kappa -> infinity limit of N(a, kappa I)
  !> (DK 5.1). With time-varying matrices inside the diffuse period, the gap
  !> between the exact smoother and the conventional one at kappa must shrink
  !> like 1/kappa. (Beyond kappa ~ 1e7 rounding in the conventional path
  !> dominates, so kappa = 1e5 and 1e6 are used.)
  subroutine test_diffuse_vs_approximate(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: exact
    type(filter_result_t) :: fe
    type(smoother_result_t) :: se
    real(dp) :: ga(2), gV(2)
    integer :: info, k
    character(len=64) :: buf

    fx = load_fixture("test/fixtures/mv_diffuse_timevarying.txt")
    exact = rep_from_fixture(fx)
    call kalman_filter(exact, fe, info)
    call state_smoother(exact, fe, se, info)
    call check(error, fe%nobs_diffuse, 2, "nobs_diffuse")
    if (allocated(error)) return
    do k = 1, 2
      call gap(10.0_dp**(4 + k), ga(k), gV(k))
    end do
    write (buf, '(4es10.2)') ga, gV
    call check(error, ga(2) < 1.0e-3_dp .and. gV(2) < 2.0e-2_dp, &
               "gap at kappa = 1e6: "//buf)
    if (allocated(error)) return
    call check(error, abs(ga(1) / ga(2) - 10.0_dp) < 0.5_dp .and. &
               abs(gV(1) / gV(2) - 10.0_dp) < 0.5_dp, "gap not O(1/kappa): "//buf)

  contains

    subroutine gap(kappa, galpha, gvar)
      real(dp), intent(in) :: kappa
      real(dp), intent(out) :: galpha, gvar
      type(ssm_rep_t) :: approx
      type(filter_result_t) :: fa
      type(smoother_result_t) :: sa

      approx = exact
      call approx%initialize_approximate_diffuse(kappa)
      call kalman_filter(approx, fa, info)
      call state_smoother(approx, fa, sa, info)
      galpha = maxval(abs(se%alphahat - sa%alphahat))
      gvar = maxval(abs(se%V - sa%V))
    end subroutine gap
  end subroutine test_diffuse_vs_approximate

  ! FILTER_UNIVARIATE after the diffuse period: the univariate smoother must
  ! hand over correctly to the diffuse recursion.
  subroutine test_uv_diffuse_trend(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/nile_lltrend_exact.txt", univariate=.true.)
  end subroutine test_uv_diffuse_trend

  subroutine test_uv_diffuse_mv(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/mv_diffuse_timevarying.txt", &
                       univariate=.true., skip_smoother_until=2)
  end subroutine test_uv_diffuse_mv

  subroutine test_uv_diffuse_seasonal(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/uc_trend_seasonal_exact.txt", &
                       univariate=.true.)
  end subroutine test_uv_diffuse_seasonal

  ! DK 5.1 general form, a1 + A delta + R0 eta0, must reproduce the block
  ! initializations; the mixed model has a rank-1 P_inf.
  subroutine test_general_mv_diffuse(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/mv_diffuse.txt", general_init=.true.)
  end subroutine test_general_mv_diffuse

  subroutine test_general_mixed(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/uc_level_ar1_mixed.txt", &
                       general_init=.true.)
  end subroutine test_general_mixed

  ! Multivariate exact initial filter and smoother (DK 5.2-5.3).
  subroutine test_dmv_nile_level(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/nile_llevel_exact.txt", diffuse_mv=.true.)
  end subroutine test_dmv_nile_level

  subroutine test_dmv_nile_trend(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/nile_lltrend_exact.txt", diffuse_mv=.true.)
  end subroutine test_dmv_nile_trend

  subroutine test_dmv_mv(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/mv_diffuse.txt", diffuse_mv=.true.)
  end subroutine test_dmv_mv

  subroutine test_dmv_mv_tv(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/mv_diffuse_timevarying.txt", &
                       diffuse_mv=.true., skip_smoother_until=2)
  end subroutine test_dmv_mv_tv

  subroutine test_dmv_seasonal(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/uc_trend_seasonal_exact.txt", &
                       diffuse_mv=.true.)
  end subroutine test_dmv_seasonal

  subroutine test_dmv_mixed(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/uc_level_ar1_mixed.txt", diffuse_mv=.true.)
  end subroutine test_dmv_mixed

  !> Bivariate series with a common diffuse level: F_inf = [1 1; 1 1] is
  !> singular but nonzero, so the multivariate filter falls back to the
  !> univariate step, and all results equal the univariate path.
  subroutine test_dmv_fallback(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: uv, mv
    type(filter_result_t) :: fu, fm
    type(smoother_result_t) :: su, sm
    integer :: info

    fx = load_fixture("test/fixtures/mv_invariant.txt")
    uv = ssm_rep(fx%get2('y'), 1, 1)
    uv%Z(:, 1, 1) = [1.0_dp, 1.0_dp]
    uv%H(:, :, 1) = reshape([1.0_dp, 0.3_dp, 0.3_dp, 2.0_dp], [2, 2])
    uv%T = 1.0_dp
    uv%Q = 0.5_dp
    call uv%initialize_diffuse()
    mv = uv
    mv%diffuse_method = DIFFUSE_MULTIVARIATE
    call kalman_filter(uv, fu, info)
    call state_smoother(uv, fu, su, info)
    call kalman_filter(mv, fm, info)
    call check(error, info, SS_OK, "filter info")
    if (allocated(error)) return
    call state_smoother(mv, fm, sm, info)
    call check(error, fm%method(1), METHOD_UNIVARIATE, "fallback at t = 1")
    if (allocated(error)) return
    call check_close(error, [fm%llf], [fu%llf], "llf")
    if (allocated(error)) return
    call check_close(error, pack(fm%P, .true.), pack(fu%P, .true.), "P")
    if (allocated(error)) return
    call check_close(error, pack(sm%alphahat, .true.), pack(su%alphahat, .true.), &
                     "alphahat")
    if (allocated(error)) return
    call check_close(error, pack(sm%V, .true.), pack(su%V, .true.), "V")
  end subroutine test_dmv_fallback

  !> With time-varying matrices in the diffuse period the statsmodels fixture
  !> is wrong, so compare the multivariate path with the univariate one (which
  !> test_diffuse_vs_approximate anchors to the kappa -> infinity limit).
  subroutine test_dmv_vs_uv_tv(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: uv, mv
    type(filter_result_t) :: fu, fm
    type(smoother_result_t) :: su, sm
    integer :: info

    fx = load_fixture("test/fixtures/mv_diffuse_timevarying.txt")
    uv = rep_from_fixture(fx)
    mv = uv
    mv%diffuse_method = DIFFUSE_MULTIVARIATE
    call kalman_filter(uv, fu, info)
    call state_smoother(uv, fu, su, info)
    call kalman_filter(mv, fm, info)
    call state_smoother(mv, fm, sm, info)
    call check(error, count(fm%method(1:fm%nobs_diffuse) == METHOD_DIFFUSE_MV) > 0, &
               "multivariate diffuse update used")
    if (allocated(error)) return
    call check_close(error, fm%llf_obs, fu%llf_obs, "llf_obs")
    if (allocated(error)) return
    call check_close(error, pack(sm%alphahat, .true.), pack(su%alphahat, .true.), &
                     "alphahat")
    if (allocated(error)) return
    call check_close(error, pack(sm%V, .true.), pack(su%V, .true.), "V")
    if (allocated(error)) return
    call check_close(error, pack(sm%etahat, .true.), pack(su%etahat, .true.), "etahat")
    if (allocated(error)) return
    call check_close(error, pack(sm%epsvar, .true.), pack(su%epsvar, .true.), "epsvar")
  end subroutine test_dmv_vs_uv_tv

  !> The fixtures above must exercise the nonsingular multivariate update and
  !> the F_inf = 0 conventional case within the diffuse period.
  subroutine test_dmv_branches(error)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), parameter :: paths(3) = [character(len=48) :: &
        "test/fixtures/mv_diffuse.txt", &
        "test/fixtures/uc_trend_seasonal_exact.txt", &
        "test/fixtures/uc_level_ar1_mixed.txt"]
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    integer :: info, k, n_mv, n_zero

    n_mv = 0
    n_zero = 0
    do k = 1, size(paths)
      fx = load_fixture(trim(paths(k)))
      rep = rep_from_fixture(fx)
      rep%diffuse_method = DIFFUSE_MULTIVARIATE
      call kalman_filter(rep, fres, info)
      n_mv = n_mv + count(fres%method(1:fres%nobs_diffuse) == METHOD_DIFFUSE_MV)
      n_zero = n_zero + count(fres%method(1:fres%nobs_diffuse) == METHOD_CONVENTIONAL)
    end do
    call check(error, n_mv > 0 .and. n_zero > 0, "both diffuse branches used")
  end subroutine test_dmv_branches

  !> The square root filter and smoother (DK 6.3) must reproduce the
  !> conventional ones.
  subroutine check_sqrt(error, path)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), intent(in) :: path
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fc, fs
    type(smoother_result_t) :: sc, ss
    real(dp), allocatable :: Pchol(:, :, :)
    integer :: info, t

    fx = load_fixture(path)
    rep = rep_from_fixture(fx)
    call kalman_filter(rep, fc, info)
    call state_smoother(rep, fc, sc, info)
    call sqrt_kalman_filter(rep, fs, info, Pchol)
    call check(error, info, SS_OK, "sqrt_kalman_filter info")
    if (allocated(error)) return
    call check_close(error, pack(fs%a, .true.), pack(fc%a, .true.), "a")
    if (allocated(error)) return
    call check_close(error, pack(fs%P, .true.), pack(fc%P, .true.), "P")
    if (allocated(error)) return
    call check_close(error, pack(fs%F, .true.), pack(fc%F, .true.), "F")
    if (allocated(error)) return
    call check_close(error, pack(fs%K, .true.), pack(fc%K, .true.), "K")
    if (allocated(error)) return
    call check_close(error, pack(fs%Ptt, .true.), pack(fc%Ptt, .true.), "Ptt")
    if (allocated(error)) return
    call check_close(error, fs%llf_obs, fc%llf_obs, "llf_obs")
    if (allocated(error)) return
    do t = 1, rep%nobs + 1
      call check_close(error, pack(matmul(Pchol(:, :, t), transpose(Pchol(:, :, t))), &
                                   .true.), pack(fc%P(:, :, t), .true.), "P~ P~'")
      if (allocated(error)) return
    end do
    call sqrt_state_smoother(rep, fs, ss, info)
    call check(error, info, SS_OK, "sqrt_state_smoother info")
    if (allocated(error)) return
    call check_close(error, pack(ss%alphahat, .true.), pack(sc%alphahat, .true.), &
                     "alphahat")
    if (allocated(error)) return
    call check_close(error, pack(ss%N, .true.), pack(sc%N, .true.), "N")
    if (allocated(error)) return
    call check_close(error, pack(ss%V, .true.), pack(sc%V, .true.), "V")
  end subroutine check_sqrt

  subroutine test_sqrt_nile_known(error)
    type(error_type), allocatable, intent(out) :: error

    call check_sqrt(error, "test/fixtures/nile_llevel_known.txt")
  end subroutine test_sqrt_nile_known

  subroutine test_sqrt_nile_missing(error)
    type(error_type), allocatable, intent(out) :: error

    call check_sqrt(error, "test/fixtures/nile_llevel_missing.txt")
  end subroutine test_sqrt_nile_missing

  subroutine test_sqrt_mv_invariant(error)
    type(error_type), allocatable, intent(out) :: error

    call check_sqrt(error, "test/fixtures/mv_invariant.txt")
  end subroutine test_sqrt_mv_invariant

  subroutine test_sqrt_mv_tv(error)
    type(error_type), allocatable, intent(out) :: error

    call check_sqrt(error, "test/fixtures/mv_timevarying.txt")
  end subroutine test_sqrt_mv_tv

  subroutine test_sqrt_mv_missing(error)
    type(error_type), allocatable, intent(out) :: error

    call check_sqrt(error, "test/fixtures/mv_missing.txt")
  end subroutine test_sqrt_mv_missing

  subroutine test_sqrt_mv_stationary(error)
    type(error_type), allocatable, intent(out) :: error

    call check_sqrt(error, "test/fixtures/mv_stationary.txt")
  end subroutine test_sqrt_mv_stationary

  subroutine test_sqrt_diffuse(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    integer :: info

    fx = load_fixture("test/fixtures/nile_llevel_exact.txt")
    rep = rep_from_fixture(fx)
    call sqrt_kalman_filter(rep, fres, info)
    call check(error, info, SS_ERR_UNSUPPORTED, "exact diffuse not supported")
  end subroutine test_sqrt_diffuse

  !> A p = 5, m = 2 local linear trend with correlated H, simulated, with
  !> partly and fully missing periods.
  function wide_trend_model() result(rep)
    type(ssm_rep_t) :: rep
    integer, parameter :: p = 5, m = 2, n = 60
    real(dp) :: y(p, n), alpha(m, n), eps(p, n), eta(m, n), u_init(m), u_eps(p, n), &
                u_eta(m, n)
    real(dp) :: A(p, p)
    integer :: info, i, j

    rep = ssm_rep(y, m, m)
    do j = 1, m
      do i = 1, p
        rep%Z(i, j, 1) = 1.0_dp + 0.3_dp * i * (j - 1) &
                         + 0.1_dp * sin(real(i + 3 * j, dp))
      end do
    end do
    do j = 1, p
      do i = 1, p
        A(i, j) = 0.2_dp * cos(real(i * j, dp))
      end do
    end do
    rep%H(:, :, 1) = matmul(A, transpose(A)) + eye(p)
    rep%T(:, :, 1) = reshape([1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [2, 2])
    rep%Q(:, :, 1) = reshape([0.5_dp, 0.0_dp, 0.0_dp, 0.05_dp], [2, 2])
    rep%d(:, 1) = [(0.1_dp * i, i=1, p)]
    call rep%initialize_known([10.0_dp, 0.5_dp], eye(m))
    do j = 1, n
      do i = 1, p
        u_eps(i, j) = sin(1.7_dp * i + 2.3_dp * j)
      end do
      u_eta(:, j) = [cos(0.9_dp * j), sin(1.3_dp * j)]
    end do
    u_init = [0.3_dp, -0.2_dp]
    call simulate(rep, u_init, u_eps, u_eta, y, alpha, eps, eta, info)
    y(2, 5) = ieee_value(1.0_dp, ieee_quiet_nan)
    y(1:3, 12) = ieee_value(1.0_dp, ieee_quiet_nan)
    y(:, 20) = ieee_value(1.0_dp, ieee_quiet_nan)
    rep%y = y
  end function wide_trend_model

  !> Collapsing (DK 6.5): the collapsed model plus the adjustment must give
  !> the same log likelihood and the same predicted and smoothed states.
  subroutine check_collapse(error, rep)
    type(error_type), allocatable, intent(out) :: error
    type(ssm_rep_t), intent(in) :: rep
    type(ssm_rep_t) :: crep
    type(filter_result_t) :: ff, fc
    type(smoother_result_t) :: sf, sc
    real(dp), allocatable :: adj(:)
    integer :: info

    allocate (adj(rep%nobs))
    call collapse_observations(rep, crep, adj, info)
    call check(error, info, SS_OK, "collapse info")
    if (allocated(error)) return
    call check(error, crep%k_endog, rep%k_states, "collapsed dimension")
    if (allocated(error)) return
    call kalman_filter(rep, ff, info)
    call state_smoother(rep, ff, sf, info)
    call kalman_filter(crep, fc, info)
    call check(error, info, SS_OK, "collapsed filter info")
    if (allocated(error)) return
    call state_smoother(crep, fc, sc, info)
    call check_close(error, [fc%llf + sum(adj)], [ff%llf], "llf")
    if (allocated(error)) return
    call check_close(error, [loglike(crep, info) + sum(adj)], [ff%llf], "loglike")
    if (allocated(error)) return
    call check_close(error, pack(fc%a, .true.), pack(ff%a, .true.), "a")
    if (allocated(error)) return
    call check_close(error, pack(fc%P, .true.), pack(ff%P, .true.), "P")
    if (allocated(error)) return
    call check_close(error, pack(sc%alphahat, .true.), pack(sf%alphahat, .true.), &
                     "alphahat")
    if (allocated(error)) return
    call check_close(error, pack(sc%V, .true.), pack(sf%V, .true.), "V")
  end subroutine check_collapse

  subroutine test_collapse_known(error)
    type(error_type), allocatable, intent(out) :: error

    call check_collapse(error, wide_trend_model())
  end subroutine test_collapse_known

  subroutine test_collapse_diffuse(error)
    type(error_type), allocatable, intent(out) :: error
    type(ssm_rep_t) :: rep

    rep = wide_trend_model()
    call rep%initialize_diffuse()
    call check_collapse(error, rep)
  end subroutine test_collapse_diffuse

  !> Linear restrictions (DK 6.6): alpha_1,t + 2 alpha_2,t = 12 in periods
  !> 10..30 only. The smoothed states must satisfy them exactly, with
  !> R* V_t = 0, and the conventional and univariate filters must agree.
  subroutine test_restrictions(error)
    type(error_type), allocatable, intent(out) :: error
    type(ssm_rep_t) :: rep, rrep, urep
    type(filter_result_t) :: fres, ures
    type(smoother_result_t) :: sres, usres
    real(dp), allocatable :: Rmat(:, :, :), rval(:, :)
    real(dp) :: worst_mean, worst_var
    integer :: info, t
    character(len=64) :: buf

    rep = wide_trend_model()
    allocate (Rmat(1, 2, 1), rval(1, rep%nobs))
    Rmat(1, :, 1) = [1.0_dp, 2.0_dp]
    rval = ieee_value(1.0_dp, ieee_quiet_nan)
    rval(1, 10:30) = 12.0_dp
    call add_state_restrictions(rep, Rmat, rval, rrep, info)
    call check(error, info, SS_OK, "add_state_restrictions info")
    if (allocated(error)) return

    call kalman_filter(rrep, fres, info)
    call check(error, info, SS_OK, "filter info")
    if (allocated(error)) return
    call state_smoother(rrep, fres, sres, info)
    worst_mean = 0.0_dp
    worst_var = 0.0_dp
    do t = 10, 30
      worst_mean = max(worst_mean, abs(dot_product(Rmat(1, :, 1), sres%alphahat(:, t)) &
                                       - 12.0_dp))
      worst_var = max(worst_var, maxval(abs(matmul(Rmat(1, :, 1), sres%V(:, :, t)))))
    end do
    write (buf, '(2es10.2)') worst_mean, worst_var
    call check(error, worst_mean < 1.0e-8_dp .and. worst_var < 1.0e-8_dp, &
               "restriction violated: "//trim(buf))
    if (allocated(error)) return
    ! Outside 10..30 the restriction is not imposed.
    call check(error, abs(dot_product(Rmat(1, :, 1), sres%alphahat(:, 45)) - 12.0_dp) &
               > 1.0e-3_dp, "restriction leaked outside its periods")
    if (allocated(error)) return

    urep = rrep
    urep%filter_method = FILTER_UNIVARIATE
    call kalman_filter(urep, ures, info)
    call state_smoother(urep, ures, usres, info)
    call check_close(error, [ures%llf], [fres%llf], "llf, univariate vs conventional")
    if (allocated(error)) return
    call check_close(error, pack(usres%alphahat, .true.), pack(sres%alphahat, .true.), &
                     "alphahat, univariate vs conventional")
  end subroutine test_restrictions

  !> The augmented filter and smoother (DK 5.7) must reproduce the exact
  !> initial filter: the diffuse log likelihood and the smoothed state.
  subroutine check_augmented(error, path)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), intent(in) :: path
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    type(augmented_result_t) :: ares
    real(dp), allocatable :: alphahat(:, :), V(:, :, :)
    integer :: info

    fx = load_fixture(path)
    rep = rep_from_fixture(fx)
    call kalman_filter(rep, fres, info)
    call state_smoother(rep, fres, sres, info)
    call augmented_filter(rep, ares, info)
    call check(error, info, SS_OK, "augmented_filter info")
    if (allocated(error)) return
    call check(error, ares%k_diffuse, fres%k_diffuse, "k_diffuse")
    if (allocated(error)) return
    call check_close(error, [ares%llf], [fres%llf], "llf")
    if (allocated(error)) return
    allocate (alphahat(rep%k_states, rep%nobs), V(rep%k_states, rep%k_states, rep%nobs))
    call augmented_smoother(rep, ares, alphahat, V, info)
    call check(error, info, SS_OK, "augmented_smoother info")
    if (allocated(error)) return
    call check_close(error, pack(alphahat, .true.), pack(sres%alphahat, .true.), &
                     "alphahat")
    if (allocated(error)) return
    call check_close(error, pack(V, .true.), pack(sres%V, .true.), "V")
  end subroutine check_augmented

  subroutine test_aug_nile_level(error)
    type(error_type), allocatable, intent(out) :: error

    call check_augmented(error, "test/fixtures/nile_llevel_exact.txt")
  end subroutine test_aug_nile_level

  subroutine test_aug_nile_trend(error)
    type(error_type), allocatable, intent(out) :: error

    call check_augmented(error, "test/fixtures/nile_lltrend_exact.txt")
  end subroutine test_aug_nile_trend

  subroutine test_aug_mv(error)
    type(error_type), allocatable, intent(out) :: error

    call check_augmented(error, "test/fixtures/mv_diffuse.txt")
  end subroutine test_aug_mv

  subroutine test_aug_mv_tv(error)
    type(error_type), allocatable, intent(out) :: error

    call check_augmented(error, "test/fixtures/mv_diffuse_timevarying.txt")
  end subroutine test_aug_mv_tv

  subroutine test_aug_seasonal(error)
    type(error_type), allocatable, intent(out) :: error

    call check_augmented(error, "test/fixtures/uc_trend_seasonal_exact.txt")
  end subroutine test_aug_seasonal

  subroutine test_aug_mixed(error)
    type(error_type), allocatable, intent(out) :: error

    call check_augmented(error, "test/fixtures/uc_level_ar1_mixed.txt")
  end subroutine test_aug_mixed

  ! The univariate treatment must reproduce the conventional results, including
  ! with correlated (and time-varying) H and missing data.
  subroutine test_uv_nile_known(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/nile_llevel_known.txt", univariate=.true.)
  end subroutine test_uv_nile_known

  subroutine test_uv_nile_missing(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/nile_llevel_missing.txt", &
                       univariate=.true.)
  end subroutine test_uv_nile_missing

  subroutine test_uv_mv_invariant(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/mv_invariant.txt", univariate=.true.)
  end subroutine test_uv_mv_invariant

  subroutine test_uv_mv_timevarying(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/mv_timevarying.txt", univariate=.true.)
  end subroutine test_uv_mv_timevarying

  subroutine test_uv_mv_missing(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/mv_missing.txt", univariate=.true.)
  end subroutine test_uv_mv_missing

  subroutine test_uv_mv_stationary(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/mv_stationary.txt", univariate=.true.)
  end subroutine test_uv_mv_stationary

  subroutine test_nile_missing(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/nile_llevel_missing.txt")
  end subroutine test_nile_missing

  subroutine test_mv_missing(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/mv_missing.txt")
  end subroutine test_mv_missing

  subroutine test_mv_stationary(error)
    type(error_type), allocatable, intent(out) :: error

    call check_fixture(error, "test/fixtures/mv_stationary.txt")
  end subroutine test_mv_stationary

  !> Out-of-sample forecasts from the observed part of mv_missing must equal
  !> statsmodels' one-step forecasts over its trailing missing periods.
  subroutine test_forecast(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: h = 8
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    type(forecast_result_t) :: fc
    real(dp), allocatable :: yhat(:, :), F(:, :, :), a(:, :), P(:, :, :)
    integer :: info, n

    fx = load_fixture("test/fixtures/mv_missing.txt")
    rep = rep_from_fixture(fx)
    n = rep%nobs - h
    rep%y = rep%y(:, 1:n)
    rep%nobs = n
    call kalman_filter(rep, fres, info)
    call check(error, info, SS_OK, "kalman_filter info")
    if (allocated(error)) return
    call forecast(rep, fres, h, fc, info)
    call check(error, info, SS_OK, "forecast info")
    if (allocated(error)) return

    yhat = fx%get2('yhat')
    F = fx%get3('F')
    a = fx%get2('a')
    P = fx%get3('P')
    call check_close(error, pack(fc%mean, .true.), pack(yhat(:, n + 1:), .true.), &
                     "mean")
    if (allocated(error)) return
    call check_close(error, pack(fc%cov, .true.), pack(F(:, :, n + 1:), .true.), "cov")
    if (allocated(error)) return
    call check_close(error, pack(fc%state, .true.), pack(a(:, n + 1:n + h), .true.), &
                     "state")
    if (allocated(error)) return
    call check_close(error, pack(fc%state_cov, .true.), &
                     pack(P(:, :, n + 1:n + h), .true.), "state_cov")
  end subroutine test_forecast

  subroutine test_nonstationary(error)
    type(error_type), allocatable, intent(out) :: error
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    real(dp) :: y(1, 5)
    integer :: info

    y = 1.0_dp
    rep = ssm_rep(y, 1, 1)
    rep%Z = 1.0_dp
    rep%T = 1.0_dp     ! random walk
    rep%Q = 1.0_dp
    call rep%initialize_stationary()
    call kalman_filter(rep, fres, info)
    call check(error, info, SS_ERR_NOT_STATIONARY, &
               "random walk has no stationary distribution")
  end subroutine test_nonstationary

  subroutine test_burn(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    real(dp), allocatable :: llf_obs(:)
    integer :: info

    fx = load_fixture("test/fixtures/nile_llevel_known.txt")
    rep = rep_from_fixture(fx)
    rep%loglikelihood_burn = 1
    llf_obs = fx%get1('llf_obs')
    call kalman_filter(rep, fres, info)
    call check_close(error, [fres%llf], [sum(llf_obs(2:))], "burned llf")
    if (allocated(error)) return
    call check_close(error, [loglike(rep, info)], [sum(llf_obs(2:))], "burned loglike")
  end subroutine test_burn

  subroutine test_validate(error)
    type(error_type), allocatable, intent(out) :: error
    type(ssm_rep_t) :: rep
    integer :: info
    real(dp) :: y(1, 5)

    y = 1.0_dp
    rep = ssm_rep(y, 1, 1)
    call rep%validate(info)
    call check(error, info, SS_ERR_INIT, "uninitialized")
    if (allocated(error)) return

    call rep%initialize_approximate_diffuse()
    call rep%validate(info)
    call check(error, info, SS_OK, "valid")
    if (allocated(error)) return

    rep%T = reshape([1.0_dp, 1.0_dp], [1, 1, 2])   ! time dimension neither 1 nor n
    call rep%validate(info)
    call check(error, info, SS_ERR_DIM, "bad time dimension")
    if (allocated(error)) return

    rep%T = reshape([1.0_dp], [1, 1, 1])
    rep%y(1, 3) = ieee_value(1.0_dp, ieee_quiet_nan)
    call rep%validate(info)
    call check(error, info, SS_OK, "missing data is valid")
  end subroutine test_validate

  !> The steady-state shortcut (DK 4.3.4) against the full recursion
  !> (tol_steady < 0): filter and smoother output and loglike.
  subroutine check_steady_pair(error, rep, t_min)
    type(error_type), allocatable, intent(out) :: error
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: t_min            !< the shortcut must start after this
    type(ssm_rep_t) :: full
    type(filter_result_t) :: f1, f2
    type(smoother_result_t) :: s1, s2
    real(dp) :: s1_scale, s2_scale
    integer :: info

    full = rep
    full%tol_steady = -1.0_dp
    call kalman_filter(rep, f1, info)
    call check(error, info, SS_OK, "filter")
    if (allocated(error)) return
    call kalman_filter(full, f2, info)
    call check(error, f1%t_steady > t_min .and. f1%t_steady < rep%nobs / 2, &
               "steady state used from a sensible period")
    if (allocated(error)) return
    call check(error, f2%t_steady == 0, "no steady state when disabled")
    if (allocated(error)) return
    call check_close(error, [f1%llf, loglike(rep, info)], &
                     [f2%llf, loglike(full, info)], "llf")
    if (allocated(error)) return
    call check_close(error, [loglike_concentrated(rep, s1_scale, info)], &
                     [loglike_concentrated(full, s2_scale, info)], "concentrated llf")
    if (allocated(error)) return
    call check_close(error, [s1_scale], [s2_scale], "concentrated scale")
    if (allocated(error)) return
    call check_close(error, pack(f1%a, .true.), pack(f2%a, .true.), "a")
    if (allocated(error)) return
    call check_close(error, pack(f1%P, .true.), pack(f2%P, .true.), "P")
    if (allocated(error)) return
    call check_close(error, pack(f1%v, .true.), pack(f2%v, .true.), "v")
    if (allocated(error)) return
    call check_close(error, pack(f1%F, .true.), pack(f2%F, .true.), "F")
    if (allocated(error)) return
    call check_close(error, pack(f1%K, .true.), pack(f2%K, .true.), "K")
    if (allocated(error)) return
    call state_smoother(rep, f1, s1, info)
    call state_smoother(full, f2, s2, info)
    call check_close(error, pack(s1%alphahat, .true.), pack(s2%alphahat, .true.), &
                     "alphahat")
    if (allocated(error)) return
    call check_close(error, pack(s1%V, .true.), pack(s2%V, .true.), "V")
  end subroutine check_steady_pair

  !> The mv_invariant model (correlated H) on a longer simulated series.
  function long_mv(n) result(rep)
    integer, intent(in) :: n
    type(ssm_rep_t) :: rep
    type(fixture_t) :: fx
    real(dp), allocatable :: e(:)

    fx = load_fixture("test/fixtures/mv_invariant.txt")
    rep = rep_from_fixture(fx)
    allocate (e(rep%k_endog * n))
    call draw_standard_normal(e)
    rep%y = reshape(e, [rep%k_endog, n])
    rep%nobs = n
  end function long_mv

  subroutine test_steady_mv(error)
    type(error_type), allocatable, intent(out) :: error

    call check_steady_pair(error, long_mv(400), 0)
  end subroutine test_steady_mv

  !> Missing observations after convergence end the steady state; the full
  !> recursion resumes and converges again.
  subroutine test_steady_missing(error)
    type(error_type), allocatable, intent(out) :: error
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    integer :: info

    rep = long_mv(400)
    rep%y(:, 250) = ieee_value(1.0_dp, ieee_quiet_nan)
    rep%y(1, 300) = ieee_value(1.0_dp, ieee_quiet_nan)
    call kalman_filter(rep, fres, info)
    call check(error, fres%t_steady < 250, "steady before the gap")
    if (allocated(error)) return
    call check_steady_pair(error, rep, 0)
  end subroutine test_steady_missing

  !> Exact diffuse local linear trend: the shortcut starts only after the
  !> diffuse period.
  subroutine test_steady_diffuse(error)
    type(error_type), allocatable, intent(out) :: error
    type(component_holder_t) :: comps(2)
    type(structural_model_t) :: mod
    type(filter_result_t) :: fres
    real(dp) :: y(1, 400), e(400)
    integer :: info, t

    call draw_standard_normal(e)
    y(1, 1) = e(1)
    do t = 2, 400
      y(1, t) = y(1, t - 1) + 0.1_dp * t / 400 + e(t)
    end do
    allocate (irregular_t :: comps(1)%c)
    allocate (trend_t :: comps(2)%c)
    mod = structural_model(y, comps, info)
    call mod%update([1.0_dp, 0.1_dp, 0.01_dp])
    call kalman_filter(mod%rep, fres, info)
    call check_steady_pair(error, mod%rep, fres%nobs_diffuse)
  end subroutine test_steady_diffuse

  !> Time-varying intercepts c_t and d_t do not affect P, so the shortcut
  !> still applies; its a and v must use each period's c and d.
  subroutine test_steady_cd(error)
    type(error_type), allocatable, intent(out) :: error
    type(ssm_rep_t) :: rep
    real(dp), allocatable :: e(:)
    integer :: n

    rep = long_mv(400)
    n = rep%nobs
    allocate (e(rep%k_endog * n))
    call draw_standard_normal(e)
    rep%d = reshape(e, [rep%k_endog, n])
    deallocate (e)
    allocate (e(rep%k_states * n))
    call draw_standard_normal(e)
    rep%c = 0.1_dp * reshape(e, [rep%k_states, n])
    call check_steady_pair(error, rep, 0)
  end subroutine test_steady_cd
end module test_filter
