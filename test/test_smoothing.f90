!> Other state smoothers and the Whittle relation (DK 4.6), cross-time smoothed
!> state covariances (DK 4.7) and weights (DK 4.8).
module test_smoothing
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use statespace
  use fixture_io, only: fixture_t, load_fixture, rep_from_fixture
  implicit none
  private

  public :: collect_smoothing

  real(dp), parameter :: rtol = 1.0e-9_dp

contains

  subroutine collect_smoothing(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("known_init", test_known), &
                new_unittest("time_varying", test_timevarying), &
                new_unittest("missing", test_missing), &
                new_unittest("stationary_init", test_stationary), &
                new_unittest("exact_diffuse", test_diffuse), &
                new_unittest("exact_diffuse_seasonal_missing", test_diffuse_seasonal), &
                new_unittest("univariate_time_varying", test_uv_timevarying), &
                new_unittest("univariate_missing", test_uv_missing) &
                ]
  end subroutine collect_smoothing

  !> Relative error over the elements where `expected` is not NaN; `actual`
  !> must not be NaN there.
  subroutine check_close(error, actual, expected, label)
    type(error_type), allocatable, intent(out) :: error
    real(dp), intent(in) :: actual(:), expected(:)
    character(len=*), intent(in) :: label
    logical, allocatable :: keep(:)
    real(dp) :: err
    character(len=32) :: buf

    keep = .not. ieee_is_nan(expected)
    call check(error, .not. any(ieee_is_nan(actual) .and. keep), &
               label//": unexpected NaN")
    if (allocated(error)) return
    if (.not. any(keep)) return
    err = maxval(abs(actual - expected), mask=keep) &
          / max(1.0_dp, maxval(abs(expected), mask=keep))
    write (buf, '(es10.3)') err
    call check(error, err <= rtol, label//": relative error "//trim(buf))
  end subroutine check_close

  subroutine check_extras(error, base, univariate)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), intent(in) :: base
    logical, intent(in), optional :: univariate
    type(fixture_t) :: fx, ex
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: acov(:, :, :), cov3(:, :, :), acov_fx(:, :, :), &
                             cov3_fx(:, :, :)
    real(dp), allocatable :: W(:, :, :, :), C(:, :, :, :), A(:, :, :), recon(:, :), &
                             yd(:, :)
    real(dp), allocatable :: c_t(:, :), a1(:), Pstar(:, :), Pinf(:, :), ahat(:, :), &
                             Vhat(:, :, :)
    integer :: info, m, p, n, nd, t, j, ic, it, ir, id, iz
    logical :: weights

    fx = load_fixture("test/fixtures/"//base//".txt")
    ex = load_fixture("test/fixtures/extras_"//base//".txt")
    rep = rep_from_fixture(fx)
    weights = .true.
    if (present(univariate)) then
      if (univariate) then
        rep%filter_method = FILTER_UNIVARIATE
        weights = .false.     ! weights do not depend on the filter method
      end if
    end if
    m = rep%k_states; p = rep%k_endog; n = rep%nobs
    call kalman_filter(rep, fres, info)
    call state_smoother(rep, fres, sres, info)
    call check(error, info, SS_OK, "smoother info")
    if (allocated(error)) return
    nd = fres%nobs_diffuse

    ! Other smoothing algorithms (DK 4.6) must reproduce the state smoother.
    allocate (ahat(m, n), Vhat(m, m, n))
    call fast_state_smoother(rep, fres, ahat, info)
    call check(error, info, SS_OK, "fast smoother info")
    if (allocated(error)) return
    call check_close(error, pack(ahat, .true.), pack(sres%alphahat, .true.), &
                     "fast smoother")
    if (allocated(error)) return
    call classical_state_smoother(rep, fres, ahat, Vhat, info)
    if (nd > 0) then
      call check(error, info, SS_ERR_UNSUPPORTED, &
                 "classical smoother with diffuse states")
      if (allocated(error)) return
    else
      call check(error, info, SS_OK, "classical smoother info")
      if (allocated(error)) return
      call check_close(error, pack(ahat, .true.), pack(sres%alphahat, .true.), &
                       "classical smoother alphahat")
      if (allocated(error)) return
      call check_close(error, pack(Vhat, .true.), pack(sres%V, .true.), &
                       "classical smoother V")
      if (allocated(error)) return
    end if

    ! Two-filter formula (DK 4.6.4), without diffuse states.
    call two_filter_smoother(rep, fres, ahat, Vhat, info)
    if (nd > 0) then
      call check(error, info, SS_ERR_UNSUPPORTED, &
                 "two-filter smoother with diffuse states")
      if (allocated(error)) return
    else
      call check(error, info, SS_OK, "two-filter smoother info")
      if (allocated(error)) return
      call check_close(error, pack(ahat, .true.), pack(sres%alphahat, .true.), &
                       "two-filter alphahat")
      if (allocated(error)) return
      call check_close(error, pack(Vhat, .true.), pack(sres%V, .true.), "two-filter V")
      if (allocated(error)) return
    end if

    ! Whittle relation (DK 4.6.3): the smoothed estimates satisfy the model
    ! equations with the disturbances replaced by their smoothed values,
    !   alphahat_t+1 = c_t + T_t alphahat_t + R_t etahat_t
    !   y_t = d_t + Z_t alphahat_t + epshat_t   (observed elements)
    do t = 1, n
      ic = tidx(size(rep%c, 2), t)
      it = tidx(size(rep%T, 3), t)
      ir = tidx(size(rep%R, 3), t)
      id = tidx(size(rep%d, 2), t)
      iz = tidx(size(rep%Z, 3), t)
      if (t < n) then
        call check_close(error, sres%alphahat(:, t + 1), &
                         rep%c(:, ic) + matmul(rep%T(:, :, it), sres%alphahat(:, t)) &
                         + matmul(rep%R(:, :, ir), sres%etahat(:, t)), &
                         "Whittle: state equation")
        if (allocated(error)) return
      end if
      call check_close(error, pack(rep%y(:, t), .not. ieee_is_nan(rep%y(:, t))), &
                       pack(rep%d(:, id) + matmul(rep%Z(:, :, iz), &
                                                  sres%alphahat(:, t)) &
                            + sres%epshat(:, t), .not. ieee_is_nan(rep%y(:, t))), &
                       "Whittle: observation equation")
      if (allocated(error)) return
    end do

    ! Cross-time covariances outside the diffuse period.
    allocate (acov(m, m, n - 1), cov3(m, m, n - 3))
    call smoothed_state_autocov(rep, fres, sres, acov, info)
    call check(error, info, SS_OK, "autocov info")
    if (allocated(error)) return
    acov_fx = ex%get3('autocov')
    acov_fx(:, :, 1:nd) = acov(:, :, 1:nd)   ! not defined in the diffuse period
    call check(error, all(ieee_is_nan(acov(:, :, 1:nd))), &
               "autocov NaN in diffuse period")
    if (allocated(error)) return
    call check_close(error, pack(acov(:, :, nd + 1:), .true.), &
                     pack(acov_fx(:, :, nd + 1:), .true.), "autocov")
    if (allocated(error)) return
    cov3_fx = ex%get3('cov_shift3')
    do t = nd + 1, n - 3
      call smoothed_state_cov_between(rep, fres, sres, t, t + 3, cov3(:, :, t), info)
      call check(error, info, SS_OK, "cov_between info")
      if (allocated(error)) return
    end do
    call check_close(error, pack(cov3(:, :, nd + 1:), .true.), &
                     pack(cov3_fx(:, :, nd + 1:), .true.), "Cov(alpha_t, alpha_t+3)")
    if (allocated(error)) return
    if (.not. weights) return

    ! Weights: statsmodels leaves diffuse periods as NaN.
    allocate (W(m, p, n, n), C(m, m, n, n), A(m, m, n))
    call smoothed_state_weights(rep, W, C, A, info)
    call check(error, info, SS_OK, "weights info")
    if (allocated(error)) return
    call check_close(error, pack(W, .true.), ex%get1('weights'), "weights")
    if (allocated(error)) return
    call check_close(error, pack(C, .true.), ex%get1('c_weights'), "c_weights")
    if (allocated(error)) return
    call check_close(error, pack(A, .true.), ex%get1('prior_weights'), "prior_weights")
    if (allocated(error)) return

    ! alphahat_t = sum_j W_tj (y_j - d_j) + sum_j C_tj c_j + A_t a1 with a1 the
    ! initial state mean, in every period including the diffuse ones.
    allocate (recon(m, n), source=0.0_dp)
    allocate (yd(p, n), c_t(m, n))
    do j = 1, n
      yd(:, j) = merge(0.0_dp, rep%y(:, j) - rep%d(:, tidx(size(rep%d, 2), j)), &
                       ieee_is_nan(rep%y(:, j)))
      c_t(:, j) = rep%c(:, tidx(size(rep%c, 2), j))
    end do
    allocate (a1(m), Pstar(m, m), Pinf(m, m))
    call rep%initial_state(a1, Pstar, Pinf, info)
    do t = 1, n
      do j = 1, n
        recon(:, t) = recon(:, t) + matmul(W(:, :, t, j), yd(:, j)) &
                      + matmul(C(:, :, t, j), c_t(:, j))
      end do
      recon(:, t) = recon(:, t) + matmul(A(:, :, t), a1)
    end do
    call check_close(error, pack(recon, .true.), pack(sres%alphahat, .true.), &
                     "alphahat from weights")
  end subroutine check_extras

  subroutine test_known(error)
    type(error_type), allocatable, intent(out) :: error

    call check_extras(error, "mv_invariant")
  end subroutine test_known

  subroutine test_timevarying(error)
    type(error_type), allocatable, intent(out) :: error

    call check_extras(error, "mv_timevarying")
  end subroutine test_timevarying

  subroutine test_missing(error)
    type(error_type), allocatable, intent(out) :: error

    call check_extras(error, "mv_missing")
  end subroutine test_missing

  subroutine test_stationary(error)
    type(error_type), allocatable, intent(out) :: error

    call check_extras(error, "mv_stationary")
  end subroutine test_stationary

  subroutine test_diffuse(error)
    type(error_type), allocatable, intent(out) :: error

    call check_extras(error, "mv_diffuse")
  end subroutine test_diffuse

  subroutine test_diffuse_seasonal(error)
    type(error_type), allocatable, intent(out) :: error

    call check_extras(error, "uc_trend_seasonal_exact")
  end subroutine test_diffuse_seasonal

  ! L_t from the per-element gains of the univariate filter.
  subroutine test_uv_timevarying(error)
    type(error_type), allocatable, intent(out) :: error

    call check_extras(error, "mv_timevarying", univariate=.true.)
  end subroutine test_uv_timevarying

  subroutine test_uv_missing(error)
    type(error_type), allocatable, intent(out) :: error

    call check_extras(error, "mv_missing", univariate=.true.)
  end subroutine test_uv_missing
end module test_smoothing
