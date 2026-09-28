!> The C interface (statespace_capi), called from Fortran as a C caller
!> would: handles, arrays through pointers, status codes.
module test_capi
  use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_null_ptr, c_int, c_double, &
                                         c_char, c_null_char
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use statespace
  use statespace_capi
  use fixture_io, only: fixture_t, load_fixture, rep_from_fixture
  implicit none
  private

  public :: collect_capi

contains

  subroutine collect_capi(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("capi_matches_fortran", test_capi_mv), &
                new_unittest("capi_errors", test_capi_errors) &
                ]
  end subroutine collect_capi

  !> mv_invariant built through the C interface: loglike, filter and
  !> smoother output equal the Fortran calls on the same model.
  subroutine test_capi_mv(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    type(c_ptr) :: hr, fh, sh
    real(dp), allocatable, target :: y(:, :), Z(:, :, :), H(:, :, :), T(:, :, :), &
                                     R(:, :, :), Q(:, :, :), c(:, :), d(:, :), a1(:), &
                                     P1(:, :), a(:, :), alphahat(:, :)
    real(c_double) :: llf, llf_f
    integer(c_int) :: nd, kd, ts
    integer :: info, kp, km, kr, kn

    fx = load_fixture("test/fixtures/mv_invariant.txt")
    rep = rep_from_fixture(fx)
    kp = rep%k_endog; km = rep%k_states; kr = rep%k_posdef; kn = rep%nobs
    y = rep%y; Z = rep%Z; H = rep%H; T = rep%T
    R = rep%R; Q = rep%Q; c = rep%c; d = rep%d
    a1 = fx%get1('a1'); P1 = fx%get2('P1')

    call check(error, ss_rep_new(kp, km, kr, kn, hr), SS_OK, "new")
    if (allocated(error)) return
    call check(error, ss_rep_set(hr, SS_ARR_Y, kn, c_loc(y)) == SS_OK .and. &
               ss_rep_set(hr, SS_ARR_Z, size(Z, 3), c_loc(Z)) == SS_OK .and. &
               ss_rep_set(hr, SS_ARR_H, size(H, 3), c_loc(H)) == SS_OK .and. &
               ss_rep_set(hr, SS_ARR_T, size(T, 3), c_loc(T)) == SS_OK .and. &
               ss_rep_set(hr, SS_ARR_R, size(R, 3), c_loc(R)) == SS_OK .and. &
               ss_rep_set(hr, SS_ARR_Q, size(Q, 3), c_loc(Q)) == SS_OK .and. &
               ss_rep_set(hr, SS_ARR_C, size(c, 2), c_loc(c)) == SS_OK .and. &
               ss_rep_set(hr, SS_ARR_D, size(d, 2), c_loc(d)) == SS_OK .and. &
               ss_rep_init_known(hr, c_loc(a1), c_loc(P1)) == SS_OK, "set")
    if (allocated(error)) return

    call check(error, ss_loglike(hr, llf), SS_OK, "loglike")
    if (allocated(error)) return
    call check(error, abs(llf - loglike(rep, info)) <= 1.0e-12_dp * abs(llf), &
               "same loglike")
    if (allocated(error)) return

    call kalman_filter(rep, fres, info)
    call state_smoother(rep, fres, sres, info)
    call check(error, ss_filter(hr, fh), SS_OK, "filter")
    if (allocated(error)) return
    call check(error, ss_filter_scalars(fh, llf_f, nd, kd, ts), SS_OK, "scalars")
    if (allocated(error)) return
    allocate (a(km, kn + 1), alphahat(km, kn))
    call check(error, ss_filter_get(fh, SS_F_A, c_loc(a)), SS_OK, "get a")
    if (allocated(error)) return
    call check(error, all(a == fres%a) .and. llf_f == fres%llf &
               .and. ts == fres%t_steady, "same filter output")
    if (allocated(error)) return
    call check(error, ss_smooth(hr, fh, sh), SS_OK, "smooth")
    if (allocated(error)) return
    call check(error, ss_smoother_get(sh, SS_S_ALPHAHAT, c_loc(alphahat)), SS_OK, &
               "get alphahat")
    if (allocated(error)) return
    call check(error, all(alphahat == sres%alphahat), "same smoother output")
    if (allocated(error)) return
    call check(error, ss_smoother_free(sh) == SS_OK .and. ss_filter_free(fh) == SS_OK &
               .and. ss_rep_free(hr) == SS_OK, "free")
  end subroutine test_capi_mv

  !> Bad handles, shapes and codes return errors instead of crashing.
  subroutine test_capi_errors(error)
    type(error_type), allocatable, intent(out) :: error
    type(c_ptr) :: hr, fh
    real(c_double) :: llf
    real(dp), target :: x(2, 2, 3)
    character(kind=c_char) :: buf(16)

    call check(error, ss_rep_new(0, 1, 1, 5, hr), SS_ERR_DIM, "p = 0")
    if (allocated(error)) return
    call check(error, ss_loglike(c_null_ptr, llf), SS_ERR_DIM, "null handle")
    if (allocated(error)) return
    call check(error, ss_rep_new(2, 2, 2, 5, hr), SS_OK, "new")
    if (allocated(error)) return
    x = 0.0_dp
    call check(error, ss_rep_set(hr, SS_ARR_T, 3, c_loc(x)), SS_ERR_DIM, &
               "nt neither 1 nor n")
    if (allocated(error)) return
    call check(error, ss_rep_set(hr, 99_c_int, 1, c_loc(x)), SS_ERR_DIM, "bad code")
    if (allocated(error)) return
    call check(error, ss_rep_init_stationary(hr), SS_OK, "init")
    if (allocated(error)) return
    ! Z = 0, H = 0: F is singular
    call check(error, ss_filter(hr, fh) /= SS_OK, "filter error reported")
    if (allocated(error)) return
    call check(error, ss_rep_free(hr), SS_OK, "free")
    if (allocated(error)) return
    call check(error, ss_version(buf, 16), SS_OK, "version")
    if (allocated(error)) return
    call check(error, buf(1) == "0" .and. any(buf == c_null_char), "version string")
  end subroutine test_capi_errors
end module test_capi
