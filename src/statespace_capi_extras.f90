!> C interface to the remaining algorithms: the smoothers and weights of DK
!> ch. 4, the square root and augmented filters, EM, collapsing, linear
!> restrictions, residual diagnostics, de Jong-Shephard simulation and the
!> effect of parameter estimation. Conventions as in statespace_capi; periods
!> passed as arguments are 1-based.
module statespace_capi_extras
  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_loc, &
                                         c_f_pointer, c_associated
  use statespace_kinds, only: dp, SS_OK, SS_ERR_DIM
  use statespace_rep, only: ssm_rep_t
  use statespace_filter, only: filter_result_t
  use statespace_smoother, only: smoother_result_t
  use statespace_capi, only: get_rep, get_filter, get_smoother, copy_out, filter_box, &
                             smoother_box
  use statespace_capi_models, only: model_box, get_model
  use statespace_smoothing, only: fast_state_smoother, classical_state_smoother, &
                                  two_filter_smoother, whittle_smoother, &
                                  fixed_point_smoother, fixed_lag_smoother, &
                                  update_smoothed, smoothed_state_autocov, &
                                  smoothed_state_cov_between, filtered_state_weights, &
                                  smoothed_state_weights, innovation_transition
  use statespace_sqrt, only: sqrt_kalman_filter, sqrt_state_smoother
  use statespace_augmented, only: augmented_result_t, augmented_filter, &
                                  augmented_smoother
  use statespace_em, only: em_variances
  use statespace_collapse, only: collapse_observations
  use statespace_restrict, only: add_state_restrictions
  use statespace_diagnostics, only: auxiliary_residuals, auxiliary_residuals_vector, &
                                    de_jong_penzer, least_squares_residuals, r2_diffuse
  use statespace_simsmooth, only: djs_measurement_disturbances, djs_state_disturbances
  use statespace_mle, only: fit_result_t, estimation_bias
  implicit none
  private

  type :: augmented_box
    type(augmented_result_t) :: res
  end type augmented_box

  public :: ss_fast_smoother, ss_classical_smoother, ss_two_filter_smoother, &
            ss_whittle_smoother
  public :: ss_fixed_point_smoother, ss_fixed_lag_smoother, ss_update_smoothed
  public :: ss_smoothed_state_autocov, ss_smoothed_state_cov_between, &
            ss_filtered_state_weights
  public :: ss_smoothed_state_weights, ss_innovation_transition
  public :: ss_sqrt_filter, ss_sqrt_smoother
  public :: ss_augmented_filter, ss_augmented_free, ss_augmented_info, ss_augmented_get
  public :: ss_augmented_filter_result, ss_augmented_smoother
  public :: ss_em, ss_collapse, ss_add_restrictions
  public :: ss_auxiliary_residuals, ss_de_jong_penzer, ss_least_squares_residuals, &
            ss_r2_diffuse
  public :: ss_djs_simulation_smoother, ss_estimation_bias

contains

  !> The representation and filter result behind two handles, or null.
  subroutine rep_filter(handle, fhandle, b, f)
    type(c_ptr), intent(in) :: handle, fhandle
    type(ssm_rep_t), pointer, intent(out) :: b
    type(filter_box), pointer, intent(out) :: f

    b => get_rep(handle)
    f => get_filter(fhandle)
    if (.not. associated(f)) b => null()
  end subroutine rep_filter

  function m2(p, n1, n2) result(x)
    type(c_ptr), intent(in) :: p
    integer, intent(in) :: n1, n2
    real(dp), pointer :: x(:, :)

    call c_f_pointer(p, x, [n1, n2])
  end function m2

  function m3(p, n1, n2, n3) result(x)
    type(c_ptr), intent(in) :: p
    integer, intent(in) :: n1, n2, n3
    real(dp), pointer :: x(:, :, :)

    call c_f_pointer(p, x, [n1, n2, n3])
  end function m3

  function m4(p, n1, n2, n3, n4) result(x)
    type(c_ptr), intent(in) :: p
    integer, intent(in) :: n1, n2, n3, n4
    real(dp), pointer :: x(:, :, :, :)

    call c_f_pointer(p, x, [n1, n2, n3, n4])
  end function m4

  ! ------------------------------------------------------- DK 4.4 - 4.8

  !> Fast state smoother (DK 4.6.2): alphahat (m, n).
  integer(c_int) function ss_fast_smoother(handle, fhandle, alphahat) &
    bind(C, name="ss_fast_smoother") result(info)
    type(c_ptr), value :: handle, fhandle, alphahat
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    integer :: stat

    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    if (.not. associated(b)) return
    call fast_state_smoother(b, f%res, m2(alphahat, b%k_states, b%nobs), stat)
    info = stat
  end function ss_fast_smoother

  !> Classical (RTS) smoother (DK 4.6.1): alphahat (m, n), V (m, m, n).
  integer(c_int) function ss_classical_smoother(handle, fhandle, alphahat, V) &
    bind(C, name="ss_classical_smoother") result(info)
    type(c_ptr), value :: handle, fhandle, alphahat, V
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    integer :: stat

    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    if (.not. associated(b)) return
    call classical_state_smoother(b, f%res, m2(alphahat, b%k_states, b%nobs), &
                                  m3(V, b%k_states, b%k_states, b%nobs), stat)
    info = stat
  end function ss_classical_smoother

  !> Two-filter smoother (DK 4.6.4): alphahat (m, n), V (m, m, n).
  integer(c_int) function ss_two_filter_smoother(handle, fhandle, alphahat, V) &
    bind(C, name="ss_two_filter_smoother") result(info)
    type(c_ptr), value :: handle, fhandle, alphahat, V
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    integer :: stat

    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    if (.not. associated(b)) return
    call two_filter_smoother(b, f%res, m2(alphahat, b%k_states, b%nobs), &
                             m3(V, b%k_states, b%k_states, b%nobs), stat)
    info = stat
  end function ss_two_filter_smoother

  !> Smoothed state by the Whittle relation (DK 4.6.3): alphahat (m, n).
  integer(c_int) function ss_whittle_smoother(handle, fhandle, alphahat) &
    bind(C, name="ss_whittle_smoother") result(info)
    type(c_ptr), value :: handle, fhandle, alphahat
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    integer :: stat

    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    if (.not. associated(b)) return
    call whittle_smoother(b, f%res, m2(alphahat, b%k_states, b%nobs), stat)
    info = stat
  end function ss_whittle_smoother

  !> Fixed-point smoother (DK 4.4.6) for period t: path_a (m, n-t+1) and
  !> path_V (m, m, n-t+1), the estimates of alpha_t given y_1..y_k,
  !> k = t..n.
  integer(c_int) function ss_fixed_point_smoother(handle, fhandle, t, path_a, path_V) &
    bind(C, name="ss_fixed_point_smoother") result(info)
    type(c_ptr), value :: handle, fhandle, path_a, path_V
    integer(c_int), value :: t
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    integer :: stat, k

    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    if (.not. associated(b) .or. t < 1 .or. t > b%nobs) return
    k = b%nobs - t + 1
    call fixed_point_smoother(b, f%res, int(t), m2(path_a, b%k_states, k), &
                              m3(path_V, b%k_states, b%k_states, k), stat)
    info = stat
  end function ss_fixed_point_smoother

  !> Fixed-lag smoother (DK 4.4.6) with lag j: alphahat (m, n) and
  !> V (m, m, n), NaN in the last j columns.
  integer(c_int) function ss_fixed_lag_smoother(handle, fhandle, j, alphahat, V) &
    bind(C, name="ss_fixed_lag_smoother") result(info)
    type(c_ptr), value :: handle, fhandle, alphahat, V
    integer(c_int), value :: j
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    integer :: stat

    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    if (.not. associated(b) .or. j < 0) return
    call fixed_lag_smoother(b, f%res, int(j), m2(alphahat, b%k_states, b%nobs), &
                            m3(V, b%k_states, b%k_states, b%nobs), stat)
    info = stat
  end function ss_fixed_lag_smoother

  !> Update smoothed estimates for new observations (DK 4.4.5): on entry
  !> alphahat (m, n) and V (m, m, n) hold, in columns 1..n0, the estimates
  !> given y_1..y_n0; on exit, all columns given y_1..y_n.
  integer(c_int) function ss_update_smoothed(handle, fhandle, n0, alphahat, V) &
    bind(C, name="ss_update_smoothed") result(info)
    type(c_ptr), value :: handle, fhandle, alphahat, V
    integer(c_int), value :: n0
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    integer :: stat

    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    if (.not. associated(b) .or. n0 < 1 .or. n0 > b%nobs) return
    call update_smoothed(b, f%res, int(n0), m2(alphahat, b%k_states, b%nobs), &
                         m3(V, b%k_states, b%k_states, b%nobs), stat)
    info = stat
  end function ss_update_smoothed

  !> acov (m, m, n-1): Cov(alpha_t+1, alpha_t | Y_n) (DK 4.7).
  integer(c_int) function ss_smoothed_state_autocov(handle, fhandle, shandle, acov) &
    bind(C, name="ss_smoothed_state_autocov") result(info)
    type(c_ptr), value :: handle, fhandle, shandle, acov
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    type(smoother_box), pointer :: s
    integer :: stat

    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    s => get_smoother(shandle)
    if (.not. associated(b) .or. .not. associated(s) .or. b%nobs < 2) return
    call smoothed_state_autocov(b, f%res, s%res, &
                                m3(acov, b%k_states, b%k_states, b%nobs - 1), stat)
    info = stat
  end function ss_smoothed_state_autocov

  !> C (m, m) = Cov(alpha_t, alpha_j | Y_n) (DK 4.7).
  integer(c_int) function ss_smoothed_state_cov_between(handle, fhandle, shandle, t, &
                                                        j, C) &
    bind(C, name="ss_smoothed_state_cov_between") result(info)
    type(c_ptr), value :: handle, fhandle, shandle, C
    integer(c_int), value :: t, j
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    type(smoother_box), pointer :: s
    integer :: stat

    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    s => get_smoother(shandle)
    if (.not. associated(b) .or. .not. associated(s)) return
    if (min(t, j) < 1 .or. max(t, j) > b%nobs) return
    call smoothed_state_cov_between(b, f%res, s%res, int(t), int(j), &
                                    m2(C, b%k_states, b%k_states), stat)
    info = stat
  end function ss_smoothed_state_cov_between

  !> Filtering weights (DK 4.8.2): Wa, Watt (m, p, n, n).
  integer(c_int) function ss_filtered_state_weights(handle, fhandle, Wa, Watt) &
    bind(C, name="ss_filtered_state_weights") result(info)
    type(c_ptr), value :: handle, fhandle, Wa, Watt
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    integer :: stat, m, p, n

    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    if (.not. associated(b)) return
    m = b%k_states; p = b%k_endog; n = b%nobs
    call filtered_state_weights(b, f%res, m4(Wa, m, p, n, n), m4(Watt, m, p, n, n), &
                                stat)
    info = stat
  end function ss_filtered_state_weights

  !> Smoothing weights (DK 4.8.3): W (m, p, n, n) on y - d, C (m, m, n, n)
  !> on c, A (m, m, n) on a1.
  integer(c_int) function ss_smoothed_state_weights(handle, W, C, A) &
    bind(C, name="ss_smoothed_state_weights") result(info)
    type(c_ptr), value :: handle, W, C, A
    type(ssm_rep_t), pointer :: b
    integer :: stat, m, p, n

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    m = b%k_states; p = b%k_endog; n = b%nobs
    call smoothed_state_weights(b, m4(W, m, p, n, n), m4(C, m, m, n, n), &
                                m3(A, m, m, n), stat)
    info = stat
  end function ss_smoothed_state_weights

  !> L_t (m, m), the transition of the state prediction error (DK 4.3).
  integer(c_int) function ss_innovation_transition(handle, fhandle, t, L) &
    bind(C, name="ss_innovation_transition") result(info)
    type(c_ptr), value :: handle, fhandle, L
    integer(c_int), value :: t
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    integer :: stat

    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    if (.not. associated(b) .or. t < 1 .or. t > b%nobs) return
    call innovation_transition(b, f%res, int(t), m2(L, b%k_states, b%k_states), stat)
    info = stat
  end function ss_innovation_transition

  ! --------------------------------------------------------- DK 6.3, 5.7

  !> Square root filter (DK 6.3); fhandle as from ss_filter.
  integer(c_int) function ss_sqrt_filter(handle, fhandle) &
    bind(C, name="ss_sqrt_filter") result(info)
    type(c_ptr), value :: handle
    type(c_ptr), intent(out) :: fhandle
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    integer :: stat

    fhandle = c_null_ptr
    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    allocate (f)
    call sqrt_kalman_filter(b, f%res, stat)
    info = stat
    if (stat /= SS_OK) then
      deallocate (f)
      return
    end if
    fhandle = c_loc(f)
  end function ss_sqrt_filter

  !> Square root smoother (DK 6.3); shandle as from ss_smooth.
  integer(c_int) function ss_sqrt_smoother(handle, fhandle, shandle) &
    bind(C, name="ss_sqrt_smoother") result(info)
    type(c_ptr), value :: handle, fhandle
    type(c_ptr), intent(out) :: shandle
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    type(smoother_box), pointer :: s
    integer :: stat

    shandle = c_null_ptr
    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    if (.not. associated(b)) return
    allocate (s)
    call sqrt_state_smoother(b, f%res, s%res, stat)
    info = stat
    if (stat /= SS_OK) then
      deallocate (s)
      return
    end if
    shandle = c_loc(s)
  end function ss_sqrt_smoother

  !> Augmented filter (DK 5.7); ahandle is freed with ss_augmented_free.
  integer(c_int) function ss_augmented_filter(handle, ahandle) &
    bind(C, name="ss_augmented_filter") result(info)
    type(c_ptr), value :: handle
    type(c_ptr), intent(out) :: ahandle
    type(ssm_rep_t), pointer :: b
    type(augmented_box), pointer :: a
    integer :: stat

    ahandle = c_null_ptr
    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    allocate (a)
    call augmented_filter(b, a%res, stat)
    info = stat
    if (stat /= SS_OK) then
      deallocate (a)
      return
    end if
    ahandle = c_loc(a)
  end function ss_augmented_filter

  function get_aug(handle) result(a)
    type(c_ptr), intent(in) :: handle
    type(augmented_box), pointer :: a

    a => null()
    if (c_associated(handle)) call c_f_pointer(handle, a)
  end function get_aug

  integer(c_int) function ss_augmented_free(ahandle) &
    bind(C, name="ss_augmented_free") result(info)
    type(c_ptr), value :: ahandle
    type(augmented_box), pointer :: a

    a => get_aug(ahandle)
    if (associated(a)) deallocate (a)
    info = SS_OK
  end function ss_augmented_free

  !> k (number of diffuse elements), the diffuse, fixed-but-unknown and
  !> marginal log likelihoods (DK 7.2.2-7.2.6).
  integer(c_int) function ss_augmented_info(ahandle, k, llf, llf_fixed, llf_marginal) &
    bind(C, name="ss_augmented_info") result(info)
    type(c_ptr), value :: ahandle
    integer(c_int), intent(out) :: k
    real(c_double), intent(out) :: llf, llf_fixed, llf_marginal
    type(augmented_box), pointer :: a

    info = SS_ERR_DIM
    a => get_aug(ahandle)
    if (.not. associated(a)) return
    k = a%res%k_diffuse
    llf = a%res%llf; llf_fixed = a%res%llf_fixed; llf_marginal = a%res%llf_marginal
    info = SS_OK
  end function ss_augmented_info

  !> delta (1, k) or delta_cov (2, k x k).
  integer(c_int) function ss_augmented_get(ahandle, code, out) &
    bind(C, name="ss_augmented_get") result(info)
    type(c_ptr), value :: ahandle, out
    integer(c_int), value :: code
    type(augmented_box), pointer :: a

    info = SS_ERR_DIM
    a => get_aug(ahandle)
    if (.not. associated(a) .or. .not. c_associated(out)) return
    select case (code)
    case (1)
      call copy_out(a%res%delta, out)
    case (2)
      call copy_out(a%res%delta_cov, out)
    case default
      return
    end select
    info = SS_OK
  end function ss_augmented_get

  !> A borrowed filter handle to the delta = 0 filter (valid while ahandle
  !> exists; do not free).
  integer(c_int) function ss_augmented_filter_result(ahandle, fhandle) &
    bind(C, name="ss_augmented_filter_result") result(info)
    type(c_ptr), value :: ahandle
    type(c_ptr), intent(out) :: fhandle
    type(augmented_box), pointer :: a

    fhandle = c_null_ptr
    info = SS_ERR_DIM
    a => get_aug(ahandle)
    if (.not. associated(a)) return
    fhandle = c_loc(a%res%filter)
    info = SS_OK
  end function ss_augmented_filter_result

  !> Smoothed state (m, n) and variance (m, m, n) from the augmented filter.
  integer(c_int) function ss_augmented_smoother(handle, ahandle, alphahat, V) &
    bind(C, name="ss_augmented_smoother") result(info)
    type(c_ptr), value :: handle, ahandle, alphahat, V
    type(ssm_rep_t), pointer :: b
    type(augmented_box), pointer :: a
    integer :: stat

    info = SS_ERR_DIM
    b => get_rep(handle)
    a => get_aug(ahandle)
    if (.not. associated(b) .or. .not. associated(a)) return
    call augmented_smoother(b, a%res, m2(alphahat, b%k_states, b%nobs), &
                            m3(V, b%k_states, b%k_states, b%nobs), stat)
    info = stat
  end function ss_augmented_smoother

  ! ------------------------------------------------ EM, collapse, restrict

  !> EM for H and Q (DK 7.3.4), updating the representation in place. path
  !> (maxiter) receives the log likelihood of every iteration, or NULL.
  integer(c_int) function ss_em(handle, maxiter, tol, diagonal_H, diagonal_Q, llf, &
                                niter, path) &
    bind(C, name="ss_em") result(info)
    type(c_ptr), value :: handle, path
    integer(c_int), value :: maxiter, diagonal_H, diagonal_Q
    real(c_double), value :: tol
    real(c_double), intent(out) :: llf
    integer(c_int), intent(out) :: niter
    type(ssm_rep_t), pointer :: b
    real(dp), allocatable :: lp(:)
    real(dp), pointer :: po(:)
    integer :: stat, it

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b) .or. maxiter < 1) return
    call em_variances(b, int(maxiter), real(tol, dp), llf, it, stat, &
                      diagonal_H=diagonal_H /= 0, diagonal_Q=diagonal_Q /= 0, &
                      llf_path=lp)
    niter = it
    info = stat
    if (c_associated(path) .and. allocated(lp)) then
      call c_f_pointer(path, po, [maxiter])
      po = ieee_nan()
      po(1:min(size(lp), int(maxiter))) = lp(1:min(size(lp), int(maxiter)))
    end if
  end function ss_em

  real(dp) function ieee_nan()
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    ieee_nan = ieee_value(1.0_dp, ieee_quiet_nan)
  end function ieee_nan

  !> Collapse the observations to the state dimension (DK 6.5): a new
  !> representation crhandle and the log likelihood adjustment (n).
  integer(c_int) function ss_collapse(handle, crhandle, llf_adjust) &
    bind(C, name="ss_collapse") result(info)
    type(c_ptr), value :: handle, llf_adjust
    type(c_ptr), intent(out) :: crhandle
    type(ssm_rep_t), pointer :: b, cr
    real(dp), pointer :: adj(:)
    integer :: stat

    crhandle = c_null_ptr
    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    call c_f_pointer(llf_adjust, adj, [b%nobs])
    allocate (cr)
    call collapse_observations(b, cr, adj, stat)
    info = stat
    if (stat /= SS_OK) then
      deallocate (cr)
      return
    end if
    crhandle = c_loc(cr)
  end function ss_collapse

  !> Linear restrictions R*_t alpha_t = r*_t (DK 6.6): Rmat (q, m, nt),
  !> rval (q, n) with NaN where inactive; a new representation rrhandle.
  integer(c_int) function ss_add_restrictions(handle, q, nt, Rmat, rval, rrhandle) &
    bind(C, name="ss_add_restrictions") result(info)
    type(c_ptr), value :: handle, Rmat, rval
    integer(c_int), value :: q, nt
    type(c_ptr), intent(out) :: rrhandle
    type(ssm_rep_t), pointer :: b, rr
    integer :: stat

    rrhandle = c_null_ptr
    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b) .or. q < 1 .or. (nt /= 1 .and. nt /= b%nobs)) return
    allocate (rr)
    call add_state_restrictions(b, m3(Rmat, int(q), b%k_states, int(nt)), &
                                m2(rval, int(q), b%nobs), rr, stat)
    info = stat
    if (stat /= SS_OK) then
      deallocate (rr)
      return
    end if
    rrhandle = c_loc(rr)
  end function ss_add_restrictions

  ! --------------------------------------------------------- diagnostics

  !> Auxiliary residuals (DK 7.5), element-wise (vector = 0) or standardized
  !> as vectors: eps_std (p, n), eta_std (r, n).
  integer(c_int) function ss_auxiliary_residuals(handle, shandle, vector, eps_std, &
                                                 eta_std) &
    bind(C, name="ss_auxiliary_residuals") result(info)
    type(c_ptr), value :: handle, shandle, eps_std, eta_std
    integer(c_int), value :: vector
    type(ssm_rep_t), pointer :: b
    type(smoother_box), pointer :: s

    info = SS_ERR_DIM
    b => get_rep(handle)
    s => get_smoother(shandle)
    if (.not. associated(b) .or. .not. associated(s)) return
    if (vector /= 0) then
      call auxiliary_residuals_vector(b, s%res, m2(eps_std, b%k_endog, b%nobs), &
                                      m2(eta_std, b%k_posdef, b%nobs))
    else
      call auxiliary_residuals(b, s%res, m2(eps_std, b%k_endog, b%nobs), &
                               m2(eta_std, b%k_posdef, b%nobs))
    end if
    info = SS_OK
  end function ss_auxiliary_residuals

  !> de Jong-Penzer statistics (DK 7.5): r_stat (m, n), e_stat (p, n).
  integer(c_int) function ss_de_jong_penzer(handle, fhandle, shandle, r_stat, e_stat) &
    bind(C, name="ss_de_jong_penzer") result(info)
    type(c_ptr), value :: handle, fhandle, shandle, r_stat, e_stat
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    type(smoother_box), pointer :: s

    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    s => get_smoother(shandle)
    if (.not. associated(b) .or. .not. associated(s)) return
    call de_jong_penzer(b, f%res, s%res, m2(r_stat, b%k_states, b%nobs), &
                        m2(e_stat, b%k_endog, b%nobs))
    info = SS_OK
  end function ss_de_jong_penzer

  !> Least squares residuals (DK 6.2.4) for regression coefficients in
  !> states first..last: vplus (p, n).
  integer(c_int) function ss_least_squares_residuals(handle, first, last, vplus) &
    bind(C, name="ss_least_squares_residuals") result(info)
    type(c_ptr), value :: handle, vplus
    integer(c_int), value :: first, last
    type(ssm_rep_t), pointer :: b
    integer :: stat

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b) .or. first < 1 .or. last < first &
        .or. last > b%k_states) return
    call least_squares_residuals(b, int(first), int(last), &
                                 m2(vplus, b%k_endog, b%nobs), stat)
    info = stat
  end function ss_least_squares_residuals

  !> R^2_D against a random walk with drift (Harvey 1989) for series i.
  integer(c_int) function ss_r2_diffuse(handle, fhandle, i, r2) &
    bind(C, name="ss_r2_diffuse") result(info)
    type(c_ptr), value :: handle, fhandle
    integer(c_int), value :: i
    real(c_double), intent(out) :: r2
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f

    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    if (.not. associated(b) .or. i < 1 .or. i > b%k_endog) return
    r2 = r2_diffuse(b, f%res, int(i))
    info = SS_OK
  end function ss_r2_diffuse

  ! ------------------------------------------------ simulation, DK 7.3.7

  !> de Jong-Shephard simulation smoother (DK 4.9.3) from standard normal
  !> variates u_eps (p, n), u_eta (r, n), u_init (m): eps (p, n),
  !> eta (r, n), alpha (m, n).
  integer(c_int) function ss_djs_simulation_smoother(handle, fhandle, u_eps, u_eta, &
                                                     u_init, eps, eta, alpha) &
    bind(C, name="ss_djs_simulation_smoother") result(info)
    type(c_ptr), value :: handle, fhandle, u_eps, u_eta, u_init, eps, eta, alpha
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    real(dp), pointer :: ui(:)
    integer :: stat

    info = SS_ERR_DIM
    call rep_filter(handle, fhandle, b, f)
    if (.not. associated(b)) return
    call djs_measurement_disturbances(b, f%res, m2(u_eps, b%k_endog, b%nobs), &
                                      m2(eps, b%k_endog, b%nobs), stat)
    if (stat /= SS_OK) then
      info = stat
      return
    end if
    call c_f_pointer(u_init, ui, [b%k_states])
    call djs_state_disturbances(b, f%res, m2(u_eta, b%k_posdef, b%nobs), ui, &
                                m2(eta, b%k_posdef, b%nobs), &
                                m2(alpha, b%k_states, b%nobs), stat)
    info = stat
  end function ss_djs_simulation_smoother

  !> Bias of the smoothed state from estimating the parameters (DK 7.3.7):
  !> params (k) and cov_params (k, k) from a fit; ndraw draws (even when
  !> antithetic); seed > 0 seeds the random numbers. bias_alpha (m, n),
  !> bias_V (m, m, n) or NULL; failed receives the number of skipped draws.
  integer(c_int) function ss_estimation_bias(mhandle, params, cov_params, ndraw, &
                                             antithetic, seed, bias_alpha, bias_V, &
                                             failed) &
    bind(C, name="ss_estimation_bias") result(info)
    type(c_ptr), value :: mhandle, params, cov_params, bias_alpha, bias_V
    integer(c_int), value :: ndraw, antithetic, seed
    integer(c_int), intent(out) :: failed
    type(model_box), pointer :: mb
    type(fit_result_t) :: fr
    real(dp), pointer :: pp(:), cp(:, :)
    integer, allocatable :: sv(:)
    integer :: k, m, n, stat, nf, ns, i

    failed = 0
    info = SS_ERR_DIM
    mb => get_model(mhandle)
    if (.not. associated(mb) .or. ndraw < 1) return
    k = mb%model%k_params; m = mb%model%rep%k_states; n = mb%model%rep%nobs
    call c_f_pointer(params, pp, [k])
    call c_f_pointer(cov_params, cp, [k, k])
    fr%params = pp
    fr%cov_params = cp
    if (seed > 0) then
      call random_seed(size=ns)
      sv = [(int(seed) + 7919 * i, i=1, ns)]
      call random_seed(put=sv)
    end if
    if (c_associated(bias_V)) then
      call estimation_bias(mb%model, fr, int(ndraw), m2(bias_alpha, m, n), stat, &
                           bias_V=m3(bias_V, m, m, n), antithetic=antithetic /= 0, &
                           failed=nf)
    else
      call estimation_bias(mb%model, fr, int(ndraw), m2(bias_alpha, m, n), stat, &
                           antithetic=antithetic /= 0, failed=nf)
    end if
    failed = nf
    info = stat
  end function ss_estimation_bias
end module statespace_capi_extras
