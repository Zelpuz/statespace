!> C interface (bind(C)) to the representation, filter and smoother, for the
!> Python package (ctypes) and other C callers.
!>
!> Conventions:
!>   * Objects are opaque handles (`void *`), created by `ss_*_new` or by the
!>     routine that computes them, and released with the matching `ss_*_free`.
!>   * Arrays are float64 in Fortran (column-major) order. Inputs are copied
!>     in. Outputs are copied into buffers the caller allocates; their sizes
!>     follow from (p, m, r, n) and the time dimension of each matrix.
!>   * Every routine returns the status code (SS_OK = 0, see statespace_kinds).
!>   * Missing observations are NaN in y.
module statespace_capi
  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, c_null_ptr, c_loc, &
                                         c_f_pointer, c_associated, c_null_char
  use statespace_kinds, only: dp, SS_OK, SS_ERR_DIM
  use statespace_rep, only: ssm_rep_t, ssm_rep
  use statespace_filter, only: filter_result_t, kalman_filter, loglike, loglike_concentrated
  use statespace_smoother, only: smoother_result_t, state_smoother
  use statespace_forecast, only: forecast_result_t, forecast
  use statespace_simsmooth, only: simsmooth_result_t, simulate, simulation_smoother
  use statespace_diagnostics, only: standardized_residuals, diagnostic_start, ljung_box, &
                                    jarque_bera, breakvar_test
  use statespace_filter, only: steady_state
  implicit none
  private

  character(len=*), parameter :: version = "0.1.1"

  !> Codes for the representation's arrays in ss_rep_set.
  integer(c_int), parameter, public :: SS_ARR_Y = 1, SS_ARR_Z = 2, SS_ARR_H = 3, SS_ARR_T = 4, &
                                       SS_ARR_R = 5, SS_ARR_Q = 6, SS_ARR_C = 7, SS_ARR_D = 8
  !> Codes for ss_rep_set_int and ss_rep_set_real.
  integer(c_int), parameter, public :: SS_OPT_FILTER_METHOD = 1, SS_OPT_DIFFUSE_METHOD = 2, &
                                       SS_OPT_LOGLIKELIHOOD_BURN = 3, SS_OPT_MARGINAL = 4
  integer(c_int), parameter, public :: SS_OPT_TOL_DIFFUSE = 1, SS_OPT_TOL_STEADY = 2
  !> Codes for ss_filter_get.
  integer(c_int), parameter, public :: SS_F_A = 1, SS_F_P = 2, SS_F_PINF = 3, SS_F_ATT = 4, &
                                       SS_F_PTT = 5, SS_F_YHAT = 6, SS_F_V = 7, SS_F_F = 8, &
                                       SS_F_FINF = 9, SS_F_FINV = 10, SS_F_K = 11, &
                                       SS_F_LLF_OBS = 12
  !> Codes for ss_smoother_get.
  integer(c_int), parameter, public :: SS_S_ALPHAHAT = 1, SS_S_V = 2, SS_S_R = 3, SS_S_N = 4, &
                                       SS_S_EPSHAT = 5, SS_S_EPSVAR = 6, SS_S_ETAHAT = 7, &
                                       SS_S_ETAVAR = 8

  type, public :: filter_box
    type(filter_result_t) :: res
  end type filter_box

  type, public :: smoother_box
    type(smoother_result_t) :: res
  end type smoother_box

  public :: ss_version, ss_rep_new, ss_rep_free, ss_rep_set, ss_rep_set_int, ss_rep_set_real
  public :: get_rep, get_filter, get_smoother, copy_out, ss_rep_info, ss_rep_get
  public :: ss_rep_init_known, ss_rep_init_diffuse, ss_rep_init_approx_diffuse
  public :: ss_rep_init_stationary, ss_rep_init_general, ss_rep_init_block
  public :: ss_loglike, ss_loglike_concentrated
  public :: ss_filter, ss_filter_free, ss_filter_get, ss_filter_scalars
  public :: ss_smooth, ss_smoother_free, ss_smoother_get
  public :: ss_forecast, ss_simulate, ss_simulation_smoother, ss_steady_state
  public :: ss_standardized_residuals, ss_diagnostic_start, ss_ljung_box, ss_jarque_bera
  public :: ss_breakvar

contains

  ! ------------------------------------------------------------- general

  !> Copy the version string (NUL-terminated) into buf(len).
  integer(c_int) function ss_version(buf, len) bind(C, name="ss_version") result(info)
    integer(c_int), value :: len
    character(kind=c_char), intent(out) :: buf(len)
    integer :: i, k

    k = min(len - 1, len_trim(version))
    do i = 1, k
      buf(i) = version(i:i)
    end do
    if (len > 0) buf(k + 1) = c_null_char
    info = SS_OK
  end function ss_version

  ! ------------------------------------------------------ representation

  !> A representation handle points at an ssm_rep_t: one owned by the
  !> caller (ss_rep_new) or one inside a model (ss_model_rep).
  function get_rep(handle) result(rp)
    type(c_ptr), intent(in) :: handle
    type(ssm_rep_t), pointer :: rp

    rp => null()
    if (c_associated(handle)) call c_f_pointer(handle, rp)
  end function get_rep

  !> A representation for p series, m states, r disturbances and n periods,
  !> with zero system matrices (time-invariant) and y = 0.
  integer(c_int) function ss_rep_new(p, m, r, n, handle) bind(C, name="ss_rep_new") result(info)
    integer(c_int), value :: p, m, r, n
    type(c_ptr), intent(out) :: handle
    type(ssm_rep_t), pointer :: b
    real(dp), allocatable :: y(:, :)

    handle = c_null_ptr
    info = SS_ERR_DIM
    if (p < 1 .or. m < 1 .or. r < 0 .or. n < 1) return
    allocate (b)
    allocate (y(p, n), source=0.0_dp)
    b = ssm_rep(y, m=m, r=r)
    handle = c_loc(b)
    info = SS_OK
  end function ss_rep_new

  integer(c_int) function ss_rep_free(handle) bind(C, name="ss_rep_free") result(info)
    type(c_ptr), value :: handle
    type(ssm_rep_t), pointer :: b

    b => get_rep(handle)
    if (associated(b)) deallocate (b)
    info = SS_OK
  end function ss_rep_free

  !> Set one of the arrays (SS_ARR_*) from data. nt is its time dimension,
  !> 1 or n (ignored for y). Shapes: y (p, n), Z (p, m, nt), H (p, p, nt),
  !> T (m, m, nt), R (m, r, nt), Q (r, r, nt), c (m, nt), d (p, nt).
  integer(c_int) function ss_rep_set(handle, code, nt, data) bind(C, name="ss_rep_set") result(info)
    type(c_ptr), value :: handle
    integer(c_int), value :: code, nt
    type(c_ptr), value :: data
    type(ssm_rep_t), pointer :: b
    real(dp), pointer :: x2(:, :), x3(:, :, :)
    integer :: p, m, r, n

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b) .or. .not. c_associated(data)) return
    p = b%k_endog; m = b%k_states; r = b%k_posdef; n = b%nobs
    if (code /= SS_ARR_Y .and. nt /= 1 .and. nt /= n) return
    select case (code)
    case (SS_ARR_Y)
      call c_f_pointer(data, x2, [p, n]); b%y = x2
    case (SS_ARR_Z)
      call c_f_pointer(data, x3, [p, m, nt]); b%Z = x3
    case (SS_ARR_H)
      call c_f_pointer(data, x3, [p, p, nt]); b%H = x3
    case (SS_ARR_T)
      call c_f_pointer(data, x3, [m, m, nt]); b%T = x3
    case (SS_ARR_R)
      call c_f_pointer(data, x3, [m, r, nt]); b%R = x3
    case (SS_ARR_Q)
      call c_f_pointer(data, x3, [r, r, nt]); b%Q = x3
    case (SS_ARR_C)
      call c_f_pointer(data, x2, [m, nt]); b%c = x2
    case (SS_ARR_D)
      call c_f_pointer(data, x2, [p, nt]); b%d = x2
    case default
      return
    end select
    info = SS_OK
  end function ss_rep_set

  function get_filter(handle) result(f)
    type(c_ptr), intent(in) :: handle
    type(filter_box), pointer :: f

    f => null()
    if (c_associated(handle)) call c_f_pointer(handle, f)
  end function get_filter

  function get_smoother(handle) result(s)
    type(c_ptr), intent(in) :: handle
    type(smoother_box), pointer :: s

    s => null()
    if (c_associated(handle)) call c_f_pointer(handle, s)
  end function get_smoother

  !> Dimensions (p, m, r, n) and the time dimension nts(i) of each array
  !> SS_ARR_i (i = 1..8; n for y).
  integer(c_int) function ss_rep_info(handle, p, m, r, n, nts) bind(C, name="ss_rep_info") &
    result(info)
    type(c_ptr), value :: handle
    integer(c_int), intent(out) :: p, m, r, n, nts(8)
    type(ssm_rep_t), pointer :: b

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    p = b%k_endog; m = b%k_states; r = b%k_posdef; n = b%nobs
    nts = [size(b%y, 2), size(b%Z, 3), size(b%H, 3), size(b%T, 3), size(b%R, 3), &
           size(b%Q, 3), size(b%c, 2), size(b%d, 2)]
    info = SS_OK
  end function ss_rep_info

  !> Copy one of the arrays (SS_ARR_*) into out, with the shape given by
  !> ss_rep_info.
  integer(c_int) function ss_rep_get(handle, code, out) bind(C, name="ss_rep_get") result(info)
    type(c_ptr), value :: handle, out
    integer(c_int), value :: code
    type(ssm_rep_t), pointer :: b

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b) .or. .not. c_associated(out)) return
    select case (code)
    case (SS_ARR_Y)
      call copy_out(b%y, out)
    case (SS_ARR_Z)
      call copy_out(b%Z, out)
    case (SS_ARR_H)
      call copy_out(b%H, out)
    case (SS_ARR_T)
      call copy_out(b%T, out)
    case (SS_ARR_R)
      call copy_out(b%R, out)
    case (SS_ARR_Q)
      call copy_out(b%Q, out)
    case (SS_ARR_C)
      call copy_out(b%c, out)
    case (SS_ARR_D)
      call copy_out(b%d, out)
    case default
      return
    end select
    info = SS_OK
  end function ss_rep_get

  !> Integer options (SS_OPT_*): filter_method, diffuse_method,
  !> loglikelihood_burn, marginal_likelihood (0/1).
  integer(c_int) function ss_rep_set_int(handle, code, value) bind(C, name="ss_rep_set_int") &
    result(info)
    type(c_ptr), value :: handle
    integer(c_int), value :: code, value
    type(ssm_rep_t), pointer :: b

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    select case (code)
    case (SS_OPT_FILTER_METHOD)
      b%filter_method = value
    case (SS_OPT_DIFFUSE_METHOD)
      b%diffuse_method = value
    case (SS_OPT_LOGLIKELIHOOD_BURN)
      b%loglikelihood_burn = value
    case (SS_OPT_MARGINAL)
      b%marginal_likelihood = value /= 0
    case default
      return
    end select
    info = SS_OK
  end function ss_rep_set_int

  !> Real options (SS_OPT_*): tol_diffuse, tol_steady.
  integer(c_int) function ss_rep_set_real(handle, code, value) bind(C, name="ss_rep_set_real") &
    result(info)
    type(c_ptr), value :: handle
    integer(c_int), value :: code
    real(c_double), value :: value
    type(ssm_rep_t), pointer :: b

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    select case (code)
    case (SS_OPT_TOL_DIFFUSE)
      b%tol_diffuse = value
    case (SS_OPT_TOL_STEADY)
      b%tol_steady = value
    case default
      return
    end select
    info = SS_OK
  end function ss_rep_set_real

  integer(c_int) function ss_rep_init_known(handle, a1, P1) bind(C, name="ss_rep_init_known") &
    result(info)
    type(c_ptr), value :: handle, a1, P1
    type(ssm_rep_t), pointer :: b
    real(dp), pointer :: a(:), P(:, :)
    integer :: m

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    m = b%k_states
    call c_f_pointer(a1, a, [m])
    call c_f_pointer(P1, P, [m, m])
    call b%initialize_known(a, P)
    info = SS_OK
  end function ss_rep_init_known

  integer(c_int) function ss_rep_init_diffuse(handle) bind(C, name="ss_rep_init_diffuse") &
    result(info)
    type(c_ptr), value :: handle
    type(ssm_rep_t), pointer :: b

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    call b%initialize_diffuse()
    info = SS_OK
  end function ss_rep_init_diffuse

  integer(c_int) function ss_rep_init_approx_diffuse(handle, kappa) &
    bind(C, name="ss_rep_init_approx_diffuse") result(info)
    type(c_ptr), value :: handle
    real(c_double), value :: kappa
    type(ssm_rep_t), pointer :: b

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    call b%initialize_approximate_diffuse(kappa=kappa)
    info = SS_OK
  end function ss_rep_init_approx_diffuse

  integer(c_int) function ss_rep_init_stationary(handle) bind(C, name="ss_rep_init_stationary") &
    result(info)
    type(c_ptr), value :: handle
    type(ssm_rep_t), pointer :: b

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    call b%initialize_stationary()
    info = SS_OK
  end function ss_rep_init_stationary

  integer(c_int) function ss_rep_init_general(handle, a1, Pstar, Pinf) &
    bind(C, name="ss_rep_init_general") result(info)
    type(c_ptr), value :: handle, a1, Pstar, Pinf
    type(ssm_rep_t), pointer :: b
    real(dp), pointer :: a(:), Ps(:, :), Pi(:, :)
    integer :: m

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    m = b%k_states
    call c_f_pointer(a1, a, [m])
    call c_f_pointer(Pstar, Ps, [m, m])
    call c_f_pointer(Pinf, Pi, [m, m])
    call b%initialize_general(a, Ps, Pi)
    info = SS_OK
  end function ss_rep_init_general

  !> Initialize states first..last (1-based) with kind (INIT_*). For
  !> INIT_KNOWN, a1 (length last-first+1) and P1 (square) are required;
  !> otherwise pass NULL. kappa is used by INIT_APPROX_DIFFUSE.
  integer(c_int) function ss_rep_init_block(handle, first, last, kind, a1, P1, kappa) &
    bind(C, name="ss_rep_init_block") result(info)
    type(c_ptr), value :: handle, a1, P1
    integer(c_int), value :: first, last, kind
    real(c_double), value :: kappa
    type(ssm_rep_t), pointer :: b
    real(dp), pointer :: a(:), P(:, :)
    integer :: nb

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    if (first < 1 .or. last < first .or. last > b%k_states) return
    nb = last - first + 1
    if (c_associated(a1) .and. c_associated(P1)) then
      call c_f_pointer(a1, a, [nb])
      call c_f_pointer(P1, P, [nb, nb])
      call b%initialize_block(first, last, kind, a1=a, P1=P, kappa=kappa)
    else
      call b%initialize_block(first, last, kind, kappa=kappa)
    end if
    info = SS_OK
  end function ss_rep_init_block

  ! ----------------------------------------------------------- filtering

  integer(c_int) function ss_loglike(handle, llf) bind(C, name="ss_loglike") result(info)
    type(c_ptr), value :: handle
    real(c_double), intent(out) :: llf
    type(ssm_rep_t), pointer :: b
    integer :: stat

    llf = 0.0_dp
    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    llf = loglike(b, stat)
    info = stat
  end function ss_loglike

  integer(c_int) function ss_loglike_concentrated(handle, llf, scale) &
    bind(C, name="ss_loglike_concentrated") result(info)
    type(c_ptr), value :: handle
    real(c_double), intent(out) :: llf, scale
    type(ssm_rep_t), pointer :: b
    integer :: stat

    llf = 0.0_dp; scale = 0.0_dp
    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    llf = loglike_concentrated(b, scale, stat)
    info = stat
  end function ss_loglike_concentrated

  !> Run the Kalman filter; fhandle receives the result (free with
  !> ss_filter_free).
  integer(c_int) function ss_filter(handle, fhandle) bind(C, name="ss_filter") result(info)
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
    call kalman_filter(b, f%res, stat)
    info = stat
    if (stat /= SS_OK) then
      deallocate (f)
      return
    end if
    fhandle = c_loc(f)
  end function ss_filter

  integer(c_int) function ss_filter_free(fhandle) bind(C, name="ss_filter_free") result(info)
    type(c_ptr), value :: fhandle
    type(filter_box), pointer :: f

    if (c_associated(fhandle)) then
      call c_f_pointer(fhandle, f)
      deallocate (f)
    end if
    info = SS_OK
  end function ss_filter_free

  !> Scalars of a filter result.
  integer(c_int) function ss_filter_scalars(fhandle, llf, nobs_diffuse, k_diffuse, t_steady) &
    bind(C, name="ss_filter_scalars") result(info)
    type(c_ptr), value :: fhandle
    real(c_double), intent(out) :: llf
    integer(c_int), intent(out) :: nobs_diffuse, k_diffuse, t_steady
    type(filter_box), pointer :: f

    info = SS_ERR_DIM
    if (.not. c_associated(fhandle)) return
    call c_f_pointer(fhandle, f)
    llf = f%res%llf
    nobs_diffuse = f%res%nobs_diffuse
    k_diffuse = f%res%k_diffuse
    t_steady = f%res%t_steady
    info = SS_OK
  end function ss_filter_scalars

  !> Copy a filter output (SS_F_*) into out. Sizes: a (m, n+1), P and Pinf
  !> (m, m, n+1), att (m, n), Ptt (m, m, n), yhat and v (p, n), F, Finf and
  !> Finv (p, p, n), K (m, p, n), llf_obs (n).
  integer(c_int) function ss_filter_get(fhandle, code, out) bind(C, name="ss_filter_get") &
    result(info)
    type(c_ptr), value :: fhandle, out
    integer(c_int), value :: code
    type(filter_box), pointer :: f

    info = SS_ERR_DIM
    if (.not. c_associated(fhandle) .or. .not. c_associated(out)) return
    call c_f_pointer(fhandle, f)
    associate (r => f%res)
      select case (code)
      case (SS_F_A)
        call copy_out(r%a, out)
      case (SS_F_P)
        call copy_out(r%P, out)
      case (SS_F_PINF)
        call copy_out(r%Pinf, out)
      case (SS_F_ATT)
        call copy_out(r%att, out)
      case (SS_F_PTT)
        call copy_out(r%Ptt, out)
      case (SS_F_YHAT)
        call copy_out(r%yhat, out)
      case (SS_F_V)
        call copy_out(r%v, out)
      case (SS_F_F)
        call copy_out(r%F, out)
      case (SS_F_FINF)
        call copy_out(r%Finf, out)
      case (SS_F_FINV)
        call copy_out(r%Finv, out)
      case (SS_F_K)
        call copy_out(r%K, out)
      case (SS_F_LLF_OBS)
        call copy_out(r%llf_obs, out)
      case default
        return
      end select
    end associate
    info = SS_OK
  end function ss_filter_get

  ! ----------------------------------------------------------- smoothing

  !> Run the state and disturbance smoother on a filter result; shandle
  !> receives the result (free with ss_smoother_free).
  integer(c_int) function ss_smooth(handle, fhandle, shandle) bind(C, name="ss_smooth") &
    result(info)
    type(c_ptr), value :: handle, fhandle
    type(c_ptr), intent(out) :: shandle
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    type(smoother_box), pointer :: s
    integer :: stat

    shandle = c_null_ptr
    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b) .or. .not. c_associated(fhandle)) return
    call c_f_pointer(fhandle, f)
    allocate (s)
    call state_smoother(b, f%res, s%res, stat)
    info = stat
    if (stat /= SS_OK) then
      deallocate (s)
      return
    end if
    shandle = c_loc(s)
  end function ss_smooth

  integer(c_int) function ss_smoother_free(shandle) bind(C, name="ss_smoother_free") result(info)
    type(c_ptr), value :: shandle
    type(smoother_box), pointer :: s

    if (c_associated(shandle)) then
      call c_f_pointer(shandle, s)
      deallocate (s)
    end if
    info = SS_OK
  end function ss_smoother_free

  !> Copy a smoother output (SS_S_*) into out. Sizes: alphahat (m, n),
  !> V (m, m, n), r (m, n+1) for t = 0..n, N (m, m, n+1), epshat (p, n),
  !> epsvar (p, p, n), etahat (r, n), etavar (r, r, n).
  integer(c_int) function ss_smoother_get(shandle, code, out) bind(C, name="ss_smoother_get") &
    result(info)
    type(c_ptr), value :: shandle, out
    integer(c_int), value :: code
    type(smoother_box), pointer :: s

    info = SS_ERR_DIM
    if (.not. c_associated(shandle) .or. .not. c_associated(out)) return
    call c_f_pointer(shandle, s)
    associate (r => s%res)
      select case (code)
      case (SS_S_ALPHAHAT)
        call copy_out(r%alphahat, out)
      case (SS_S_V)
        call copy_out(r%V, out)
      case (SS_S_R)
        call copy_out(r%r, out)
      case (SS_S_N)
        call copy_out(r%N, out)
      case (SS_S_EPSHAT)
        call copy_out(r%epshat, out)
      case (SS_S_EPSVAR)
        call copy_out(r%epsvar, out)
      case (SS_S_ETAHAT)
        call copy_out(r%etahat, out)
      case (SS_S_ETAVAR)
        call copy_out(r%etavar, out)
      case default
        return
      end select
    end associate
    info = SS_OK
  end function ss_smoother_get

  ! ------------------------------------------ forecasting and simulation

  !> Forecasts h steps past the end of the sample from a filter result
  !> (time-invariant models, DK 4.11): mean (p, h) and cov (p, p, h).
  integer(c_int) function ss_forecast(handle, fhandle, h, mean, cov) bind(C, name="ss_forecast") &
    result(info)
    type(c_ptr), value :: handle, fhandle, mean, cov
    integer(c_int), value :: h
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f
    type(forecast_result_t) :: fc
    integer :: stat

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b) .or. .not. c_associated(fhandle)) return
    call c_f_pointer(fhandle, f)
    call forecast(b, f%res, int(h), fc, stat)
    info = stat
    if (stat /= SS_OK) return
    call copy_out(fc%mean, mean)
    call copy_out(fc%cov, cov)
  end function ss_forecast

  !> Simulate from the model with standard normal variates u_init (m),
  !> u_eps (p, n) and u_eta (r, n): y (p, n), alpha (m, n), eps, eta.
  integer(c_int) function ss_simulate(handle, u_init, u_eps, u_eta, y, alpha, eps, eta) &
    bind(C, name="ss_simulate") result(info)
    type(c_ptr), value :: handle, u_init, u_eps, u_eta, y, alpha, eps, eta
    type(ssm_rep_t), pointer :: b
    real(dp), pointer :: ui(:), ue(:, :), un(:, :), yo(:, :), ao(:, :), eo(:, :), no(:, :)
    integer :: p, m, r, n, stat

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    p = b%k_endog; m = b%k_states; r = b%k_posdef; n = b%nobs
    call c_f_pointer(u_init, ui, [m]); call c_f_pointer(u_eps, ue, [p, n])
    call c_f_pointer(u_eta, un, [r, n]); call c_f_pointer(y, yo, [p, n])
    call c_f_pointer(alpha, ao, [m, n]); call c_f_pointer(eps, eo, [p, n])
    call c_f_pointer(eta, no, [r, n])
    call simulate(b, ui, ue, un, yo, ao, eo, no, stat)
    info = stat
  end function ss_simulate

  !> One draw of alpha, eps and eta given the data (DK 4.9.2) from standard
  !> normal variates; outputs state (m, n), eps (p, n), eta (r, n).
  integer(c_int) function ss_simulation_smoother(handle, u_init, u_eps, u_eta, state, eps, eta) &
    bind(C, name="ss_simulation_smoother") result(info)
    type(c_ptr), value :: handle, u_init, u_eps, u_eta, state, eps, eta
    type(ssm_rep_t), pointer :: b
    type(simsmooth_result_t) :: sim
    real(dp), pointer :: ui(:), ue(:, :), un(:, :)
    integer :: p, m, r, n, stat

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    p = b%k_endog; m = b%k_states; r = b%k_posdef; n = b%nobs
    call c_f_pointer(u_init, ui, [m]); call c_f_pointer(u_eps, ue, [p, n])
    call c_f_pointer(u_eta, un, [r, n])
    call simulation_smoother(b, ui, ue, un, sim, stat)
    info = stat
    if (stat /= SS_OK) return
    call copy_out(sim%state, state)
    call copy_out(sim%eps, eps)
    call copy_out(sim%eta, eta)
  end function ss_simulation_smoother

  !> Steady state P (m, m) and F (p, p) of a time-invariant model (DK 4.3.4).
  integer(c_int) function ss_steady_state(handle, P, F) bind(C, name="ss_steady_state") &
    result(info)
    type(c_ptr), value :: handle, P, F
    type(ssm_rep_t), pointer :: b
    real(dp), pointer :: Po(:, :), Fo(:, :)
    integer :: stat

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b)) return
    call c_f_pointer(P, Po, [b%k_states, b%k_states])
    call c_f_pointer(F, Fo, [b%k_endog, b%k_endog])
    call steady_state(b, Po, Fo, stat)
    info = stat
  end function ss_steady_state

  ! --------------------------------------------------------- diagnostics

  !> Standardized one-step residuals (p, n) of a filter result (DK 7.5).
  integer(c_int) function ss_standardized_residuals(fhandle, out) &
    bind(C, name="ss_standardized_residuals") result(info)
    type(c_ptr), value :: fhandle, out
    type(filter_box), pointer :: f
    real(dp), pointer :: e(:, :)

    info = SS_ERR_DIM
    if (.not. c_associated(fhandle) .or. .not. c_associated(out)) return
    call c_f_pointer(fhandle, f)
    call c_f_pointer(out, e, [f%res%k_endog, f%res%nobs])
    call standardized_residuals(f%res, e)
    info = SS_OK
  end function ss_standardized_residuals

  !> First period (1-based) for the diagnostics: after the burn-in and the
  !> diffuse period.
  integer(c_int) function ss_diagnostic_start(handle, fhandle, t0) &
    bind(C, name="ss_diagnostic_start") result(info)
    type(c_ptr), value :: handle, fhandle
    integer(c_int), intent(out) :: t0
    type(ssm_rep_t), pointer :: b
    type(filter_box), pointer :: f

    info = SS_ERR_DIM
    b => get_rep(handle)
    if (.not. associated(b) .or. .not. c_associated(fhandle)) return
    call c_f_pointer(fhandle, f)
    t0 = diagnostic_start(b, f%res)
    info = SS_OK
  end function ss_diagnostic_start

  !> Ljung-Box Q(k), k = 1..lags, and p-values for x (n); NaNs skipped.
  integer(c_int) function ss_ljung_box(n, x, lags, model_df, stat, pvalue) &
    bind(C, name="ss_ljung_box") result(info)
    integer(c_int), value :: n, lags, model_df
    type(c_ptr), value :: x, stat, pvalue
    real(dp), pointer :: xi(:), so(:), po(:)

    info = SS_ERR_DIM
    if (n < 2 .or. lags < 1) return
    call c_f_pointer(x, xi, [n]); call c_f_pointer(stat, so, [lags])
    call c_f_pointer(pvalue, po, [lags])
    call ljung_box(xi, so, po, int(model_df))
    info = SS_OK
  end function ss_ljung_box

  !> Jarque-Bera normality test for x (n); NaNs skipped.
  integer(c_int) function ss_jarque_bera(n, x, jb, pvalue, skew, kurtosis) &
    bind(C, name="ss_jarque_bera") result(info)
    integer(c_int), value :: n
    type(c_ptr), value :: x
    real(c_double), intent(out) :: jb, pvalue, skew, kurtosis
    real(dp), pointer :: xi(:)

    info = SS_ERR_DIM
    if (n < 3) return
    call c_f_pointer(x, xi, [n])
    call jarque_bera(xi, jb, pvalue, skew, kurtosis)
    info = SS_OK
  end function ss_jarque_bera

  !> Heteroskedasticity test H(h) for x (n); h <= 0 for the default.
  integer(c_int) function ss_breakvar(n, x, h, stat, pvalue) bind(C, name="ss_breakvar") &
    result(info)
    integer(c_int), value :: n, h
    type(c_ptr), value :: x
    real(c_double), intent(out) :: stat, pvalue
    real(dp), pointer :: xi(:)

    info = SS_ERR_DIM
    if (n < 3) return
    call c_f_pointer(x, xi, [n])
    if (h > 0) then
      call breakvar_test(xi, stat, pvalue, int(h))
    else
      call breakvar_test(xi, stat, pvalue)
    end if
    info = SS_OK
  end function ss_breakvar

  !> Copy an array of any rank into the buffer at out, in array element order.
  subroutine copy_out(x, out)
    real(dp), intent(in) :: x(..)
    type(c_ptr), intent(in) :: out
    real(dp), pointer :: buf(:)

    call c_f_pointer(out, buf, [size(x)])
    select rank (x)
    rank (1)
      buf = x
    rank (2)
      buf = reshape(x, [size(x)])
    rank (3)
      buf = reshape(x, [size(x)])
    end select
  end subroutine copy_out
end module statespace_capi
