!> Further smoothing results: updating smoothed estimates and fixed-point /
!> fixed-lag smoothing (DK 4.4.5-4.4.6), other state smoothing algorithms
!> (DK 4.6), covariances of the smoothed state across time (DK 4.7) and
!> filtering and smoothing weights (DK 4.8).
!>
!> The DK recursions here are written in terms of the conventional K_t and
!> F_t^-1; `conventional_gain` supplies them in any period after the diffuse
!> period, whichever filter was used.
module statespace_smoothing
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use statespace_kinds, only: dp, SS_OK, SS_ERR_DIM, SS_ERR_UNSUPPORTED
  use statespace_linalg, only: gemm, gemv, eye, psd_solve, symmetrize, chol_inv, solve
  use statespace_rep, only: ssm_rep_t, tidx
  use statespace_filter, only: filter_result_t, kalman_filter, METHOD_UNIVARIATE, &
                               METHOD_DIFFUSE_MV
  use statespace_smoother, only: smoother_result_t, state_smoother, diffuse_mv_gains
  implicit none
  private

  public :: innovation_transition, smoothed_state_cov_between, smoothed_state_autocov
  public :: smoothed_state_weights, fast_state_smoother, classical_state_smoother
  public :: two_filter_smoother, conventional_gain, update_smoothed, fixed_point_smoother
  public :: fixed_lag_smoother, filtered_state_weights, whittle_smoother

contains

  !> L_t, the transition of the one-step state prediction error
  !> a_t+1 - alpha_t+1 = L_t (a_t - alpha_t) + ...  (DK 4.3):
  !> T_t - K_t Z_t in conventional periods, and T_t L_t,p ... L_t,1 with
  !> L_t,i = I - K_t,i Z*_i in univariate ones. Not defined in diffuse periods.
  subroutine innovation_transition(rep, fres, t, L, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    integer, intent(in) :: t
    real(dp), intent(out), contiguous :: L(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: A(:, :), K(:)
    integer :: m, i, j, it, iz

    info = SS_ERR_UNSUPPORTED
    if (t <= fres%nobs_diffuse) return
    info = SS_OK
    m = rep%k_states
    it = tidx(size(rep%T, 3), t)
    iz = tidx(size(rep%Z, 3), t)
    if (fres%method(t) /= METHOD_UNIVARIATE) then
      L = rep%T(:, :, it)
      call gemm('N', 'N', -1.0_dp, fres%K(:, :, t), rep%Z(:, :, iz), 1.0_dp, L)
      return
    end if

    A = eye(m)
    do i = 1, fres%uv_n(t)
      if (.not. fres%uv_Fstar(i, t) > 0.0_dp) cycle
      K = fres%uv_Mstar(:, i, t) / fres%uv_Fstar(i, t)
      ! A <- (I - K z') A
      associate (z => fres%uv_Z(i, :, t))
        do j = 1, m
          A(:, j) = A(:, j) - K * dot_product(z, A(:, j))
        end do
      end associate
    end do
    call gemm('N', 'N', 1.0_dp, rep%T(:, :, it), A, 0.0_dp, L)
  end subroutine innovation_transition

  !> Conventional K_t = T P Z' F^-1 and F_t^-1 (over the observed elements,
  !> zero-padded) for a period after the diffuse period. Univariate periods
  !> have none stored, so they are formed from the stored P_t and F_t.
  subroutine conventional_gain(rep, fres, t, K, Finv, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    integer, intent(in) :: t
    real(dp), intent(out), contiguous :: K(:, :), Finv(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: Fo(:, :)
    integer, allocatable :: idx(:)
    real(dp) :: logdet
    integer :: i, it, iz

    info = SS_ERR_UNSUPPORTED
    if (t <= fres%nobs_diffuse) return
    info = SS_OK
    if (fres%method(t) /= METHOD_UNIVARIATE) then
      K = fres%K(:, :, t)
      Finv = fres%Finv(:, :, t)
      return
    end if
    it = tidx(size(rep%T, 3), t)
    iz = tidx(size(rep%Z, 3), t)
    idx = pack([(i, i=1, rep%k_endog)], .not. ieee_is_nan(rep%y(:, t)))
    Fo = fres%F(idx, idx, t)
    call chol_inv(Fo, logdet, info)
    if (info /= SS_OK) return
    Finv = 0.0_dp
    Finv(idx, idx) = Fo
    K = matmul(rep%T(:, :, it), matmul(fres%P(:, :, t), matmul(transpose(rep%Z(:, :, iz)), Finv)))
  end subroutine conventional_gain

  !> Z_t' F_t^-1 v_t (missing elements dropped), Z_t' F_t^-1 Z_t and
  !> L_t = T_t - K_t Z_t for a period after the diffuse period.
  subroutine period_terms(rep, fres, t, zfv, ZFZ, L, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    integer, intent(in) :: t
    real(dp), allocatable, intent(out) :: zfv(:), ZFZ(:, :), L(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: K(:, :), Finv(:, :), vz(:), Zt(:, :)
    integer :: it, iz

    allocate (K(rep%k_states, rep%k_endog), Finv(rep%k_endog, rep%k_endog))
    call conventional_gain(rep, fres, t, K, Finv, info)
    if (info /= SS_OK) return
    it = tidx(size(rep%T, 3), t)
    iz = tidx(size(rep%Z, 3), t)
    Zt = rep%Z(:, :, iz)
    vz = merge(0.0_dp, fres%v(:, t), ieee_is_nan(fres%v(:, t)))
    zfv = matmul(transpose(Zt), matmul(Finv, vz))
    ZFZ = matmul(transpose(Zt), matmul(Finv, Zt))
    L = rep%T(:, :, it) - matmul(K, Zt)
  end subroutine period_terms

  !> Updating smoothed estimates as observations arrive (DK 4.4.5). On
  !> entry alphahat(:, 1:n0) and V(:, :, 1:n0) are smoothed given
  !> y_1..y_n0; `fres` is the filter over all n observations. On exit
  !> alphahat and V are smoothed given y_1..y_n. For each new y_k+1:
  !>
  !>   b_t|k+1 = b_t|k L_k',  b_k|k = I
  !>   alphahat_t|k+1 = alphahat_t|k + P_t b_t|k+1 Z_k+1' F_k+1^-1 v_k+1    (4.49)
  !>   V_t|k+1 = V_t|k - P_t b_t|k+1 Z_k+1' F_k+1^-1 Z_k+1 b_t|k+1' P_t    (4.50)
  !>
  !> for t <= k, and alphahat_k+1|k+1 = a_k+1|k+1, V_k+1|k+1 = P_k+1|k+1.
  !> Not defined with diffuse states.
  subroutine update_smoothed(rep, fres, n0, alphahat, V, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    integer, intent(in) :: n0
    real(dp), intent(inout), contiguous :: alphahat(:, :)   !< (m, n)
    real(dp), intent(inout), contiguous :: V(:, :, :)       !< (m, m, n)
    integer, intent(out) :: info
    real(dp), allocatable :: B(:, :, :), zfv(:), ZFZ(:, :), L(:, :), PB(:, :)
    integer :: m, n, t, k

    m = rep%k_states; n = rep%nobs
    info = SS_ERR_UNSUPPORTED
    if (fres%nobs_diffuse > 0) return
    info = SS_ERR_DIM
    if (n0 < 1 .or. n0 > n) return
    allocate (B(m, m, n))
    ! b_t|n0 = L_t' L_t+1' ... L_n0-1'
    B(:, :, n0) = eye(m)
    do t = n0 - 1, 1, -1
      call period_terms(rep, fres, t, zfv, ZFZ, L, info)
      if (info /= SS_OK) return
      B(:, :, t) = matmul(transpose(L), B(:, :, t + 1))
    end do
    do k = n0, n - 1
      call period_terms(rep, fres, k, zfv, ZFZ, L, info)
      if (info /= SS_OK) return
      do t = 1, k
        B(:, :, t) = matmul(B(:, :, t), transpose(L))
      end do
      call period_terms(rep, fres, k + 1, zfv, ZFZ, L, info)
      if (info /= SS_OK) return
      do t = 1, k
        PB = matmul(fres%P(:, :, t), B(:, :, t))
        alphahat(:, t) = alphahat(:, t) + matmul(PB, zfv)
        V(:, :, t) = V(:, :, t) - matmul(PB, matmul(ZFZ, transpose(PB)))
        call symmetrize(V(:, :, t))
      end do
      alphahat(:, k + 1) = fres%att(:, k + 1)
      V(:, :, k + 1) = fres%Ptt(:, :, k + 1)
      B(:, :, k + 1) = eye(m)
    end do
    info = SS_OK
  end subroutine update_smoothed

  !> Fixed-point smoother (DK 4.4.6): alphahat_t|k and V_t|k for fixed t and
  !> k = t..n, returned in path_a(:, k - t + 1) and path_V(:, :, k - t + 1),
  !> by the updating recursions (4.49)-(4.50) started from a_t|t, P_t|t.
  subroutine fixed_point_smoother(rep, fres, t, path_a, path_V, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    integer, intent(in) :: t
    real(dp), intent(out), contiguous :: path_a(:, :)       !< (m, n-t+1)
    real(dp), intent(out), contiguous :: path_V(:, :, :)    !< (m, m, n-t+1)
    integer, intent(out) :: info
    real(dp), allocatable :: b(:, :), zfv(:), ZFZ(:, :), L(:, :), PB(:, :)
    integer :: k, j

    info = SS_ERR_UNSUPPORTED
    if (t <= fres%nobs_diffuse) return
    info = SS_ERR_DIM
    if (t < 1 .or. t > rep%nobs .or. size(path_a, 2) /= rep%nobs - t + 1) return
    path_a(:, 1) = fres%att(:, t)
    path_V(:, :, 1) = fres%Ptt(:, :, t)
    b = eye(rep%k_states)
    do k = t, rep%nobs - 1
      j = k - t + 1
      call period_terms(rep, fres, k, zfv, ZFZ, L, info)
      if (info /= SS_OK) return
      b = matmul(b, transpose(L))
      call period_terms(rep, fres, k + 1, zfv, ZFZ, L, info)
      if (info /= SS_OK) return
      PB = matmul(fres%P(:, :, t), b)
      path_a(:, j + 1) = path_a(:, j) + matmul(PB, zfv)
      path_V(:, :, j + 1) = path_V(:, :, j) - matmul(PB, matmul(ZFZ, transpose(PB)))
      call symmetrize(path_V(:, :, j + 1))
    end do
    info = SS_OK
  end subroutine fixed_point_smoother

  !> Fixed-lag smoother (DK 4.4.6): alphahat_s|s+j and V_s|s+j for a fixed lag
  !> j, stored in column s = 1..n-j (later columns are NaN), from
  !>
  !>   alphahat_s|s+j = a_s + P_s r_s-1,   V_s|s+j = P_s - P_s N_s-1 P_s   (4.51)
  !>   r_t-1 = Z_t' F_t^-1 v_t + L_t' r_t,  N_t-1 = Z_t' F_t^-1 Z_t + L_t' N_t L_t
  !>
  !> run backwards from t = s+j with r = 0, N = 0.
  subroutine fixed_lag_smoother(rep, fres, j, alphahat, V, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    integer, intent(in) :: j
    real(dp), intent(out), contiguous :: alphahat(:, :)     !< (m, n)
    real(dp), intent(out), contiguous :: V(:, :, :)         !< (m, m, n)
    integer, intent(out) :: info
    real(dp), allocatable :: r(:), N(:, :), zfv(:), ZFZ(:, :), L(:, :)
    integer :: s, t, m

    m = rep%k_states
    info = SS_ERR_DIM
    if (j < 0 .or. j >= rep%nobs) return
    alphahat = ieee_value(1.0_dp, ieee_quiet_nan)
    V = ieee_value(1.0_dp, ieee_quiet_nan)
    info = SS_ERR_UNSUPPORTED
    if (fres%nobs_diffuse > 0) return
    do s = 1, rep%nobs - j
      allocate (r(m), N(m, m), source=0.0_dp)
      do t = s + j, s, -1
        call period_terms(rep, fres, t, zfv, ZFZ, L, info)
        if (info /= SS_OK) return
        r = zfv + matmul(transpose(L), r)
        N = ZFZ + matmul(transpose(L), matmul(N, L))
      end do
      alphahat(:, s) = fres%a(:, s) + matmul(fres%P(:, :, s), r)
      V(:, :, s) = fres%P(:, :, s) - matmul(fres%P(:, :, s), matmul(N, fres%P(:, :, s)))
      call symmetrize(V(:, :, s))
      deallocate (r, N)
    end do
    info = SS_OK
  end subroutine fixed_lag_smoother

  !> Filtering weights (DK 4.8.2): the predicted and filtered states as
  !> weighted sums of past observations (intercepts and the initial mean
  !> aside). Wa(:, :, t, j) is the weight of y_j in a_t and Watt(:, :, t, j)
  !> that in a_t|t (Table 4.5, dropping H_j):
  !>
  !>   a_t:    omega_jt = L_t-1 ... L_j+1 K_j                   (j < t)
  !>   a_t|t:  (I - P_t Z_t' F_t^-1 Z_t) omega_jt (j < t),  P_t Z_t' F_t^-1 (j = t)
  !>
  !> and zero otherwise. Not defined with diffuse states.
  subroutine filtered_state_weights(rep, fres, Wa, Watt, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    real(dp), intent(out), contiguous :: Wa(:, :, :, :)     !< (m, p, n, n)
    real(dp), intent(out), contiguous :: Watt(:, :, :, :)   !< (m, p, n, n)
    integer, intent(out) :: info
    real(dp), allocatable :: K(:, :), Finv(:, :), zfv(:), ZFZ(:, :), L(:, :), A(:, :), PZF(:, :)
    integer :: m, p, n, t, j, iz

    m = rep%k_states; p = rep%k_endog; n = rep%nobs
    info = SS_ERR_UNSUPPORTED
    if (fres%nobs_diffuse > 0) return
    allocate (K(m, p), Finv(p, p))
    Wa = 0.0_dp
    Watt = 0.0_dp
    ! omega_j,t+1 = L_t omega_jt, omega_j,j+1 = K_j
    do t = 1, n
      call conventional_gain(rep, fres, t, K, Finv, info)
      if (info /= SS_OK) return
      call period_terms(rep, fres, t, zfv, ZFZ, L, info)
      iz = tidx(size(rep%Z, 3), t)
      PZF = matmul(fres%P(:, :, t), matmul(transpose(rep%Z(:, :, iz)), Finv))
      A = eye(m) - matmul(PZF, rep%Z(:, :, iz))
      do j = 1, t - 1
        Watt(:, :, t, j) = matmul(A, Wa(:, :, t, j))
      end do
      Watt(:, :, t, t) = PZF
      if (t < n) then
        do j = 1, t - 1
          Wa(:, :, t + 1, j) = matmul(L, Wa(:, :, t, j))
        end do
        Wa(:, :, t + 1, t) = K
      end if
    end do
  end subroutine filtered_state_weights

  !> The Whittle relation (DK 4.6.3): the smoothed states solve the first-
  !> order conditions of the joint density of alpha and Y, which for
  !> t = n, ..., 2 give the backward recursion
  !>
  !>   T_t-1 alphahat_t-1 = alphahat_t - c_t-1 - W_t-1 [ Z_o' H_oo^-1 (y_o - d_o - Z_o alphahat_t)
  !>                        + T_t' W_t^-1 (alphahat_t+1 - c_t - T_t alphahat_t) ]
  !>
  !> with W = R Q R', started from alphahat_n = a_n|n and
  !> alphahat_n+1 = c_n + T_n alphahat_n (for the local level model this is
  !> alphahat_t-1 = 2 alphahat_t - alphahat_t+1 - q (y_t - alphahat_t)). It
  !> needs no stored filter output besides a_n|n, but needs T_t, W_t and H_t
  !> nonsingular, and as DK note it is numerically unstable: errors grow
  !> geometrically going back, so it is only usable for short series.
  subroutine whittle_smoother(rep, fres, alphahat, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    real(dp), intent(out), contiguous :: alphahat(:, :)     !< (m, n)
    integer, intent(out) :: info
    real(dp), allocatable :: anext(:), Winv(:, :), Wprev(:, :), Hinv(:, :), rhs(:, :), e(:), Tt(:, :)
    integer, allocatable :: idx(:)
    real(dp) :: logdet
    integer :: m, n, t, i, iz, ih, id, it, ic, itp, icp

    m = rep%k_states; n = rep%nobs
    alphahat(:, n) = fres%att(:, n)
    it = tidx(size(rep%T, 3), n)
    ic = tidx(size(rep%c, 2), n)
    anext = rep%c(:, ic) + matmul(rep%T(:, :, it), alphahat(:, n))
    do t = n, 2, -1
      iz = tidx(size(rep%Z, 3), t)
      ih = tidx(size(rep%H, 3), t)
      id = tidx(size(rep%d, 2), t)
      it = tidx(size(rep%T, 3), t)
      ic = tidx(size(rep%c, 2), t)
      itp = tidx(size(rep%T, 3), t - 1)
      icp = tidx(size(rep%c, 2), t - 1)
      ! T_t' W_t^-1 (alphahat_t+1 - c_t - T_t alphahat_t)
      Winv = state_noise(rep, t)
      call chol_inv(Winv, logdet, info)
      if (info /= SS_OK) return
      e = matmul(transpose(rep%T(:, :, it)), matmul(Winv, anext - rep%c(:, ic) &
                                                    - matmul(rep%T(:, :, it), alphahat(:, t))))
      ! + Z_o' H_oo^-1 (y_o - d_o - Z_o alphahat_t)
      idx = pack([(i, i=1, rep%k_endog)], .not. ieee_is_nan(rep%y(:, t)))
      if (size(idx) > 0) then
        Hinv = rep%H(idx, idx, ih)
        call chol_inv(Hinv, logdet, info)
        if (info /= SS_OK) return
        e = e + matmul(transpose(rep%Z(idx, :, iz)), matmul(Hinv, rep%y(idx, t) - rep%d(idx, id) &
                                                              - matmul(rep%Z(idx, :, iz), alphahat(:, t))))
      end if
      ! T_t-1 alphahat_t-1 = alphahat_t - c_t-1 - W_t-1 e
      Wprev = state_noise(rep, t - 1)
      rhs = reshape(alphahat(:, t) - rep%c(:, icp) - matmul(Wprev, e), [m, 1])
      Tt = rep%T(:, :, itp)
      call solve(Tt, rhs, info)
      if (info /= SS_OK) return
      anext = alphahat(:, t)
      alphahat(:, t - 1) = rhs(:, 1)
    end do
    info = SS_OK

  contains

    function state_noise(rep, t) result(W)
      type(ssm_rep_t), intent(in) :: rep
      integer, intent(in) :: t
      real(dp), allocatable :: W(:, :)
      integer :: ir, iq

      ir = tidx(size(rep%R, 3), t)
      iq = tidx(size(rep%Q, 3), t)
      W = matmul(rep%R(:, :, ir), matmul(rep%Q(:, :, iq), transpose(rep%R(:, :, ir))))
    end function state_noise
  end subroutine whittle_smoother

  !> Fast state smoother (DK 4.6): a backward pass for r_t only, then
  !>     alphahat_1 = a_1 + P_1 r_0  (+ P_inf,1 r^(1)_0 with diffuse states),
  !>     alphahat_t+1 = c_t + T_t alphahat_t + R_t Q_t R_t' r_t.
  !> Same alphahat as `state_smoother` without computing N_t or V_t.
  subroutine fast_state_smoother(rep, fres, alphahat, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    real(dp), intent(out), contiguous :: alphahat(:, :)    !< (m, n)
    integer, intent(out) :: info
    real(dp), allocatable :: r(:, :), r0(:), r1(:), u(:), vz(:), K0(:), K1(:), tmp(:), RQR(:, :)
    real(dp), allocatable :: RQ(:, :), Zo(:, :), vo(:), F1(:, :), F2(:, :), L0(:, :), L1(:, :)
    real(dp) :: v, Fs, Fi, k0r, k1r
    logical :: diffuse
    integer :: m, n, t, i, it, iz, ir, iq, ic

    info = SS_ERR_DIM
    if (fres%nobs /= rep%nobs .or. any(shape(alphahat) /= [rep%k_states, rep%nobs])) return
    info = SS_OK
    m = rep%k_states; n = rep%nobs
    allocate (r(m, 0:n), r0(m), r1(m), tmp(m), RQR(m, m), RQ(m, rep%k_posdef), source=0.0_dp)

    do t = n, 1, -1
      it = tidx(size(rep%T, 3), t)
      if (fres%method(t) == METHOD_DIFFUSE_MV) then
        ! r1 <- Z_o' F1 v_o + L0' r1 + L1' r0,  r0 <- L0' r0
        call diffuse_mv_gains(rep, fres, t, Zo, vo, F1, F2, L0, L1)
        r1 = matmul(transpose(Zo), matmul(F1, vo)) + matmul(transpose(L0), r1) &
             + matmul(transpose(L1), r0)
        r0 = matmul(transpose(L0), r0)
      else if (fres%method(t) /= METHOD_UNIVARIATE) then
        ! u = F^-1 v - K' r,  r_t-1 = Z' u + T' r_t  (and r1 <- T' r1 if diffuse)
        iz = tidx(size(rep%Z, 3), t)
        vz = merge(0.0_dp, fres%v(:, t), ieee_is_nan(fres%v(:, t)))
        u = matmul(fres%Finv(:, :, t), vz) - matmul(transpose(fres%K(:, :, t)), r0)
        call gemv('T', 1.0_dp, rep%T(:, :, it), r0, 0.0_dp, tmp)
        call gemv('T', 1.0_dp, rep%Z(:, :, iz), u, 1.0_dp, tmp)
        r0 = tmp
        if (t <= fres%nobs_diffuse) r1 = matmul(transpose(rep%T(:, :, it)), r1)
      else
        diffuse = t <= fres%nobs_diffuse
        call gemv('T', 1.0_dp, rep%T(:, :, it), r0, 0.0_dp, tmp)
        r0 = tmp
        if (diffuse) then
          call gemv('T', 1.0_dp, rep%T(:, :, it), r1, 0.0_dp, tmp)
          r1 = tmp
        end if
        do i = fres%uv_n(t), 1, -1
          associate (z => fres%uv_Z(i, :, t))
            v = fres%uv_v(i, t)
            Fs = fres%uv_Fstar(i, t)
            Fi = fres%uv_Finf(i, t)
            if (diffuse .and. Fi > rep%tol_diffuse) then
              ! r1 <- z v / Fi + L0' r1 + L1' r0,  r0 <- L0' r0
              K0 = fres%uv_Minf(:, i, t) / Fi
              K1 = (fres%uv_Mstar(:, i, t) - K0 * Fs) / Fi
              k0r = dot_product(K0, r1)
              k1r = dot_product(K1, r0)
              r1 = r1 + z * (v / Fi - k0r - k1r)
              r0 = r0 - z * dot_product(K0, r0)
            else if (Fs > merge(rep%tol_diffuse, 0.0_dp, diffuse)) then
              K0 = fres%uv_Mstar(:, i, t) / Fs
              r0 = r0 + z * (v / Fs - dot_product(K0, r0))
            end if
          end associate
        end do
        if (.not. diffuse) r1 = 0.0_dp
      end if
      r(:, t - 1) = r0
    end do

    alphahat(:, 1) = fres%a(:, 1)
    call gemv('N', 1.0_dp, fres%P(:, :, 1), r(:, 0), 1.0_dp, alphahat(:, 1))
    call gemv('N', 1.0_dp, fres%Pinf(:, :, 1), r1, 1.0_dp, alphahat(:, 1))
    do t = 1, n - 1
      it = tidx(size(rep%T, 3), t)
      ir = tidx(size(rep%R, 3), t)
      iq = tidx(size(rep%Q, 3), t)
      ic = tidx(size(rep%c, 2), t)
      call gemm('N', 'N', 1.0_dp, rep%R(:, :, ir), rep%Q(:, :, iq), 0.0_dp, RQ)
      call gemm('N', 'T', 1.0_dp, RQ, rep%R(:, :, ir), 0.0_dp, RQR)
      alphahat(:, t + 1) = rep%c(:, ic)
      call gemv('N', 1.0_dp, rep%T(:, :, it), alphahat(:, t), 1.0_dp, alphahat(:, t + 1))
      call gemv('N', 1.0_dp, RQR, r(:, t), 1.0_dp, alphahat(:, t + 1))
    end do
  end subroutine fast_state_smoother

  !> Classical fixed-interval (Rauch-Tung-Striebel) smoother (DK 4.6), from
  !> the filtered states:
  !>     J_t = P_t|t T_t' P_t+1^-1
  !>     alphahat_t = a_t|t + J_t (alphahat_t+1 - a_t+1)
  !>     V_t = P_t|t + J_t (V_t+1 - P_t+1) J_t'
  !> P_t+1^-1 is applied as a pseudo-inverse, which is exact when P_t+1 is
  !> singular. Not defined with diffuse states.
  subroutine classical_state_smoother(rep, fres, alphahat, V, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    real(dp), intent(out), contiguous :: alphahat(:, :)    !< (m, n)
    real(dp), intent(out), contiguous :: V(:, :, :)        !< (m, m, n)
    integer, intent(out) :: info
    real(dp), allocatable :: TP(:, :), J(:, :), D(:, :), JD(:, :)
    integer :: m, n, t, it

    info = SS_ERR_DIM
    if (fres%nobs /= rep%nobs .or. any(shape(alphahat) /= [rep%k_states, rep%nobs]) .or. &
        any(shape(V) /= [rep%k_states, rep%k_states, rep%nobs])) return
    info = SS_ERR_UNSUPPORTED
    if (fres%nobs_diffuse > 0) return
    info = SS_OK
    m = rep%k_states; n = rep%nobs
    allocate (TP(m, m), JD(m, m))

    alphahat(:, n) = fres%att(:, n)
    V(:, :, n) = fres%Ptt(:, :, n)
    do t = n - 1, 1, -1
      it = tidx(size(rep%T, 3), t)
      ! J' = P_t+1^+ T P_t|t
      call gemm('N', 'N', 1.0_dp, rep%T(:, :, it), fres%Ptt(:, :, t), 0.0_dp, TP)
      J = transpose(psd_solve(fres%P(:, :, t + 1), TP))
      alphahat(:, t) = fres%att(:, t)
      call gemv('N', 1.0_dp, J, alphahat(:, t + 1) - fres%a(:, t + 1), 1.0_dp, alphahat(:, t))
      D = V(:, :, t + 1) - fres%P(:, :, t + 1)
      call gemm('N', 'N', 1.0_dp, J, D, 0.0_dp, JD)
      V(:, :, t) = fres%Ptt(:, :, t)
      call gemm('N', 'T', 1.0_dp, JD, J, 1.0_dp, V(:, :, t))
      call symmetrize(V(:, :, t))
    end do
  end subroutine classical_state_smoother

  !> Two-filter formula for smoothing (DK 4.6.4). A backward information
  !> filter summarizes y_t..y_n alone by (i_t, I_t), with I_t the information
  !> matrix and i_t the information vector about alpha_t:
  !>
  !>   update:   I_t|t = I_t|t+1 + Z' H^-1 Z,  i_t|t = i_t|t+1 + Z' H^-1 (y_t - d_t)
  !>   predict:  G = I - I_t+1|t+1 R Q (I + R' I_t+1|t+1 R Q)^-1 R'
  !>             I_t|t+1 = T' G I_t+1|t+1 T,
  !>             i_t|t+1 = T' G (i_t+1|t+1 - I_t+1|t+1 c_t)
  !>
  !> (observed elements only, starting from I_n|n+1 = 0). It is combined with
  !> the forward prediction N(a_t, P_t) without inverting P_t:
  !>
  !>   alphahat_t = a_t + P_t (I + I_t|t P_t)^-1 (i_t|t - I_t|t a_t)
  !>   V_t = P_t (I + I_t|t P_t)^-1
  !>
  !> Needs H_t nonsingular over the observed elements; not defined with
  !> diffuse states.
  subroutine two_filter_smoother(rep, fres, alphahat, V, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    real(dp), intent(out), contiguous :: alphahat(:, :)    !< (m, n)
    real(dp), intent(out), contiguous :: V(:, :, :)        !< (m, m, n)
    integer, intent(out) :: info
    real(dp), allocatable :: Iinfo(:, :), ivec(:), Hinv(:, :), Zo(:, :), ZtHinv(:, :), A(:, :)
    real(dp), allocatable :: B(:, :), RQ(:, :), G(:, :), W(:, :), rhs(:, :)
    integer, allocatable :: idx(:)
    real(dp) :: logdet
    integer :: m, n, t, i, iz, ih, id, it, ir, iq, ic

    info = SS_ERR_DIM
    if (fres%nobs /= rep%nobs .or. any(shape(alphahat) /= [rep%k_states, rep%nobs]) .or. &
        any(shape(V) /= [rep%k_states, rep%k_states, rep%nobs])) return
    info = SS_ERR_UNSUPPORTED
    if (fres%nobs_diffuse > 0) return
    info = SS_OK
    m = rep%k_states; n = rep%nobs
    allocate (Iinfo(m, m), ivec(m), source=0.0_dp)

    do t = n, 1, -1
      iz = tidx(size(rep%Z, 3), t)
      ih = tidx(size(rep%H, 3), t)
      id = tidx(size(rep%d, 2), t)
      if (t < n) then
        ! Backward prediction from t+1 to t with the system matrices of t.
        it = tidx(size(rep%T, 3), t)
        ir = tidx(size(rep%R, 3), t)
        iq = tidx(size(rep%Q, 3), t)
        ic = tidx(size(rep%c, 2), t)
        RQ = matmul(rep%R(:, :, ir), rep%Q(:, :, iq))
        A = eye(rep%k_posdef) + matmul(transpose(rep%R(:, :, ir)), matmul(Iinfo, RQ))
        rhs = transpose(rep%R(:, :, ir))
        call solve(A, rhs, info)
        if (info /= SS_OK) return
        G = eye(m) - matmul(matmul(Iinfo, RQ), rhs)
        ivec = matmul(transpose(rep%T(:, :, it)), matmul(G, ivec - matmul(Iinfo, rep%c(:, ic))))
        Iinfo = matmul(transpose(rep%T(:, :, it)), matmul(matmul(G, Iinfo), rep%T(:, :, it)))
        call symmetrize(Iinfo)
      end if

      ! Update with the observed elements of y_t.
      idx = pack([(i, i=1, rep%k_endog)], .not. ieee_is_nan(rep%y(:, t)))
      if (size(idx) > 0) then
        Zo = rep%Z(idx, :, iz)
        Hinv = rep%H(idx, idx, ih)
        call chol_inv(Hinv, logdet, info)
        if (info /= SS_OK) return
        ZtHinv = matmul(transpose(Zo), Hinv)
        Iinfo = Iinfo + matmul(ZtHinv, Zo)
        ivec = ivec + matmul(ZtHinv, rep%y(idx, t) - rep%d(idx, id))
        call symmetrize(Iinfo)
      end if

      ! Combine with the forward prediction: W = (I + I P)^-1.
      B = eye(m) + matmul(Iinfo, fres%P(:, :, t))
      W = eye(m)
      call solve(B, W, info)
      if (info /= SS_OK) return
      alphahat(:, t) = fres%a(:, t) + matmul(fres%P(:, :, t), &
                                             matmul(W, ivec - matmul(Iinfo, fres%a(:, t))))
      V(:, :, t) = matmul(fres%P(:, :, t), W)
      call symmetrize(V(:, :, t))
    end do
  end subroutine two_filter_smoother

  !> Cov(alpha_t, alpha_j | Y_n) (DK 4.7). For j > t:
  !>     P_t L_t' L_t+1' ... L_j-1' (I - N_j-1 P_j),
  !> for j = t it is V_t, and for j < t the transpose of the (j, t) result.
  !> Periods t..j-1 must be past the diffuse period.
  subroutine smoothed_state_cov_between(rep, fres, sres, t, j, C, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    type(smoother_result_t), intent(in) :: sres
    integer, intent(in) :: t, j
    real(dp), intent(out), contiguous :: C(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: G(:, :), L(:, :), tmp(:, :), B(:, :)
    integer :: m, s, t0, t1

    info = SS_ERR_DIM
    if (min(t, j) < 1 .or. max(t, j) > rep%nobs) return
    info = SS_OK
    if (t == j) then
      C = sres%V(:, :, t)
      return
    end if
    t0 = min(t, j)
    t1 = max(t, j)
    m = rep%k_states
    allocate (L(m, m), tmp(m, m))

    ! G = P_t0 L_t0' ... L_t1-1'
    G = fres%P(:, :, t0)
    do s = t0, t1 - 1
      call innovation_transition(rep, fres, s, L, info)
      if (info /= SS_OK) return
      call gemm('N', 'T', 1.0_dp, G, L, 0.0_dp, tmp)
      G = tmp
    end do
    ! C = G (I - N_t1-1 P_t1)
    B = eye(m)
    call gemm('N', 'N', -1.0_dp, sres%N(:, :, t1 - 1), fres%P(:, :, t1), 1.0_dp, B)
    call gemm('N', 'N', 1.0_dp, G, B, 0.0_dp, tmp)
    if (t < j) then
      C = tmp
    else
      C = transpose(tmp)
    end if
  end subroutine smoothed_state_cov_between

  !> acov(:, :, t) = Cov(alpha_t+1, alpha_t | Y_n) for t = 1..n-1, NaN where
  !> it involves the diffuse period.
  subroutine smoothed_state_autocov(rep, fres, sres, acov, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    type(smoother_result_t), intent(in) :: sres
    real(dp), intent(out), contiguous :: acov(:, :, :)     !< (m, m, n-1)
    integer, intent(out) :: info
    integer :: t, stat

    info = SS_ERR_DIM
    if (size(acov, 3) /= rep%nobs - 1) return
    info = SS_OK
    do t = 1, rep%nobs - 1
      call smoothed_state_cov_between(rep, fres, sres, t + 1, t, acov(:, :, t), stat)
      if (stat /= SS_OK) acov(:, :, t) = ieee_value(1.0_dp, ieee_quiet_nan)
    end do
  end subroutine smoothed_state_autocov

  !> Weights of the smoothed state on the data (DK 4.8), extended to the state
  !> intercepts and the prior mean:
  !>
  !>     alphahat_t = sum_j W_tj (y_j - d_j) + sum_j C_tj c_j + A_t a1
  !>
  !> with W(k, i, t, j) the weight of y_j,i in alphahat_t,k, C(k, l, t, j) that
  !> of c_j,l and A(k, l, t) that of a1_l. Here a1 is the mean of alpha_1 from
  !> `rep%initial_state`, taken as given (for a stationary block it is itself
  !> a function of c_1). alphahat is affine in (y, c, a1), so
  !> each weight is the smoothed state of a unit input in the model with
  !> c = d = a1 = 0 and the same missing pattern. This holds in diffuse periods
  !> too (the diffuse directions of a1 get zero weight). It costs one smoother
  !> run per observation element, per state intercept element and per state,
  !> so O(n^2) in total; `j_list` restricts the observation and intercept
  !> weights to the periods listed (others are NaN). Weights of missing
  !> observations are zero.
  subroutine smoothed_state_weights(rep, W, C, A, info, j_list)
    type(ssm_rep_t), intent(in) :: rep
    real(dp), intent(out), contiguous :: W(:, :, :, :)   !< (m, p, n, n)
    real(dp), intent(out), contiguous :: C(:, :, :, :)   !< (m, m, n, n)
    real(dp), intent(out), contiguous :: A(:, :, :)      !< (m, m, n)
    integer, intent(out) :: info
    integer, intent(in), optional :: j_list(:)
    type(ssm_rep_t) :: base, unit
    real(dp), allocatable :: alphahat(:, :), a1(:), Pstar(:, :), Pinf(:, :)
    integer, allocatable :: js(:)
    integer :: p, m, n, i, j, k, l

    p = rep%k_endog; m = rep%k_states; n = rep%nobs
    info = SS_ERR_DIM
    if (any(shape(W) /= [m, p, n, n]) .or. any(shape(C) /= [m, m, n, n]) .or. &
        any(shape(A) /= [m, m, n])) return
    call rep%validate(info)
    if (info /= SS_OK) return

    if (present(j_list)) then
      js = j_list
    else
      js = [(j, j=1, n)]
    end if
    W = ieee_value(1.0_dp, ieee_quiet_nan)
    C = ieee_value(1.0_dp, ieee_quiet_nan)

    ! Base model: zero inputs, same missing pattern, and the initialization
    ! frozen as (0, P_star, P_inf) so the prior mean is an input of its own;
    ! c is time-varying so that one period's intercept can be set.
    allocate (a1(m), Pstar(m, m), Pinf(m, m))
    call rep%initial_state(a1, Pstar, Pinf, info)
    if (info /= SS_OK) return
    base = rep
    base%y = merge(rep%y, 0.0_dp, ieee_is_nan(rep%y))
    base%d = spread(spread(0.0_dp, 1, p), 2, 1)
    base%c = spread(spread(0.0_dp, 1, m), 2, n)
    call base%initialize_general(0.0_dp * a1, Pstar, Pinf)

    do k = 1, size(js)
      j = js(k)
      do i = 1, p
        if (ieee_is_nan(rep%y(i, j))) then
          W(:, i, :, j) = 0.0_dp
          cycle
        end if
        unit = base
        unit%y(i, j) = 1.0_dp
        call smoothed(unit, alphahat, info)
        if (info /= SS_OK) return
        W(:, i, :, j) = alphahat
      end do
      do l = 1, m
        unit = base
        unit%c(l, j) = 1.0_dp
        call smoothed(unit, alphahat, info)
        if (info /= SS_OK) return
        C(:, l, :, j) = alphahat
      end do
    end do

    do l = 1, m
      unit = base
      unit%a1(l) = 1.0_dp
      call smoothed(unit, alphahat, info)
      if (info /= SS_OK) return
      A(:, l, :) = alphahat
    end do

  contains

    subroutine smoothed(model, alphahat, stat)
      type(ssm_rep_t), intent(in) :: model
      real(dp), allocatable, intent(out) :: alphahat(:, :)
      integer, intent(out) :: stat
      type(filter_result_t) :: fres
      type(smoother_result_t) :: sres

      call kalman_filter(model, fres, stat)
      if (stat /= SS_OK) return
      call state_smoother(model, fres, sres, stat)
      if (stat /= SS_OK) return
      alphahat = sres%alphahat
    end subroutine smoothed
  end subroutine smoothed_state_weights
end module statespace_smoothing
