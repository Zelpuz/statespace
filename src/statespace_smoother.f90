!> State smoother and disturbance smoother.
!>
!> Conventional periods use DK 4.4.4 and 4.5.3. Univariate periods use the
!> univariate exact initial state smoother (DK 5.3 with 6.4.4), processing
!> the observed elements of y_t in reverse; periods filtered with the
!> multivariate exact initial update use its smoother (DK 5.3). In diffuse
!> periods eps moments come from the smoothed state (see
!> `measurement_disturbance`).
module statespace_smoother
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use statespace_kinds, only: dp, SS_OK, SS_ERR_DIM
  use statespace_linalg, only: gemm, gemv, symmetrize, eye, psd_solve, chol_inv
  use statespace_rep, only: ssm_rep_t, tidx
  use statespace_filter, only: filter_result_t, METHOD_UNIVARIATE, METHOD_DIFFUSE_MV
  implicit none
  private

  public :: smoother_result_t, state_smoother, diffuse_mv_gains

  !> Smoother output. `r` and `N` use DK's indexing, 0..n with r_n = 0 and
  !> N_n = 0; in diffuse periods they hold r^(0) and N^(0).
  !>
  !> For a missing element of y_t, `epshat` and `epsvar` are the conditional
  !> moments of that element of eps_t given the observed data, which are
  !> nonzero when H_t has correlations. statsmodels reports 0 and H_t there.
  !> All eps moments are in the original coordinates of y_t, also in
  !> univariate periods (statsmodels reports transformed coordinates there).
  type :: smoother_result_t
    integer :: k_endog = 0, k_states = 0, k_posdef = 0, nobs = 0
    real(dp), allocatable :: alphahat(:, :)  !< (m, n) E[alpha_t | Y_n]
    real(dp), allocatable :: V(:, :, :)      !< (m, m, n) Var[alpha_t | Y_n]
    real(dp), allocatable :: r(:, :)         !< (m, 0:n)
    real(dp), allocatable :: N(:, :, :)      !< (m, m, 0:n)
    real(dp), allocatable :: epshat(:, :)    !< (p, n) E[eps_t | Y_n]
    real(dp), allocatable :: epsvar(:, :, :) !< (p, p, n) Var[eps_t | Y_n]
    real(dp), allocatable :: etahat(:, :)    !< (r, n) E[eta_t | Y_n]
    real(dp), allocatable :: etavar(:, :, :) !< (r, r, n) Var[eta_t | Y_n]
  end type smoother_result_t

  !> Backward recursion state. r0, N0 are DK's r_t, N_t; r1, N1, N2 are the
  !> diffuse parts r^(1), N^(1), N^(2), zero outside the diffuse period.
  type :: smoother_ws_t
    real(dp), allocatable :: r0(:), r1(:), N0(:, :), N1(:, :), N2(:, :)
    !> Scratch for the conventional step and the state disturbances,
    !> allocated once per run.
    real(dp), allocatable :: u(:), vz(:), r_prev(:), NL(:, :), PN(:, :), NK(:, :), HD(:, :), &
                             ZtFinv(:, :), N_prev(:, :), D(:, :), L(:, :), NRQ(:, :)
    !> R Q for the slices ir, iq of R and Q.
    real(dp), allocatable :: RQ(:, :)
    integer :: ir = 0, iq = 0
  end type smoother_ws_t

contains

  !> Run the smoother over the output of `kalman_filter`.
  subroutine state_smoother(rep, fres, sres, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    type(smoother_result_t), intent(out) :: sres
    integer, intent(out) :: info
    type(smoother_ws_t) :: ws
    integer :: p, m, r, n, t

    info = SS_ERR_DIM
    if (fres%nobs /= rep%nobs .or. fres%k_states /= rep%k_states .or. &
        fres%k_endog /= rep%k_endog) return
    info = SS_OK

    p = rep%k_endog; m = rep%k_states; r = rep%k_posdef; n = rep%nobs
    sres%k_endog = p; sres%k_states = m; sres%k_posdef = r; sres%nobs = n
    allocate (sres%alphahat(m, n), sres%V(m, m, n), sres%r(m, 0:n), sres%N(m, m, 0:n), &
              sres%epshat(p, n), sres%epsvar(p, p, n), sres%etahat(r, n), sres%etavar(r, r, n))
    allocate (ws%r0(m), ws%r1(m), ws%N0(m, m), ws%N1(m, m), ws%N2(m, m), source=0.0_dp)
    allocate (ws%u(p), ws%vz(p), ws%r_prev(m), ws%NL(m, m), ws%PN(m, m), ws%NK(m, p), &
              ws%HD(p, p), ws%ZtFinv(m, p), ws%N_prev(m, m), ws%D(p, p), ws%L(m, m), &
              ws%NRQ(m, r), ws%RQ(m, r))

    sres%r(:, n) = 0.0_dp
    sres%N(:, :, n) = 0.0_dp
    do t = n, 1, -1
      call state_disturbance(rep, t, ws, sres%etahat(:, t), sres%etavar(:, :, t))
      select case (fres%method(t))
      case (METHOD_UNIVARIATE)
        call univariate_step(rep, fres, t, ws, sres)
      case (METHOD_DIFFUSE_MV)
        call diffuse_mv_step(rep, fres, t, ws, sres)
      case default
        call conventional_step(rep, fres, t, ws, sres)
      end select
      sres%r(:, t - 1) = ws%r0
      sres%N(:, :, t - 1) = ws%N0
    end do
  end subroutine state_smoother

  !> eta_hat_t = Q R' r_t,  Var = Q - Q R' N_t R Q  (DK 4.5.3, also in the
  !> diffuse period, DK 5.3), with r_t, N_t before the step at t.
  subroutine state_disturbance(rep, t, ws, etahat, etavar)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: t
    type(smoother_ws_t), intent(inout) :: ws
    real(dp), intent(out), contiguous :: etahat(:), etavar(:, :)
    integer :: ir, iq

    ir = tidx(size(rep%R, 3), t)
    iq = tidx(size(rep%Q, 3), t)
    if (ir /= ws%ir .or. iq /= ws%iq) then
      call gemm('N', 'N', 1.0_dp, rep%R(:, :, ir), rep%Q(:, :, iq), 0.0_dp, ws%RQ)
      ws%ir = ir
      ws%iq = iq
    end if
    call gemv('T', 1.0_dp, ws%RQ, ws%r0, 0.0_dp, etahat)
    call gemm('N', 'N', 1.0_dp, ws%N0, ws%RQ, 0.0_dp, ws%NRQ)
    etavar = rep%Q(:, :, iq)
    call gemm('T', 'N', -1.0_dp, ws%RQ, ws%NRQ, 1.0_dp, etavar)
    call symmetrize(etavar)
  end subroutine state_disturbance

  !> Conventional backward step at t (DK 4.4.4, 4.5.3):
  !>
  !>     u_t     = F_t^-1 v_t - K_t' r_t,          D_t = F_t^-1 + K_t' N_t K_t
  !>     r_t-1   = Z_t' u_t + T_t' r_t              (= Z' F^-1 v + L' r)
  !>     N_t-1   = Z_t' F_t^-1 Z_t + L_t' N_t L_t,   L_t = T_t - K_t Z_t
  !>     alphahat_t = a_t + P_t r_t-1,  V_t = P_t - P_t N_t-1 P_t
  !>     epshat_t = H_t u_t,             Var = H_t - H_t D_t H_t
  !>
  !> In a diffuse period this is the F_inf = 0 case of DK 5.3, which adds
  !>     r^(1)_t-1 = T' r^(1)_t,  N^(1)_t-1 = T' N^(1)_t L,  N^(2)_t-1 = T' N^(2)_t T
  !> and the P_inf terms of alphahat and V.
  subroutine conventional_step(rep, fres, t, ws, sres)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    integer, intent(in) :: t
    type(smoother_ws_t), intent(inout) :: ws
    type(smoother_result_t), intent(inout) :: sres
    integer :: p, m, iz, ih, it

    p = rep%k_endog; m = rep%k_states
    iz = tidx(size(rep%Z, 3), t)
    ih = tidx(size(rep%H, 3), t)
    it = tidx(size(rep%T, 3), t)

    ! u = F^-1 v - K' r  (missing elements of v contribute nothing)
    ws%vz = merge(0.0_dp, fres%v(:, t), ieee_is_nan(fres%v(:, t)))
    call gemv('N', 1.0_dp, fres%Finv(:, :, t), ws%vz, 0.0_dp, ws%u)
    call gemv('T', -1.0_dp, fres%K(:, :, t), ws%r0, 1.0_dp, ws%u)
    call gemv('N', 1.0_dp, rep%H(:, :, ih), ws%u, 0.0_dp, sres%epshat(:, t))

    ! Var[eps_t] = H - H D H,  D = F^-1 + K' N_t K
    call gemm('N', 'N', 1.0_dp, ws%N0, fres%K(:, :, t), 0.0_dp, ws%NK)
    ws%D = fres%Finv(:, :, t)
    call gemm('T', 'N', 1.0_dp, fres%K(:, :, t), ws%NK, 1.0_dp, ws%D)
    call gemm('N', 'N', 1.0_dp, rep%H(:, :, ih), ws%D, 0.0_dp, ws%HD)
    sres%epsvar(:, :, t) = rep%H(:, :, ih)
    call gemm('N', 'N', -1.0_dp, ws%HD, rep%H(:, :, ih), 1.0_dp, sres%epsvar(:, :, t))
    call symmetrize(sres%epsvar(:, :, t))

    ! r_t-1 = Z' u + T' r_t
    call gemv('T', 1.0_dp, rep%T(:, :, it), ws%r0, 0.0_dp, ws%r_prev)
    call gemv('T', 1.0_dp, rep%Z(:, :, iz), ws%u, 1.0_dp, ws%r_prev)

    ! N_t-1 = Z' F^-1 Z + L' N_t L
    ws%L = rep%T(:, :, it)
    call gemm('N', 'N', -1.0_dp, fres%K(:, :, t), rep%Z(:, :, iz), 1.0_dp, ws%L)
    call gemm('T', 'N', 1.0_dp, rep%Z(:, :, iz), fres%Finv(:, :, t), 0.0_dp, ws%ZtFinv)
    call gemm('N', 'N', 1.0_dp, ws%ZtFinv, rep%Z(:, :, iz), 0.0_dp, ws%N_prev)
    call gemm('N', 'N', 1.0_dp, ws%N0, ws%L, 0.0_dp, ws%NL)
    call gemm('T', 'N', 1.0_dp, ws%L, ws%NL, 1.0_dp, ws%N_prev)
    call symmetrize(ws%N_prev)
    ws%r0 = ws%r_prev
    ws%N0 = ws%N_prev

    if (t > fres%nobs_diffuse) then
      ! alphahat_t = a_t + P_t r_t-1,  V_t = P_t - P_t N_t-1 P_t
      sres%alphahat(:, t) = fres%a(:, t)
      call gemv('N', 1.0_dp, fres%P(:, :, t), ws%r0, 1.0_dp, sres%alphahat(:, t))
      call gemm('N', 'N', 1.0_dp, fres%P(:, :, t), ws%N0, 0.0_dp, ws%PN)
      sres%V(:, :, t) = fres%P(:, :, t)
      call gemm('N', 'N', -1.0_dp, ws%PN, fres%P(:, :, t), 1.0_dp, sres%V(:, :, t))
      call symmetrize(sres%V(:, :, t))
      return
    end if

    ! Diffuse period with F_inf = 0.
    ws%r1 = matmul(transpose(rep%T(:, :, it)), ws%r1)
    ws%N1 = matmul(transpose(rep%T(:, :, it)), matmul(ws%N1, ws%L))
    ws%N2 = matmul(transpose(rep%T(:, :, it)), matmul(ws%N2, rep%T(:, :, it)))
    call symmetrize(ws%N2)
    call diffuse_state(fres, t, ws, sres)
    call measurement_disturbance(rep, t, sres%alphahat(:, t), sres%V(:, :, t), &
                                 sres%epshat(:, t), sres%epsvar(:, :, t))
  end subroutine conventional_step

  !> Multivariate exact initial smoothing step for a period with F_inf
  !> nonsingular over the observed elements o (DK 5.3):
  !>
  !>   K0 = T M_inf F1,  K1 = T (M_star F1 + M_inf F2),  L0 = T - K0 Z_o,  L1 = -K1 Z_o
  !>   r0 <- L0' r0,   r1 <- Z_o' F1 v_o + L0' r1 + L1' r0
  !>   N0 <- L0' N0 L0
  !>   N1 <- Z_o' F1 Z_o + L0' N1 L0 + L1' N0 L0
  !>   N2 <- Z_o' F2 Z_o + L0' N2 L0 + L0' N1 L1 + L1' N1' L0 + L1' N0 L1
  !>
  !> with F1 = F_inf^-1, F2 = -F1 F_star F1, M = P Z_o'.
  subroutine diffuse_mv_step(rep, fres, t, ws, sres)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    integer, intent(in) :: t
    type(smoother_ws_t), intent(inout) :: ws
    type(smoother_result_t), intent(inout) :: sres
    real(dp), allocatable :: Zo(:, :), vo(:), F1(:, :), F2(:, :), L0(:, :), L1(:, :)
    real(dp), allocatable :: r0(:), r1(:), N0(:, :), N1(:, :), N2(:, :)

    call diffuse_mv_gains(rep, fres, t, Zo, vo, F1, F2, L0, L1)
    r0 = ws%r0; r1 = ws%r1; N0 = ws%N0; N1 = ws%N1; N2 = ws%N2
    ws%r0 = matmul(transpose(L0), r0)
    ws%r1 = matmul(transpose(Zo), matmul(F1, vo)) + matmul(transpose(L0), r1) &
            + matmul(transpose(L1), r0)
    ws%N0 = matmul(transpose(L0), matmul(N0, L0))
    ws%N1 = matmul(transpose(Zo), matmul(F1, Zo)) + matmul(transpose(L0), matmul(N1, L0)) &
            + matmul(transpose(L1), matmul(N0, L0))
    ws%N2 = matmul(transpose(Zo), matmul(F2, Zo)) + matmul(transpose(L0), matmul(N2, L0)) &
            + matmul(transpose(L0), matmul(N1, L1)) + matmul(transpose(L1), matmul(transpose(N1), L0)) &
            + matmul(transpose(L1), matmul(N0, L1))
    call symmetrize(ws%N0)
    call symmetrize(ws%N2)

    call diffuse_state(fres, t, ws, sres)
    call measurement_disturbance(rep, t, sres%alphahat(:, t), sres%V(:, :, t), &
                                 sres%epshat(:, t), sres%epsvar(:, :, t))
  end subroutine diffuse_mv_step

  !> Quantities of the multivariate exact initial smoother at t (F_inf
  !> nonsingular over the observed elements o): Z_o, v_o, F1 = F_inf^-1,
  !> F2 = -F1 F_star F1, L0 = T - K0 Z_o and L1 = -K1 Z_o with
  !> K0 = T P_inf Z_o' F1, K1 = T (P_star Z_o' F1 + P_inf Z_o' F2).
  subroutine diffuse_mv_gains(rep, fres, t, Zo, vo, F1, F2, L0, L1)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    integer, intent(in) :: t
    real(dp), allocatable, intent(out) :: Zo(:, :), vo(:), F1(:, :), F2(:, :), L0(:, :), L1(:, :)
    real(dp), allocatable :: Tt(:, :), K0(:, :), K1(:, :)
    integer, allocatable :: idx(:)
    real(dp) :: logdet
    integer :: i, info, it, iz

    idx = pack([(i, i=1, rep%k_endog)], .not. ieee_is_nan(rep%y(:, t)))
    it = tidx(size(rep%T, 3), t)
    iz = tidx(size(rep%Z, 3), t)
    Tt = rep%T(:, :, it)
    Zo = rep%Z(idx, :, iz)
    vo = fres%v(idx, t)
    F1 = fres%Finf(idx, idx, t)
    call chol_inv(F1, logdet, info)   ! nonsingular: checked by the filter
    F2 = -matmul(F1, matmul(fres%F(idx, idx, t), F1))
    K0 = matmul(Tt, matmul(fres%Pinf(:, :, t), matmul(transpose(Zo), F1)))
    K1 = matmul(Tt, matmul(fres%P(:, :, t), matmul(transpose(Zo), F1)) &
                + matmul(fres%Pinf(:, :, t), matmul(transpose(Zo), F2)))
    L0 = Tt - matmul(K0, Zo)
    L1 = -matmul(K1, Zo)
  end subroutine diffuse_mv_gains

  !> alphahat = a + P_star r0 + P_inf r1 and
  !> V = P_star - P_star N0 P_star - P_inf N1 P_star - (P_inf N1 P_star)' - P_inf N2 P_inf.
  subroutine diffuse_state(fres, t, ws, sres)
    type(filter_result_t), intent(in) :: fres
    integer, intent(in) :: t
    type(smoother_ws_t), intent(in) :: ws
    type(smoother_result_t), intent(inout) :: sres
    real(dp), allocatable :: B(:, :)

    associate (Ps => fres%P(:, :, t), Pi => fres%Pinf(:, :, t))
      sres%alphahat(:, t) = fres%a(:, t) + matmul(Ps, ws%r0) + matmul(Pi, ws%r1)
      B = matmul(Pi, matmul(ws%N1, Ps))
      sres%V(:, :, t) = Ps - matmul(Ps, matmul(ws%N0, Ps)) - B - transpose(B) &
                        - matmul(Pi, matmul(ws%N2, Pi))
      call symmetrize(sres%V(:, :, t))
    end associate
  end subroutine diffuse_state

  !> Univariate backward step at t. First r and N move from alpha_t+1 to
  !> alpha_t (r <- T_t' r, N <- T_t' N T_t), then the elements i = n_o..1 are
  !> processed with L0 = I - K0 z', L1 = -K1 z':
  !>
  !>   F_inf > 0:  r1 <- z v / F_inf + L0' r1 + L1' r0,   r0 <- L0' r0
  !>               N2 <- z z' F2 + L0' N2 L0 + L0' N1 L1 + L1' N1' L0 + L1' N0 L1
  !>               N1 <- z z' / F_inf + L0' N1 L0 + L1' N0 L0,   N0 <- L0' N0 L0
  !>               (K0 = M_inf / F_inf, K1 = (M_star - K0 F_star) / F_inf,
  !>                F2 = -F_star / F_inf^2)
  !>   otherwise:  r0 <- z v / F_star + L0' r0,  N0 <- z z' / F_star + L0' N0 L0,
  !>               N1 <- N1 L0   (K0 = M_star / F_star)
  !>
  !> Then alphahat = a + P_star r0 + P_inf r1 and
  !> V = P_star - P_star N0 P_star - P_inf N1 P_star - (P_inf N1 P_star)' - P_inf N2 P_inf.
  subroutine univariate_step(rep, fres, t, ws, sres)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    integer, intent(in) :: t
    type(smoother_ws_t), intent(inout) :: ws
    type(smoother_result_t), intent(inout) :: sres
    real(dp), allocatable :: z(:), K0(:), K1(:), L0(:, :), L1(:, :), zz(:, :), A(:, :), B(:, :)
    real(dp), allocatable :: N0(:, :), N1(:, :), N2(:, :), r0(:), r1(:), PN(:, :)
    real(dp) :: v, Fs, Fi
    logical :: diffuse
    integer :: m, i, j, it

    m = rep%k_states
    it = tidx(size(rep%T, 3), t)
    diffuse = t <= fres%nobs_diffuse
    allocate (K0(m), K1(m), L1(m, m), zz(m, m), A(m, m), B(m, m), PN(m, m))

    ! Move from alpha_t+1 to alpha_t.
    call gemv('T', 1.0_dp, rep%T(:, :, it), ws%r0, 0.0_dp, K0)
    ws%r0 = K0
    call congruence(rep%T(:, :, it), ws%N0, A)
    ws%N0 = A
    if (diffuse) then
      call gemv('T', 1.0_dp, rep%T(:, :, it), ws%r1, 0.0_dp, K0)
      ws%r1 = K0
      ! N1 is not symmetric: T' N1 T
      call gemm('T', 'N', 1.0_dp, rep%T(:, :, it), ws%N1, 0.0_dp, A)
      call gemm('N', 'N', 1.0_dp, A, rep%T(:, :, it), 0.0_dp, ws%N1)
      call congruence(rep%T(:, :, it), ws%N2, A)
      ws%N2 = A
    end if

    do i = fres%uv_n(t), 1, -1
      z = fres%uv_Z(i, :, t)
      v = fres%uv_v(i, t)
      Fs = fres%uv_Fstar(i, t)
      Fi = fres%uv_Finf(i, t)
      do j = 1, m
        zz(:, j) = z * z(j)
      end do
      if (diffuse .and. Fi > rep%tol_diffuse) then
        K0 = fres%uv_Minf(:, i, t) / Fi
        K1 = (fres%uv_Mstar(:, i, t) - K0 * Fs) / Fi
        L0 = eye(m)
        do j = 1, m
          L0(:, j) = L0(:, j) - K0 * z(j)
          L1(:, j) = -K1 * z(j)
        end do
        r0 = ws%r0; r1 = ws%r1; N0 = ws%N0; N1 = ws%N1; N2 = ws%N2

        ! r1 <- z v / F_inf + L0' r1 + L1' r0,  r0 <- L0' r0
        ws%r1 = z * (v / Fi)
        call gemv('T', 1.0_dp, L0, r1, 1.0_dp, ws%r1)
        call gemv('T', 1.0_dp, L1, r0, 1.0_dp, ws%r1)
        call gemv('T', 1.0_dp, L0, r0, 0.0_dp, ws%r0)

        ! N2 <- z z' F2 + L0' N2 L0 + L0' N1 L1 + L1' N1' L0 + L1' N0 L1
        ws%N2 = zz * (-Fs / Fi**2)
        call add_lt_a_l(L0, N2, L0, ws%N2)
        call add_lt_a_l(L0, N1, L1, ws%N2)
        call add_lt_a_l(L1, transpose(N1), L0, ws%N2)
        call add_lt_a_l(L1, N0, L1, ws%N2)
        call symmetrize(ws%N2)

        ! N1 <- z z' / F_inf + L0' N1 L0 + L1' N0 L0
        ws%N1 = zz / Fi
        call add_lt_a_l(L0, N1, L0, ws%N1)
        call add_lt_a_l(L1, N0, L0, ws%N1)

        ! N0 <- L0' N0 L0
        ws%N0 = 0.0_dp
        call add_lt_a_l(L0, N0, L0, ws%N0)
        call symmetrize(ws%N0)
      else if (Fs > merge(rep%tol_diffuse, 0.0_dp, diffuse)) then
        K0 = fres%uv_Mstar(:, i, t) / Fs
        L0 = eye(m)
        do j = 1, m
          L0(:, j) = L0(:, j) - K0 * z(j)
        end do
        r0 = ws%r0
        ws%r0 = z * (v / Fs)
        call gemv('T', 1.0_dp, L0, r0, 1.0_dp, ws%r0)
        N0 = ws%N0
        ws%N0 = zz / Fs
        call add_lt_a_l(L0, N0, L0, ws%N0)
        call symmetrize(ws%N0)
        if (diffuse) then
          call gemm('N', 'N', 1.0_dp, ws%N1, L0, 0.0_dp, A)
          ws%N1 = A
        end if
      end if
    end do

    ! alphahat = a + P_star r0 + P_inf r1
    associate (Ps => fres%P(:, :, t), Pi => fres%Pinf(:, :, t))
      sres%alphahat(:, t) = fres%a(:, t)
      call gemv('N', 1.0_dp, Ps, ws%r0, 1.0_dp, sres%alphahat(:, t))
      call gemm('N', 'N', 1.0_dp, Ps, ws%N0, 0.0_dp, PN)
      sres%V(:, :, t) = Ps
      call gemm('N', 'N', -1.0_dp, PN, Ps, 1.0_dp, sres%V(:, :, t))
      if (diffuse) then
        call gemv('N', 1.0_dp, Pi, ws%r1, 1.0_dp, sres%alphahat(:, t))
        ! B = P_inf N1 P_star;  V -= B + B' + P_inf N2 P_inf
        call gemm('N', 'N', 1.0_dp, Pi, ws%N1, 0.0_dp, A)
        call gemm('N', 'N', 1.0_dp, A, Ps, 0.0_dp, B)
        sres%V(:, :, t) = sres%V(:, :, t) - B - transpose(B)
        call gemm('N', 'N', 1.0_dp, Pi, ws%N2, 0.0_dp, A)
        call gemm('N', 'N', -1.0_dp, A, Pi, 1.0_dp, sres%V(:, :, t))
      end if
      call symmetrize(sres%V(:, :, t))
    end associate
    if (.not. diffuse) then
      ws%r1 = 0.0_dp
      ws%N1 = 0.0_dp
      ws%N2 = 0.0_dp
    end if

    call measurement_disturbance(rep, t, sres%alphahat(:, t), sres%V(:, :, t), &
                                 sres%epshat(:, t), sres%epsvar(:, :, t))
  end subroutine univariate_step

  !> eps moments from the smoothed state. For the observed elements o,
  !> eps_o = y_o - d_o - Z_o alpha_t exactly, so
  !>     epshat_o = y_o - d_o - Z_o alphahat_t,   Var = Z_o V_t Z_o'.
  !> Missing elements m depend on the data only through eps_o:
  !>     epshat_m = G epshat_o,   G = H_mo H_oo^+,
  !>     Var_mm = H_mm - G H_om + G Var_oo G',   Var_mo = G Var_oo.
  subroutine measurement_disturbance(rep, t, alphahat, V, epshat, epsvar)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: t
    real(dp), intent(in), contiguous :: alphahat(:), V(:, :)
    real(dp), intent(out), contiguous :: epshat(:), epsvar(:, :)
    real(dp), allocatable :: Zo(:, :), ZV(:, :), Voo(:, :), G(:, :), GV(:, :), eo(:)
    logical, allocatable :: obs(:)
    integer, allocatable :: io(:), im(:)
    integer :: p, i, iz, ih, id

    p = rep%k_endog
    iz = tidx(size(rep%Z, 3), t)
    ih = tidx(size(rep%H, 3), t)
    id = tidx(size(rep%d, 2), t)
    obs = .not. ieee_is_nan(rep%y(:, t))
    io = pack([(i, i=1, p)], obs)
    im = pack([(i, i=1, p)], .not. obs)

    epshat = 0.0_dp
    epsvar = rep%H(:, :, ih)
    if (size(io) == 0) return

    Zo = rep%Z(io, :, iz)
    eo = rep%y(io, t) - rep%d(io, id)
    call gemv('N', -1.0_dp, Zo, alphahat, 1.0_dp, eo)
    epshat(io) = eo
    allocate (ZV(size(io), rep%k_states), Voo(size(io), size(io)))
    call gemm('N', 'N', 1.0_dp, Zo, V, 0.0_dp, ZV)
    call gemm('N', 'T', 1.0_dp, ZV, Zo, 0.0_dp, Voo)
    call symmetrize(Voo)
    epsvar(io, io) = Voo
    if (size(im) == 0) return

    ! G' = H_oo^+ H_om
    G = transpose(psd_solve(rep%H(io, io, ih), rep%H(io, im, ih)))
    epshat(im) = matmul(G, eo)
    GV = matmul(G, Voo)
    epsvar(im, io) = GV
    epsvar(io, im) = transpose(GV)
    epsvar(im, im) = rep%H(im, im, ih) - matmul(G, rep%H(io, im, ih)) + matmul(GV, transpose(G))
  end subroutine measurement_disturbance

  !> S <- T' S T for symmetric S.
  subroutine congruence(T, S, out)
    real(dp), intent(in), contiguous :: T(:, :), S(:, :)
    real(dp), intent(out), contiguous :: out(:, :)
    real(dp), allocatable :: ST(:, :)

    allocate (ST(size(S, 1), size(T, 2)))
    call gemm('N', 'N', 1.0_dp, S, T, 0.0_dp, ST)
    call gemm('T', 'N', 1.0_dp, T, ST, 0.0_dp, out)
    call symmetrize(out)
  end subroutine congruence

  !> C <- C + L1' A L2.
  subroutine add_lt_a_l(L1, A, L2, C)
    real(dp), intent(in), contiguous :: L1(:, :), A(:, :), L2(:, :)
    real(dp), intent(inout), contiguous :: C(:, :)
    real(dp), allocatable :: AL(:, :)

    allocate (AL(size(A, 1), size(L2, 2)))
    call gemm('N', 'N', 1.0_dp, A, L2, 0.0_dp, AL)
    call gemm('T', 'N', 1.0_dp, L1, AL, 1.0_dp, C)
  end subroutine add_lt_a_l
end module statespace_smoother
