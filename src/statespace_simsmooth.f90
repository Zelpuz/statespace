!> Simulation from the model and the simulation smoother (DK 4.9, 5.5).
!>
!> `simulation_smoother` draws (alpha, eps, eta) from their distribution given
!> the data with the mean-correction method of Durbin and Koopman (2002):
!>
!>   1. Simulate (y+, alpha+, eps+, eta+) from the model with `simulate`.
!>   2. Smooth y* = y - y+ in the model with c = d = 0 and a1 = 0. Smoothing is
!>      affine, so this gives alphahat(y) - alphahat(y+) exactly.
!>   3. alpha~ = alpha+ + alphahat*, and likewise for eps and eta.
!>
!> The random input is standard normal variates, so draws are reproducible:
!>   u_init (m), u_eps (p, n), u_eta (r, n)
!> They are scaled by lower triangular square roots of P_star, H_t and Q_t
!> (the Cholesky factor when the matrix is positive definite), as in
!> statsmodels, so the same variates give the same draws. Diffuse directions
!> of alpha_1 get no draw: the exact diffuse smoother reproduces any shift in
!> them, so they cancel in step 3 (DK 5.5).
module statespace_simsmooth
  use statespace_kinds, only: dp, SS_OK, SS_ERR_DIM
  use statespace_linalg, only: gemv, ldl_psd, psd_solve, symmetrize
  use statespace_rep, only: ssm_rep_t, tidx
  use statespace_filter, only: filter_result_t, kalman_filter
  use statespace_smoother, only: smoother_result_t, state_smoother
  use statespace_smoothing, only: conventional_gain
  use statespace_kinds, only: SS_ERR_UNSUPPORTED
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  private

  public :: simsmooth_result_t, simulate, simulation_smoother, draw_standard_normal, &
            psd_sqrt
  public :: djs_measurement_disturbances, djs_state_disturbances

  type :: simsmooth_result_t
    real(dp), allocatable :: state(:, :)  !< (m, n) draw of alpha_t | Y_n
    real(dp), allocatable :: eps(:, :)    !< (p, n) draw of eps_t | Y_n
    real(dp), allocatable :: eta(:, :)    !< (r, n) draw of eta_t | Y_n
  end type simsmooth_result_t

contains

  !> Simulate from the model: alpha_1 = a1 + S(P_star) u_init,
  !> y_t = d + Z alpha_t + S(H_t) u_eps_t,
  !> alpha_t+1 = c + T alpha_t + R S(Q_t) u_eta_t,
  !> where S(A) is a lower triangular square root of A. Missing values in
  !> rep%y are ignored; y is fully simulated.
  subroutine simulate(rep, u_init, u_eps, u_eta, y, alpha, eps, eta, info)
    type(ssm_rep_t), intent(in) :: rep
    real(dp), intent(in), contiguous :: u_init(:), u_eps(:, :), u_eta(:, :)
    real(dp), intent(out), contiguous :: y(:, :)       !< (p, n)
    real(dp), intent(out), contiguous :: alpha(:, :)   !< (m, n)
    real(dp), intent(out), contiguous :: eps(:, :)     !< (p, n)
    real(dp), intent(out), contiguous :: eta(:, :)     !< (r, n)
    integer, intent(out) :: info
    real(dp), allocatable :: a1(:), Pstar(:, :), Pinf(:, :), S(:, :), SH(:, :), SQ(:, :)
    real(dp), allocatable :: anext(:)
    integer :: p, m, r, n, t, ih, iq, iz, it, ir, ic, id, last_ih, last_iq

    p = rep%k_endog; m = rep%k_states; r = rep%k_posdef; n = rep%nobs
    info = SS_ERR_DIM
    if (size(u_init) /= m .or. any(shape(u_eps) /= [p, n]) &
        .or. any(shape(u_eta) /= [r, n])) return
    call rep%validate(info)
    if (info /= SS_OK) return

    allocate (a1(m), Pstar(m, m), Pinf(m, m), anext(m))
    call rep%initial_state(a1, Pstar, Pinf, info)
    if (info /= SS_OK) return
    S = psd_sqrt(Pstar)
    alpha(:, 1) = a1
    call gemv('N', 1.0_dp, S, u_init, 1.0_dp, alpha(:, 1))

    last_ih = 0
    last_iq = 0
    do t = 1, n
      ih = tidx(size(rep%H, 3), t)
      iq = tidx(size(rep%Q, 3), t)
      if (ih /= last_ih) SH = psd_sqrt(rep%H(:, :, ih))
      if (iq /= last_iq) SQ = psd_sqrt(rep%Q(:, :, iq))
      last_ih = ih
      last_iq = iq

      iz = tidx(size(rep%Z, 3), t)
      it = tidx(size(rep%T, 3), t)
      ir = tidx(size(rep%R, 3), t)
      ic = tidx(size(rep%c, 2), t)
      id = tidx(size(rep%d, 2), t)

      call gemv('N', 1.0_dp, SH, u_eps(:, t), 0.0_dp, eps(:, t))
      y(:, t) = rep%d(:, id) + eps(:, t)
      call gemv('N', 1.0_dp, rep%Z(:, :, iz), alpha(:, t), 1.0_dp, y(:, t))

      call gemv('N', 1.0_dp, SQ, u_eta(:, t), 0.0_dp, eta(:, t))
      if (t < n) then
        anext = rep%c(:, ic)
        call gemv('N', 1.0_dp, rep%T(:, :, it), alpha(:, t), 1.0_dp, anext)
        call gemv('N', 1.0_dp, rep%R(:, :, ir), eta(:, t), 1.0_dp, anext)
        alpha(:, t + 1) = anext
      end if
    end do
  end subroutine simulate

  !> Draw alpha, eps and eta given the data (DK 4.9.2, mean correction).
  subroutine simulation_smoother(rep, u_init, u_eps, u_eta, sim, info)
    type(ssm_rep_t), intent(in) :: rep
    real(dp), intent(in), contiguous :: u_init(:), u_eps(:, :), u_eta(:, :)
    type(simsmooth_result_t), intent(out) :: sim
    integer, intent(out) :: info
    type(ssm_rep_t) :: star
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: yplus(:, :)

    allocate (yplus(rep%k_endog, rep%nobs), sim%state(rep%k_states, rep%nobs), &
              sim%eps(rep%k_endog, rep%nobs), sim%eta(rep%k_posdef, rep%nobs))
    call simulate(rep, u_init, u_eps, u_eta, yplus, sim%state, sim%eps, sim%eta, info)
    if (info /= SS_OK) return

    ! y* = y - y+ keeps the NaNs of y, so missing data carries over.
    star = rep
    star%y = rep%y - yplus
    star%c = 0.0_dp
    star%d = 0.0_dp
    star%a1 = 0.0_dp
    call kalman_filter(star, fres, info)
    if (info /= SS_OK) return
    call state_smoother(star, fres, sres, info)
    if (info /= SS_OK) return

    sim%state = sim%state + sres%alphahat
    sim%eps = sim%eps + sres%epshat
    sim%eta = sim%eta + sres%etahat
  end subroutine simulation_smoother

  !> de Jong-Shephard simulation of the measurement disturbances (DK 4.9.3):
  !> eps_t is drawn from N(epsbar_t, C_t), its distribution given the data
  !> and eps_t+1..eps_n, going backwards from t = n (4.83-4.88):
  !>
  !>   epsbar_t = H (F^-1 v - K' r~_t),     C_t = H - H (F^-1 + K' N~_t K) H
  !>   W~_t = H (F^-1 Z - K' N~_t L),       d_t = eps_t - epsbar_t
  !>   r~_t-1 = Z' F^-1 v - W~' C^-1 d_t + L' r~_t
  !>   N~_t-1 = Z' F^-1 Z + W~' C^-1 W~ + L' N~_t L
  !>
  !> with u_eps the N(0, 1) variates (C^-1 applied as a pseudo-inverse).
  !> Not defined with diffuse states.
  subroutine djs_measurement_disturbances(rep, fres, u_eps, eps, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    real(dp), intent(in), contiguous :: u_eps(:, :)    !< (p, n)
    real(dp), intent(out), contiguous :: eps(:, :)     !< (p, n)
    integer, intent(out) :: info
    real(dp), allocatable :: r(:), N(:, :), K(:, :), Finv(:, :), vz(:), L(:, :), &
                             C(:, :), W(:, :)
    real(dp), allocatable :: Zt(:, :), Ht(:, :), dvec(:), CW(:, :)
    integer :: t, p, m, iz, ih, it

    p = rep%k_endog; m = rep%k_states
    info = SS_ERR_UNSUPPORTED
    if (fres%nobs_diffuse > 0) return
    allocate (r(m), N(m, m), source=0.0_dp)
    allocate (K(m, p), Finv(p, p))
    do t = rep%nobs, 1, -1
      call conventional_gain(rep, fres, t, K, Finv, info)
      if (info /= SS_OK) return
      iz = tidx(size(rep%Z, 3), t)
      ih = tidx(size(rep%H, 3), t)
      it = tidx(size(rep%T, 3), t)
      Zt = rep%Z(:, :, iz)
      Ht = rep%H(:, :, ih)
      L = rep%T(:, :, it) - matmul(K, Zt)
      vz = merge(0.0_dp, fres%v(:, t), ieee_is_nan(fres%v(:, t)))
      eps(:, t) = matmul(Ht, matmul(Finv, vz) - matmul(transpose(K), r))
      C = Ht - matmul(Ht, matmul(Finv + matmul(transpose(K), matmul(N, K)), Ht))
      call symmetrize(C)
      dvec = matmul(psd_sqrt(C), u_eps(:, t))
      eps(:, t) = eps(:, t) + dvec
      W = matmul(Ht, matmul(Finv, Zt) - matmul(transpose(K), matmul(N, L)))
      CW = psd_solve(C, W)
      r = matmul(transpose(Zt), matmul(Finv, vz)) - matmul(transpose(CW), dvec) &
          + matmul(transpose(L), r)
      N = matmul(transpose(Zt), matmul(Finv, Zt)) + matmul(transpose(W), CW) &
          + matmul(transpose(L), matmul(N, L))
      call symmetrize(N)
    end do
  end subroutine djs_measurement_disturbances

  !> de Jong-Shephard simulation of the state disturbances and the state
  !> (DK 4.9.3): eta_t from N(etabar_t, Cbar_t) going backwards (4.89-4.91),
  !>
  !>   etabar_t = Q R' r~_t,      Cbar_t = Q - Q R' N~_t R Q,   W-_t = Q R' N~_t L
  !>   r~_t-1 = Z' F^-1 v - W-' Cbar^-1 d_t + L' r~_t
  !>   N~_t-1 = Z' F^-1 Z + W-' Cbar^-1 W- + L' N~_t L
  !>
  !> then the state forwards (4.92): alpha_t+1 = c_t + T_t alpha_t + R_t eta_t.
  !> DK's (4.92) starts from alpha_1 = a_1 + P_1 r~_0, the mean of alpha_1
  !> given Y and the eta draws; its variance given those is
  !> P_1 - P_1 N~_0 P_1, so alpha_1 is drawn from that distribution with the
  !> variates u_init (without it the state draws are too concentrated).
  !> Not defined with diffuse states.
  subroutine djs_state_disturbances(rep, fres, u_eta, u_init, eta, alpha, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    real(dp), intent(in), contiguous :: u_eta(:, :)    !< (r, n)
    real(dp), intent(in), contiguous :: u_init(:)      !< (m)
    real(dp), intent(out), contiguous :: eta(:, :)     !< (r, n)
    real(dp), intent(out), contiguous :: alpha(:, :)   !< (m, n)
    integer, intent(out) :: info
    real(dp), allocatable :: rtil(:), N(:, :), K(:, :), Finv(:, :), vz(:), L(:, :), &
                             C(:, :), W(:, :)
    real(dp), allocatable :: Zt(:, :), Qt(:, :), Rmat(:, :), dvec(:), CW(:, :)
    integer :: t, p, m, it, ir, ic, iz, iq

    p = rep%k_endog; m = rep%k_states
    info = SS_ERR_UNSUPPORTED
    if (fres%nobs_diffuse > 0) return
    allocate (rtil(m), N(m, m), source=0.0_dp)
    allocate (K(m, p), Finv(p, p))
    do t = rep%nobs, 1, -1
      call conventional_gain(rep, fres, t, K, Finv, info)
      if (info /= SS_OK) return
      iz = tidx(size(rep%Z, 3), t)
      ir = tidx(size(rep%R, 3), t)
      iq = tidx(size(rep%Q, 3), t)
      it = tidx(size(rep%T, 3), t)
      Zt = rep%Z(:, :, iz)
      Rmat = rep%R(:, :, ir)
      Qt = rep%Q(:, :, iq)
      L = rep%T(:, :, it) - matmul(K, Zt)
      vz = merge(0.0_dp, fres%v(:, t), ieee_is_nan(fres%v(:, t)))
      eta(:, t) = matmul(Qt, matmul(transpose(Rmat), rtil))
      C = Qt - matmul(Qt, matmul(transpose(Rmat), matmul(N, matmul(Rmat, Qt))))
      call symmetrize(C)
      dvec = matmul(psd_sqrt(C), u_eta(:, t))
      eta(:, t) = eta(:, t) + dvec
      W = matmul(Qt, matmul(transpose(Rmat), matmul(N, L)))
      CW = psd_solve(C, W)
      rtil = matmul(transpose(Zt), matmul(Finv, vz)) - matmul(transpose(CW), dvec) &
           + matmul(transpose(L), rtil)
      N = matmul(transpose(Zt), matmul(Finv, Zt)) + matmul(transpose(W), CW) &
          + matmul(transpose(L), matmul(N, L))
      call symmetrize(N)
    end do
    C = fres%P(:, :, 1) - matmul(fres%P(:, :, 1), matmul(N, fres%P(:, :, 1)))
    call symmetrize(C)
    alpha(:, 1) = fres%a(:, 1) + matmul(fres%P(:, :, 1), rtil) &
                  + matmul(psd_sqrt(C), u_init)
    do t = 1, rep%nobs - 1
      it = tidx(size(rep%T, 3), t)
      ir = tidx(size(rep%R, 3), t)
      ic = tidx(size(rep%c, 2), t)
      alpha(:, t + 1) = rep%c(:, ic) + matmul(rep%T(:, :, it), alpha(:, t)) &
                        + matmul(rep%R(:, :, ir), eta(:, t))
    end do
  end subroutine djs_state_disturbances

  !> Lower triangular S with S S' = A for symmetric positive semi-definite A:
  !> S = L D^(1/2) from A = L D L'. For positive definite A this is the
  !> Cholesky factor.
  function psd_sqrt(A) result(S)
    real(dp), intent(in), contiguous :: A(:, :)
    real(dp), allocatable :: S(:, :)
    real(dp), allocatable :: D(:)
    integer :: j

    allocate (S(size(A, 1), size(A, 1)), D(size(A, 1)))
    call ldl_psd(A, S, D)
    do j = 1, size(D)
      S(:, j) = S(:, j) * sqrt(D(j))
    end do
  end function psd_sqrt

  !> Fill x with independent N(0, 1) draws (Box-Muller on random_number).
  subroutine draw_standard_normal(x)
    real(dp), intent(out), contiguous :: x(:)
    real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp
    real(dp) :: u(2)
    integer :: i

    do i = 1, size(x), 2
      call random_number(u)
      u(1) = 1.0_dp - u(1)      ! in (0, 1], so the log is finite
      x(i) = sqrt(-2.0_dp * log(u(1))) * cos(twopi * u(2))
      if (i < size(x)) x(i + 1) = sqrt(-2.0_dp * log(u(1))) * sin(twopi * u(2))
    end do
  end subroutine draw_standard_normal
end module statespace_simsmooth
