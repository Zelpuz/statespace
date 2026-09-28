!> Kalman filter for the linear Gaussian state space model.
!>
!> Two entry points share the step routines:
!>   * `kalman_filter` stores the full filter output needed by the smoother.
!>   * `loglike` computes only the log likelihood and stores no history. MLE
!>     should call this one.
!>
!> Each period is processed in one of two ways:
!>   * Conventional (DK 4.3, eq. 4.24): y_t as a vector.
!>   * Univariate treatment (DK 6.4): the observed elements of y_t one at a
!>     time. When H_t is not diagonal, the observed block is first transformed
!>     with H_oo = L D L' (DK 6.4.3), y* = L^-1 (y - d), Z* = L^-1 Z. Periods
!>     with a diffuse state use this with the exact initial Kalman filter
!>     (DK 5.2) and contribute the diffuse log likelihood (DK 7.2.2). With
!>     FILTER_UNIVARIATE every period uses it.
!>
!> Missing observations (NaN in y) are handled as in DK 4.10: only the
!> observed elements of y_t enter the update. In conventional periods F_t^-1
!> is stored zero-padded over the missing rows and columns, which makes the
!> full-size update and smoother formulas exact.
module statespace_filter
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use statespace_kinds, only: dp, log2pi, SS_OK, SS_ERR_UNSUPPORTED, &
                              SS_ERR_NOT_CONVERGED, SS_ERR_NOT_PD
  use statespace_linalg, only: gemm, gemv, chol_inv, symmetrize, ldl_psd, &
                               solve_unit_lower, is_diagonal, solve, eye
  use statespace_rep, only: ssm_rep_t, tidx, FILTER_UNIVARIATE, DIFFUSE_MULTIVARIATE
  implicit none
  private

  public :: filter_result_t, kalman_filter, loglike, loglike_concentrated, &
            marginal_correction
  public :: steady_state

  !> How each period was filtered (filter_result_t%method).
  integer, parameter, public :: METHOD_CONVENTIONAL = 0  !< vector update (DK 4.3);
                                                        !! with diffuse states, the
                                                        !! F_inf = 0 case
  integer, parameter, public :: METHOD_UNIVARIATE = 1    !< univariate treatment
                                                        !! (DK 6.4)
  integer, parameter, public :: METHOD_DIFFUSE_MV = 2    !< multivariate exact
                                                        !! initial update, F_inf
                                                        !! nonsingular (DK 5.2)

  !> Off-diagonal elements of H_t below this are treated as zero, as in
  !> statsmodels.
  real(dp), parameter :: tol_diagonal = 1.0e-9_dp

  !> Kalman filter output. `a` and `P` hold the predicted a_t and P_t (the
  !> P_star part in diffuse periods) for t = 1..n+1, `Pinf` the diffuse part
  !> P_inf,t; `att`, `Ptt` the filtered a_t|t, P_t|t.
  !>
  !> All per-period quantities are in the original coordinates of y_t:
  !> yhat = d + Z a, v = y - yhat, F = Z P Z' + H, Finf = Z P_inf Z'.
  !> For missing elements v is NaN and F is still the full forecast variance.
  !>
  !> Conventional periods (method(t) = METHOD_CONVENTIONAL) also store Finv
  !> (F^-1 over the observed block, zero-padded) and K = T P Z' F^-1 (zero
  !> columns for missing elements). Other periods have no conventional gain,
  !> so Finv and K are NaN there. Univariate periods store the per-element
  !> quantities of the univariate filter in the `uv_` arrays, in transformed
  !> coordinates and in the order the observed elements were processed.
  type :: filter_result_t
    integer :: k_endog = 0, k_states = 0, nobs = 0
    integer :: nobs_diffuse = 0       !< periods with a diffuse state (DK's d)
    integer :: t_steady = 0           !< first period in the steady state
                                      !! (0: none; see tol_steady)
    integer, allocatable :: method(:)  !< (n) METHOD_* used for each period
    integer :: k_diffuse = 0          !< rank of P_inf,1
    real(dp), allocatable :: a(:, :)        !< (m, n+1)
    real(dp), allocatable :: P(:, :, :)     !< (m, m, n+1)
    real(dp), allocatable :: Pinf(:, :, :)  !< (m, m, n+1), zero after the
                                            !! diffuse period
    real(dp), allocatable :: att(:, :)      !< (m, n)
    real(dp), allocatable :: Ptt(:, :, :)   !< (m, m, n)
    real(dp), allocatable :: yhat(:, :)     !< (p, n) one-step forecasts d + Z a_t
    real(dp), allocatable :: v(:, :)        !< (p, n) one-step forecast errors
    real(dp), allocatable :: F(:, :, :)     !< (p, p, n) forecast error variances
    real(dp), allocatable :: Finf(:, :, :)  !< (p, p, n) diffuse part Z P_inf Z'
    real(dp), allocatable :: Finv(:, :, :)  !< (p, p, n)
    real(dp), allocatable :: K(:, :, :)     !< (m, p, n)
    real(dp), allocatable :: llf_obs(:)     !< (n) log likelihood contributions
    real(dp) :: llf = 0.0_dp                !< sum of llf_obs after the burn-in

    ! Univariate periods: element i = 1..uv_n(t) is observation uv_idx(i, t).
    ! The time dimension is 0 unless the model is diffuse or FILTER_UNIVARIATE.
    integer, allocatable :: uv_n(:)               !< (n)
    integer, allocatable :: uv_idx(:, :)          !< (p, n)
    real(dp), allocatable :: uv_Z(:, :, :)        !< (p, m, n) rows of Z*
    real(dp), allocatable :: uv_sig2(:, :)        !< (p, n) diagonal of D
    real(dp), allocatable :: uv_v(:, :)           !< (p, n) v_t,i
    real(dp), allocatable :: uv_Fstar(:, :)       !< (p, n) F_star,t,i
    real(dp), allocatable :: uv_Finf(:, :)        !< (p, n) F_inf,t,i
    real(dp), allocatable :: uv_Mstar(:, :, :)    !< (m, p, n) P_star,t,i Z*_i'
    real(dp), allocatable :: uv_Minf(:, :, :)     !< (m, p, n) P_inf,t,i Z*_i'
  end type filter_result_t

  !> Scratch arrays for one filter step, allocated once per filter run. The
  !> e_ arrays hold the per-element output of the last univariate step.
  type :: filter_ws_t
    real(dp), allocatable :: PZt(:, :), M(:, :), vz(:), Finv_v(:), TP(:, :), RQ(:, :), &
                             RQR(:, :)
    logical, allocatable :: obs(:)
    integer :: ir = 0, iq = 0   !< slices of R and Q that RQR was computed from
    !> Set by each step: the sum of v^2 / F (v' F^-1 v) terms in llf_t and how
    !> many observed elements they cover (excluding diffuse F_inf > 0 updates).
    real(dp) :: quad_t = 0.0_dp
    integer :: nstar_t = 0
    !> llf_t without its quadratic terms: llf_t = llf_nq_t - quad_t / 2.
    real(dp) :: llf_nq_t = 0.0_dp
    integer :: e_n = 0
    integer, allocatable :: e_idx(:)
    !> e_Zt(:, i) is row i of Z* (stored by columns so each is contiguous)
    real(dp), allocatable :: e_Zt(:, :), e_sig2(:), e_v(:), e_Fstar(:), e_Finf(:)
    real(dp), allocatable :: e_Mstar(:, :), e_Minf(:, :)
    !> Cache of the transformed observation equation: e_Zt and e_sig2 and the
    !> LDL factor L of H_oo stay valid while the slices of Z and H and the
    !> missing pattern are unchanged; then only y is transformed.
    logical :: c_valid = .false., c_diag = .false.
    integer :: c_iz = 0, c_ih = 0
    logical, allocatable :: c_obs(:)
    real(dp), allocatable :: c_L(:, :), c_y(:, :)
    real(dp), allocatable :: K0(:), K1(:)   !< univariate gains
    !> Steady state: whether the model allows it, whether it is in effect,
    !> the period t_ss whose F^-1, M = P Z' F^-1 and llf_nq it reuses.
    logical :: can_steady = .false., steady = .false.
    integer :: t_ss = 0
    real(dp), allocatable :: ss_Finv(:, :), ss_M(:, :)
    real(dp) :: ss_llf_nq = 0.0_dp
  end type filter_ws_t

contains

  !> Run the Kalman filter and store everything needed for smoothing.
  subroutine kalman_filter(rep, res, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(out) :: res
    integer, intent(out) :: info
    type(filter_ws_t) :: ws
    integer :: p, m, n, t, nu, it
    logical :: diffuse
    real(dp) :: nan

    call rep%validate(info)
    if (info /= SS_OK) return
    p = rep%k_endog; m = rep%k_states; n = rep%nobs
    nan = ieee_value(1.0_dp, ieee_quiet_nan)
    res%k_endog = p; res%k_states = m; res%nobs = n
    res%k_diffuse = rep%k_diffuse()
    allocate (res%a(m, n + 1), res%P(m, m, n + 1), res%Pinf(m, m, n + 1), &
              res%att(m, n), res%Ptt(m, m, n), res%yhat(p, n), res%v(p, n), &
              res%F(p, p, n), res%Finf(p, p, n), res%Finv(p, p, n), res%K(m, p, n), &
              res%llf_obs(n), res%method(n))
    call ws_init(ws, rep)

    call rep%initial_state(res%a(:, 1), res%P(:, :, 1), res%Pinf(:, :, 1), info)
    if (info /= SS_OK) return
    diffuse = sum(res%Pinf(:, :, 1)**2) > rep%tol_diffuse

    ! Per-element storage only if some period uses the univariate treatment.
    nu = 0
    if (diffuse .or. rep%filter_method == FILTER_UNIVARIATE) nu = n
    allocate (res%uv_n(n), res%uv_idx(p, nu), res%uv_Z(p, m, nu), res%uv_sig2(p, nu), &
              res%uv_v(p, nu), res%uv_Fstar(p, nu), res%uv_Finf(p, nu), &
              res%uv_Mstar(m, p, nu), res%uv_Minf(m, p, nu))
    res%uv_n = 0
    do t = 1, n
      res%method(t) = period_method(rep, t, diffuse, res%Pinf(:, :, t))
      if (diffuse) res%nobs_diffuse = t
      if (res%method(t) == METHOD_CONVENTIONAL .and. .not. diffuse) then
        if (use_steady(rep, t, ws)) then
          associate (ts => ws%t_ss)
            call steady_step(rep, t, ws, res%a(:, t), res%yhat(:, t), res%v(:, t), &
                             res%att(:, t), res%a(:, t + 1), res%llf_obs(t))
            res%F(:, :, t) = res%F(:, :, ts)
            res%Finv(:, :, t) = res%Finv(:, :, ts)
            res%K(:, :, t) = res%K(:, :, ts)
            res%Ptt(:, :, t) = res%Ptt(:, :, ts)
            res%P(:, :, t + 1) = res%P(:, :, ts + 1)
            res%Finf(:, :, t) = 0.0_dp
            res%Pinf(:, :, t + 1) = 0.0_dp
          end associate
          if (res%t_steady == 0) res%t_steady = t
          cycle
        end if
      end if
      select case (res%method(t))
      case (METHOD_UNIVARIATE)
        call observation_moments(rep, t, res%a(:, t), res%P(:, :, t), &
                                 res%Pinf(:, :, t), res%yhat(:, t), res%v(:, t), &
                                 res%F(:, :, t), res%Finf(:, :, t))
        res%Finv(:, :, t) = nan
        res%K(:, :, t) = nan
        call univariate_step(rep, t, ws, diffuse, res%a(:, t), res%P(:, :, t), &
                             res%Pinf(:, :, t), res%att(:, t), res%Ptt(:, :, t), &
                             res%a(:, t + 1), res%P(:, :, t + 1), &
                             res%Pinf(:, :, t + 1), res%llf_obs(t), info)
        if (info /= SS_OK) return
        call store_elements(ws, t, res)
      case (METHOD_DIFFUSE_MV)
        call observation_moments(rep, t, res%a(:, t), res%P(:, :, t), &
                                 res%Pinf(:, :, t), res%yhat(:, t), res%v(:, t), &
                                 res%F(:, :, t), res%Finf(:, :, t))
        res%Finv(:, :, t) = nan
        res%K(:, :, t) = nan
        call diffuse_mv_step(rep, t, ws, res%a(:, t), res%P(:, :, t), &
                             res%Pinf(:, :, t), res%att(:, t), res%Ptt(:, :, t), &
                             res%a(:, t + 1), res%P(:, :, t + 1), &
                             res%Pinf(:, :, t + 1), res%llf_obs(t), info)
        if (info /= SS_OK) return
      case default
        ! Conventional; in a diffuse period this is the F_inf = 0 case, where
        ! P_inf only propagates: P_inf,t+1 = T P_inf T'.
        call filter_step(rep, t, ws, res%a(:, t), res%P(:, :, t), res%yhat(:, t), &
                         res%v(:, t), res%F(:, :, t), res%Finv(:, :, t), &
                         res%K(:, :, t), res%att(:, t), res%Ptt(:, :, t), &
                         res%a(:, t + 1), res%P(:, :, t + 1), res%llf_obs(t), info)
        if (info /= SS_OK) return
        res%Finf(:, :, t) = 0.0_dp
        res%Pinf(:, :, t + 1) = 0.0_dp
        if (.not. diffuse) then
          call check_steady(rep, t, ws, res%P(:, :, t), res%P(:, :, t + 1), &
                            res%Finv(:, :, t))
        end if
        if (diffuse) then
          call observation_moments(rep, t, res%a(:, t), res%P(:, :, t), &
                                   res%Pinf(:, :, t), res%yhat(:, t), res%v(:, t), &
                                   res%F(:, :, t), res%Finf(:, :, t))
          it = tidx(size(rep%T, 3), t)
          call gemm('N', 'N', 1.0_dp, rep%T(:, :, it), res%Pinf(:, :, t), 0.0_dp, ws%TP)
          call gemm('N', 'T', 1.0_dp, ws%TP, rep%T(:, :, it), 0.0_dp, &
                    res%Pinf(:, :, t + 1))
        end if
      end select
      if (diffuse) then
        diffuse = sum(res%Pinf(:, :, t + 1)**2) > rep%tol_diffuse
        if (.not. diffuse) res%Pinf(:, :, t + 1) = 0.0_dp
      end if
    end do
    res%llf = sum(res%llf_obs(rep%loglikelihood_burn + 1:))
    if (rep%marginal_likelihood) res%llf = res%llf + marginal_correction(rep, info)
  end subroutine kalman_filter

  !> Log likelihood only (DK 7.2), with no stored history.
  function loglike(rep, info) result(llf)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(out) :: info
    real(dp) :: llf
    real(dp) :: quad
    integer :: nstar

    call loglike_core(rep, llf, quad, nstar, info)
    if (info == SS_OK .and. rep%marginal_likelihood) then
      llf = llf + marginal_correction(rep, info)
    end if
  end function loglike

  !> Log likelihood with a scale factor sigma^2 concentrated out (DK 2.10.2,
  !> 7.3): H, Q and P_star are taken relative to sigma^2 (P_inf is not
  !> affected). With S the sum of the v' F^-1 v terms at sigma^2 = 1 over the
  !> n* observed elements outside the diffuse (F_inf > 0) updates,
  !>
  !>     sigma2_hat = S / n*,   log L_c = log L_1 + S/2 - n* (log sigma2_hat + 1) / 2,
  !>
  !> which is the ordinary log likelihood of the model with H, Q and P_star
  !> multiplied by sigma2_hat (see `ssm_rep_t%scale_by`).
  function loglike_concentrated(rep, scale, info) result(llf)
    type(ssm_rep_t), intent(in) :: rep
    real(dp), intent(out) :: scale
    integer, intent(out) :: info
    real(dp) :: llf
    real(dp) :: quad, llf_nq
    integer :: nstar

    ! Built from the non-quadratic part, avoiding the cancellation of
    ! log L_1 + S/2 (S can be large when the scale is).
    call loglike_core(rep, llf, quad, nstar, info, llf_nq)
    scale = 0.0_dp
    if (info /= SS_OK .or. nstar == 0) return
    scale = quad / nstar
    llf = llf_nq - 0.5_dp * nstar * (log(scale) + 1.0_dp)
    if (rep%marginal_likelihood) llf = llf + marginal_correction(rep, info)
  end function loglike_concentrated

  !> Marginal minus diffuse log likelihood (Francke, Koopman and de Vos 2010,
  !> eqs. 16 and 21; DK 7.2.6):
  !>
  !>     log L_M = log L_d + k/2 log 2 pi + 1/2 log|S*|,
  !>     S* = X' X = sum_t V*_t' V*_t,  V*_t = Z_t,o A*_t,  A*_t+1 = T_t A*_t,
  !>
  !> with A*_1 = A the k diffuse directions of alpha_1 (A A' = P_inf) and Z_t,o
  !> the rows of Z_t for the observed elements. (The k/2 log 2 pi arises because
  !> the marginal likelihood is the density of n - k transformed observations.)
  !> The correction depends on Z, T and A only. Zero without diffuse states.
  real(dp) function marginal_correction(rep, info) result(corr)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(out) :: info
    real(dp), allocatable :: a1(:), Pstar(:, :), Pinf(:, :), L(:, :), D(:), A(:, :), &
                             V(:, :)
    real(dp), allocatable :: S(:, :), TA(:, :)
    integer, allocatable :: idx(:)
    real(dp) :: logdet
    integer :: m, k, j, t, i, iz, it

    corr = 0.0_dp
    m = rep%k_states
    allocate (a1(m), Pstar(m, m), Pinf(m, m), L(m, m), D(m))
    call rep%initial_state(a1, Pstar, Pinf, info)
    if (info /= SS_OK) return
    call ldl_psd(Pinf, L, D)
    k = count(D > 0.0_dp)
    if (k == 0) return
    allocate (A(m, k), S(k, k), source=0.0_dp)
    k = 0
    do j = 1, m
      if (D(j) > 0.0_dp) then
        k = k + 1
        A(:, k) = L(:, j) * sqrt(D(j))
      end if
    end do
    do t = 1, rep%nobs
      iz = tidx(size(rep%Z, 3), t)
      it = tidx(size(rep%T, 3), t)
      idx = pack([(i, i=1, rep%k_endog)], .not. ieee_is_nan(rep%y(:, t)))
      if (size(idx) > 0) then
        V = matmul(rep%Z(idx, :, iz), A)
        S = S + matmul(transpose(V), V)
      end if
      TA = matmul(rep%T(:, :, it), A)
      A = TA
    end do
    call symmetrize(S)
    call chol_inv(S, logdet, info)
    if (info /= SS_OK) return
    corr = 0.5_dp * (k * log2pi + logdet)
  end function marginal_correction

  !> Filter pass without stored history: the log likelihood, the sum of the
  !> quadratic terms and their count, and optionally the log likelihood
  !> without the quadratic terms (all after the burn-in).
  subroutine loglike_core(rep, llf, quad, nstar, info, llf_nq)
    type(ssm_rep_t), intent(in) :: rep
    real(dp), intent(out) :: llf, quad
    real(dp), intent(out), optional :: llf_nq
    integer, intent(out) :: nstar
    integer, intent(out) :: info
    type(filter_ws_t) :: ws
    real(dp), allocatable :: at(:), Pt(:, :), Pinft(:, :), anext(:), Pnext(:, :), &
                             Pinfnext(:, :)
    real(dp), allocatable :: att(:), Ptt(:, :), yhat(:), v(:), F(:, :), Finv(:, :), &
                             K(:, :)
    real(dp) :: llf_t
    integer :: p, m, t, it
    logical :: diffuse

    llf = 0.0_dp
    quad = 0.0_dp
    nstar = 0
    if (present(llf_nq)) llf_nq = 0.0_dp
    call rep%validate(info)
    if (info /= SS_OK) return
    p = rep%k_endog; m = rep%k_states
    allocate (at(m), Pt(m, m), Pinft(m, m), anext(m), Pnext(m, m), Pinfnext(m, m), &
              att(m), Ptt(m, m), yhat(p), v(p), F(p, p), Finv(p, p), K(m, p))
    call ws_init(ws, rep)

    call rep%initial_state(at, Pt, Pinft, info)
    if (info /= SS_OK) return
    diffuse = sum(Pinft**2) > rep%tol_diffuse
    do t = 1, rep%nobs
      if (.not. diffuse) then
        if (use_steady(rep, t, ws)) then
          ! P (Pt) stays at its steady-state value
          call steady_step(rep, t, ws, at, yhat, v, att, anext, llf_t)
          if (t > rep%loglikelihood_burn) then
            llf = llf + llf_t
            quad = quad + ws%quad_t
            nstar = nstar + ws%nstar_t
            if (present(llf_nq)) llf_nq = llf_nq + ws%llf_nq_t
          end if
          at = anext
          cycle
        end if
      end if
      select case (period_method(rep, t, diffuse, Pinft))
      case (METHOD_UNIVARIATE)
        call univariate_step(rep, t, ws, diffuse, at, Pt, Pinft, att, Ptt, anext, &
                             Pnext, Pinfnext, llf_t, info)
      case (METHOD_DIFFUSE_MV)
        call diffuse_mv_step(rep, t, ws, at, Pt, Pinft, att, Ptt, anext, Pnext, &
                             Pinfnext, llf_t, info)
      case default
        call filter_step(rep, t, ws, at, Pt, yhat, v, F, Finv, K, att, Ptt, anext, &
                         Pnext, llf_t, info)
        if (info == SS_OK .and. .not. diffuse) then
          call check_steady(rep, t, ws, Pt, Pnext, Finv)
        end if
        if (diffuse) then
          it = tidx(size(rep%T, 3), t)
          call gemm('N', 'N', 1.0_dp, rep%T(:, :, it), Pinft, 0.0_dp, ws%TP)
          call gemm('N', 'T', 1.0_dp, ws%TP, rep%T(:, :, it), 0.0_dp, Pinfnext)
        end if
      end select
      if (info /= SS_OK) return
      if (t > rep%loglikelihood_burn) then
        llf = llf + llf_t
        quad = quad + ws%quad_t
        nstar = nstar + ws%nstar_t
        if (present(llf_nq)) llf_nq = llf_nq + ws%llf_nq_t
      end if
      at = anext
      Pt = Pnext
      if (diffuse) then
        Pinft = Pinfnext
        diffuse = sum(Pinft**2) > rep%tol_diffuse
      end if
    end do
  end subroutine loglike_core

  !> The method for period t (see METHOD_*). With diffuse states and
  !> DIFFUSE_MULTIVARIATE, it depends on F_inf over the observed elements:
  !> zero (conventional), nonsingular (multivariate exact initial update), or
  !> singular but nonzero (univariate step, as DK recommend).
  integer function period_method(rep, t, diffuse, Pinf) result(method)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: t
    logical, intent(in) :: diffuse
    real(dp), intent(in), contiguous :: Pinf(:, :)
    real(dp), allocatable :: Zo(:, :), Fi(:, :), L(:, :), D(:)
    logical, allocatable :: obs(:)
    integer :: n_o, iz

    if (rep%filter_method == FILTER_UNIVARIATE) then
      method = METHOD_UNIVARIATE
    else if (.not. diffuse) then
      method = METHOD_CONVENTIONAL
    else if (rep%diffuse_method /= DIFFUSE_MULTIVARIATE) then
      method = METHOD_UNIVARIATE
    else
      obs = .not. ieee_is_nan(rep%y(:, t))
      n_o = count(obs)
      iz = tidx(size(rep%Z, 3), t)
      Zo = reshape(pack(rep%Z(:, :, iz), spread(obs, 2, rep%k_states)), &
                   [n_o, rep%k_states])
      Fi = matmul(Zo, matmul(Pinf, transpose(Zo)))
      if (maxval(abs(Fi)) <= rep%tol_diffuse) then
        method = METHOD_CONVENTIONAL
      else
        allocate (L(n_o, n_o), D(n_o))
        call ldl_psd(Fi, L, D)
        if (all(D > rep%tol_diffuse)) then
          method = METHOD_DIFFUSE_MV
        else
          method = METHOD_UNIVARIATE
        end if
      end if
    end if
  end function period_method

  !> Multivariate exact initial update for a period with F_inf nonsingular
  !> over the observed elements o (DK 5.2.1), in filtered form:
  !>
  !>   F1 = F_inf^-1,  F2 = -F1 F_star F1,  M_star = P_star Z_o', M_inf = P_inf Z_o'
  !>   a_t|t      = a + M_inf F1 v
  !>   P_inf,t|t  = P_inf - M_inf F1 M_inf'
  !>   P_star,t|t = P_star - M_inf F1 M_star' - M_star F1 M_inf' - M_inf F2 M_inf'
  !>
  !> followed by the usual prediction, with P_inf,t+1 = T P_inf,t|t T'. The
  !> period contributes -(n_o log 2 pi + log|F_inf|) / 2 (DK 7.2.2).
  subroutine diffuse_mv_step(rep, t, ws, a, Pstar, Pinf, att, Ptt, anext, Pnext, &
                             Pinfnext, llf_t, info)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: t
    type(filter_ws_t), intent(inout) :: ws
    real(dp), intent(in), contiguous :: a(:), Pstar(:, :), Pinf(:, :)
    real(dp), intent(out), contiguous :: att(:), Ptt(:, :), anext(:), Pnext(:, :), &
                                         Pinfnext(:, :)
    real(dp), intent(out) :: llf_t
    integer, intent(out) :: info
    real(dp), allocatable :: Zo(:, :), vo(:), Fs(:, :), F1(:, :), F2(:, :), Ms(:, :), &
                             Mi(:, :)
    real(dp), allocatable :: MiF1(:, :), Pinftt(:, :)
    integer, allocatable :: idx(:)
    real(dp) :: logdet
    integer :: i, n_o, iz, ih, id, it

    iz = tidx(size(rep%Z, 3), t)
    ih = tidx(size(rep%H, 3), t)
    id = tidx(size(rep%d, 2), t)
    it = tidx(size(rep%T, 3), t)
    idx = pack([(i, i=1, rep%k_endog)], .not. ieee_is_nan(rep%y(:, t)))
    n_o = size(idx)

    Zo = rep%Z(idx, :, iz)
    vo = rep%y(idx, t) - rep%d(idx, id) - matmul(Zo, a)
    Ms = matmul(Pstar, transpose(Zo))
    Mi = matmul(Pinf, transpose(Zo))
    Fs = matmul(Zo, Ms) + rep%H(idx, idx, ih)
    F1 = matmul(Zo, Mi)
    call symmetrize(F1)
    call chol_inv(F1, logdet, info)
    if (info /= SS_OK) return
    F2 = -matmul(F1, matmul(Fs, F1))
    MiF1 = matmul(Mi, F1)

    att = a + matmul(MiF1, vo)
    Pinftt = Pinf - matmul(MiF1, transpose(Mi))
    Ptt = Pstar - matmul(MiF1, transpose(Ms)) - matmul(Ms, transpose(MiF1)) &
          - matmul(Mi, matmul(F2, transpose(Mi)))
    call symmetrize(Pinftt)
    call symmetrize(Ptt)

    call predict(rep, t, ws, att, Ptt, anext, Pnext)
    call gemm('N', 'N', 1.0_dp, rep%T(:, :, it), Pinftt, 0.0_dp, ws%TP)
    call gemm('N', 'T', 1.0_dp, ws%TP, rep%T(:, :, it), 0.0_dp, Pinfnext)
    call symmetrize(Pinfnext)
    llf_t = -0.5_dp * (n_o * log2pi + logdet)
    ws%quad_t = 0.0_dp
    ws%nstar_t = 0
    ws%llf_nq_t = llf_t
  end subroutine diffuse_mv_step

  subroutine ws_init(ws, rep)
    type(filter_ws_t), intent(out) :: ws
    type(ssm_rep_t), intent(in) :: rep
    integer :: p, m, r

    p = rep%k_endog; m = rep%k_states; r = rep%k_posdef
    allocate (ws%PZt(m, p), ws%M(m, p), ws%vz(p), ws%Finv_v(p), ws%TP(m, m), &
              ws%RQ(m, r), ws%RQR(m, m), ws%obs(p))
    allocate (ws%e_idx(p), ws%e_Zt(m, p), ws%e_sig2(p), ws%e_v(p), ws%e_Fstar(p), &
              ws%e_Finf(p), ws%e_Mstar(m, p), ws%e_Minf(m, p))
    ws%can_steady = rep%tol_steady >= 0.0_dp .and. rep%filter_method &
                    /= FILTER_UNIVARIATE .and. size(rep%Z, 3) == 1 &
                    .and. size(rep%H, 3) == 1 .and. size(rep%T, 3) == 1 &
                    .and. size(rep%R, 3) == 1 .and. size(rep%Q, 3) == 1
    if (ws%can_steady) allocate (ws%ss_Finv(p, p), ws%ss_M(m, p))
    allocate (ws%c_obs(p), ws%c_L(p, p), ws%c_y(p, 1), ws%K0(m), ws%K1(m))
  end subroutine ws_init

  !> After a conventional, fully observed, non-diffuse step at t from P to
  !> Pnext: enter the steady state if P has converged (see tol_steady).
  subroutine check_steady(rep, t, ws, P, Pnext, Finv)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: t
    type(filter_ws_t), intent(inout) :: ws
    real(dp), intent(in), contiguous :: P(:, :), Pnext(:, :), Finv(:, :)

    if (.not. ws%can_steady .or. ws%steady) return
    if (.not. all(ws%obs)) return
    if (sum((Pnext - P)**2) > rep%tol_steady**2 * sum(Pnext**2)) return
    ws%steady = .true.
    ws%t_ss = t
    ws%ss_Finv = Finv
    ws%ss_M = ws%M
    ws%ss_llf_nq = ws%llf_nq_t
  end subroutine check_steady

  !> Whether period t can use the steady state: it is in effect and y_t is
  !> fully observed. A missing element ends it.
  logical function use_steady(rep, t, ws)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: t
    type(filter_ws_t), intent(inout) :: ws

    use_steady = .false.
    if (.not. ws%steady) return
    if (any(ieee_is_nan(rep%y(:, t)))) then
      ws%steady = .false.
      return
    end if
    use_steady = .true.
  end function use_steady

  !> One steady-state step (DK 4.3.4): P, F and K are fixed, so only
  !>   v = y - d - Z a,  a_t|t = a + M v,  a_t+1 = c + T a_t|t
  !> and the log likelihood are computed, with M = P Z' F^-1.
  subroutine steady_step(rep, t, ws, a, yhat, v, att, anext, llf_t)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: t
    type(filter_ws_t), intent(inout) :: ws
    real(dp), intent(in), contiguous :: a(:)
    real(dp), intent(out), contiguous :: yhat(:), v(:), att(:), anext(:)
    real(dp), intent(out) :: llf_t
    integer :: id, ic

    id = tidx(size(rep%d, 2), t)
    ic = tidx(size(rep%c, 2), t)
    yhat = rep%d(:, id)
    call gemv('N', 1.0_dp, rep%Z(:, :, 1), a, 1.0_dp, yhat)
    v = rep%y(:, t) - yhat
    call gemv('N', 1.0_dp, ws%ss_Finv, v, 0.0_dp, ws%Finv_v)
    att = a
    call gemv('N', 1.0_dp, ws%ss_M, v, 1.0_dp, att)
    anext = rep%c(:, ic)
    call gemv('N', 1.0_dp, rep%T(:, :, 1), att, 1.0_dp, anext)
    ws%quad_t = dot_product(v, ws%Finv_v)
    ws%nstar_t = rep%k_endog
    ws%llf_nq_t = ws%ss_llf_nq
    llf_t = ws%ss_llf_nq - 0.5_dp * ws%quad_t
  end subroutine steady_step

  !> Cache R Q R' for the slices of R and Q at time t.
  subroutine update_RQR(rep, t, ws)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: t
    type(filter_ws_t), intent(inout) :: ws
    integer :: ir, iq

    ir = tidx(size(rep%R, 3), t)
    iq = tidx(size(rep%Q, 3), t)
    if (ir == ws%ir .and. iq == ws%iq) return
    call gemm('N', 'N', 1.0_dp, rep%R(:, :, ir), rep%Q(:, :, iq), 0.0_dp, ws%RQ)
    call gemm('N', 'T', 1.0_dp, ws%RQ, rep%R(:, :, ir), 0.0_dp, ws%RQR)
    call symmetrize(ws%RQR)
    ws%ir = ir
    ws%iq = iq
  end subroutine update_RQR

  !> Predicted: a_t+1 = c + T a_t|t,  P_t+1 = T P_t|t T' + R Q R'.
  subroutine predict(rep, t, ws, att, Ptt, anext, Pnext)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: t
    type(filter_ws_t), intent(inout) :: ws
    real(dp), intent(in), contiguous :: att(:), Ptt(:, :)
    real(dp), intent(out), contiguous :: anext(:), Pnext(:, :)
    integer :: it, ic

    it = tidx(size(rep%T, 3), t)
    ic = tidx(size(rep%c, 2), t)
    anext = rep%c(:, ic)
    call gemv('N', 1.0_dp, rep%T(:, :, it), att, 1.0_dp, anext)
    call update_RQR(rep, t, ws)
    call gemm('N', 'N', 1.0_dp, rep%T(:, :, it), Ptt, 0.0_dp, ws%TP)
    Pnext = ws%RQR
    call gemm('N', 'T', 1.0_dp, ws%TP, rep%T(:, :, it), 1.0_dp, Pnext)
    call symmetrize(Pnext)
  end subroutine predict

  !> One conventional step: from the predicted (a, P) at time t, compute the
  !> forecast error quantities, the filtered state, and the next prediction.
  subroutine filter_step(rep, t, ws, a, P, yhat, v, F, Finv, K, att, Ptt, anext, &
                         Pnext, llf_t, info)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: t
    type(filter_ws_t), intent(inout) :: ws
    real(dp), intent(in), contiguous :: a(:), P(:, :)
    real(dp), intent(out), contiguous :: yhat(:), v(:), F(:, :), Finv(:, :), K(:, :)
    real(dp), intent(out), contiguous :: att(:), Ptt(:, :), anext(:), Pnext(:, :)
    real(dp), intent(out) :: llf_t
    integer, intent(out) :: info
    integer :: iz, ih, it, id, nobs_t
    real(dp) :: logdetF

    iz = tidx(size(rep%Z, 3), t)
    ih = tidx(size(rep%H, 3), t)
    it = tidx(size(rep%T, 3), t)
    id = tidx(size(rep%d, 2), t)

    ws%obs = .not. ieee_is_nan(rep%y(:, t))
    nobs_t = count(ws%obs)

    ! yhat = d + Z a,  v = y - yhat (NaN where missing; vz has zeros there)
    yhat = rep%d(:, id)
    call gemv('N', 1.0_dp, rep%Z(:, :, iz), a, 1.0_dp, yhat)
    v = rep%y(:, t) - yhat
    ws%vz = merge(v, 0.0_dp, ws%obs)

    ! F = Z P Z' + H, and F^-1 over the observed elements
    call gemm('N', 'T', 1.0_dp, P, rep%Z(:, :, iz), 0.0_dp, ws%PZt)
    F = rep%H(:, :, ih)
    call gemm('N', 'N', 1.0_dp, rep%Z(:, :, iz), ws%PZt, 1.0_dp, F)
    call symmetrize(F)
    call observed_inverse(F, ws%obs, nobs_t, Finv, logdetF, info)
    if (info /= SS_OK) return
    call gemv('N', 1.0_dp, Finv, ws%vz, 0.0_dp, ws%Finv_v)

    ! Filtered: a_t|t = a + P Z' F^-1 v,  P_t|t = P - P Z' F^-1 Z P
    call gemm('N', 'N', 1.0_dp, ws%PZt, Finv, 0.0_dp, ws%M)
    att = a
    call gemv('N', 1.0_dp, ws%PZt, ws%Finv_v, 1.0_dp, att)
    Ptt = P
    call gemm('N', 'T', -1.0_dp, ws%M, ws%PZt, 1.0_dp, Ptt)
    call symmetrize(Ptt)

    ! K = T P Z' F^-1
    call gemm('N', 'N', 1.0_dp, rep%T(:, :, it), ws%M, 0.0_dp, K)

    call predict(rep, t, ws, att, Ptt, anext, Pnext)

    if (nobs_t == 0) then
      llf_t = 0.0_dp
    else
      llf_t = -0.5_dp * (nobs_t * log2pi + logdetF + dot_product(ws%vz, ws%Finv_v))
    end if
    ws%quad_t = dot_product(ws%vz, ws%Finv_v)
    ws%nstar_t = nobs_t
    ws%llf_nq_t = -0.5_dp * (nobs_t * log2pi + logdetF)
    if (nobs_t == 0) ws%llf_nq_t = 0.0_dp
  end subroutine filter_step

  !> One univariate step (DK 6.4, with DK 5.2 when `diffuse`). The observed
  !> elements of y_t are processed in turn, each updating (a, P_star, P_inf):
  !>
  !>   F_inf > 0:  K0 = M_inf / F_inf,  K1 = (M_star - K0 F_star) / F_inf
  !>               a += K0 v,  P_star -= M_star K0' + M_inf K1',
  !>               P_inf -= M_inf K0',  llf += -(log 2pi + log F_inf) / 2
  !>   otherwise:  K = M_star / F_star,  a += K v,  P_star -= M_star K',
  !>               llf += -(log 2pi + log F_star + v^2 / F_star) / 2
  !>
  !> with M = P Z*_i' and F_star = Z*_i P_star Z*_i' + sigma2_i,
  !> F_inf = Z*_i P_inf Z*_i'. Elements with F_star = 0 (or <= tol_diffuse
  !> in diffuse periods) carry no information and are skipped.
  subroutine univariate_step(rep, t, ws, diffuse, a, Pstar, Pinf, att, Ptt, anext, &
                             Pnext, Pinfnext, llf_t, info)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: t
    type(filter_ws_t), intent(inout) :: ws
    logical, intent(in) :: diffuse
    real(dp), intent(in), contiguous :: a(:), Pstar(:, :), Pinf(:, :)
    real(dp), intent(out), contiguous :: att(:), Ptt(:, :), anext(:), Pnext(:, :), &
                                         Pinfnext(:, :)
    real(dp), intent(out) :: llf_t
    integer, intent(out) :: info
    real(dp), allocatable :: Pinftt(:, :)
    real(dp) :: v, Fs, Fi, tol_Fs
    integer :: i, j, n_o, it

    info = SS_OK
    ! As in statsmodels, diffuse periods skip elements with F_star <= tol.
    tol_Fs = merge(rep%tol_diffuse, 0.0_dp, diffuse)
    call transform_observation(rep, t, ws)
    n_o = ws%e_n
    att = a
    Ptt = Pstar
    if (diffuse) Pinftt = Pinf
    llf_t = 0.0_dp
    ws%quad_t = 0.0_dp
    ws%nstar_t = 0
    ws%llf_nq_t = 0.0_dp
    do i = 1, n_o
      associate (z => ws%e_Zt(:, i), Ms => ws%e_Mstar(:, i), Mi => ws%e_Minf(:, i))
        v = ws%e_v(i) - dot_product(z, att)
        call gemv('N', 1.0_dp, Ptt, z, 0.0_dp, Ms)
        Fs = max(dot_product(z, Ms) + ws%e_sig2(i), 0.0_dp)
        if (diffuse) then
          call gemv('N', 1.0_dp, Pinftt, z, 0.0_dp, Mi)
          Fi = max(dot_product(z, Mi), 0.0_dp)
        else
          Mi = 0.0_dp
          Fi = 0.0_dp
        end if
        ws%e_v(i) = v
        ws%e_Fstar(i) = Fs
        ws%e_Finf(i) = Fi

        if (diffuse .and. Fi > rep%tol_diffuse) then
          ws%K0 = Mi / Fi
          ws%K1 = (Ms - ws%K0 * Fs) / Fi
          att = att + ws%K0 * v
          do j = 1, rep%k_states
            Ptt(:, j) = Ptt(:, j) - Ms * ws%K0(j) - Mi * ws%K1(j)
            Pinftt(:, j) = Pinftt(:, j) - Mi * ws%K0(j)
          end do
          llf_t = llf_t - 0.5_dp * (log2pi + log(Fi))
          ws%llf_nq_t = ws%llf_nq_t - 0.5_dp * (log2pi + log(Fi))
        else if (Fs > tol_Fs) then
          ws%K0 = Ms / Fs
          att = att + ws%K0 * v
          do j = 1, rep%k_states
            Ptt(:, j) = Ptt(:, j) - Ms * ws%K0(j)
          end do
          llf_t = llf_t - 0.5_dp * (log2pi + log(Fs) + v**2 / Fs)
          ws%quad_t = ws%quad_t + v**2 / Fs
          ws%nstar_t = ws%nstar_t + 1
          ws%llf_nq_t = ws%llf_nq_t - 0.5_dp * (log2pi + log(Fs))
        end if
      end associate
    end do
    call symmetrize(Ptt)

    call predict(rep, t, ws, att, Ptt, anext, Pnext)
    if (diffuse) then
      call symmetrize(Pinftt)
      it = tidx(size(rep%T, 3), t)
      call gemm('N', 'N', 1.0_dp, rep%T(:, :, it), Pinftt, 0.0_dp, ws%TP)
      call gemm('N', 'T', 1.0_dp, ws%TP, rep%T(:, :, it), 0.0_dp, Pinfnext)
      call symmetrize(Pinfnext)
    else
      Pinfnext = 0.0_dp
    end if
  end subroutine univariate_step

  !> Select the observed elements of y_t and, if their block of H_t is not
  !> diagonal, transform them (DK 6.4.3): H_oo = L D L', y* = L^-1 (y - d),
  !> Z* = L^-1 Z, with variances D. Results go to ws%e_n, e_idx, e_Zt,
  !> e_sig2, and e_v (holding y* until the step computes v).
  subroutine transform_observation(rep, t, ws)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: t
    type(filter_ws_t), intent(inout) :: ws
    real(dp), allocatable :: Ho(:, :), L(:, :), D(:), Zo(:, :), yo(:, :)
    integer :: i, n_o, iz, ih, id

    iz = tidx(size(rep%Z, 3), t)
    ih = tidx(size(rep%H, 3), t)
    id = tidx(size(rep%d, 2), t)
    ws%obs = .not. ieee_is_nan(rep%y(:, t))
    n_o = count(ws%obs)
    ws%e_n = n_o
    if (n_o == 0) return

    if (ws%c_valid .and. iz == ws%c_iz .and. ih == ws%c_ih) then
      if (all(ws%obs .eqv. ws%c_obs)) then
        ! Same observation equation as last time: transform y only
        associate (idx => ws%e_idx(1:n_o))
          ws%c_y(1:n_o, 1) = rep%y(idx, t) - rep%d(idx, id)
          if (.not. ws%c_diag) then
            call solve_unit_lower(ws%c_L(1:n_o, 1:n_o), ws%c_y(1:n_o, :))
          end if
          ws%e_v(1:n_o) = ws%c_y(1:n_o, 1)
        end associate
        return
      end if
    end if

    ws%e_idx(1:n_o) = pack([(i, i=1, rep%k_endog)], ws%obs)
    associate (idx => ws%e_idx(1:n_o))
      Ho = rep%H(idx, idx, ih)
      Zo = rep%Z(idx, :, iz)
      yo = reshape(rep%y(idx, t) - rep%d(idx, id), [n_o, 1])
      allocate (L(n_o, n_o), D(n_o))
      if (is_diagonal(Ho, tol_diagonal)) then
        L = 0.0_dp
        do i = 1, n_o
          L(i, i) = 1.0_dp
          D(i) = Ho(i, i)
        end do
      else
        call ldl_psd(Ho, L, D)
        call solve_unit_lower(L, Zo)
        call solve_unit_lower(L, yo)
      end if
    end associate
    ws%e_Zt(:, 1:n_o) = transpose(Zo)
    ws%e_sig2(1:n_o) = D
    ws%e_v(1:n_o) = yo(:, 1)
    ws%c_valid = .true.
    ws%c_iz = iz
    ws%c_ih = ih
    ws%c_obs = ws%obs
    ws%c_diag = all(L == eye(n_o))
    ws%c_L(1:n_o, 1:n_o) = L
  end subroutine transform_observation

  !> Copy the per-element output of the last univariate step into `res`.
  subroutine store_elements(ws, t, res)
    type(filter_ws_t), intent(in) :: ws
    integer, intent(in) :: t
    type(filter_result_t), intent(inout) :: res
    integer :: n_o

    n_o = ws%e_n
    res%uv_n(t) = n_o
    res%uv_idx(:, t) = 0
    res%uv_Z(:, :, t) = 0.0_dp
    res%uv_sig2(:, t) = 0.0_dp
    res%uv_v(:, t) = 0.0_dp
    res%uv_Fstar(:, t) = 0.0_dp
    res%uv_Finf(:, t) = 0.0_dp
    res%uv_Mstar(:, :, t) = 0.0_dp
    res%uv_Minf(:, :, t) = 0.0_dp
    if (n_o == 0) return
    res%uv_idx(1:n_o, t) = ws%e_idx(1:n_o)
    res%uv_Z(1:n_o, :, t) = transpose(ws%e_Zt(:, 1:n_o))
    res%uv_sig2(1:n_o, t) = ws%e_sig2(1:n_o)
    res%uv_v(1:n_o, t) = ws%e_v(1:n_o)
    res%uv_Fstar(1:n_o, t) = ws%e_Fstar(1:n_o)
    res%uv_Finf(1:n_o, t) = ws%e_Finf(1:n_o)
    res%uv_Mstar(:, 1:n_o, t) = ws%e_Mstar(:, 1:n_o)
    res%uv_Minf(:, 1:n_o, t) = ws%e_Minf(:, 1:n_o)
  end subroutine store_elements

  !> Original-coordinate forecast quantities at the predicted (a, P, P_inf):
  !> yhat = d + Z a, v = y - yhat, F = Z P Z' + H, Finf = Z P_inf Z'.
  subroutine observation_moments(rep, t, a, P, Pinf, yhat, v, F, Finf)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: t
    real(dp), intent(in), contiguous :: a(:), P(:, :), Pinf(:, :)
    real(dp), intent(out), contiguous :: yhat(:), v(:), F(:, :), Finf(:, :)
    real(dp), allocatable :: PZt(:, :)
    integer :: iz, ih, id

    iz = tidx(size(rep%Z, 3), t)
    ih = tidx(size(rep%H, 3), t)
    id = tidx(size(rep%d, 2), t)
    yhat = rep%d(:, id)
    call gemv('N', 1.0_dp, rep%Z(:, :, iz), a, 1.0_dp, yhat)
    v = rep%y(:, t) - yhat
    allocate (PZt(rep%k_states, rep%k_endog))
    call gemm('N', 'T', 1.0_dp, P, rep%Z(:, :, iz), 0.0_dp, PZt)
    F = rep%H(:, :, ih)
    call gemm('N', 'N', 1.0_dp, rep%Z(:, :, iz), PZt, 1.0_dp, F)
    call symmetrize(F)
    call gemm('N', 'T', 1.0_dp, Pinf, rep%Z(:, :, iz), 0.0_dp, PZt)
    call gemm('N', 'N', 1.0_dp, rep%Z(:, :, iz), PZt, 0.0_dp, Finf)
    call symmetrize(Finf)
  end subroutine observation_moments

  !> Inverse and log determinant of the observed submatrix of F, returned
  !> zero-padded to full size in Finv.
  subroutine observed_inverse(F, obs, nobs_t, Finv, logdet, info)
    real(dp), intent(in), contiguous :: F(:, :)
    logical, intent(in) :: obs(:)
    integer, intent(in) :: nobs_t
    real(dp), intent(out), contiguous :: Finv(:, :)
    real(dp), intent(out) :: logdet
    integer, intent(out) :: info
    real(dp), allocatable :: Fo(:, :)
    integer, allocatable :: idx(:)
    integer :: i

    info = SS_OK
    logdet = 0.0_dp
    if (nobs_t == size(F, 1)) then
      Finv = F
      call chol_inv(Finv, logdet, info)
    else if (nobs_t == 0) then
      Finv = 0.0_dp
    else
      idx = pack([(i, i=1, size(F, 1))], obs)
      Fo = F(idx, idx)
      call chol_inv(Fo, logdet, info)
      Finv = 0.0_dp
      Finv(idx, idx) = Fo
    end if
  end subroutine observed_inverse

  !> Steady state of the Kalman filter for a time-invariant model (DK 2.11,
  !> 4.3.4): the limit P of
  !>
  !>     P_t+1 = T P_t T' - T P_t Z' F_t^-1 Z P_t T' + R Q R',  F_t = Z P_t Z' + H,
  !>
  !> started from P_1 = 0, and F = Z P Z' + H. With H positive definite the
  !> structure-preserving doubling algorithm (Chu, Fan, Lin and Wang 2004) is
  !> used: writing the recursion
  !> as X = A' X (I + G X)^-1 A + W with A = T', G = Z' H^-1 Z, W = R Q R',
  !>
  !>     A_k+1 = A_k (I + G_k X_k)^-1 A_k
  !>     G_k+1 = G_k + A_k (I + G_k X_k)^-1 G_k A_k'
  !>     X_k+1 = X_k + A_k' X_k (I + G_k X_k)^-1 A_k,
  !>
  !> X_0 = W, so that X_k = P_(2^k + 1). Convergence like 1/t (a fixed
  !> regression effect, a nearly fixed seasonal) then takes about 40 steps.
  !> Otherwise the recursion itself is iterated, at most `maxiter` times
  !> (default 100000), and a singular F_t gives SS_ERR_NOT_PD. Converged when the
  !> largest change in P is at most `tol` (default 1e-12) times max(1, max|P|);
  !> else SS_ERR_NOT_CONVERGED, with the last iterate returned. Time-varying
  !> system matrices give SS_ERR_UNSUPPORTED. Missing observations and the
  !> initialization play no role.
  subroutine steady_state(rep, P, F, info, tol, maxiter, niter)
    type(ssm_rep_t), intent(in) :: rep
    real(dp), intent(out), contiguous :: P(:, :)          !< (m, m) steady-state P
    real(dp), intent(out), contiguous :: F(:, :)          !< (p, p) steady-state F
    integer, intent(out) :: info
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    integer, intent(out), optional :: niter   !< doubling steps or iterations taken
    real(dp), allocatable :: Z(:, :), H(:, :), T(:, :), W(:, :), A(:, :), G(:, :), &
                             M(:, :), X(:, :), Xn(:, :), Hinv(:, :), Fi(:, :)
    real(dp) :: tol_, logdet
    integer :: maxit, it, mm, stat

    tol_ = 1.0e-12_dp
    if (present(tol)) tol_ = tol
    maxit = 100000
    if (present(maxiter)) maxit = maxiter
    info = SS_ERR_UNSUPPORTED
    if (size(rep%Z, 3) > 1 .or. size(rep%H, 3) > 1 .or. size(rep%T, 3) > 1 .or. &
        size(rep%R, 3) > 1 .or. size(rep%Q, 3) > 1) return

    mm = rep%k_states
    Z = rep%Z(:, :, 1); H = rep%H(:, :, 1); T = rep%T(:, :, 1)
    W = matmul(rep%R(:, :, 1), matmul(rep%Q(:, :, 1), transpose(rep%R(:, :, 1))))
    Hinv = H
    call chol_inv(Hinv, logdet, stat)
    info = SS_ERR_NOT_CONVERGED
    if (stat == SS_OK) then
      ! Doubling
      A = transpose(T)
      G = matmul(transpose(Z), matmul(Hinv, Z))
      X = W
      do it = 1, 100
        M = eye(mm) + matmul(G, X)
        Fi = A                       ! (I + G_k X_k)^-1 A_k
        call solve(M, Fi, stat)
        if (stat /= SS_OK) exit
        Xn = X + matmul(transpose(A), matmul(X, Fi))
        G = G + matmul(A, matmul(inverse_times(eye(mm) + matmul(G, X), G), &
                                 transpose(A)))
        A = matmul(A, Fi)
        call symmetrize(Xn); call symmetrize(G)
        if (converged(X, Xn)) then
          X = Xn
          info = SS_OK
          exit
        end if
        X = Xn
      end do
    else
      ! The recursion itself (singular H), from P_2 = R Q R'
      X = W
      do it = 1, maxit
        Fi = matmul(Z, matmul(X, transpose(Z))) + H
        call chol_inv(Fi, logdet, stat)
        if (stat /= SS_OK) then
          info = SS_ERR_NOT_PD
          exit
        end if
        M = matmul(T, matmul(X, transpose(Z)))          ! T P Z'
        Xn = matmul(T, matmul(X, transpose(T))) - matmul(M, matmul(Fi, transpose(M))) &
             + W
        call symmetrize(Xn)
        if (converged(X, Xn)) then
          X = Xn
          info = SS_OK
          exit
        end if
        X = Xn
      end do
    end if
    if (present(niter)) niter = it
    P = X
    F = matmul(Z, matmul(X, transpose(Z))) + H

  contains

    logical function converged(Xo, Xnew)
      real(dp), intent(in), contiguous :: Xo(:, :), Xnew(:, :)

      converged = maxval(abs(Xnew - Xo)) <= tol_ * max(1.0_dp, maxval(abs(Xnew)))
    end function converged

    !> B^-1 C for square B.
    function inverse_times(B, C) result(D)
      real(dp), intent(in), contiguous :: B(:, :), C(:, :)
      real(dp), allocatable :: D(:, :)

      D = C
      call solve(B, D, stat)
    end function inverse_times
  end subroutine steady_state
end module statespace_filter
