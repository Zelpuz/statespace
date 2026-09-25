!> Analytic score vector (DK 7.3.3) from one smoother pass.
!>
!> By the Fisher identity, d log L / d psi = E[d log p(Y, alpha; psi) / d psi | Y]
!> (DK 7.13-7.15). For parameters in H_t, R_t, Q_t and the P_star part of the
!> initialization this gives (DK 7.16, Koopman and Shephard 1992)
!>
!>   d log L / d psi_i = 1/2 sum_t tr(dH_t G^eps_t)
!>                       + 1/2 sum_t<n tr(d(R_t Q_t R_t') (r_t r_t' - N_t))
!>                       + 1/2 tr(dP_star G^1),
!>
!>   G^eps_t = H_t^-1 (epshat epshat' + Var(eps_t|Y) - H_t) H_t^-1  (= u u' - D),
!>   G^1     = P_star^+ (V_1 + (alphahat_1 - a_1)(...)' - P_star) P_star^+.
!>
!> The derivatives dH, d(RQR'), dP_star are central differences of what the
!> model's `update` sets, so no extra user code is needed. DK recommend
!> numerical scores for parameters in Z_t and T_t: if a parameter moves Z, T,
!> c, d, a_1 or P_inf, SS_ERR_UNSUPPORTED is returned (use numerical
!> derivatives). With missing
!> data the eps moments of the missing elements are the conditional ones
!> (see `smoother_result_t`), which this identity needs. It also holds for the
!> diffuse log likelihood (DK 7.3.5), and for a concentrated scale (envelope
!> theorem), but not with a log likelihood burn-in.
module statespace_score
  use statespace_kinds, only: dp, SS_OK, SS_ERR_UNSUPPORTED, SS_ERR_NOT_PD
  use statespace_linalg, only: chol_inv, psd_solve, symmetrize, ldl_psd
  use statespace_rep, only: ssm_rep_t
  use statespace_filter, only: filter_result_t, kalman_filter
  use statespace_smoother, only: smoother_result_t, state_smoother
  use statespace_model, only: ssm_model_t
  implicit none
  private

  public :: analytic_score

contains

  !> Score of the log likelihood at constrained `params`. `llf` receives the
  !> log likelihood at `params` from the same filter run.
  subroutine analytic_score(model, params, score, llf, info)
    class(ssm_model_t), intent(inout) :: model
    real(dp), intent(in) :: params(:)
    real(dp), intent(out) :: score(:)
    real(dp), intent(out) :: llf
    integer, intent(out) :: info
    type(ssm_rep_t) :: rep, plus, minus
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: Gh(:, :, :), Gw(:, :, :), G1(:, :), dH(:, :, :), dW(:, :, :)
    real(dp), allocatable :: dPs(:, :), p(:), a1(:), Pstar(:, :), Pinf(:, :)
    real(dp), allocatable :: a1p(:), Psp(:, :), Pip(:, :), a1m(:), Psm(:, :), Pim(:, :)
    logical, allocatable :: h_singular(:)
    real(dp) :: h, scale
    integer :: i, k, m

    score = 0.0_dp
    llf = 0.0_dp
    k = size(params)
    m = model%rep%k_states
    call model%rep_at(params, rep, info)
    if (info /= SS_OK) return
    scale = 1.0_dp
    if (model%concentrate_scale) scale = model%scale
    info = SS_ERR_UNSUPPORTED
    if (rep%loglikelihood_burn > 0) return

    call kalman_filter(rep, fres, info)
    if (info /= SS_OK) return
    call state_smoother(rep, fres, sres, info)
    if (info /= SS_OK) return
    llf = fres%llf

    allocate (a1(m), Pstar(m, m), Pinf(m, m), a1p(m), Psp(m, m), Pip(m, m), a1m(m), &
              Psm(m, m), Pim(m, m))
    call rep%initial_state(a1, Pstar, Pinf, info)
    if (info /= SS_OK) return
    call score_weights(rep, sres, a1, Pstar, Gh, Gw, G1, h_singular, info)
    if (info /= SS_OK) return

    ! dH, d(RQR'), dP_star by central differences through `update`.
    do i = 1, k
      h = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * max(abs(params(i)), 1.0_dp)
      p = params
      p(i) = params(i) + h
      call model%update(p)
      plus = model%rep
      call plus%initial_state(a1p, Psp, Pip, info)
      if (info /= SS_OK) exit
      p(i) = params(i) - h
      call model%update(p)
      minus = model%rep
      call minus%initial_state(a1m, Psm, Pim, info)
      if (info /= SS_OK) exit

      info = SS_ERR_UNSUPPORTED
      if (any(plus%Z /= minus%Z) .or. any(plus%T /= minus%T) .or. &
          any(plus%c /= minus%c) .or. any(plus%d /= minus%d) .or. any(a1p /= a1m) .or. &
          any(Pip /= Pim)) exit
      info = SS_OK

      dH = scale * (plus%H - minus%H) / (2.0_dp * h)
      dW = scale * (rqr(plus) - rqr(minus)) / (2.0_dp * h)
      dPs = scale * (Psp - Psm) / (2.0_dp * h)
      ! A parameter moving a singular H_t or Q_t is outside this formula.
      info = SS_ERR_UNSUPPORTED
      if (moves_singular(dH, h_singular)) exit
      info = SS_OK
      score(i) = 0.5_dp * (contract(dH, Gh) + contract(dW, Gw) + sum(dPs * G1))
    end do
    call model%update(params)
    if (info /= SS_OK) score = 0.0_dp
  end subroutine analytic_score

  !> The weights of the score formula, one slice per period: Gh as above and
  !> Gw_t = r_t r_t' - N_t (zero for t = n: eta_n does not affect the data).
  subroutine score_weights(rep, sres, a1, Pstar, Gh, Gw, G1, h_singular, info)
    type(ssm_rep_t), intent(in) :: rep
    type(smoother_result_t), intent(in) :: sres
    real(dp), intent(in) :: a1(:), Pstar(:, :)
    real(dp), allocatable, intent(out) :: Gh(:, :, :), Gw(:, :, :), G1(:, :)
    logical, allocatable, intent(out) :: h_singular(:)   !< (n)
    integer, intent(out) :: info
    real(dp), allocatable :: Hinv(:, :), M(:, :), e1(:)
    real(dp) :: logdet
    integer :: t, n, ih, iq

    n = rep%nobs
    allocate (Gh(rep%k_endog, rep%k_endog, n), Gw(rep%k_states, rep%k_states, n))
    allocate (h_singular(n), source=.false.)
    info = SS_OK
    do t = 1, n
      ih = min(t, size(rep%H, 3))
      iq = min(t, size(rep%Q, 3))
      ! G^eps_t = H^-1 (epshat epshat' + Var - H) H^-1. Singular H_t is
      ! flagged; the caller rejects parameters that move it.
      Hinv = rep%H(:, :, ih)
      call chol_inv(Hinv, logdet, info)
      if (info /= SS_OK) then
        if (.not. is_psd(rep%H(:, :, ih))) then
          info = SS_ERR_NOT_PD        ! invalid parameters, not a scope limit
          return
        end if
        Gh(:, :, t) = 0.0_dp
        h_singular(t) = .true.
      else
        M = outer(sres%epshat(:, t)) + sres%epsvar(:, :, t) - rep%H(:, :, ih)
        Gh(:, :, t) = matmul(Hinv, matmul(M, Hinv))
      end if
      ! An indefinite Q means invalid parameters.
      if (.not. is_psd(rep%Q(:, :, iq))) then
        info = SS_ERR_NOT_PD
        return
      end if
      if (t < n) then
        Gw(:, :, t) = outer(sres%r(:, t)) - sres%N(:, :, t)
      else
        Gw(:, :, t) = 0.0_dp
      end if
    end do
    info = SS_OK

    ! G^1 = P_star^+ (V_1 + e e' - P_star) P_star^+, e = alphahat_1 - a1
    e1 = sres%alphahat(:, 1) - a1
    M = sres%V(:, :, 1) + outer(e1) - Pstar
    G1 = psd_solve(Pstar, transpose(psd_solve(Pstar, M)))
    call symmetrize(G1)
  end subroutine score_weights

  !> R_t Q_t R_t', one slice per distinct slice of R and Q.
  function rqr(rep) result(W)
    type(ssm_rep_t), intent(in) :: rep
    real(dp), allocatable :: W(:, :, :)
    integer :: t, nw, ir, iq

    nw = max(size(rep%R, 3), size(rep%Q, 3))
    allocate (W(rep%k_states, rep%k_states, nw))
    do t = 1, nw
      ir = min(t, size(rep%R, 3))
      iq = min(t, size(rep%Q, 3))
      W(:, :, t) = matmul(rep%R(:, :, ir), matmul(rep%Q(:, :, iq), transpose(rep%R(:, :, ir))))
    end do
  end function rqr

  !> sum_t tr(dX_t G_t) for symmetric G_t, with dX time-invariant (one slice)
  !> or time-varying.
  real(dp) function contract(dX, G)
    real(dp), intent(in) :: dX(:, :, :), G(:, :, :)
    integer :: t

    contract = 0.0_dp
    if (size(dX, 3) == 1) then
      contract = sum(dX(:, :, 1) * sum(G, dim=3))
    else
      do t = 1, size(G, 3)
        contract = contract + sum(dX(:, :, t) * G(:, :, t))
      end do
    end if
  end function contract

  !> True if A is positive semi-definite (LDL' pivots >= 0).
  logical function is_psd(A)
    real(dp), intent(in) :: A(:, :)
    real(dp), allocatable :: L(:, :), D(:)
    integer :: i

    ! ldl_psd zeroes small or negative pivots, so check them on its input
    ! diagonal recursion directly: any negative diagonal element, or a
    ! reconstruction L D L' that differs from A, means A is indefinite.
    allocate (L(size(A, 1), size(A, 1)), D(size(A, 1)))
    is_psd = all([(A(i, i) >= 0.0_dp, i=1, size(A, 1))])
    if (.not. is_psd) return
    call ldl_psd(A, L, D)
    is_psd = maxval(abs(matmul(L, matmul(diag(D), transpose(L))) - A)) &
             <= 1.0e-8_dp * max(maxval(abs(A)), tiny(1.0_dp))
  end function is_psd

  pure function diag(x) result(A)
    real(dp), intent(in) :: x(:)
    real(dp) :: A(size(x), size(x))
    integer :: i

    A = 0.0_dp
    do i = 1, size(x)
      A(i, i) = x(i)
    end do
  end function diag

  !> True if dX is nonzero in a period flagged singular.
  logical function moves_singular(dX, singular)
    real(dp), intent(in) :: dX(:, :, :)
    logical, intent(in) :: singular(:)
    integer :: t

    moves_singular = .false.
    do t = 1, size(singular)
      if (singular(t) .and. any(dX(:, :, min(t, size(dX, 3))) /= 0.0_dp)) then
        moves_singular = .true.
        return
      end if
    end do
  end function moves_singular

  pure function outer(x) result(A)
    real(dp), intent(in) :: x(:)
    real(dp) :: A(size(x), size(x))
    integer :: j

    do j = 1, size(x)
      A(:, j) = x * x(j)
    end do
  end function outer
end module statespace_score
