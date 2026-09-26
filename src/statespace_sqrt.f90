!> Square root filter and smoother (DK 6.3).
!>
!> The filter propagates a square root P~_t of P_t (P = P~ P~') rather than
!> P_t itself, so P_t stays positive semi-definite by construction. Each step
!> triangularizes the pre-array
!>
!>     U = [ Z_o P~   H~_oo   0     ]      [ F~   0       0 ]
!>         [ T P~     0       R Q~  ]  ->  [ K~   P~_t+1  0 ]  = U G
!>
!> for orthogonal G (LQ via Householder QR of U'), where H~ H~' = H,
!> Q~ Q~' = Q and o are the observed elements. Then F = F~ F~',
!> K = K~ F~^-1 and P_t+1 = P~_t+1 P~_t+1'. The smoother does the same for
!> N_t-1 = Z' F^-1 Z + L' N_t L = [Z' F~^-T, L' N~_t] [...]', giving N~_t-1.
!>
!> The output uses the same types as `kalman_filter` and `state_smoother`
!> (with P, F, N formed from their square roots), so the rest of the library
!> works on it. Diffuse initialization is not supported; use approximate
!> diffuse or the exact initial filter.
module statespace_sqrt
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use statespace_kinds, only: dp, log2pi, SS_OK, SS_ERR_UNSUPPORTED, SS_ERR_NOT_PD
  use statespace_linalg, only: tria, chol_inv, symmetrize
  use statespace_rep, only: ssm_rep_t, tidx
  use statespace_filter, only: filter_result_t, METHOD_CONVENTIONAL
  use statespace_smoother, only: smoother_result_t, state_smoother
  use statespace_simsmooth, only: psd_sqrt
  implicit none
  private

  public :: sqrt_kalman_filter, sqrt_state_smoother

contains

  !> Square root Kalman filter. `res` holds the same quantities as from
  !> `kalman_filter`; `Pchol`, if present, receives the square roots P~_t.
  subroutine sqrt_kalman_filter(rep, res, info, Pchol)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(out) :: res
    integer, intent(out) :: info
    real(dp), allocatable, intent(out), optional :: Pchol(:, :, :)   !< (m, m, n+1)
    real(dp), allocatable :: a1(:), Pstar(:, :), Pinf(:, :), Ps(:, :, :), Hs(:, :), Qs(:, :)
    real(dp), allocatable :: U(:, :), L(:, :), Zo(:, :), Fs(:, :), Kt(:, :), Finv(:, :)
    real(dp), allocatable :: vo(:), Fv(:), PZt(:, :)
    integer, allocatable :: idx(:)
    real(dp) :: logdet
    integer :: p, m, r, n, t, i, n_o, iz, ih, it, ir, iq, ic, id

    call rep%validate(info)
    if (info /= SS_OK) return
    p = rep%k_endog; m = rep%k_states; r = rep%k_posdef; n = rep%nobs
    allocate (a1(m), Pstar(m, m), Pinf(m, m))
    call rep%initial_state(a1, Pstar, Pinf, info)
    if (info /= SS_OK) return
    info = SS_ERR_UNSUPPORTED
    if (any(Pinf /= 0.0_dp)) return
    info = SS_OK

    res%k_endog = p; res%k_states = m; res%nobs = n
    allocate (res%a(m, n + 1), res%P(m, m, n + 1), res%Pinf(m, m, n + 1), res%att(m, n), &
              res%Ptt(m, m, n), res%yhat(p, n), res%v(p, n), res%F(p, p, n), res%Finf(p, p, n), &
              res%Finv(p, p, n), res%K(m, p, n), res%llf_obs(n), res%method(n), &
              res%uv_n(n), res%uv_idx(p, 0), res%uv_Z(p, m, 0), res%uv_sig2(p, 0), &
              res%uv_v(p, 0), res%uv_Fstar(p, 0), res%uv_Finf(p, 0), res%uv_Mstar(m, p, 0), &
              res%uv_Minf(m, p, 0))
    res%Pinf = 0.0_dp
    res%Finf = 0.0_dp
    res%method = METHOD_CONVENTIONAL
    res%uv_n = 0
    allocate (Ps(m, m, n + 1))

    res%a(:, 1) = a1
    Ps(:, :, 1) = psd_sqrt(Pstar)
    do t = 1, n
      iz = tidx(size(rep%Z, 3), t)
      ih = tidx(size(rep%H, 3), t)
      it = tidx(size(rep%T, 3), t)
      ir = tidx(size(rep%R, 3), t)
      iq = tidx(size(rep%Q, 3), t)
      ic = tidx(size(rep%c, 2), t)
      id = tidx(size(rep%d, 2), t)
      idx = pack([(i, i=1, p)], .not. ieee_is_nan(rep%y(:, t)))
      n_o = size(idx)

      associate (P => Ps(:, :, t), a => res%a(:, t))
        res%P(:, :, t) = matmul(P, transpose(P))
        res%yhat(:, t) = rep%d(:, id) + matmul(rep%Z(:, :, iz), a)
        res%v(:, t) = rep%y(:, t) - res%yhat(:, t)
        res%F(:, :, t) = matmul(rep%Z(:, :, iz), matmul(res%P(:, :, t), &
                                transpose(rep%Z(:, :, iz)))) + rep%H(:, :, ih)
        call symmetrize(res%F(:, :, t))

        ! Pre-array [Z_o P~, H~, 0; T P~, 0, R Q~] and its triangularization.
        Zo = rep%Z(idx, :, iz)
        Hs = psd_sqrt(rep%H(idx, idx, ih))
        Qs = psd_sqrt(rep%Q(:, :, iq))
        allocate (U(n_o + m, m + n_o + r), source=0.0_dp)
        U(1:n_o, 1:m) = matmul(Zo, P)
        U(1:n_o, m + 1:m + n_o) = Hs
        U(n_o + 1:, 1:m) = matmul(rep%T(:, :, it), P)
        U(n_o + 1:, m + n_o + 1:) = matmul(rep%R(:, :, ir), Qs)
        L = tria(U)
        deallocate (U)
        Fs = L(1:n_o, 1:n_o)
        Kt = L(n_o + 1:, 1:n_o)
        Ps(:, :, t + 1) = L(n_o + 1:, n_o + 1:)

        ! F^-1 over the observed block from F = F~ F~', K = K~ F~^-1.
        Finv = matmul(Fs, transpose(Fs))
        call chol_inv(Finv, logdet, info)
        if (info /= SS_OK) then
          info = SS_ERR_NOT_PD
          return
        end if
        res%Finv(:, :, t) = 0.0_dp
        res%Finv(idx, idx, t) = Finv
        res%K(:, :, t) = 0.0_dp
        res%K(:, idx, t) = matmul(Kt, matmul(transpose(Fs), Finv))

        vo = res%v(idx, t)
        Fv = matmul(Finv, vo)
        PZt = matmul(res%P(:, :, t), transpose(Zo))
        res%att(:, t) = a + matmul(PZt, Fv)
        res%Ptt(:, :, t) = res%P(:, :, t) - matmul(PZt, matmul(Finv, transpose(PZt)))
        call symmetrize(res%Ptt(:, :, t))
        res%a(:, t + 1) = rep%c(:, ic) + matmul(rep%T(:, :, it), a) &
                          + matmul(res%K(:, idx, t), vo)
        if (n_o == 0) then
          res%llf_obs(t) = 0.0_dp
        else
          res%llf_obs(t) = -0.5_dp * (n_o * log2pi + logdet + dot_product(vo, Fv))
        end if
      end associate
    end do
    res%P(:, :, n + 1) = matmul(Ps(:, :, n + 1), transpose(Ps(:, :, n + 1)))
    res%llf = sum(res%llf_obs(rep%loglikelihood_burn + 1:))
    if (present(Pchol)) Pchol = Ps
  end subroutine sqrt_kalman_filter

  !> Square root state smoother: r_t and alphahat as in `state_smoother`, with
  !> N_t-1 propagated through its square root N~_t-1 = tria([Z_o' F~^-T, L' N~_t])
  !> and V_t = P_t - (P_t N~_t-1)(P_t N~_t-1)'.
  subroutine sqrt_state_smoother(rep, fres, sres, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    type(smoother_result_t), intent(out) :: sres
    integer, intent(out) :: info
    real(dp), allocatable :: Ns(:, :), U(:, :), Lt(:, :), Zo(:, :), Fchol(:, :), PN(:, :)
    integer, allocatable :: idx(:)
    integer :: m, n, t, i, n_o, iz, it

    info = SS_ERR_UNSUPPORTED
    if (fres%nobs_diffuse > 0 .or. any(fres%method /= METHOD_CONVENTIONAL)) return
    ! The means (r_t, alphahat, disturbances) do not involve N.
    call state_smoother(rep, fres, sres, info)
    if (info /= SS_OK) return

    m = rep%k_states; n = rep%nobs
    allocate (Ns(m, m), source=0.0_dp)
    do t = n, 1, -1
      iz = tidx(size(rep%Z, 3), t)
      it = tidx(size(rep%T, 3), t)
      idx = pack([(i, i=1, rep%k_endog)], .not. ieee_is_nan(rep%y(:, t)))
      n_o = size(idx)
      Zo = rep%Z(idx, :, iz)
      ! Z_o' F~^-T with F~ the lower Cholesky factor of F_oo
      Fchol = psd_sqrt(fres%F(idx, idx, t))
      Lt = rep%T(:, :, it) - matmul(fres%K(:, :, t), rep%Z(:, :, iz))
      allocate (U(m, n_o + m))
      U(:, 1:n_o) = matmul(transpose(Zo), transpose(invert_lower(Fchol)))
      U(:, n_o + 1:) = matmul(transpose(Lt), Ns)
      Ns = tria(U)
      deallocate (U)
      sres%N(:, :, t - 1) = matmul(Ns, transpose(Ns))
      PN = matmul(fres%P(:, :, t), Ns)
      sres%V(:, :, t) = fres%P(:, :, t) - matmul(PN, transpose(PN))
      call symmetrize(sres%V(:, :, t))
    end do
    info = SS_OK
  end subroutine sqrt_state_smoother

  !> Inverse of a nonsingular lower triangular matrix by forward substitution.
  function invert_lower(L) result(X)
    real(dp), intent(in), contiguous :: L(:, :)
    real(dp), allocatable :: X(:, :)
    integer :: i, j, n

    n = size(L, 1)
    allocate (X(n, n), source=0.0_dp)
    do j = 1, n
      X(j, j) = 1.0_dp / L(j, j)
      do i = j + 1, n
        X(i, j) = -dot_product(L(i, j:i - 1), X(j:i - 1, j)) / L(i, i)
      end do
    end do
  end function invert_lower
end module statespace_sqrt
