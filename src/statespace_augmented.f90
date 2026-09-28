!> Augmented Kalman filter and smoother (DK 5.7), an alternative to the exact
!> initial Kalman filter for diffuse initialization, and the basis of
!> regression estimation in DK 6.2.
!>
!> With alpha_1 = a + A delta + R_0 eta_0 and delta unknown (diffuse), the
!> filter is run once with delta = 0 (giving v*_t, F_t, K_t, P_t), and the
!> effect of delta is carried along as extra columns:
!>
!>     A_1 = A,   V^A_t = -Z_t A_t,   A_t+1 = T_t A_t + K_t V^A_t,
!>     s = sum_t V^A_t' F_t^-1 v*_t,   S = sum_t V^A_t' F_t^-1 V^A_t,
!>
!> so that v_t(delta) = v*_t + V^A_t delta. The GLS estimate of delta is
!> delta_hat = -S^-1 s with variance S^-1, and the diffuse log likelihood
!> (DK 7.2.2-7.2.3) is
!>
!>     log L_d = log L* + s' S^-1 s / 2 - log|S| / 2,
!>
!> where log L* is the log likelihood of the delta = 0 filter. Smoothing uses
!> R^A_t-1 = Z_t' F_t^-1 V^A_t + L_t' R^A_t alongside r*_t and N_t:
!>
!>     alphahat_t = a*_t + A_t delta_hat + P_t (r*_t-1 + R^A_t-1 delta_hat)
!>     V_t = P_t - P_t N_t-1 P_t + B_t S^-1 B_t',   B_t = A_t + P_t R^A_t-1,
!>
!> the second from the law of total variance over delta | Y_n ~
!> N(delta_hat, S^-1). The filter needs F_t = Z P Z' + H nonsingular over the
!> observed elements, which the exact initial filter does not.
module statespace_augmented
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use statespace_kinds, only: dp, SS_OK, SS_ERR_NOT_PD
  use statespace_linalg, only: gemm, gemv, chol_inv, ldl_psd, symmetrize
  use statespace_rep, only: ssm_rep_t, tidx
  use statespace_filter, only: filter_result_t, kalman_filter, marginal_correction
  use statespace_smoother, only: smoother_result_t, state_smoother
  implicit none
  private

  public :: augmented_result_t, augmented_filter, augmented_smoother

  type :: augmented_result_t
    integer :: k_diffuse = 0
    type(filter_result_t) :: filter        !< the delta = 0 filter
    real(dp), allocatable :: A(:, :, :)    !< (m, k, n+1) A_t
    real(dp), allocatable :: VA(:, :, :)   !< (p, k, n) V^A_t
    real(dp), allocatable :: s(:)          !< (k)
    real(dp), allocatable :: S_mat(:, :)   !< (k, k)
    real(dp), allocatable :: delta(:)      !< (k) GLS estimate of delta
    real(dp), allocatable :: delta_cov(:, :) !< (k, k) its variance, S^-1
    real(dp) :: llf = 0.0_dp               !< diffuse log likelihood (DK 7.7)
    !> Log likelihood with delta fixed but unknown, concentrated at
    !> delta_hat (DK 7.2.4, eq. 7.9): llf + log|S| / 2.
    real(dp) :: llf_fixed = 0.0_dp
    !> Marginal log likelihood (Francke, Koopman and de Vos 2010; DK 7.2.6):
    !> llf + k/2 log 2 pi + log|S*| / 2, see `marginal_correction`.
    real(dp) :: llf_marginal = 0.0_dp
  end type augmented_result_t

contains

  !> Run the augmented Kalman filter.
  subroutine augmented_filter(rep, res, info)
    type(ssm_rep_t), intent(in) :: rep
    type(augmented_result_t), intent(out) :: res
    integer, intent(out) :: info
    type(ssm_rep_t) :: star
    real(dp), allocatable :: a1(:), Pstar(:, :), Pinf(:, :), L(:, :), D(:), FV(:, :), &
                             vz(:)
    real(dp) :: logdetS
    integer :: p, m, n, k, t, j, iz, it

    p = rep%k_endog; m = rep%k_states; n = rep%nobs
    call rep%validate(info)
    if (info /= SS_OK) return
    allocate (a1(m), Pstar(m, m), Pinf(m, m), L(m, m), D(m))
    call rep%initial_state(a1, Pstar, Pinf, info)
    if (info /= SS_OK) return

    ! A with A A' = P_inf: the nonzero columns of L D^(1/2).
    call ldl_psd(Pinf, L, D)
    k = count(D > 0.0_dp)
    res%k_diffuse = k
    allocate (res%A(m, k, n + 1), res%VA(p, k, n), res%s(k), res%S_mat(k, k), FV(p, k))
    res%A(:, :, 1) = 0.0_dp
    k = 0
    do j = 1, m
      if (D(j) > 0.0_dp) then
        k = k + 1
        res%A(:, k, 1) = L(:, j) * sqrt(D(j))
      end if
    end do

    ! The delta = 0 filter.
    star = rep
    call star%initialize_general(a1, Pstar, 0.0_dp * Pinf)
    call kalman_filter(star, res%filter, info)
    if (info /= SS_OK) return

    res%s = 0.0_dp
    res%S_mat = 0.0_dp
    associate (f => res%filter)
      do t = 1, n
        iz = tidx(size(rep%Z, 3), t)
        it = tidx(size(rep%T, 3), t)
        ! V^A = -Z A,  A_t+1 = T A + K V^A
        call gemm('N', 'N', -1.0_dp, rep%Z(:, :, iz), res%A(:, :, t), 0.0_dp, &
                  res%VA(:, :, t))
        call gemm('N', 'N', 1.0_dp, rep%T(:, :, it), res%A(:, :, t), 0.0_dp, &
                  res%A(:, :, t + 1))
        call gemm('N', 'N', 1.0_dp, f%K(:, :, t), res%VA(:, :, t), 1.0_dp, &
                  res%A(:, :, t + 1))
        ! s += V^A' F^-1 v*,  S += V^A' F^-1 V^A  (F^-1 is zero-padded over missing)
        vz = merge(0.0_dp, f%v(:, t), ieee_is_nan(f%v(:, t)))
        call gemm('N', 'N', 1.0_dp, f%Finv(:, :, t), res%VA(:, :, t), 0.0_dp, FV)
        call gemv('T', 1.0_dp, FV, vz, 1.0_dp, res%s)
        call gemm('T', 'N', 1.0_dp, res%VA(:, :, t), FV, 1.0_dp, res%S_mat)
      end do
    end associate
    call symmetrize(res%S_mat)

    res%delta_cov = res%S_mat
    call chol_inv(res%delta_cov, logdetS, info)
    if (info /= SS_OK) then
      info = SS_ERR_NOT_PD      ! the data do not identify delta
      return
    end if
    res%delta = -matmul(res%delta_cov, res%s)
    res%llf = res%filter%llf + 0.5_dp * dot_product(res%s, -res%delta) &
              - 0.5_dp * logdetS
    res%llf_fixed = res%llf + 0.5_dp * logdetS
    res%llf_marginal = res%llf + marginal_correction(rep, info)
  end subroutine augmented_filter

  !> Smoothed state and its variance from the augmented filter.
  subroutine augmented_smoother(rep, res, alphahat, V, info)
    type(ssm_rep_t), intent(in) :: rep
    type(augmented_result_t), intent(in) :: res
    real(dp), intent(out), contiguous :: alphahat(:, :)    !< (m, n)
    real(dp), intent(out), contiguous :: V(:, :, :)        !< (m, m, n)
    integer, intent(out) :: info
    type(ssm_rep_t) :: star
    type(smoother_result_t) :: sres
    real(dp), allocatable :: RA(:, :), RAprev(:, :), Lt(:, :), B(:, :), BS(:, :), &
                             ZtFinv(:, :)
    real(dp), allocatable :: a1(:), Pstar(:, :), Pinf(:, :)
    integer :: m, k, n, t, iz, it

    m = rep%k_states; n = rep%nobs; k = res%k_diffuse
    allocate (a1(m), Pstar(m, m), Pinf(m, m))
    call rep%initial_state(a1, Pstar, Pinf, info)
    if (info /= SS_OK) return
    star = rep
    call star%initialize_general(a1, Pstar, 0.0_dp * Pinf)
    call state_smoother(star, res%filter, sres, info)
    if (info /= SS_OK) return

    allocate (RA(m, k), RAprev(m, k), ZtFinv(m, rep%k_endog), B(m, k), BS(m, k), &
              source=0.0_dp)
    associate (f => res%filter)
      do t = n, 1, -1
        iz = tidx(size(rep%Z, 3), t)
        it = tidx(size(rep%T, 3), t)
        ! R^A_t-1 = Z' F^-1 V^A + L' R^A_t,  L = T - K Z
        Lt = rep%T(:, :, it)
        call gemm('N', 'N', -1.0_dp, f%K(:, :, t), rep%Z(:, :, iz), 1.0_dp, Lt)
        call gemm('T', 'N', 1.0_dp, rep%Z(:, :, iz), f%Finv(:, :, t), 0.0_dp, ZtFinv)
        call gemm('N', 'N', 1.0_dp, ZtFinv, res%VA(:, :, t), 0.0_dp, RAprev)
        call gemm('T', 'N', 1.0_dp, Lt, RA, 1.0_dp, RAprev)
        RA = RAprev

        ! B = A_t + P_t R^A_t-1
        B = res%A(:, :, t)
        call gemm('N', 'N', 1.0_dp, f%P(:, :, t), RA, 1.0_dp, B)
        alphahat(:, t) = sres%alphahat(:, t)
        call gemv('N', 1.0_dp, B, res%delta, 1.0_dp, alphahat(:, t))
        V(:, :, t) = sres%V(:, :, t)
        call gemm('N', 'N', 1.0_dp, B, res%delta_cov, 0.0_dp, BS)
        call gemm('N', 'T', 1.0_dp, BS, B, 1.0_dp, V(:, :, t))
        call symmetrize(V(:, :, t))
      end do
    end associate
  end subroutine augmented_smoother
end module statespace_augmented
