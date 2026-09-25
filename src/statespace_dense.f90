!> Matrix formulation of the model (DK 4.13): the whole sample as one
!> Gaussian vector. With alpha = (alpha_1', ..., alpha_n')' and Y_n the
!> observed elements of y_1..y_n stacked,
!>
!>   E(alpha) = a*,  Var(alpha) = V*   (DK 4.100, built recursively:
!>             Var(alpha_t+1) = T Var(alpha_t) T' + R Q R',
!>             Cov(alpha_t+1, alpha_s) = T Cov(alpha_t, alpha_s), s <= t)
!>   E(Y_n) = mu = d + Z a*,  Var(Y_n) = Omega = Z V* Z' + H   (DK 4.102)
!>
!>   log L = -(N log 2 pi + log|Omega| + (Y_n - mu)' Omega^-1 (Y_n - mu)) / 2
!>   E(alpha | Y_n) = a* + V* Z' Omega^-1 (Y_n - mu)
!>   Var(alpha | Y_n) = V* - V* Z' Omega^-1 Z V*.
!>
!> This costs O((n m)^3) and is meant as a reference for small problems and
!> for checking the recursive algorithms. Known or approximate diffuse
!> initialization only (no P_inf).
module statespace_dense
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use statespace_kinds, only: dp, log2pi, SS_OK, SS_ERR_UNSUPPORTED
  use statespace_linalg, only: chol_inv
  use statespace_rep, only: ssm_rep_t, tidx
  implicit none
  private

  public :: dense_loglike_smooth

contains

  subroutine dense_loglike_smooth(rep, llf, alphahat, V, info)
    type(ssm_rep_t), intent(in) :: rep
    real(dp), intent(out) :: llf
    real(dp), intent(out) :: alphahat(:, :)     !< (m, n)
    real(dp), intent(out) :: V(:, :, :)         !< (m, m, n)
    integer, intent(out) :: info
    real(dp), allocatable :: a1(:), Pstar(:, :), Pinf(:, :), mu_a(:), Va(:, :), W(:, :)
    real(dp), allocatable :: Zb(:, :), Omega(:, :), mu_y(:), yv(:), e(:), VZ(:, :), post(:, :)
    integer, allocatable :: row_t(:), row_i(:)
    real(dp) :: logdet
    integer :: m, n, t, s, i, k, nobs_total, it, ir, iq, ic, iz, ih, id
    integer :: bt, bs

    m = rep%k_states; n = rep%nobs
    call rep%validate(info)
    if (info /= SS_OK) return
    allocate (a1(m), Pstar(m, m), Pinf(m, m))
    call rep%initial_state(a1, Pstar, Pinf, info)
    if (info /= SS_OK) return
    info = SS_ERR_UNSUPPORTED
    if (any(Pinf /= 0.0_dp)) return

    ! E(alpha) and Var(alpha), block by block.
    allocate (mu_a(n * m), Va(n * m, n * m))
    mu_a(1:m) = a1
    Va(1:m, 1:m) = Pstar
    do t = 1, n - 1
      it = tidx(size(rep%T, 3), t)
      ir = tidx(size(rep%R, 3), t)
      iq = tidx(size(rep%Q, 3), t)
      ic = tidx(size(rep%c, 2), t)
      bt = t * m        ! offset of block t+1
      mu_a(bt + 1:bt + m) = rep%c(:, ic) + matmul(rep%T(:, :, it), mu_a(bt - m + 1:bt))
      do s = 1, t
        bs = (s - 1) * m
        ! Cov(alpha_t+1, alpha_s) = T Cov(alpha_t, alpha_s)
        Va(bt + 1:bt + m, bs + 1:bs + m) = matmul(rep%T(:, :, it), Va(bt - m + 1:bt, bs + 1:bs + m))
        Va(bs + 1:bs + m, bt + 1:bt + m) = transpose(Va(bt + 1:bt + m, bs + 1:bs + m))
      end do
      W = matmul(rep%R(:, :, ir), matmul(rep%Q(:, :, iq), transpose(rep%R(:, :, ir))))
      Va(bt + 1:bt + m, bt + 1:bt + m) = matmul(rep%T(:, :, it), &
                                                matmul(Va(bt - m + 1:bt, bt - m + 1:bt), &
                                                       transpose(rep%T(:, :, it)))) + W
    end do

    ! Observed elements of Y_n, and Z, mu, Omega over them.
    nobs_total = count(.not. ieee_is_nan(rep%y))
    allocate (row_t(nobs_total), row_i(nobs_total))
    allocate (Zb(nobs_total, n * m), mu_y(nobs_total), yv(nobs_total), source=0.0_dp)
    k = 0
    do t = 1, n
      iz = tidx(size(rep%Z, 3), t)
      id = tidx(size(rep%d, 2), t)
      do i = 1, rep%k_endog
        if (ieee_is_nan(rep%y(i, t))) cycle
        k = k + 1
        row_t(k) = t
        row_i(k) = i
        Zb(k, (t - 1) * m + 1:t * m) = rep%Z(i, :, iz)
        mu_y(k) = rep%d(i, id) + dot_product(rep%Z(i, :, iz), mu_a((t - 1) * m + 1:t * m))
        yv(k) = rep%y(i, t)
      end do
    end do
    VZ = matmul(Va, transpose(Zb))
    Omega = matmul(Zb, VZ)
    do k = 1, nobs_total
      do i = 1, nobs_total
        if (row_t(i) == row_t(k)) then
          ih = tidx(size(rep%H, 3), row_t(k))
          Omega(i, k) = Omega(i, k) + rep%H(row_i(i), row_i(k), ih)
        end if
      end do
    end do

    call chol_inv(Omega, logdet, info)   ! Omega <- Omega^-1
    if (info /= SS_OK) return
    e = yv - mu_y
    llf = -0.5_dp * (nobs_total * log2pi + logdet + dot_product(e, matmul(Omega, e)))

    mu_a = mu_a + matmul(VZ, matmul(Omega, e))
    post = Va - matmul(VZ, matmul(Omega, transpose(VZ)))
    do t = 1, n
      bt = (t - 1) * m
      alphahat(:, t) = mu_a(bt + 1:bt + m)
      V(:, :, t) = post(bt + 1:bt + m, bt + 1:bt + m)
    end do
    info = SS_OK
  end subroutine dense_loglike_smooth
end module statespace_dense
