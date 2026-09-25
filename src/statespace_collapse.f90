!> Collapsing large observation vectors (DK 6.5; Jungbacker and Koopman 2008).
!>
!> When p is large relative to m, each y_t is replaced by its GLS projection
!> on the state:
!>
!>     ybar_t = (Z' H^-1 Z)^-1 Z' H^-1 (y_t - d_t) = alpha_t + epsbar_t,
!>     Var(epsbar_t) = Hbar_t = (Z' H^-1 Z)^-1,
!>
!> over the observed elements of y_t. The GLS residual
!> e_t = (y_t - d_t) - Z ybar_t is independent of ybar_t and of the states,
!> so the collapsed model (observation y*_t = ybar_t, Z* = I, H* = Hbar_t,
!> d* = 0) gives the same filtered and smoothed states, and
!>
!>     log L = log L* + sum_t [ -(p_t - m)/2 log 2 pi - log(|H_t| / |Hbar_t|) / 2
!>                              - e_t' H_t^-1 e_t / 2 ].
!>
!> Every period with observations needs H_t (observed block) nonsingular and
!> Z_t (observed rows) of full column rank m.
module statespace_collapse
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use statespace_kinds, only: dp, log2pi, SS_OK, SS_ERR_NOT_PD
  use statespace_linalg, only: chol_inv, eye, symmetrize
  use statespace_rep, only: ssm_rep_t, tidx
  implicit none
  private

  public :: collapse_observations

contains

  !> Build the collapsed representation `crep` and the per-period log
  !> likelihood adjustment, so that log L = loglike(crep) + sum(llf_adjust).
  subroutine collapse_observations(rep, crep, llf_adjust, info)
    type(ssm_rep_t), intent(in) :: rep
    type(ssm_rep_t), intent(out) :: crep
    real(dp), intent(out) :: llf_adjust(:)   !< (n)
    integer, intent(out) :: info
    real(dp), allocatable :: Zo(:, :), Hinv(:, :), ZtHinv(:, :), Hbar(:, :), yo(:), ybar(:), e(:)
    integer, allocatable :: idx(:)
    real(dp) :: logdetH, logdetHbarinv, nan
    integer :: m, n, t, i, n_o, iz, ih, id

    call rep%validate(info)
    if (info /= SS_OK) return
    m = rep%k_states; n = rep%nobs
    nan = ieee_value(1.0_dp, ieee_quiet_nan)

    crep = rep
    crep%k_endog = m
    deallocate (crep%y, crep%Z, crep%H, crep%d)
    allocate (crep%y(m, n), crep%H(m, m, n), crep%d(m, 1), source=0.0_dp)
    crep%Z = reshape(eye(m), [m, m, 1])

    do t = 1, n
      iz = tidx(size(rep%Z, 3), t)
      ih = tidx(size(rep%H, 3), t)
      id = tidx(size(rep%d, 2), t)
      idx = pack([(i, i=1, rep%k_endog)], .not. ieee_is_nan(rep%y(:, t)))
      n_o = size(idx)
      llf_adjust(t) = 0.0_dp
      if (n_o == 0) then
        crep%y(:, t) = nan
        crep%H(:, :, t) = eye(m)
        cycle
      end if

      Zo = rep%Z(idx, :, iz)
      yo = rep%y(idx, t) - rep%d(idx, id)
      Hinv = rep%H(idx, idx, ih)
      call chol_inv(Hinv, logdetH, info)
      if (info /= SS_OK) return
      ZtHinv = matmul(transpose(Zo), Hinv)
      Hbar = matmul(ZtHinv, Zo)             ! Hbar^-1 = Z' H^-1 Z
      call symmetrize(Hbar)
      call chol_inv(Hbar, logdetHbarinv, info)
      if (info /= SS_OK) then
        info = SS_ERR_NOT_PD                ! Z_o not of full column rank
        return
      end if
      ybar = matmul(Hbar, matmul(ZtHinv, yo))
      e = yo - matmul(Zo, ybar)

      crep%y(:, t) = ybar
      crep%H(:, :, t) = Hbar
      ! log|Hbar| = -log|Hbar^-1|
      llf_adjust(t) = -0.5_dp * ((n_o - m) * log2pi + logdetH + logdetHbarinv &
                                 + dot_product(e, matmul(Hinv, e)))
    end do
  end subroutine collapse_observations
end module statespace_collapse
