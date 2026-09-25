!> Filtering and smoothing under linear restrictions on the state (DK 6.6).
!>
!> Restrictions R*_t alpha_t = r*_t (q of them) are imposed by appending them
!> to the observation equation as exact observations:
!>
!>     y+_t = [ y_t  ],   Z+_t = [ Z_t  ],   H+_t = [ H_t  0 ],   d+_t = [ d_t ]
!>            [ r*_t ]           [ R*_t ]           [ 0    0 ]           [ 0   ]
!>
!> A NaN in r*_t leaves that restriction inactive in period t (it is treated
!> as a missing observation). The filtered and smoothed states of the result
!> satisfy the active restrictions exactly. Its log likelihood includes the
!> restriction rows, so estimate parameters with the unrestricted model.
!> Exact restrictions give zero-variance observations, which the univariate
!> treatment (FILTER_UNIVARIATE) handles most robustly.
module statespace_restrict
  use statespace_kinds, only: dp, SS_OK, SS_ERR_DIM
  use statespace_rep, only: ssm_rep_t, tidx
  implicit none
  private

  public :: add_state_restrictions

contains

  subroutine add_state_restrictions(rep, Rmat, rval, rrep, info)
    type(ssm_rep_t), intent(in) :: rep
    real(dp), intent(in) :: Rmat(:, :, :)     !< (q, m, 1|n) R*_t
    real(dp), intent(in) :: rval(:, :)        !< (q, n) r*_t, NaN where inactive
    type(ssm_rep_t), intent(out) :: rrep
    integer, intent(out) :: info
    integer :: p, q, m, n, nz, nh, t, iz, ir

    call rep%validate(info)
    if (info /= SS_OK) return
    p = rep%k_endog; m = rep%k_states; n = rep%nobs; q = size(Rmat, 1)
    info = SS_ERR_DIM
    if (size(Rmat, 2) /= m .or. .not. (size(Rmat, 3) == 1 .or. size(Rmat, 3) == n)) return
    if (any(shape(rval) /= [q, n])) return
    info = SS_OK

    rrep = rep
    rrep%k_endog = p + q
    nz = max(size(rep%Z, 3), size(Rmat, 3))
    nh = size(rep%H, 3)
    deallocate (rrep%y, rrep%Z, rrep%H, rrep%d)
    allocate (rrep%y(p + q, n), rrep%Z(p + q, m, nz), rrep%H(p + q, p + q, nh), &
              rrep%d(p + q, size(rep%d, 2)), source=0.0_dp)
    rrep%y(1:p, :) = rep%y
    rrep%y(p + 1:, :) = rval
    do t = 1, nz
      iz = tidx(size(rep%Z, 3), t)
      ir = tidx(size(Rmat, 3), t)
      rrep%Z(1:p, :, t) = rep%Z(:, :, iz)
      rrep%Z(p + 1:, :, t) = Rmat(:, :, ir)
    end do
    rrep%H(1:p, 1:p, :) = rep%H
    rrep%d(1:p, :) = rep%d
  end subroutine add_state_restrictions
end module statespace_restrict
