!> EM algorithm for the variance matrices H and Q (DK 7.3.4).
!>
!> With the disturbances as complete data, the M-step for time-invariant,
!> otherwise unrestricted H and Q is closed form:
!>
!>     H_new = 1/n sum_t=1..n (epshat_t epshat_t' + Var(eps_t | Y))
!>     Q_new = 1/(n-1) sum_t=1..n-1 (etahat_t etahat_t' + Var(eta_t | Y))
!>
!> (eta_n does not affect the data). Missing elements of y_t enter through
!> their conditional moments, which is why the smoother reports those rather
!> than 0 and H. Each iteration does not decrease the log likelihood, and a
!> fixed point is a stationary point of it. EM needs an initialization that
!> does not depend on H or Q (not stationary blocks) and no log likelihood
!> burn-in.
module statespace_em
  use statespace_kinds, only: dp, SS_OK, SS_ERR_UNSUPPORTED
  use statespace_linalg, only: symmetrize
  use statespace_rep, only: ssm_rep_t, INIT_STATIONARY
  use statespace_filter, only: filter_result_t, kalman_filter
  use statespace_smoother, only: smoother_result_t, state_smoother
  implicit none
  private

  public :: em_step, em_variances

contains

  !> One EM iteration: E-step at the current H and Q of `rep`, M-step
  !> replacing them. `llf` is the log likelihood before the update. With
  !> `diagonal_H` / `diagonal_Q`, the M-step keeps only the diagonal
  !> (restricted maximization for diagonal matrices).
  subroutine em_step(rep, llf, info, diagonal_H, diagonal_Q)
    type(ssm_rep_t), intent(inout) :: rep
    real(dp), intent(out) :: llf
    integer, intent(out) :: info
    logical, intent(in), optional :: diagonal_H, diagonal_Q
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: Hn(:, :), Qn(:, :)
    integer :: t, n, j

    llf = 0.0_dp
    info = SS_ERR_UNSUPPORTED
    if (size(rep%H, 3) > 1 .or. size(rep%Q, 3) > 1 .or. rep%loglikelihood_burn > 0) then
      return
    end if
    if (allocated(rep%blk_kind)) then
      if (any(rep%blk_kind == INIT_STATIONARY)) return
    end if
    n = rep%nobs
    if (n < 2) return

    call kalman_filter(rep, fres, info)
    if (info /= SS_OK) return
    call state_smoother(rep, fres, sres, info)
    if (info /= SS_OK) return
    llf = fres%llf

    allocate (Hn(rep%k_endog, rep%k_endog), Qn(rep%k_posdef, rep%k_posdef), &
              source=0.0_dp)
    do t = 1, n
      do j = 1, rep%k_endog
        Hn(:, j) = Hn(:, j) + sres%epshat(:, t) * sres%epshat(j, t)
      end do
      Hn = Hn + sres%epsvar(:, :, t)
      if (t < n) then
        do j = 1, rep%k_posdef
          Qn(:, j) = Qn(:, j) + sres%etahat(:, t) * sres%etahat(j, t)
        end do
        Qn = Qn + sres%etavar(:, :, t)
      end if
    end do
    Hn = Hn / n
    Qn = Qn / (n - 1)
    call symmetrize(Hn)
    call symmetrize(Qn)
    if (present(diagonal_H)) then
      if (diagonal_H) Hn = diagonal_part(Hn)
    end if
    if (present(diagonal_Q)) then
      if (diagonal_Q) Qn = diagonal_part(Qn)
    end if
    rep%H(:, :, 1) = Hn
    rep%Q(:, :, 1) = Qn
    info = SS_OK
  end subroutine em_step

  !> Iterate `em_step` until the log likelihood improves by less than `tol`
  !> or `maxiter` iterations. `llf` is the log likelihood at the start of the
  !> last iteration (rep holds the matrices after it); `llf_path` (optional)
  !> receives it for every iteration.
  subroutine em_variances(rep, maxiter, tol, llf, niter, info, diagonal_H, diagonal_Q, &
                          llf_path)
    type(ssm_rep_t), intent(inout) :: rep
    integer, intent(in) :: maxiter
    real(dp), intent(in) :: tol
    real(dp), intent(out) :: llf
    integer, intent(out) :: niter, info
    logical, intent(in), optional :: diagonal_H, diagonal_Q
    real(dp), allocatable, intent(out), optional :: llf_path(:)
    real(dp), allocatable :: path(:)
    real(dp) :: llf_prev
    integer :: k

    allocate (path(maxiter))
    llf_prev = -huge(1.0_dp)
    niter = 0
    do k = 1, maxiter
      call em_step(rep, llf, info, diagonal_H, diagonal_Q)
      if (info /= SS_OK) return
      niter = k
      path(k) = llf
      if (abs(llf - llf_prev) < tol) exit
      llf_prev = llf
    end do
    if (present(llf_path)) llf_path = path(1:niter)
  end subroutine em_variances

  pure function diagonal_part(A) result(D)
    real(dp), intent(in), contiguous :: A(:, :)
    real(dp) :: D(size(A, 1), size(A, 2))
    integer :: i

    D = 0.0_dp
    do i = 1, min(size(A, 1), size(A, 2))
      D(i, i) = A(i, i)
    end do
  end function diagonal_part
end module statespace_em
