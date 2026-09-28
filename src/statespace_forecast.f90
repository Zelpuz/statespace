!> Out-of-sample forecasts (DK 4.11).
!>
!> Forecasting is filtering with missing observations. `forecast` runs the
!> prediction recursion directly from the end of a filter run, which must be
!> past the diffuse period. For models with
!> time-varying matrices, append NaN observations to y instead, supply the
!> matrices for the forecast horizon, and read `yhat` and `F` from the filter.
module statespace_forecast
  use statespace_kinds, only: dp, SS_OK, SS_ERR_DIM, SS_ERR_UNSUPPORTED
  use statespace_linalg, only: gemm, gemv, symmetrize
  use statespace_rep, only: ssm_rep_t
  use statespace_filter, only: filter_result_t
  implicit none
  private

  public :: forecast_result_t, forecast

  type :: forecast_result_t
    integer :: horizon = 0
    real(dp), allocatable :: mean(:, :)         !< (p, h) E[y_n+j | Y_n]
    real(dp), allocatable :: cov(:, :, :)       !< (p, p, h) Var[y_n+j | Y_n]
    real(dp), allocatable :: state(:, :)        !< (m, h) E[alpha_n+j | Y_n]
    real(dp), allocatable :: state_cov(:, :, :) !< (m, m, h) Var[alpha_n+j | Y_n]
  end type forecast_result_t

contains

  !> Forecast j = 1..h steps past the end of the sample. Requires
  !> time-invariant system matrices.
  subroutine forecast(rep, fres, h, fc, info)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    integer, intent(in) :: h
    type(forecast_result_t), intent(out) :: fc
    integer, intent(out) :: info
    real(dp), allocatable :: PZt(:, :), RQ(:, :), RQR(:, :), TP(:, :)
    integer :: p, m, j

    info = SS_ERR_UNSUPPORTED
    if (size(rep%Z, 3) > 1 .or. size(rep%H, 3) > 1 .or. size(rep%T, 3) > 1 .or. &
        size(rep%R, 3) > 1 .or. size(rep%Q, 3) > 1 .or. size(rep%c, 2) > 1 .or. &
        size(rep%d, 2) > 1) return
    info = SS_ERR_DIM
    if (fres%nobs /= rep%nobs .or. fres%k_states /= rep%k_states .or. h < 0) return
    ! Still diffuse at the end of the sample: forecasts have infinite variance.
    info = SS_ERR_UNSUPPORTED
    if (any(fres%Pinf(:, :, rep%nobs + 1) /= 0.0_dp)) return
    info = SS_OK

    p = rep%k_endog; m = rep%k_states
    fc%horizon = h
    allocate (fc%mean(p, h), fc%cov(p, p, h), fc%state(m, h), fc%state_cov(m, m, h))
    if (h == 0) return
    allocate (PZt(m, p), RQ(m, rep%k_posdef), RQR(m, m), TP(m, m))
    call gemm('N', 'N', 1.0_dp, rep%R(:, :, 1), rep%Q(:, :, 1), 0.0_dp, RQ)
    call gemm('N', 'T', 1.0_dp, RQ, rep%R(:, :, 1), 0.0_dp, RQR)
    fc%state(:, 1) = fres%a(:, rep%nobs + 1)
    fc%state_cov(:, :, 1) = fres%P(:, :, rep%nobs + 1)
    do j = 1, h
      if (j > 1) then
        ! a_j = c + T a_j-1,  P_j = T P_j-1 T' + R Q R'
        fc%state(:, j) = rep%c(:, 1)
        call gemv('N', 1.0_dp, rep%T(:, :, 1), fc%state(:, j - 1), 1.0_dp, &
                  fc%state(:, j))
        call gemm('N', 'N', 1.0_dp, rep%T(:, :, 1), fc%state_cov(:, :, j - 1), 0.0_dp, &
                  TP)
        fc%state_cov(:, :, j) = RQR
        call gemm('N', 'T', 1.0_dp, TP, rep%T(:, :, 1), 1.0_dp, fc%state_cov(:, :, j))
        call symmetrize(fc%state_cov(:, :, j))
      end if
      ! y: mean d + Z a, variance Z P Z' + H
      fc%mean(:, j) = rep%d(:, 1)
      call gemv('N', 1.0_dp, rep%Z(:, :, 1), fc%state(:, j), 1.0_dp, fc%mean(:, j))
      call gemm('N', 'T', 1.0_dp, fc%state_cov(:, :, j), rep%Z(:, :, 1), 0.0_dp, PZt)
      fc%cov(:, :, j) = rep%H(:, :, 1)
      call gemm('N', 'N', 1.0_dp, rep%Z(:, :, 1), PZt, 1.0_dp, fc%cov(:, :, j))
      call symmetrize(fc%cov(:, :, j))
    end do
  end subroutine forecast
end module statespace_forecast
