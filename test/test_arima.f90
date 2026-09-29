!> ARMA/ARIMA (DK 3.4, 3.6.2) against SARIMAX, DK's worked initialization
!> examples (5.6.1, 5.6.2, 5.6.4) and exponential smoothing (3.5).
module test_arima
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use statespace
  use fixture_io, only: fixture_t, load_fixture
  implicit none
  private

  public :: collect_arima

contains

  subroutine collect_arima(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("arma_2_1", test_201), &
                new_unittest("arima_2_1_1", test_211), &
                new_unittest("arima_2_2_1", test_221), &
                new_unittest("seasonal_arima_missing", test_seasonal), &
                new_unittest("regression_with_arma_errors", test_regression), &
                new_unittest("fit_arima_2_1_1", test_fit), &
                new_unittest("dk_5_6_1_local_linear_trend", test_dk_561), &
                new_unittest("dk_5_6_2_arma_1_1_variance", test_dk_562), &
                new_unittest("dk_5_6_4_ar1_constant", test_dk_564), &
                new_unittest("dk_3_5_exponential_smoothing", test_ewma) &
                ]
  end subroutine collect_arima

  subroutine check_rel(error, actual, expected, tol, label)
    type(error_type), allocatable, intent(out) :: error
    real(dp), intent(in) :: actual(:), expected(:), tol
    character(len=*), intent(in) :: label
    real(dp) :: err
    character(len=32) :: buf

    err = maxval(abs(actual - expected)) / max(1.0_dp, maxval(abs(expected)))
    write (buf, '(es10.3)') err
    call check(error, err <= tol, label//": relative error "//trim(buf))
  end subroutine check_rel

  !> ARIMA model for a fixture, orders from the fixture.
  function arima_model(fx, info) result(model)
    type(fixture_t), intent(in) :: fx
    integer, intent(out) :: info
    type(structural_model_t) :: model
    type(component_holder_t) :: comps(1)
    integer, allocatable :: o(:)

    o = nint(fx%get1('order'))
    comps(1)%c = arima_t(ar=o(1), d=o(2), ma=o(3), sar=o(4), sd=o(5), sma=o(6), s=o(7))
    model = structural_model(fx%get2('y'), comps, info)
  end function arima_model

  subroutine check_arima(error, path)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), intent(in) :: path
    type(fixture_t) :: fx
    type(structural_model_t) :: model
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    integer :: info

    fx = load_fixture(path)
    model = arima_model(fx, info)
    call check(error, info, SS_OK, "model info")
    if (allocated(error)) return
    call model%smooth(fx%get1('params'), fres, sres, info)
    call check(error, info, SS_OK, "smooth info")
    if (allocated(error)) return
    call check_rel(error, fres%llf_obs, fx%get1('llf_obs'), 1.0e-9_dp, "llf_obs")
    if (allocated(error)) return
    ! Same state layout as SARIMAX (DK 3.4)
    call check_rel(error, pack(sres%alphahat, .true.), fx%get1('alphahat'), 1.0e-8_dp, &
                   "alphahat")
  end subroutine check_arima

  subroutine test_201(error)
    type(error_type), allocatable, intent(out) :: error

    call check_arima(error, "test/fixtures/arima_201.txt")
  end subroutine test_201

  subroutine test_211(error)
    type(error_type), allocatable, intent(out) :: error

    call check_arima(error, "test/fixtures/arima_211.txt")
  end subroutine test_211

  subroutine test_221(error)
    type(error_type), allocatable, intent(out) :: error

    call check_arima(error, "test/fixtures/arima_221.txt")
  end subroutine test_221

  subroutine test_seasonal(error)
    type(error_type), allocatable, intent(out) :: error

    call check_arima(error, "test/fixtures/arima_111_011_4.txt")
  end subroutine test_seasonal

  subroutine test_regression(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: comps(2)
    type(structural_model_t) :: model
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: z(:)
    integer :: info

    fx = load_fixture("test/fixtures/arima_regression.txt")
    z = fx%get1('z')
    comps(1)%c = regression_t(x=reshape([spread(1.0_dp, 1, size(z)), z], [size(z), 2]))
    comps(2)%c = arima_t(ar=1, ma=1)
    model = structural_model(fx%get2('y'), comps, info)
    call model%smooth(fx%get1('params'), fres, sres, info)
    call check(error, info, SS_OK, "smooth info")
    if (allocated(error)) return
    call check_rel(error, fres%llf_obs, fx%get1('llf_obs'), 1.0e-9_dp, "llf_obs")
    if (allocated(error)) return
    call check_rel(error, sres%alphahat(1:2, model%rep%nobs), fx%get1('beta'), &
                   1.0e-9_dp, "beta")
  end subroutine test_regression

  subroutine test_fit(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(structural_model_t) :: model
    type(fit_result_t) :: res
    type(fit_options_t) :: opts
    integer :: info

    fx = load_fixture("test/fixtures/arima_211.txt")
    model = arima_model(fx, info)
    opts%factr = 10.0_dp
    opts%pgtol = 1.0e-9_dp
    call fit(model, res, options=opts, info=info)
    call check(error, info, SS_OK, "fit info")
    if (allocated(error)) return
    call check(error, abs(res%llf - sum(fx%get1('llf_tight'))) < 1.0e-7_dp, "llf")
    if (allocated(error)) return
    call check_rel(error, res%params, fx%get1('params_tight'), 1.0e-4_dp, "params")
  end subroutine test_fit

  !> DK 5.6.1: exact initial filter for the local linear trend, reaching
  !> a_3 = (2 y_2 - y_1, y_2 - y_1)', P_star,3 first column
  !> sigma2_eps (5 + 2 q_xi + q_zeta, 3 + q_xi + q_zeta)', P_inf,3 = 0.
  subroutine test_dk_561(error)
    type(error_type), allocatable, intent(out) :: error
    type(component_holder_t) :: comps(2)
    type(structural_model_t) :: model
    type(filter_result_t) :: fres
    real(dp) :: y(1, 10), s2, qx, qz
    integer :: info, t

    y(1, :) = [(sin(1.3_dp * t) + 0.2_dp * t, t=1, 10)]
    s2 = 2.0_dp; qx = 0.3_dp; qz = 0.1_dp
    allocate (irregular_t :: comps(1)%c)
    allocate (trend_t :: comps(2)%c)
    model = structural_model(y, comps, info)
    call model%filter([s2, qx * s2, qz * s2], fres, info)
    call check(error, fres%nobs_diffuse, 2, "d = 2")
    if (allocated(error)) return
    call check_rel(error, fres%a(:, 3), [2 * y(1, 2) - y(1, 1), y(1, 2) - y(1, 1)], &
                   1.0e-12_dp, "a_3")
    if (allocated(error)) return
    call check_rel(error, fres%P(:, 1, 3), s2 * [5 + 2 * qx + qz, 3 + qx + qz], &
                   1.0e-12_dp, "P_star,3")
    if (allocated(error)) return
    call check(error, all(fres%Pinf(:, :, 3) == 0.0_dp), "P_inf,3 = 0")
  end subroutine test_dk_561

  !> DK 5.6.2: ARMA(1, 1) stationary variance
  !> Q0 = [(1 + theta^2 + 2 phi theta)/(1 - phi^2), theta; theta, theta^2].
  subroutine test_dk_562(error)
    type(error_type), allocatable, intent(out) :: error
    type(component_holder_t) :: comps(1)
    type(structural_model_t) :: model
    type(filter_result_t) :: fres
    real(dp) :: y(1, 5), phi, theta
    integer :: info

    y = 1.0_dp
    phi = 0.6_dp; theta = 0.3_dp
    comps(1)%c = arima_t(ar=1, ma=1)
    model = structural_model(y, comps, info)
    call model%filter([phi, theta, 1.0_dp], fres, info)
    call check_rel(error, pack(fres%P(:, :, 1), .true.), &
                   [(1 + theta**2 + 2 * phi * theta) / (1 - phi**2), theta, theta, &
                    theta**2], 1.0e-12_dp, "Q0")
  end subroutine test_dk_562

  !> DK 5.6.4: AR(1) with a constant, alpha = (mu, xi)'. After one step,
  !> a_2 = (y_1, 0)', P_star,2 = sigma2 / (1 - phi^2) [1, -phi; -phi, 1],
  !> P_inf,2 = 0.
  subroutine test_dk_564(error)
    type(error_type), allocatable, intent(out) :: error
    type(component_holder_t) :: comps(2)
    type(structural_model_t) :: model
    type(filter_result_t) :: fres
    real(dp) :: y(1, 6), phi, s2
    integer :: info, t

    y(1, :) = [(1.0_dp + 0.5_dp * cos(real(t, dp)), t=1, 6)]
    phi = 0.7_dp; s2 = 1.5_dp
    comps(1)%c = regression_t(x=reshape([spread(1.0_dp, 1, 6)], [6, 1]))
    comps(2)%c = arima_t(ar=1)
    model = structural_model(y, comps, info)
    call model%filter([phi, s2], fres, info)
    call check(error, fres%nobs_diffuse, 1, "d = 1")
    if (allocated(error)) return
    call check_rel(error, fres%a(:, 2), [y(1, 1), 0.0_dp], 1.0e-12_dp, "a_2")
    if (allocated(error)) return
    call check_rel(error, pack(fres%P(:, :, 2), .true.), &
                   s2 / (1 - phi**2) * [1.0_dp, -phi, -phi, 1.0_dp], 1.0e-12_dp, &
                   "P_star,2")
    if (allocated(error)) return
    call check(error, all(fres%Pinf(:, :, 2) == 0.0_dp), "P_inf,2 = 0")
  end subroutine test_dk_564

  !> DK 3.5: the steady-state Kalman filter of the local level model is the
  !> EWMA yhat_t+1 = (1 - lambda) y_t + lambda yhat_t with lambda = 1 - K.
  subroutine test_ewma(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(component_holder_t) :: comps(2)
    type(structural_model_t) :: model
    type(filter_result_t) :: fres
    real(dp) :: lambda, yhat
    integer :: info, t

    fx = load_fixture("test/fixtures/nile_llevel_known.txt")
    allocate (irregular_t :: comps(1)%c)
    allocate (level_t :: comps(2)%c)
    model = structural_model(fx%get2('y'), comps, info)
    call model%filter([15099.0_dp, 1469.1_dp], fres, info)
    lambda = 1.0_dp - fres%K(1, 1, model%rep%nobs)
    yhat = fres%a(1, 40)
    do t = 40, model%rep%nobs
      yhat = (1.0_dp - lambda) * model%rep%y(1, t) + lambda * yhat
    end do
    call check_rel(error, [yhat], [fres%a(1, model%rep%nobs + 1)], 1.0e-9_dp, &
                   "EWMA forecast = steady-state Kalman filter")
  end subroutine test_ewma
end module test_arima
