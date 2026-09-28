!> DK 8.6: dynamic Nelson-Siegel model for the yield curve (Nelson and
!> Siegel 1987; Diebold and Li 2006), a dynamic factor model (DK 3.7):
!>
!>     y_t = Lambda f_t + eps_t,        eps_t ~ N(0, sigma2 I_N),
!>     f_t+1 = Phi f_t + eta_t,         eta_t ~ N(0, Sigma_eta),
!>
!> with f_t = (level, slope, curvature) and the Nelson-Siegel loadings for
!> maturity tau_i (months): x_i1 = 1, x_i2 = (1 - z_i) / (lambda tau_i),
!> x_i3 = x_i2 - z_i, z_i = exp(-lambda tau_i). The parameters lambda,
!> sigma2, Phi and Sigma_eta are estimated by maximum likelihood. The
!> factors are initialized as diffuse. As in DK, the observation vector can
!> be collapsed to the 3 factors (DK 6.5); the example checks that this
!> gives the same log likelihood and smoothed factors.
!>
!> The model is written here by extending `ssm_model_t`, the way to define
!> a model that the built-in components don't cover.
!>
!> DK use the Diebold-Li data (17 maturities from the CRSP bond files,
!> which are licensed) and report lambda = 0.078, sigma2 = 0.014,
!> Phi = [0.994 0.029 -0.022; -0.029 0.939 0.040; 0.025 0.023 0.841],
!> Sigma_eta = [0.095 -0.014 0.044; 0.383 0.009; 0.799]. This example uses the
!> public-domain FRED constant-maturity yields for the same months (8
!> maturities, see data/README.md), so its estimates differ.
!>
!> Run from the repo root:  fpm run --example dk_8_6_yield_curve
module dynamic_nelson_siegel
  use statespace
  implicit none
  private

  public :: dns_t, dns_model

  !> params = [lambda, sigma2, Phi (column-major, 9), lower triangle of
  !> Sigma_eta by columns (6)].
  type, extends(ssm_model_t) :: dns_t
    real(dp), allocatable :: tau(:)
  contains
    procedure :: update
    procedure :: start_params
    procedure :: transform_params
    procedure :: untransform_params
    procedure :: param_names
  end type dns_t

contains

  function dns_model(y, tau) result(mod)
    real(dp), intent(in) :: y(:, :), tau(:)
    type(dns_t) :: mod

    mod%tau = tau
    mod%k_params = 17
    mod%rep = ssm_rep(y, m=3, r=3)
    mod%rep%R(:, :, 1) = eye(3)
    call mod%rep%initialize_diffuse()
    call mod%update(mod%start_params())
  end function dns_model

  !> Nelson-Siegel loadings (N x 3).
  pure function loadings(tau, lambda) result(Z)
    real(dp), intent(in) :: tau(:), lambda
    real(dp) :: Z(size(tau), 3)
    real(dp) :: zi(size(tau))

    zi = exp(-lambda * tau)
    Z(:, 1) = 1.0_dp
    Z(:, 2) = (1.0_dp - zi) / (lambda * tau)
    Z(:, 3) = Z(:, 2) - zi
  end function loadings

  subroutine update(self, params)
    class(dns_t), intent(inout) :: self
    real(dp), intent(in) :: params(:)

    self%rep%Z(:, :, 1) = loadings(self%tau, params(1))
    self%rep%H(:, :, 1) = params(2) * eye(size(self%tau))
    self%rep%T(:, :, 1) = reshape(params(3:11), [3, 3])
    self%rep%Q(:, :, 1) = lower_to_sym(params(12:17))
  end subroutine update

  function start_params(self) result(params)
    class(dns_t), intent(in) :: self
    real(dp), allocatable :: params(:)

    ! Diebold and Li's lambda = 0.0609 per month
    params = [0.0609_dp, 0.01_dp, pack(0.95_dp * eye(3), .true.), &
              0.1_dp, 0.0_dp, 0.0_dp, 0.1_dp, 0.0_dp, 0.1_dp]
  end function start_params

  !> lambda = exp(x), sigma2 = exp(2 x) (DK 7.3.2), Phi unrestricted,
  !> Sigma_eta = L L' with L lower triangular and positive diagonal exp(x).
  function transform_params(self, unconstrained) result(constrained)
    class(dns_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)
    real(dp) :: L(3, 3)

    L = lower_to_mat(unconstrained(12:17))
    L(1, 1) = exp(L(1, 1)); L(2, 2) = exp(L(2, 2)); L(3, 3) = exp(L(3, 3))
    constrained = [exp(unconstrained(1)), constrain_positive(unconstrained(2)), &
                   unconstrained(3:11), sym_to_lower(matmul(L, transpose(L)))]
  end function transform_params

  function untransform_params(self, constrained) result(unconstrained)
    class(dns_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)
    real(dp) :: L(3, 3)

    L = cholesky3(lower_to_sym(constrained(12:17)))
    L(1, 1) = log(L(1, 1)); L(2, 2) = log(L(2, 2)); L(3, 3) = log(L(3, 3))
    unconstrained = [log(constrained(1)), unconstrain_positive(constrained(2)), &
                     constrained(3:11), sym_to_lower(L)]
  end function untransform_params

  function param_names(self) result(names)
    class(dns_t), intent(in) :: self
    character(len=32), allocatable :: names(:)
    integer :: i, j

    names = [character(len=32) :: "lambda", "sigma2"]
    do j = 1, 3
      do i = 1, 3
        names = [character(len=32) :: names, "phi."//achar(48 + i)//achar(48 + j)]
      end do
    end do
    do j = 1, 3
      do i = j, 3
        names = [character(len=32) :: names, "sigma_eta."//achar(48 + i)//achar(48 + j)]
      end do
    end do
  end function param_names

  !> The 3 x 3 lower triangular matrix with lower triangle v (by columns).
  pure function lower_to_mat(v) result(L)
    real(dp), intent(in) :: v(6)
    real(dp) :: L(3, 3)

    L = 0.0_dp
    L(1:3, 1) = v(1:3)
    L(2:3, 2) = v(4:5)
    L(3, 3) = v(6)
  end function lower_to_mat

  pure function lower_to_sym(v) result(S)
    real(dp), intent(in) :: v(6)
    real(dp) :: S(3, 3)

    S = lower_to_mat(v)
    S = S + transpose(S)
    S(1, 1) = v(1); S(2, 2) = v(4); S(3, 3) = v(6)
  end function lower_to_sym

  pure function sym_to_lower(S) result(v)
    real(dp), intent(in) :: S(3, 3)
    real(dp) :: v(6)

    v = [S(1:3, 1), S(2:3, 2), S(3, 3)]
  end function sym_to_lower

  pure function cholesky3(S) result(L)
    real(dp), intent(in) :: S(3, 3)
    real(dp) :: L(3, 3)
    integer :: i, j

    L = 0.0_dp
    do j = 1, 3
      L(j, j) = sqrt(S(j, j) - sum(L(j, 1:j - 1)**2))
      do i = j + 1, 3
        L(i, j) = (S(i, j) - sum(L(i, 1:j - 1) * L(j, 1:j - 1))) / L(j, j)
      end do
    end do
  end function cholesky3
end module dynamic_nelson_siegel

program dk_8_6_yield_curve
  use statespace
  use csv_io, only: read_csv, column
  use dynamic_nelson_siegel, only: dns_t, dns_model
  implicit none

  real(dp), parameter :: tau(8) = [3, 6, 12, 24, 36, 60, 84, 120]
  character(len=3), parameter :: col(8) = ["3  ", "6  ", "12 ", "24 ", "36 ", "60 ", &
                                           "84 ", "120"]
  character(len=32), allocatable :: names(:)
  real(dp), allocatable :: data(:, :), y(:, :), adjust(:), Phi(:, :)
  type(dns_t) :: mod
  type(fit_result_t) :: res
  type(fit_options_t) :: opts
  type(filter_result_t) :: fres, cfres
  type(smoother_result_t) :: sres, csres
  type(ssm_rep_t) :: crep
  integer :: info, i, t, n

  call read_csv("data/us_yields.csv", names, data)
  n = size(data, 1)
  allocate (y(size(tau), n))
  do i = 1, size(tau)
    y(i, :) = column(names, data, "m"//trim(col(i)))
  end do

  mod = dns_model(y, tau)
  opts%factr = 1.0e3_dp
  opts%pgtol = 1.0e-7_dp
  opts%maxiter = 2000
  call fit(mod, res, options=opts, info=info)
  if (info /= SS_OK .and. info /= SS_ERR_NOT_PD) error stop "fit"
  print '(a)', trim(res%message)
  print '(a, f10.4, a, f10.5)', "lambda: ", res%params(1), "   sigma2: ", res%params(2)
  Phi = reshape(res%params(3:11), [3, 3])
  print '(a)', "Phi:"
  print '(3f9.3)', (Phi(i, :), i=1, 3)
  print '(a)', "Sigma_eta (lower triangle):"
  print '(f9.3)', res%params(12)
  print '(2f9.3)', res%params(13), res%params(15)
  print '(3f9.3)', res%params(14), res%params(16), res%params(17)
  print '(a, f12.3)', "log likelihood: ", res%llf

  ! The collapsed model (DK 6.5) gives the same likelihood and factors
  call mod%smooth(res%params, fres, sres, info)
  allocate (adjust(n))
  call collapse_observations(mod%rep, crep, adjust, info)
  if (info /= SS_OK) error stop "collapse"
  call kalman_filter(crep, cfres, info)
  call state_smoother(crep, cfres, csres, info)
  print '(a, es10.2)', "collapsed log likelihood - full: ", &
      cfres%llf + sum(adjust) - fres%llf
  print '(a, es10.2)', "largest smoothed factor difference: ", &
      maxval(abs(csres%alphahat - sres%alphahat))

  ! Smoothed factors against the data proxies of DK Fig. 8.12
  print '(/, a6, 6a10)', "month", "level", "10y", "slope", "3m-10y", "curv", "proxy"
  do t = 1, n, 24
    print '(i6, 6f10.3)', t, sres%alphahat(1, t), y(8, t), sres%alphahat(2, t), &
        y(1, t) - y(8, t), sres%alphahat(3, t), 2 * y(4, t) - y(1, t) - y(8, t)
  end do
end program dk_8_6_yield_curve
