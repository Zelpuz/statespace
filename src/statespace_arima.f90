!> ARMA and ARIMA components (DK 3.4, 3.6.2, 5.6.2-5.6.4).
!>
!> The stationary part y*_t = Delta^d Delta_s^D y_t follows the ARMA model
!> phi(L) Phi(L^s) y*_t = theta(L) Theta(L^s) zeta_t, put in DK's form
!> (3.19)-(3.20): with r = max(p + sP, q + sQ + 1),
!>
!>     T = [phi* | I; 0],  R = (1, theta*_1, ..., theta*_r-1)',  Z = e_1',
!>
!> phi* and theta* the coefficients of the multiplied polynomials. For d > 0
!> the state is extended as in DK 3.4 by the cumulative differences
!> (y_t-1, Delta y_t-1, ..., Delta^d-1 y_t-1), and for D = 1 by the s lags
!> of Delta^d y (Delta^d y_t = Delta^d y_t-s + y*_t). The differenced states
!> are diffuse and the ARMA block stationary (DK 5.6.3). A constant or
!> regressors (regression with ARMA errors, DK 3.6.2, 5.6.4) are added as
!> a `regression_t` component.
!>
!> Parameters: AR (p), seasonal AR (P), MA (q), seasonal MA (Q), then the
!> variance of zeta. AR polynomials are kept stationary and MA polynomials
!> invertible (Monahan's transform, as statsmodels).
module statespace_arima
  use statespace_kinds, only: dp
  use statespace_rep, only: INIT_DIFFUSE, INIT_STATIONARY
  use statespace_model, only: constrain_positive, unconstrain_positive, constrain_stationary, &
                              unconstrain_stationary
  use statespace_components, only: component_t, diff_variance
  implicit none
  private

  public :: arima_t

  type, extends(component_t) :: arima_t
    integer :: ar = 0, d = 0, ma = 0            !< p, d, q
    integer :: sar = 0, sd = 0, sma = 0         !< P, D, Q
    integer :: s = 0                            !< seasonal period
    integer :: series = 1                       !< the series (or signal) it drives
    logical :: enforce_stationarity = .true.
    logical :: enforce_invertibility = .true.
  contains
    procedure :: setup => arima_setup
    procedure :: fill => arima_fill
    procedure :: transform_params => arima_transform
    procedure :: untransform_params => arima_untransform
    procedure :: start_params => arima_start
    procedure :: param_names => arima_names
    procedure :: init_blocks => arima_init_blocks
  end type arima_t

contains

  pure integer function n_diff(self)
    class(arima_t), intent(in) :: self

    n_diff = self%d + self%s * self%sd
  end function n_diff

  pure integer function n_arma(self)
    class(arima_t), intent(in) :: self

    n_arma = max(self%ar + self%s * self%sar, self%ma + self%s * self%sma + 1)
  end function n_arma

  subroutine arima_setup(self, y)
    class(arima_t), intent(inout) :: self
    real(dp), intent(in) :: y(:, :)

    self%p = size(y, 1); self%n = size(y, 2)
    if (self%sd > 1) error stop "arima_t: seasonal differencing of order > 1 is not supported"
    self%m = n_diff(self) + n_arma(self)
    self%r = 1
    self%k = self%ar + self%sar + self%ma + self%sma + 1
  end subroutine arima_setup

  !> Coefficients c(1:) of (1 + sign sum_i a_i L^i)(1 + sign sum_j b_j L^(s j)),
  !> without the leading 1.
  pure function polymul(a, b, s, sign) result(c)
    real(dp), intent(in) :: a(:), b(:), sign
    integer, intent(in) :: s
    real(dp), allocatable :: c(:)
    real(dp), allocatable :: pa(:), pb(:), pc(:)
    integer :: i, j

    allocate (pa(0:size(a)), pb(0:s * size(b)), source=0.0_dp)
    pa(0) = 1.0_dp
    pa(1:) = sign * a
    pb(0) = 1.0_dp
    do j = 1, size(b)
      pb(s * j) = sign * b(j)
    end do
    allocate (pc(0:size(a) + s * size(b)), source=0.0_dp)
    do i = 0, size(a)
      do j = 0, s * size(b)
        pc(i + j) = pc(i + j) + pa(i) * pb(j)
      end do
    end do
    c = sign * pc(1:)
  end function polymul

  subroutine arima_fill(self, params, Z, H, T, R, Q)
    class(arima_t), intent(in) :: self
    real(dp), intent(in) :: params(:)
    real(dp), intent(inout) :: Z(:, :, :), H(:, :, :), T(:, :, :), R(:, :, :), Q(:, :, :)
    real(dp), allocatable :: phi(:), theta(:)
    integer :: nd, na, i, j, a0, ia, isa, ima, isma

    nd = n_diff(self); na = n_arma(self)
    ia = 0; isa = ia + self%ar; ima = isa + self%sar; isma = ima + self%ma
    ! phi(L) Phi(L^s) = 1 - sum phi*_k L^k,  theta(L) Theta(L^s) = 1 + sum theta*_k L^k
    phi = polymul(params(ia + 1:ia + self%ar), params(isa + 1:isa + self%sar), max(self%s, 1), -1.0_dp)
    theta = polymul(params(ima + 1:ima + self%ma), params(isma + 1:isma + self%sma), &
                    max(self%s, 1), 1.0_dp)

    a0 = nd       ! offset of the ARMA block
    ! ARMA block (DK 3.20)
    do i = 1, size(phi)
      T(a0 + i, a0 + 1, 1) = phi(i)
    end do
    do i = 1, na - 1
      T(a0 + i, a0 + i + 1, :) = 1.0_dp
    end do
    R(a0 + 1, 1, :) = 1.0_dp
    do i = 1, size(theta)
      R(a0 + 1 + i, 1, :) = theta(i)
    end do
    Q(1, 1, 1) = params(self%k)

    ! y_t = sum of the cumulative differences + Delta^d y_t, and
    ! Delta^d y_t = (last seasonal lag, if D = 1) + y*_t.
    do j = 1, self%d
      Z(self%series, j, :) = 1.0_dp
    end do
    if (self%sd == 1) Z(self%series, self%d + self%s, :) = 1.0_dp
    Z(self%series, a0 + 1, :) = 1.0_dp
    ! Cumulative difference rows: Delta^(i-1) y_t = sum_(j >= i) Delta^(j-1) y_t-1 + Delta^d y_t
    do i = 1, self%d
      do j = i, self%d
        T(i, j, :) = 1.0_dp
      end do
      if (self%sd == 1) T(i, self%d + self%s, :) = 1.0_dp
      T(i, a0 + 1, :) = 1.0_dp
    end do
    ! Seasonal lags of Delta^d y: new first = last lag + y*_t, then shift
    if (self%sd == 1) then
      T(self%d + 1, self%d + self%s, :) = 1.0_dp
      T(self%d + 1, a0 + 1, :) = 1.0_dp
      do i = 2, self%s
        T(self%d + i, self%d + i - 1, :) = 1.0_dp
      end do
    end if
  end subroutine arima_fill

  function arima_transform(self, unconstrained) result(constrained)
    class(arima_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)

    constrained = unconstrained
    call apply(self%ar, 0, .true., self%enforce_stationarity)
    call apply(self%sar, self%ar, .true., self%enforce_stationarity)
    call apply(self%ma, self%ar + self%sar, .false., self%enforce_invertibility)
    call apply(self%sma, self%ar + self%sar + self%ma, .false., self%enforce_invertibility)
    constrained(self%k) = constrain_positive(unconstrained(self%k))

  contains

    !> AR: phi = constrain_stationary(x); MA: theta = -constrain_stationary(x).
    subroutine apply(n, off, is_ar, enforce)
      integer, intent(in) :: n, off
      logical, intent(in) :: is_ar, enforce

      if (n == 0 .or. .not. enforce) return
      constrained(off + 1:off + n) = constrain_stationary(unconstrained(off + 1:off + n))
      if (.not. is_ar) constrained(off + 1:off + n) = -constrained(off + 1:off + n)
    end subroutine apply
  end function arima_transform

  function arima_untransform(self, constrained) result(unconstrained)
    class(arima_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)

    unconstrained = constrained
    call apply(self%ar, 0, .true., self%enforce_stationarity)
    call apply(self%sar, self%ar, .true., self%enforce_stationarity)
    call apply(self%ma, self%ar + self%sar, .false., self%enforce_invertibility)
    call apply(self%sma, self%ar + self%sar + self%ma, .false., self%enforce_invertibility)
    unconstrained(self%k) = unconstrain_positive(constrained(self%k))

  contains

    subroutine apply(n, off, is_ar, enforce)
      integer, intent(in) :: n, off
      logical, intent(in) :: is_ar, enforce

      if (n == 0 .or. .not. enforce) return
      if (is_ar) then
        unconstrained(off + 1:off + n) = unconstrain_stationary(constrained(off + 1:off + n))
      else
        unconstrained(off + 1:off + n) = unconstrain_stationary(-constrained(off + 1:off + n))
      end if
    end subroutine apply
  end function arima_untransform

  function arima_start(self, y) result(params)
    class(arima_t), intent(in) :: self
    real(dp), intent(in) :: y(:, :)
    real(dp), allocatable :: params(:)

    allocate (params(self%k), source=0.0_dp)
    params(self%k) = diff_variance(y)
  end function arima_start

  function arima_names(self) result(names)
    class(arima_t), intent(in) :: self
    character(len=32), allocatable :: names(:)
    integer :: i, k

    allocate (names(self%k))
    k = 0
    do i = 1, self%ar
      k = k + 1; write (names(k), '("ar.L", i0)') i
    end do
    do i = 1, self%sar
      k = k + 1; write (names(k), '("ar.S.L", i0)') i * self%s
    end do
    do i = 1, self%ma
      k = k + 1; write (names(k), '("ma.L", i0)') i
    end do
    do i = 1, self%sma
      k = k + 1; write (names(k), '("ma.S.L", i0)') i * self%s
    end do
    names(self%k) = "sigma2"
  end function arima_names

  !> Differenced states diffuse, ARMA block stationary (DK 5.6.3).
  function arima_init_blocks(self) result(blocks)
    class(arima_t), intent(in) :: self
    integer, allocatable :: blocks(:, :)
    integer :: nd

    nd = n_diff(self)
    if (nd == 0) then
      blocks = reshape([1, self%m, INIT_STATIONARY], [1, 3])
    else
      blocks = reshape([1, nd + 1, nd, self%m, INIT_DIFFUSE, INIT_STATIONARY], [2, 3])
    end if
  end function arima_init_blocks
end module statespace_arima
