!> Structural time series components (DK 3.2-3.3).
!>
!> Each component can be used for p series at once, as a seemingly unrelated
!> (SUTSE) model (DK 3.3): the series have their own states and the
!> disturbances are correlated across series through a p x p variance
!> matrix, which is either diagonal (COV_DIAGONAL), full (COV_FULL, via its
!> Cholesky factor) or fixed at zero (COV_NONE, a deterministic component).
!> For p = 1 the covariance types reduce to a single variance or none.
!>
!> State order within a component is series by series, DK's order within
!> each series: level-slope for the trend, (gamma_t, ..., gamma_t-s+2) for
!> the dummy seasonal, (gamma_1, gamma_1*, gamma_2, ...) for the
!> trigonometric seasonal, (c, c*) for the cycle.
module statespace_structural
  use statespace_kinds, only: dp
  use statespace_rep, only: INIT_DIFFUSE, INIT_STATIONARY
  use statespace_model, only: constrain_positive, unconstrain_positive
  use statespace_components, only: component_t, constrain_interval, unconstrain_interval, &
                                   diff_variance
  implicit none
  private

  public :: irregular_t, level_t, trend_t, seasonal_t, cycle_t, regression_t
  public :: continuous_level_t, continuous_trend_t
  public :: step_intervention, pulse_intervention, slope_intervention

  integer, parameter, public :: COV_NONE = 0, COV_DIAGONAL = 1, COV_FULL = 2
  integer, parameter, public :: SEASONAL_DUMMY = 1, SEASONAL_TRIG = 2, SEASONAL_HS = 3

  real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

  !> The irregular eps_t ~ N(0, Sigma_eps): contributes H, no states.
  !> With `weights` (n), H_t = w_t Sigma_eps: unequal observation variances,
  !> as for irregularly spaced data (DK 3.8).
  type, extends(component_t) :: irregular_t
    integer :: cov = COV_DIAGONAL
    real(dp), allocatable :: weights(:)
  contains
    procedure :: setup => irregular_setup
    procedure :: fill => irregular_fill
    procedure :: transform_params => cov_transform_irr
    procedure :: untransform_params => cov_untransform_irr
    procedure :: start_params => irregular_start
    procedure :: observation_level => irregular_observation_level
  end type irregular_t

  !> Local level (random walk) mu_t+1 = mu_t + xi_t (DK 2, 3.3).
  type, extends(component_t) :: level_t
    integer :: cov = COV_DIAGONAL
  contains
    procedure :: setup => level_setup
    procedure :: fill => level_fill
    procedure :: transform_params => cov_transform_level
    procedure :: untransform_params => cov_untransform_level
    procedure :: start_params => level_start
  end type level_t

  !> Local linear trend (DK 3.2): mu_t+1 = mu_t + nu_t + xi_t,
  !> nu_t+1 = nu_t + zeta_t. cov_level = COV_NONE gives the smooth trend
  !> (integrated random walk); both COV_NONE a deterministic linear trend.
  type, extends(component_t) :: trend_t
    integer :: cov_level = COV_DIAGONAL
    integer :: cov_slope = COV_DIAGONAL
  contains
    procedure :: setup => trend_setup
    procedure :: fill => trend_fill
    procedure :: transform_params => cov_transform_trend
    procedure :: untransform_params => cov_untransform_trend
    procedure :: start_params => trend_start
  end type trend_t

  !> Seasonal with period s (DK 3.2.2): SEASONAL_DUMMY (3.3),
  !> SEASONAL_TRIG (3.7-3.8), SEASONAL_HS (Harrison-Stevens, 3.4).
  type, extends(component_t) :: seasonal_t
    integer :: period = 4
    integer :: form = SEASONAL_DUMMY
    integer :: cov = COV_DIAGONAL
  contains
    procedure :: setup => seasonal_setup
    procedure :: fill => seasonal_fill
    procedure :: transform_params => cov_transform_seas
    procedure :: untransform_params => cov_untransform_seas
    procedure :: start_params => seasonal_start
  end type seasonal_t

  !> Cycle (DK 3.2.4, 3.13): (c, c*)_t+1 = rho C(lambda) (c, c*)_t + w_t,
  !> Var(w_t) = Sigma (x) I_2. With `damped` (rho in (0, 1)) the cycle is
  !> stationary and so initialized; lambda is bounded by the period bounds
  !> (in observations), e.g. 1.5 to 12 years of data.
  type, extends(component_t) :: cycle_t
    integer :: cov = COV_DIAGONAL
    logical :: damped = .true.
    real(dp) :: period_min = 2.0_dp
    real(dp) :: period_max = 0.0_dp      !< 0: the number of observations
  contains
    procedure :: setup => cycle_setup
    procedure :: fill => cycle_fill
    procedure :: transform_params => cycle_transform
    procedure :: untransform_params => cycle_untransform
    procedure :: start_params => cycle_start
    procedure :: init_blocks => cycle_init_blocks
  end type cycle_t

  !> Continuous-time local level (DK 3.8.1) observed at times t_1 < ... < t_n:
  !> alpha(t) = alpha(0) + sigma w(t), so between observations
  !> alpha_i+1 = alpha_i + eta_i with Var(eta_i) = sigma^2 delta_i,
  !> delta_i = t_i+1 - t_i. Univariate. Parameter: sigma^2 per unit time.
  type, extends(component_t) :: continuous_level_t
    real(dp), allocatable :: times(:)
  contains
    procedure :: setup => clevel_setup
    procedure :: fill => clevel_fill
  end type continuous_level_t

  !> Continuous-time smooth trend (DK 3.8.2, 3.9.2): d nu = sigma dw,
  !> d mu = nu dt, giving (3.41)-(3.42)
  !>     T_i = [1 delta_i; 0 1],
  !>     Q_i = sigma^2 delta_i [delta_i^2/3 delta_i/2; delta_i/2 1].
  !> With diffuse (mu(0), nu(0)) and an irregular of variance sigma_eps^2 the
  !> smoothed mu is the cubic smoothing spline with lambda = sigma_eps^2 /
  !> sigma^2 (Wahba 1978). Univariate. Parameter: sigma^2 per unit time.
  type, extends(component_t) :: continuous_trend_t
    real(dp), allocatable :: times(:)
  contains
    procedure :: setup => ctrend_setup
    procedure :: fill => ctrend_fill
  end type continuous_trend_t

  !> Regression effects for series `series` (DK 3.2.5, 3.6): Z_t = x_t,
  !> coefficients as diffuse states, fixed (T = I, no disturbance) or random
  !> walks where `random_walk` is true (3.15), with variance parameters.
  !> With signal loadings, set `at_observations` for effects on an observed
  !> series (e.g. the intercepts a* of the common levels model, DK 3.3.2).
  type, extends(component_t) :: regression_t
    real(dp), allocatable :: x(:, :)        !< (n, k_x) regressors
    logical, allocatable :: random_walk(:)  !< (k_x), default all false
    integer :: series = 1
  contains
    procedure :: setup => regression_setup
    procedure :: fill => regression_fill
  end type regression_t

contains

  ! -------------------------------------------------------------- helpers

  !> Number of parameters of a p x p covariance of the given type.
  pure integer function cov_nparams(cov, p)
    integer, intent(in) :: cov, p

    select case (cov)
    case (COV_DIAGONAL)
      cov_nparams = p
    case (COV_FULL)
      cov_nparams = p * (p + 1) / 2
    case default
      cov_nparams = 0
    end select
  end function cov_nparams

  !> Start values for a covariance: v on the diagonal, zero elsewhere.
  pure function cov_start(cov, p, v) result(c)
    integer, intent(in) :: cov, p
    real(dp), intent(in) :: v
    real(dp), allocatable :: c(:)
    integer :: i, j, k

    allocate (c(cov_nparams(cov, p)))
    select case (cov)
    case (COV_DIAGONAL)
      c = v
    case (COV_FULL)
      k = 0
      do j = 1, p
        do i = j, p
          k = k + 1
          c(k) = merge(v, 0.0_dp, i == j)
        end do
      end do
    end select
  end function cov_start

  !> Covariance matrix from constrained parameters: variances (diagonal) or
  !> the lower triangle of Sigma column by column (full).
  pure function cov_matrix(cov, p, params) result(S)
    integer, intent(in) :: cov, p
    real(dp), intent(in) :: params(:)
    real(dp) :: S(p, p)
    integer :: i, j, k

    S = 0.0_dp
    select case (cov)
    case (COV_DIAGONAL)
      do i = 1, p
        S(i, i) = params(i)
      end do
    case (COV_FULL)
      k = 0
      do j = 1, p
        do i = j, p
          k = k + 1
          S(i, j) = params(k)
          S(j, i) = params(k)
        end do
      end do
    end select
  end function cov_matrix

  !> Unconstrained -> constrained for a covariance: log-sd for variances
  !> (DK 7.3.2); for a full matrix, Sigma = L L' with L lower triangular,
  !> diagonal exp(psi), and the lower triangle of Sigma returned.
  pure function cov_constrain(cov, p, x) result(c)
    integer, intent(in) :: cov, p
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: c(:)
    real(dp) :: L(p, p), S(p, p)
    integer :: i, j, k

    select case (cov)
    case (COV_DIAGONAL)
      c = constrain_positive(x)
    case (COV_FULL)
      L = 0.0_dp
      k = 0
      do j = 1, p
        do i = j, p
          k = k + 1
          if (i == j) then
            L(i, j) = exp(x(k))
          else
            L(i, j) = x(k)
          end if
        end do
      end do
      S = matmul(L, transpose(L))
      allocate (c(size(x)))
      k = 0
      do j = 1, p
        do i = j, p
          k = k + 1
          c(k) = S(i, j)
        end do
      end do
    case default
      allocate (c(0))
    end select
  end function cov_constrain

  pure function cov_unconstrain(cov, p, c) result(x)
    integer, intent(in) :: cov, p
    real(dp), intent(in) :: c(:)
    real(dp), allocatable :: x(:)
    real(dp) :: S(p, p), L(p, p)
    integer :: i, j, k

    select case (cov)
    case (COV_DIAGONAL)
      x = unconstrain_positive(c)
    case (COV_FULL)
      S = cov_matrix(cov, p, c)
      L = 0.0_dp
      do j = 1, p
        L(j, j) = sqrt(S(j, j) - sum(L(j, 1:j - 1)**2))
        do i = j + 1, p
          L(i, j) = (S(i, j) - sum(L(i, 1:j - 1) * L(j, 1:j - 1))) / L(j, j)
        end do
      end do
      allocate (x(size(c)))
      k = 0
      do j = 1, p
        do i = j, p
          k = k + 1
          if (i == j) then
            x(k) = log(L(i, j))
          else
            x(k) = L(i, j)
          end if
        end do
      end do
    case default
      allocate (x(0))
    end select
  end function cov_unconstrain

  ! ------------------------------------------------------------ irregular

  subroutine irregular_setup(self, y)
    class(irregular_t), intent(inout) :: self
    real(dp), intent(in) :: y(:, :)

    self%p = size(y, 1); self%n = size(y, 2)
    self%m = 0; self%r = 0
    self%k = cov_nparams(self%cov, self%p)
    self%tv_H = allocated(self%weights)
  end subroutine irregular_setup

  subroutine irregular_fill(self, params, Z, H, T, R, Q)
    class(irregular_t), intent(in) :: self
    real(dp), intent(in) :: params(:)
    real(dp), intent(inout) :: Z(:, :, :), H(:, :, :), T(:, :, :), R(:, :, :), Q(:, :, :)
    real(dp) :: S(self%p, self%p)
    integer :: j

    S = cov_matrix(self%cov, self%p, params)
    do j = 1, size(H, 3)
      if (allocated(self%weights)) then
        H(:, :, j) = H(:, :, j) + self%weights(j) * S
      else
        H(:, :, j) = H(:, :, j) + S
      end if
    end do
  end subroutine irregular_fill

  !> The irregular always acts on the observations.
  logical function irregular_observation_level(self)
    class(irregular_t), intent(in) :: self

    irregular_observation_level = .true.
  end function irregular_observation_level

  function cov_transform_irr(self, unconstrained) result(constrained)
    class(irregular_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)

    constrained = cov_constrain(self%cov, self%p, unconstrained)
  end function cov_transform_irr

  function cov_untransform_irr(self, constrained) result(unconstrained)
    class(irregular_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)

    unconstrained = cov_unconstrain(self%cov, self%p, constrained)
  end function cov_untransform_irr

  ! ---------------------------------------------------------------- level

  subroutine level_setup(self, y)
    class(level_t), intent(inout) :: self
    real(dp), intent(in) :: y(:, :)

    self%p = size(y, 1); self%n = size(y, 2)
    self%m = self%p
    self%r = merge(0, self%p, self%cov == COV_NONE)
    self%k = cov_nparams(self%cov, self%p)
  end subroutine level_setup

  subroutine level_fill(self, params, Z, H, T, R, Q)
    class(level_t), intent(in) :: self
    real(dp), intent(in) :: params(:)
    real(dp), intent(inout) :: Z(:, :, :), H(:, :, :), T(:, :, :), R(:, :, :), Q(:, :, :)
    integer :: i

    do i = 1, self%p
      Z(i, i, :) = 1.0_dp
      T(i, i, :) = 1.0_dp
      if (self%r > 0) R(i, i, :) = 1.0_dp
    end do
    if (self%r > 0) Q(:, :, 1) = cov_matrix(self%cov, self%p, params)
  end subroutine level_fill

  function cov_transform_level(self, unconstrained) result(constrained)
    class(level_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)

    constrained = cov_constrain(self%cov, self%p, unconstrained)
  end function cov_transform_level

  function cov_untransform_level(self, constrained) result(unconstrained)
    class(level_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)

    unconstrained = cov_unconstrain(self%cov, self%p, constrained)
  end function cov_untransform_level

  ! ---------------------------------------------------------------- trend

  subroutine trend_setup(self, y)
    class(trend_t), intent(inout) :: self
    real(dp), intent(in) :: y(:, :)

    self%p = size(y, 1); self%n = size(y, 2)
    self%m = 2 * self%p
    self%r = merge(0, self%p, self%cov_level == COV_NONE) + merge(0, self%p, self%cov_slope == COV_NONE)
    self%k = cov_nparams(self%cov_level, self%p) + cov_nparams(self%cov_slope, self%p)
  end subroutine trend_setup

  !> States (mu_i, nu_i) for series i = 1..p; disturbances xi (if any) then zeta.
  subroutine trend_fill(self, params, Z, H, T, R, Q)
    class(trend_t), intent(in) :: self
    real(dp), intent(in) :: params(:)
    real(dp), intent(inout) :: Z(:, :, :), H(:, :, :), T(:, :, :), R(:, :, :), Q(:, :, :)
    integer :: i, p, kl, e

    p = self%p
    kl = cov_nparams(self%cov_level, p)
    do i = 1, p
      Z(i, 2 * i - 1, :) = 1.0_dp
      T(2 * i - 1, 2 * i - 1, :) = 1.0_dp
      T(2 * i - 1, 2 * i, :) = 1.0_dp
      T(2 * i, 2 * i, :) = 1.0_dp
    end do
    e = 0
    if (self%cov_level /= COV_NONE) then
      do i = 1, p
        R(2 * i - 1, i, :) = 1.0_dp
      end do
      Q(1:p, 1:p, 1) = cov_matrix(self%cov_level, p, params(1:kl))
      e = p
    end if
    if (self%cov_slope /= COV_NONE) then
      do i = 1, p
        R(2 * i, e + i, :) = 1.0_dp
      end do
      Q(e + 1:e + p, e + 1:e + p, 1) = cov_matrix(self%cov_slope, p, params(kl + 1:))
    end if
  end subroutine trend_fill

  function cov_transform_trend(self, unconstrained) result(constrained)
    class(trend_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)
    integer :: kl

    kl = cov_nparams(self%cov_level, self%p)
    constrained = [cov_constrain(self%cov_level, self%p, unconstrained(1:kl)), &
                   cov_constrain(self%cov_slope, self%p, unconstrained(kl + 1:))]
  end function cov_transform_trend

  function cov_untransform_trend(self, constrained) result(unconstrained)
    class(trend_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)
    integer :: kl

    kl = cov_nparams(self%cov_level, self%p)
    unconstrained = [cov_unconstrain(self%cov_level, self%p, constrained(1:kl)), &
                     cov_unconstrain(self%cov_slope, self%p, constrained(kl + 1:))]
  end function cov_untransform_trend

  ! ------------------------------------------------------------- seasonal

  !> States per series: s - 1 (dummy, trigonometric) or s (Harrison-Stevens);
  !> disturbances per series: 1 (dummy), s - 1 (trigonometric), s (H-S).
  subroutine seasonal_setup(self, y)
    class(seasonal_t), intent(inout) :: self
    real(dp), intent(in) :: y(:, :)
    integer :: ms, rs

    self%p = size(y, 1); self%n = size(y, 2)
    select case (self%form)
    case (SEASONAL_TRIG)
      ms = self%period - 1; rs = self%period - 1
    case (SEASONAL_HS)
      ms = self%period; rs = self%period
    case default
      ms = self%period - 1; rs = 1
    end select
    self%m = self%p * ms
    self%r = merge(0, self%p * rs, self%cov == COV_NONE)
    self%k = cov_nparams(self%cov, self%p)
  end subroutine seasonal_setup

  subroutine seasonal_fill(self, params, Z, H, T, R, Q)
    class(seasonal_t), intent(in) :: self
    real(dp), intent(in) :: params(:)
    real(dp), intent(inout) :: Z(:, :, :), H(:, :, :), T(:, :, :), R(:, :, :), Q(:, :, :)
    real(dp), allocatable :: Zs(:), Ts(:, :), Qs(:, :), S(:, :)
    integer :: i, j, ms, rs, sp, a, c, d

    sp = self%period
    call seasonal_block(self%form, sp, Zs, Ts, Qs)
    ms = size(Ts, 1); rs = size(Qs, 1)
    S = cov_matrix(self%cov, self%p, params)
    do i = 1, self%p
      a = (i - 1) * ms
      do j = 1, ms
        Z(i, a + j, :) = Zs(j)
      end do
      T(a + 1:a + ms, a + 1:a + ms, 1) = Ts
      if (self%r > 0) then
        c = (i - 1) * rs
        if (self%form == SEASONAL_DUMMY) then
          R(a + 1, c + 1, :) = 1.0_dp
        else
          do j = 1, rs
            R(a + j, c + j, :) = 1.0_dp
          end do
        end if
        ! Q = Sigma (x) Qs
        do j = 1, self%p
          d = (j - 1) * rs
          Q(c + 1:c + rs, d + 1:d + rs, 1) = S(i, j) * Qs
        end do
      end if
    end do
  end subroutine seasonal_fill

  !> One series' seasonal block: Z row, T and the disturbance correlation
  !> pattern Qs (DK 3.2.2-3.2.3).
  pure subroutine seasonal_block(form, s, Zs, Ts, Qs)
    integer, intent(in) :: form, s
    real(dp), allocatable, intent(out) :: Zs(:), Ts(:, :), Qs(:, :)
    real(dp) :: lambda
    integer :: j, ns, k

    select case (form)
    case (SEASONAL_TRIG)
      ! (gamma_1, gamma_1*, ..., ) with C_j blocks; for even s the last
      ! frequency pi has the single state with T = -1 (DK 3.2.3).
      allocate (Zs(s - 1), Ts(s - 1, s - 1), Qs(s - 1, s - 1), source=0.0_dp)
      ns = s / 2
      k = 0
      do j = 1, ns
        lambda = 2.0_dp * pi * j / s
        if (2 * j == s) then
          k = k + 1
          Zs(k) = 1.0_dp
          Ts(k, k) = -1.0_dp
        else
          Zs(k + 1) = 1.0_dp
          Ts(k + 1, k + 1) = cos(lambda); Ts(k + 1, k + 2) = sin(lambda)
          Ts(k + 2, k + 1) = -sin(lambda); Ts(k + 2, k + 2) = cos(lambda)
          k = k + 2
        end if
      end do
      do j = 1, s - 1
        Qs(j, j) = 1.0_dp
      end do
    case (SEASONAL_HS)
      ! (gamma_t, ..., gamma_t-s+1), T = [0 I; 1 0], Q = I - 1 1' / s
      allocate (Zs(s), Ts(s, s), Qs(s, s), source=0.0_dp)
      Zs(1) = 1.0_dp
      do j = 1, s - 1
        Ts(j, j + 1) = 1.0_dp
      end do
      Ts(s, 1) = 1.0_dp
      Qs = -1.0_dp / s
      do j = 1, s
        Qs(j, j) = Qs(j, j) + 1.0_dp
      end do
    case default
      ! Dummy (3.3): gamma_t+1 = -sum_j gamma_t+1-j + omega_t
      allocate (Zs(s - 1), Ts(s - 1, s - 1), Qs(1, 1), source=0.0_dp)
      Zs(1) = 1.0_dp
      Ts(1, :) = -1.0_dp
      do j = 2, s - 1
        Ts(j, j - 1) = 1.0_dp
      end do
      Qs(1, 1) = 1.0_dp
    end select
  end subroutine seasonal_block

  function cov_transform_seas(self, unconstrained) result(constrained)
    class(seasonal_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)

    constrained = cov_constrain(self%cov, self%p, unconstrained)
  end function cov_transform_seas

  function cov_untransform_seas(self, constrained) result(unconstrained)
    class(seasonal_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)

    unconstrained = cov_unconstrain(self%cov, self%p, constrained)
  end function cov_untransform_seas

  ! ---------------------------------------------------------------- cycle

  !> Parameters: the covariance of w, then lambda, then rho (if damped).
  subroutine cycle_setup(self, y)
    class(cycle_t), intent(inout) :: self
    real(dp), intent(in) :: y(:, :)

    self%p = size(y, 1); self%n = size(y, 2)
    if (self%period_max <= 0.0_dp) self%period_max = max(real(self%n, dp), self%period_min + 1.0_dp)
    self%m = 2 * self%p
    self%r = merge(0, 2 * self%p, self%cov == COV_NONE)
    self%k = cov_nparams(self%cov, self%p) + 1 + merge(1, 0, self%damped)
  end subroutine cycle_setup

  subroutine cycle_fill(self, params, Z, H, T, R, Q)
    class(cycle_t), intent(in) :: self
    real(dp), intent(in) :: params(:)
    real(dp), intent(inout) :: Z(:, :, :), H(:, :, :), T(:, :, :), R(:, :, :), Q(:, :, :)
    real(dp) :: lambda, rho, C(2, 2), S(self%p, self%p)
    integer :: i, j, kc

    kc = cov_nparams(self%cov, self%p)
    lambda = params(kc + 1)
    rho = 1.0_dp
    if (self%damped) rho = params(kc + 2)
    C = rho * reshape([cos(lambda), -sin(lambda), sin(lambda), cos(lambda)], [2, 2])
    S = cov_matrix(self%cov, self%p, params(1:kc))
    do i = 1, self%p
      Z(i, 2 * i - 1, :) = 1.0_dp
      T(2 * i - 1:2 * i, 2 * i - 1:2 * i, 1) = C
      if (self%r > 0) then
        R(2 * i - 1, 2 * i - 1, :) = 1.0_dp
        R(2 * i, 2 * i, :) = 1.0_dp
        do j = 1, self%p
          Q(2 * i - 1, 2 * j - 1, 1) = S(i, j)
          Q(2 * i, 2 * j, 1) = S(i, j)
        end do
      end if
    end do
  end subroutine cycle_fill

  function cycle_transform(self, unconstrained) result(constrained)
    class(cycle_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)
    integer :: kc

    kc = cov_nparams(self%cov, self%p)
    constrained = [cov_constrain(self%cov, self%p, unconstrained(1:kc)), &
                   constrain_interval(unconstrained(kc + 1), 2.0_dp * pi / self%period_max, &
                                      2.0_dp * pi / self%period_min)]
    if (self%damped) constrained = [constrained, constrain_interval(unconstrained(kc + 2), 0.0_dp, 1.0_dp)]
  end function cycle_transform

  function cycle_untransform(self, constrained) result(unconstrained)
    class(cycle_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)
    integer :: kc

    kc = cov_nparams(self%cov, self%p)
    unconstrained = [cov_unconstrain(self%cov, self%p, constrained(1:kc)), &
                     unconstrain_interval(constrained(kc + 1), 2.0_dp * pi / self%period_max, &
                                          2.0_dp * pi / self%period_min)]
    if (self%damped) unconstrained = [unconstrained, unconstrain_interval(constrained(kc + 2), 0.0_dp, 1.0_dp)]
  end function cycle_untransform

  function cycle_start(self, y) result(params)
    class(cycle_t), intent(in) :: self
    real(dp), intent(in) :: y(:, :)
    real(dp), allocatable :: params(:)
    integer :: kc

    kc = cov_nparams(self%cov, self%p)
    allocate (params(self%k))
    params(1:kc) = cov_start(self%cov, self%p, diff_variance(y) / 4.0_dp)
    ! frequency of the geometric mid-range period
    params(kc + 1) = 2.0_dp * pi / sqrt(self%period_min * self%period_max)
    if (self%damped) params(kc + 2) = 0.9_dp
  end function cycle_start

  !> Stationary when damped (DK 3.2.4), otherwise diffuse.
  function cycle_init_blocks(self) result(blocks)
    class(cycle_t), intent(in) :: self
    integer, allocatable :: blocks(:, :)

    blocks = reshape([1, self%m, merge(INIT_STATIONARY, INIT_DIFFUSE, self%damped)], [1, 3])
  end function cycle_init_blocks

  ! ----------------------------------------------------------- regression

  subroutine regression_setup(self, y)
    class(regression_t), intent(inout) :: self
    real(dp), intent(in) :: y(:, :)

    self%p = size(y, 1); self%n = size(y, 2)
    self%m = size(self%x, 2)
    if (.not. allocated(self%random_walk)) allocate (self%random_walk(self%m), source=.false.)
    self%r = count(self%random_walk)
    self%k = self%r
    self%tv_Z = .true.
  end subroutine regression_setup

  subroutine regression_fill(self, params, Z, H, T, R, Q)
    class(regression_t), intent(in) :: self
    real(dp), intent(in) :: params(:)
    real(dp), intent(inout) :: Z(:, :, :), H(:, :, :), T(:, :, :), R(:, :, :), Q(:, :, :)
    integer :: j, e

    Z(self%series, :, :) = transpose(self%x)
    e = 0
    do j = 1, self%m
      T(j, j, :) = 1.0_dp
      if (self%random_walk(j)) then
        e = e + 1
        R(j, e, :) = 1.0_dp
        Q(e, e, 1) = params(e)
      end if
    end do
  end subroutine regression_fill

  ! ---------------------------------------------------------- start values

  function irregular_start(self, y) result(params)
    class(irregular_t), intent(in) :: self
    real(dp), intent(in) :: y(:, :)
    real(dp), allocatable :: params(:)

    params = cov_start(self%cov, self%p, diff_variance(y) / 4.0_dp)
  end function irregular_start

  function level_start(self, y) result(params)
    class(level_t), intent(in) :: self
    real(dp), intent(in) :: y(:, :)
    real(dp), allocatable :: params(:)

    params = cov_start(self%cov, self%p, diff_variance(y) / 4.0_dp)
  end function level_start

  function trend_start(self, y) result(params)
    class(trend_t), intent(in) :: self
    real(dp), intent(in) :: y(:, :)
    real(dp), allocatable :: params(:)

    params = [cov_start(self%cov_level, self%p, diff_variance(y) / 4.0_dp), &
              cov_start(self%cov_slope, self%p, diff_variance(y) / 100.0_dp)]
  end function trend_start

  function seasonal_start(self, y) result(params)
    class(seasonal_t), intent(in) :: self
    real(dp), intent(in) :: y(:, :)
    real(dp), allocatable :: params(:)

    params = cov_start(self%cov, self%p, diff_variance(y) / 100.0_dp)
  end function seasonal_start

  ! ------------------------------------------------------- continuous time

  subroutine clevel_setup(self, y)
    class(continuous_level_t), intent(inout) :: self
    real(dp), intent(in) :: y(:, :)

    self%p = size(y, 1); self%n = size(y, 2)
    if (size(self%times) /= self%n) error stop "continuous_level_t: times must have one entry per observation"
    self%m = 1; self%r = 1; self%k = 1
    self%tv_Q = .true.
  end subroutine clevel_setup

  subroutine clevel_fill(self, params, Z, H, T, R, Q)
    class(continuous_level_t), intent(in) :: self
    real(dp), intent(in) :: params(:)
    real(dp), intent(inout) :: Z(:, :, :), H(:, :, :), T(:, :, :), R(:, :, :), Q(:, :, :)
    integer :: i

    Z(1, 1, :) = 1.0_dp
    T(1, 1, :) = 1.0_dp
    R(1, 1, :) = 1.0_dp
    do i = 1, self%n - 1
      Q(1, 1, i) = params(1) * (self%times(i + 1) - self%times(i))
    end do
    Q(1, 1, self%n) = 0.0_dp
  end subroutine clevel_fill

  subroutine ctrend_setup(self, y)
    class(continuous_trend_t), intent(inout) :: self
    real(dp), intent(in) :: y(:, :)

    self%p = size(y, 1); self%n = size(y, 2)
    if (size(self%times) /= self%n) error stop "continuous_trend_t: times must have one entry per observation"
    self%m = 2; self%r = 2; self%k = 1
    self%tv_T = .true.
    self%tv_Q = .true.
  end subroutine ctrend_setup

  subroutine ctrend_fill(self, params, Z, H, T, R, Q)
    class(continuous_trend_t), intent(in) :: self
    real(dp), intent(in) :: params(:)
    real(dp), intent(inout) :: Z(:, :, :), H(:, :, :), T(:, :, :), R(:, :, :), Q(:, :, :)
    real(dp) :: delta
    integer :: i

    Z(1, 1, :) = 1.0_dp
    R(1, 1, :) = 1.0_dp
    R(2, 2, :) = 1.0_dp
    do i = 1, self%n
      delta = 0.0_dp
      if (i < self%n) delta = self%times(i + 1) - self%times(i)
      T(:, :, i) = reshape([1.0_dp, 0.0_dp, delta, 1.0_dp], [2, 2])
      Q(:, :, i) = params(1) * delta * reshape([delta**2 / 3, delta / 2, delta / 2, 1.0_dp], [2, 2])
    end do
  end subroutine ctrend_fill

  ! --------------------------------------------------------- interventions

  !> Level shift at tau: w_t = 0 for t < tau, 1 for t >= tau (DK 3.2.5).
  pure function step_intervention(n, tau) result(w)
    integer, intent(in) :: n, tau
    real(dp) :: w(n)
    integer :: t

    w = [(merge(1.0_dp, 0.0_dp, t >= tau), t=1, n)]
  end function step_intervention

  !> Pulse at tau: w_t = 1 for t = tau, 0 otherwise.
  pure function pulse_intervention(n, tau) result(w)
    integer, intent(in) :: n, tau
    real(dp) :: w(n)
    integer :: t

    w = [(merge(1.0_dp, 0.0_dp, t == tau), t=1, n)]
  end function pulse_intervention

  !> Slope change at tau: w_t = 0 for t < tau, 1 + t - tau for t >= tau.
  pure function slope_intervention(n, tau) result(w)
    integer, intent(in) :: n, tau
    real(dp) :: w(n)
    integer :: t

    w = [(merge(1.0_dp + t - tau, 0.0_dp, t >= tau), t=1, n)]
  end function slope_intervention
end module statespace_structural
