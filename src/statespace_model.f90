!> Abstract parameterized state space model, the Fortran counterpart of
!> statsmodels' `MLEModel`.
!>
!> A user model extends `ssm_model_t`, sets up `rep` (data, fixed matrices,
!> initialization) in its constructor, sets `k_params`, and implements:
!>
!>   * `update(params)`: write constrained `params` into the system matrices.
!>   * `start_params()`: constrained starting values for the optimizer.
!>
!> and optionally overrides `transform_params` (unconstrained -> constrained),
!> `untransform_params` (its inverse), and `param_names`. The optimizer works
!> on unconstrained parameters; `update` always receives constrained ones.
module statespace_model
  use statespace_kinds, only: dp, SS_OK
  use statespace_rep, only: ssm_rep_t
  use statespace_filter, only: filter_result_t, kalman_filter, loglike, loglike_concentrated
  use statespace_smoother, only: smoother_result_t, state_smoother
  implicit none
  private

  public :: ssm_model_t, constrain_positive, unconstrain_positive
  public :: constrain_stationary, unconstrain_stationary

  type, abstract :: ssm_model_t
    type(ssm_rep_t) :: rep
    integer :: k_params = 0
    !> Concentrate a scale sigma^2 out of the likelihood (DK 2.10.2): H, Q
    !> and P_star set by `update` are relative to sigma^2, which is not among
    !> the parameters. `scale` holds its estimate at the last evaluation.
    logical :: concentrate_scale = .false.
    real(dp) :: scale = 1.0_dp
  contains
    procedure(update_iface), deferred :: update
    procedure(start_params_iface), deferred :: start_params
    procedure :: transform_params
    procedure :: untransform_params
    procedure :: param_names
    procedure :: loglike => model_loglike
    procedure :: rep_at
    procedure :: filter => model_filter
    procedure :: smooth => model_smooth
  end type ssm_model_t

  abstract interface
    subroutine update_iface(self, params)
      import :: ssm_model_t, dp
      class(ssm_model_t), intent(inout) :: self
      real(dp), intent(in) :: params(:)
    end subroutine update_iface

    function start_params_iface(self) result(params)
      import :: ssm_model_t, dp
      class(ssm_model_t), intent(in) :: self
      real(dp), allocatable :: params(:)
    end function start_params_iface
  end interface

contains

  !> Map unconstrained optimizer parameters to constrained model parameters.
  !> The default is the identity.
  function transform_params(self, unconstrained) result(constrained)
    class(ssm_model_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)

    constrained = unconstrained
  end function transform_params

  !> Inverse of `transform_params`. The default is the identity.
  function untransform_params(self, constrained) result(unconstrained)
    class(ssm_model_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)

    unconstrained = constrained
  end function untransform_params

  !> Parameter names; the default is param1, param2, ...
  function param_names(self) result(names)
    class(ssm_model_t), intent(in) :: self
    character(len=32), allocatable :: names(:)
    integer :: i

    allocate (names(self%k_params))
    do i = 1, self%k_params
      write (names(i), '("param", i0)') i
    end do
  end function param_names

  !> Log likelihood at `params`. Pass `transformed = .false.` if `params` are
  !> unconstrained.
  function model_loglike(self, params, info, transformed) result(llf)
    class(ssm_model_t), intent(inout) :: self
    real(dp), intent(in) :: params(:)
    integer, intent(out) :: info
    logical, intent(in), optional :: transformed
    real(dp) :: llf

    call update_params(self, params, transformed)
    if (self%concentrate_scale) then
      llf = loglike_concentrated(self%rep, self%scale, info)
    else
      llf = loglike(self%rep, info)
    end if
  end function model_loglike

  !> The representation at `params`, with H, Q and P_star multiplied by the
  !> estimated scale if it is concentrated out. Filtering and smoothing this
  !> gives results on the scale of the data.
  subroutine rep_at(self, params, rep, info, transformed)
    class(ssm_model_t), intent(inout) :: self
    real(dp), intent(in) :: params(:)
    type(ssm_rep_t), intent(out) :: rep
    integer, intent(out) :: info
    logical, intent(in), optional :: transformed
    real(dp) :: llf

    call update_params(self, params, transformed)
    rep = self%rep
    info = SS_OK
    if (self%concentrate_scale) then
      llf = loglike_concentrated(self%rep, self%scale, info)
      if (info /= SS_OK) return
      call rep%scale_by(self%scale)
    end if
  end subroutine rep_at

  !> Kalman filter at `params`.
  subroutine model_filter(self, params, fres, info, transformed)
    class(ssm_model_t), intent(inout) :: self
    real(dp), intent(in) :: params(:)
    type(filter_result_t), intent(out) :: fres
    integer, intent(out) :: info
    logical, intent(in), optional :: transformed
    type(ssm_rep_t) :: rep

    call self%rep_at(params, rep, info, transformed)
    if (info /= SS_OK) return
    call kalman_filter(rep, fres, info)
  end subroutine model_filter

  !> Kalman filter and state smoother at `params`.
  subroutine model_smooth(self, params, fres, sres, info, transformed)
    class(ssm_model_t), intent(inout) :: self
    real(dp), intent(in) :: params(:)
    type(filter_result_t), intent(out) :: fres
    type(smoother_result_t), intent(out) :: sres
    integer, intent(out) :: info
    logical, intent(in), optional :: transformed

    type(ssm_rep_t) :: rep

    call self%rep_at(params, rep, info, transformed)
    if (info /= SS_OK) return
    call kalman_filter(rep, fres, info)
    if (info /= SS_OK) return
    call state_smoother(rep, fres, sres, info)
  end subroutine model_smooth

  subroutine update_params(self, params, transformed)
    class(ssm_model_t), intent(inout) :: self
    real(dp), intent(in) :: params(:)
    logical, intent(in), optional :: transformed
    logical :: constrained

    constrained = .true.
    if (present(transformed)) constrained = transformed
    if (constrained) then
      call self%update(params)
    else
      call self%update(self%transform_params(params))
    end if
  end subroutine update_params

  !> Transform for positive parameters such as variances (DK 7.3.2):
  !> sigma^2 = exp(2 psi), i.e. psi = log(sigma^2) / 2.
  elemental real(dp) function constrain_positive(x)
    real(dp), intent(in) :: x

    constrain_positive = exp(2.0_dp * x)
  end function constrain_positive

  !> Inverse of `constrain_positive`: psi = log(sigma^2) / 2.
  elemental real(dp) function unconstrain_positive(x)
    real(dp), intent(in) :: x

    unconstrain_positive = 0.5_dp * log(x)
  end function unconstrain_positive

  !> Map unconstrained values to the coefficients phi of a stationary
  !> AR(n) polynomial 1 - phi_1 L - ... - phi_n L^n (Monahan 1984), matching
  !> statsmodels' `constrain_stationary_univariate`. The partial
  !> autocorrelations are x / sqrt(1 + x^2), DK's bounded transform
  !> (7.3.2) with a = 1. Also used for invertible MA coefficients.
  pure function constrain_stationary(unconstrained) result(phi)
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: phi(:)
    real(dp), allocatable :: y(:, :), r(:)
    integer :: n, k, i

    n = size(unconstrained)
    allocate (y(n, n), source=0.0_dp)
    r = unconstrained / sqrt(1.0_dp + unconstrained**2)   ! partial autocorrelations
    do k = 1, n
      do i = 1, k - 1
        y(k, i) = y(k - 1, i) + r(k) * y(k - 1, k - i)
      end do
      y(k, k) = r(k)
    end do
    phi = -y(n, :)
  end function constrain_stationary

  !> Inverse of `constrain_stationary`.
  pure function unconstrain_stationary(phi) result(unconstrained)
    real(dp), intent(in) :: phi(:)
    real(dp), allocatable :: unconstrained(:)
    real(dp), allocatable :: y(:, :), r(:)
    integer :: n, k, i

    n = size(phi)
    allocate (y(n, n), source=0.0_dp, r(n))
    y(n, :) = -phi
    do k = n, 2, -1
      do i = 1, k - 1
        y(k - 1, i) = (y(k, i) - y(k, k) * y(k, k - i)) / (1.0_dp - y(k, k)**2)
      end do
    end do
    do k = 1, n
      r(k) = y(k, k)
    end do
    unconstrained = r / sqrt(1.0_dp - r**2)
  end function unconstrain_stationary
end module statespace_model
