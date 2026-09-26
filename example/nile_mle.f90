!> Maximum likelihood estimation of the Nile local level model with exact
!> diffuse initialization (DK 2.10; estimates 15099 and 1469.1).
!> Shows how to define a model by extending `ssm_model_t`.
!>
!> Run from the repo root:  fpm run --example nile_mle
module nile_local_level_model
  use statespace
  implicit none
  private

  public :: local_level_t, local_level

  !> y_t = alpha_t + eps_t,  alpha_t+1 = alpha_t + eta_t,
  !> params = [sigma2_eps, sigma2_eta].
  type, extends(ssm_model_t) :: local_level_t
  contains
    procedure :: update
    procedure :: start_params
    procedure :: transform_params
    procedure :: untransform_params
    procedure :: param_names
  end type local_level_t

contains

  function local_level(y) result(mod)
    real(dp), intent(in) :: y(:, :)
    type(local_level_t) :: mod

    mod%k_params = 2
    mod%rep = ssm_rep(y, m=1, r=1)
    mod%rep%Z = 1.0_dp
    mod%rep%T = 1.0_dp
    call mod%rep%initialize_diffuse()
  end function local_level

  subroutine update(self, params)
    class(local_level_t), intent(inout) :: self
    real(dp), intent(in) :: params(:)

    self%rep%H(1, 1, 1) = params(1)
    self%rep%Q(1, 1, 1) = params(2)
  end subroutine update

  function start_params(self) result(params)
    class(local_level_t), intent(in) :: self
    real(dp), allocatable :: params(:)
    real(dp) :: v

    associate (y => self%rep%y(1, :))
      v = sum((y - sum(y) / size(y))**2) / size(y)
    end associate
    params = [v / 2, v / 2]
  end function start_params

  function transform_params(self, unconstrained) result(constrained)
    class(local_level_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)

    constrained = constrain_positive(unconstrained)
  end function transform_params

  function untransform_params(self, constrained) result(unconstrained)
    class(local_level_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)

    unconstrained = unconstrain_positive(constrained)
  end function untransform_params

  function param_names(self) result(names)
    class(local_level_t), intent(in) :: self
    character(len=32), allocatable :: names(:)

    names = [character(len=32) :: "sigma2.irregular", "sigma2.level"]
  end function param_names
end module nile_local_level_model

program nile_mle
  use statespace
  use nile_local_level_model, only: local_level_t, local_level
  implicit none

  type(local_level_t) :: mod
  type(fit_result_t) :: res
  type(fit_options_t) :: opts
  character(len=32), allocatable :: names(:)
  real(dp), allocatable :: y(:, :)
  integer :: info, i

  y = read_nile("data/nile.csv")
  mod = local_level(y)

  opts%factr = 10.0_dp       ! tight convergence; the likelihood is flat
  opts%pgtol = 1.0e-9_dp
  call fit(mod, res, options=opts, info=info)
  if (info /= SS_OK) error stop "fit failed"

  names = mod%param_names()
  print '(a)', trim(res%message)
  print '(a, i0, a, i0)', "iterations: ", res%niter, "  likelihood evaluations: ", res%nfev
  print '(a, f14.6)', "log likelihood: ", res%llf
  print '(a, f10.4, a, f10.4)', "AIC: ", res%aic, "  BIC: ", res%bic
  print '(a20, 2a14)', "", "estimate", "std err"
  do i = 1, mod%k_params
    print '(a20, 2f14.4)', trim(names(i)), res%params(i), res%bse(i)
  end do

contains

  function read_nile(path) result(y)
    character(len=*), intent(in) :: path
    real(dp), allocatable :: y(:, :)
    real(dp) :: year, vol
    integer :: unit, ios, n

    open (newunit=unit, file=path, status="old", action="read")
    read (unit, *)
    n = 0
    do
      read (unit, *, iostat=ios) year, vol
      if (ios /= 0) exit
      n = n + 1
    end do
    rewind (unit)
    read (unit, *)
    allocate (y(1, n))
    do n = 1, size(y, 2)
      read (unit, *) year, y(1, n)
    end do
    close (unit)
  end function read_nile
end program nile_mle
