!> A model whose `update` (and optionally the parameter transforms) are C
!> function pointers, for models defined in another language (the Python
!> package's MLEModel).
!>
!> update_cb(k, params, rep) must fill the model's representation, given as
!> a handle usable with the ss_rep_* routines (statespace_capi), from the k
!> constrained params, and return 0 on success. On failure the
!> representation's H is set to NaN, so the likelihood cannot be evaluated
!> (and the optimizer backs off); H is restored before the next update. transform_cb(k, x, out) and
!> untransform_cb(k, x, out) map between the unconstrained and constrained
!> parameters; without them the parameters are unconstrained.
module statespace_callback
  use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_funptr, c_f_procpointer, c_loc, &
                                         c_associated, c_null_funptr, c_null_ptr
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use statespace_kinds, only: dp
  use statespace_model, only: ssm_model_t
  implicit none
  private

  abstract interface
    integer(c_int) function update_fn(k, params, rep) bind(C)
      import :: c_int, c_ptr
      integer(c_int), value :: k
      type(c_ptr), value :: params, rep
    end function update_fn

    integer(c_int) function transform_fn(k, x, out) bind(C)
      import :: c_int, c_ptr
      integer(c_int), value :: k
      type(c_ptr), value :: x, out
    end function transform_fn
  end interface

  type, extends(ssm_model_t), public :: callback_model_t
    type(c_funptr) :: update_cb = c_null_funptr
    type(c_funptr) :: transform_cb = c_null_funptr
    type(c_funptr) :: untransform_cb = c_null_funptr
    !> Handle to this model's representation, passed to update_cb. Set by
    !> the owner once the model has its final address.
    type(c_ptr) :: rep_handle = c_null_ptr
    real(dp), allocatable :: start(:)
    !> Number of failed callbacks so far.
    integer :: failures = 0
    !> H before it was poisoned by a failed callback.
    real(dp), allocatable :: saved_H(:, :, :)
  contains
    procedure :: update => cb_update
    procedure :: start_params => cb_start
    procedure :: transform_params => cb_transform
    procedure :: untransform_params => cb_untransform
  end type callback_model_t

contains

  subroutine cb_update(self, params)
    class(callback_model_t), intent(inout) :: self
    real(dp), intent(in) :: params(:)
    real(dp), allocatable, target :: p(:)
    procedure(update_fn), pointer :: f
    integer(c_int) :: stat

    if (allocated(self%saved_H)) then
      call move_alloc(self%saved_H, self%rep%H)
    end if
    p = params
    call c_f_procpointer(self%update_cb, f)
    stat = f(int(size(p), c_int), c_loc(p), self%rep_handle)
    if (stat /= 0) then
      self%failures = self%failures + 1
      self%saved_H = self%rep%H
      self%rep%H = ieee_value(1.0_dp, ieee_quiet_nan)
    end if
  end subroutine cb_update

  function cb_start(self) result(params)
    class(callback_model_t), intent(in) :: self
    real(dp), allocatable :: params(:)

    params = self%start
  end function cb_start

  function cb_transform(self, unconstrained) result(constrained)
    class(callback_model_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)

    constrained = unconstrained
    if (c_associated(self%transform_cb)) call apply(self%transform_cb, unconstrained, constrained)
  end function cb_transform

  function cb_untransform(self, constrained) result(unconstrained)
    class(callback_model_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)

    unconstrained = constrained
    if (c_associated(self%untransform_cb)) call apply(self%untransform_cb, constrained, unconstrained)
  end function cb_untransform

  subroutine apply(fp, x, out)
    type(c_funptr), intent(in) :: fp
    real(dp), intent(in) :: x(:)
    real(dp), intent(inout) :: out(:)
    real(dp), allocatable, target :: xc(:), oc(:)
    procedure(transform_fn), pointer :: f
    integer(c_int) :: stat

    xc = x
    oc = out
    call c_f_procpointer(fp, f)
    stat = f(int(size(xc), c_int), c_loc(xc), c_loc(oc))
    out = oc
    if (stat /= 0) out = ieee_value(1.0_dp, ieee_quiet_nan)
  end subroutine apply
end module statespace_callback
