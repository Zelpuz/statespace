!> A model defined by a map from parameters to entries of the system
!> matrices, for models written without Fortran code (e.g. from Python).
!>
!> The representation holds every fixed entry and the initialization; the
!> map says which entries the parameters set:
!>
!>   * entries: parameter k sets X(i, j, t) = coef * psi_k, for X one of Z,
!>     H, T, R, Q, c, d (codes as SS_ARR_* in statespace_capi: 2..8) and
!>     t = 0 for every time slice;
!>   * covariance blocks: parameters k..k+q-1, q = dim (dim + 1) / 2, are the
!>     lower triangle (by columns) of a dim x dim block of H or Q starting
!>     at row and column `offset`.
!>
!> Parameter groups set the transforms to the optimizer's unconstrained
!> scale: MAP_POSITIVE (sigma^2 = exp(2 x), DK 7.3.2), MAP_INTERVAL
!> (lo, hi), MAP_STATIONARY (AR coefficients, Monahan), MAP_INVERTIBLE (MA
!> coefficients), and covariance blocks (Cholesky factor with log
!> diagonal). Parameters in no group are unconstrained.
module statespace_mapped
  use statespace_kinds, only: dp
  use statespace_model, only: ssm_model_t, constrain_positive, unconstrain_positive, &
                              constrain_stationary, unconstrain_stationary
  use statespace_components, only: constrain_interval, unconstrain_interval
  use statespace_structural, only: cov_constrain, cov_unconstrain, COV_FULL
  implicit none
  private

  integer, parameter, public :: MAP_NONE = 0, MAP_POSITIVE = 1, MAP_INTERVAL = 2, &
                                MAP_STATIONARY = 3, MAP_INVERTIBLE = 4, MAP_COV = 5

  type, public :: map_entry_t
    integer :: param = 0, array = 0, i = 0, j = 0, t = 0
    real(dp) :: coef = 1.0_dp
  end type map_entry_t

  type, public :: map_block_t
    integer :: array = 0, offset = 0, dim = 0, first = 0
  end type map_block_t

  type, public :: map_group_t
    integer :: kind = MAP_NONE, first = 0, last = 0
    real(dp) :: lo = 0.0_dp, hi = 0.0_dp
  end type map_group_t

  type, extends(ssm_model_t), public :: mapped_model_t
    type(map_entry_t), allocatable :: entries(:)
    type(map_block_t), allocatable :: blocks(:)
    type(map_group_t), allocatable :: groups(:)
    real(dp), allocatable :: start(:)
    character(len=32), allocatable :: names(:)
  contains
    procedure :: update => mapped_update
    procedure :: start_params => mapped_start
    procedure :: transform_params => mapped_transform
    procedure :: untransform_params => mapped_untransform
    procedure :: param_names => mapped_names
    procedure :: add_entry, add_block, add_group
  end type mapped_model_t

  public :: mapped_model

contains

  !> A mapped model with k parameters over `rep` (copied), with start values
  !> 0.1 until set.
  function mapped_model(rep, k) result(model)
    use statespace_rep, only: ssm_rep_t
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: k
    type(mapped_model_t) :: model
    integer :: i

    model%rep = rep
    model%k_params = k
    allocate (model%entries(0), model%blocks(0), model%groups(0))
    allocate (model%start(k), source=0.1_dp)
    allocate (model%names(k))
    do i = 1, k
      write (model%names(i), '("param", i0)') i
    end do
  end function mapped_model

  subroutine add_entry(self, param, array, i, j, t, coef)
    class(mapped_model_t), intent(inout) :: self
    integer, intent(in) :: param, array, i, j, t
    real(dp), intent(in) :: coef

    self%entries = [self%entries, map_entry_t(param, array, i, j, t, coef)]
  end subroutine add_entry

  !> A covariance block, which also becomes a MAP_COV parameter group.
  subroutine add_block(self, array, offset, dim, first)
    class(mapped_model_t), intent(inout) :: self
    integer, intent(in) :: array, offset, dim, first

    self%blocks = [self%blocks, map_block_t(array, offset, dim, first)]
    call self%add_group(MAP_COV, first, first + dim * (dim + 1) / 2 - 1, &
                        real(dim, dp), 0.0_dp)
  end subroutine add_block

  subroutine add_group(self, kind, first, last, lo, hi)
    class(mapped_model_t), intent(inout) :: self
    integer, intent(in) :: kind, first, last
    real(dp), intent(in) :: lo, hi

    self%groups = [self%groups, map_group_t(kind, first, last, lo, hi)]
  end subroutine add_group

  subroutine mapped_update(self, params)
    class(mapped_model_t), intent(inout) :: self
    real(dp), intent(in) :: params(:)
    integer :: e, b, i, j, k, o

    do e = 1, size(self%entries)
      associate (en => self%entries(e))
        call set(en%array, en%i, en%j, en%t, en%coef * params(en%param))
      end associate
    end do
    do b = 1, size(self%blocks)
      associate (bl => self%blocks(b))
        k = bl%first
        o = bl%offset - 1
        do j = 1, bl%dim
          do i = j, bl%dim
            call set(bl%array, o + i, o + j, 0, params(k))
            if (i /= j) call set(bl%array, o + j, o + i, 0, params(k))
            k = k + 1
          end do
        end do
      end associate
    end do

  contains

    !> X(i, j, t) = value; t = 0: every time slice. c and d use i only.
    subroutine set(array, i, j, t, value)
      integer, intent(in) :: array, i, j, t
      real(dp), intent(in) :: value

      select case (array)
      case (2)
        call put3(self%rep%Z, i, j, t, value)
      case (3)
        call put3(self%rep%H, i, j, t, value)
      case (4)
        call put3(self%rep%T, i, j, t, value)
      case (5)
        call put3(self%rep%R, i, j, t, value)
      case (6)
        call put3(self%rep%Q, i, j, t, value)
      case (7)
        call put2(self%rep%c, i, t, value)
      case (8)
        call put2(self%rep%d, i, t, value)
      end select
    end subroutine set
  end subroutine mapped_update

  subroutine put3(X, i, j, t, value)
    real(dp), intent(inout) :: X(:, :, :)
    integer, intent(in) :: i, j, t
    real(dp), intent(in) :: value

    if (t == 0 .or. size(X, 3) == 1) then
      X(i, j, :) = value
    else
      X(i, j, t) = value
    end if
  end subroutine put3

  subroutine put2(X, i, t, value)
    real(dp), intent(inout) :: X(:, :)
    integer, intent(in) :: i, t
    real(dp), intent(in) :: value

    if (t == 0 .or. size(X, 2) == 1) then
      X(i, :) = value
    else
      X(i, t) = value
    end if
  end subroutine put2

  function mapped_start(self) result(params)
    class(mapped_model_t), intent(in) :: self
    real(dp), allocatable :: params(:)

    params = self%start
  end function mapped_start

  function mapped_names(self) result(names)
    class(mapped_model_t), intent(in) :: self
    character(len=32), allocatable :: names(:)

    names = self%names
  end function mapped_names

  function mapped_transform(self, unconstrained) result(constrained)
    class(mapped_model_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)
    integer :: g

    constrained = unconstrained
    do g = 1, size(self%groups)
      associate (gr => self%groups(g), &
                 x => unconstrained(self%groups(g)%first:self%groups(g)%last))
        select case (gr%kind)
        case (MAP_POSITIVE)
          constrained(gr%first:gr%last) = constrain_positive(x)
        case (MAP_INTERVAL)
          constrained(gr%first:gr%last) = constrain_interval(x, gr%lo, gr%hi)
        case (MAP_STATIONARY)
          constrained(gr%first:gr%last) = constrain_stationary(x)
        case (MAP_INVERTIBLE)
          constrained(gr%first:gr%last) = -constrain_stationary(x)
        case (MAP_COV)
          constrained(gr%first:gr%last) = cov_constrain(COV_FULL, nint(gr%lo), x)
        end select
      end associate
    end do
  end function mapped_transform

  function mapped_untransform(self, constrained) result(unconstrained)
    class(mapped_model_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)
    integer :: g

    unconstrained = constrained
    do g = 1, size(self%groups)
      associate (gr => self%groups(g), &
                 c => constrained(self%groups(g)%first:self%groups(g)%last))
        select case (gr%kind)
        case (MAP_POSITIVE)
          unconstrained(gr%first:gr%last) = unconstrain_positive(c)
        case (MAP_INTERVAL)
          unconstrained(gr%first:gr%last) = unconstrain_interval(c, gr%lo, gr%hi)
        case (MAP_STATIONARY)
          unconstrained(gr%first:gr%last) = unconstrain_stationary(c)
        case (MAP_INVERTIBLE)
          unconstrained(gr%first:gr%last) = unconstrain_stationary(-c)
        case (MAP_COV)
          unconstrained(gr%first:gr%last) = cov_unconstrain(COV_FULL, nint(gr%lo), c)
        end select
      end associate
    end do
  end function mapped_untransform
end module statespace_mapped
