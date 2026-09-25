!> Reader for the named-array fixtures written by test/fixtures/make_fixtures.py.
module fixture_io
  use statespace_kinds, only: dp
  use statespace_rep, only: ssm_rep_t, ssm_rep
  implicit none
  private

  public :: fixture_t, load_fixture, rep_from_fixture

  type :: fixture_array_t
    character(len=:), allocatable :: name
    integer, allocatable :: dims(:)
    real(dp), allocatable :: values(:)
  end type fixture_array_t

  type :: fixture_t
    type(fixture_array_t), allocatable :: arrays(:)
  contains
    procedure :: has
    procedure :: get1
    procedure :: get2
    procedure :: get3
  end type fixture_t

contains

  function load_fixture(path) result(fx)
    character(len=*), intent(in) :: path
    type(fixture_t) :: fx
    type(fixture_array_t) :: arr
    character(len=256) :: line, name
    integer :: unit, ios, ndim, i

    allocate (fx%arrays(0))
    open (newunit=unit, file=path, status='old', action='read')
    do
      read (unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      read (line, *) name, ndim
      allocate (arr%dims(ndim))
      read (line, *) name, ndim, arr%dims
      arr%name = trim(name)
      allocate (arr%values(product(arr%dims)))
      do i = 1, size(arr%values)
        read (unit, *) arr%values(i)
      end do
      fx%arrays = [fx%arrays, arr]
      deallocate (arr%name, arr%dims, arr%values)
    end do
    close (unit)
  end function load_fixture

  integer function find(self, name)
    class(fixture_t), intent(in) :: self
    character(len=*), intent(in) :: name

    do find = 1, size(self%arrays)
      if (self%arrays(find)%name == name) return
    end do
    error stop 'fixture array not found: '//name
  end function find

  logical function has(self, name)
    class(fixture_t), intent(in) :: self
    character(len=*), intent(in) :: name
    integer :: i

    has = .false.
    do i = 1, size(self%arrays)
      if (self%arrays(i)%name == name) has = .true.
    end do
  end function has

  !> Flattened values of the named array.
  function get1(self, name) result(x)
    class(fixture_t), intent(in) :: self
    character(len=*), intent(in) :: name
    real(dp), allocatable :: x(:)

    x = self%arrays(find(self, name))%values
  end function get1

  !> Named array as rank 2; lower-rank arrays get trailing dimensions of 1.
  function get2(self, name) result(x)
    class(fixture_t), intent(in) :: self
    character(len=*), intent(in) :: name
    real(dp), allocatable :: x(:, :)
    integer :: dims(2)

    associate (arr => self%arrays(find(self, name)))
      dims = 1
      dims(1:size(arr%dims)) = arr%dims
      x = reshape(arr%values, dims)
    end associate
  end function get2

  !> Named array as rank 3; lower-rank arrays get trailing dimensions of 1.
  function get3(self, name) result(x)
    class(fixture_t), intent(in) :: self
    character(len=*), intent(in) :: name
    real(dp), allocatable :: x(:, :, :)
    integer :: dims(3)

    associate (arr => self%arrays(find(self, name)))
      dims = 1
      dims(1:size(arr%dims)) = arr%dims
      x = reshape(arr%values, dims)
    end associate
  end function get3

  !> Build a representation from the model inputs stored in a fixture.
  function rep_from_fixture(fx) result(rep)
    type(fixture_t), intent(in) :: fx
    type(ssm_rep_t) :: rep
    real(dp), allocatable :: y(:, :)
    integer, allocatable :: blocks(:, :)
    integer :: k

    y = fx%get2('y')
    rep = ssm_rep(y, size(fx%get3('T'), 1), size(fx%get3('Q'), 1))
    rep%Z = fx%get3('Z')
    rep%H = fx%get3('H')
    rep%T = fx%get3('T')
    rep%R = fx%get3('R')
    rep%Q = fx%get3('Q')
    rep%c = fx%get2('c')
    rep%d = fx%get2('d')
    if (fx%has('init_blocks')) then
      blocks = nint(fx%get2('init_blocks'))
      do k = 1, size(blocks, 1)
        call rep%initialize_block(blocks(k, 1), blocks(k, 2), blocks(k, 3))
      end do
    else if (fx%has('a1')) then
      call rep%initialize_known(fx%get1('a1'), fx%get2('P1'))
    else
      call rep%initialize_stationary()
    end if
  end function rep_from_fixture
end module fixture_io
