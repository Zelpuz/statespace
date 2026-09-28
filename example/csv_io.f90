!> Minimal CSV reader for the examples: a header line, then numeric fields.
!> Fields that are not numbers (row names, dates, "NA") read as NaN.
module csv_io
  use statespace, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  private

  public :: read_csv, column

contains

  !> Read `path` into names(ncol) and data(nrow, ncol).
  subroutine read_csv(path, names, data)
    character(len=*), intent(in) :: path
    character(len=32), allocatable, intent(out) :: names(:)
    real(dp), allocatable, intent(out) :: data(:, :)
    character(len=4096) :: line
    integer :: unit, ios, nrow, i

    open (newunit=unit, file=path, status="old", action="read")
    read (unit, '(a)') line
    names = split(line)
    nrow = 0
    do
      read (unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) nrow = nrow + 1
    end do
    rewind (unit)
    read (unit, '(a)') line
    allocate (data(nrow, size(names)))
    i = 0
    do
      read (unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      i = i + 1
      data(i, :) = to_real(split(line), size(names))
    end do
    close (unit)
  end subroutine read_csv

  !> The column called `name`.
  function column(names, data, name) result(x)
    character(len=*), intent(in) :: names(:), name
    real(dp), intent(in) :: data(:, :)
    real(dp), allocatable :: x(:)
    integer :: j

    do j = 1, size(names)
      if (trim(names(j)) == name) then
        x = data(:, j)
        return
      end if
    end do
    error stop "csv_io: no column "//name
  end function column

  function split(line) result(fields)
    character(len=*), intent(in) :: line
    character(len=32), allocatable :: fields(:)
    integer :: i0, i

    allocate (fields(0))
    i0 = 1
    do i = 1, len_trim(line) + 1
      if (i > len_trim(line)) then
        fields = [fields, unquote(line(i0:i - 1))]
      else if (line(i:i) == ",") then
        fields = [fields, unquote(line(i0:i - 1))]
        i0 = i + 1
      end if
    end do
  end function split

  function unquote(s) result(t)
    character(len=*), intent(in) :: s
    character(len=32) :: t
    integer :: i, k

    t = ""
    k = 0
    do i = 1, len(s)
      if (s(i:i) /= '"' .and. k < len(t)) then
        k = k + 1
        t(k:k) = s(i:i)
      end if
    end do
  end function unquote

  function to_real(fields, n) result(x)
    character(len=32), intent(in) :: fields(:)
    integer, intent(in) :: n
    real(dp) :: x(n)
    integer :: j, ios

    x = ieee_value(1.0_dp, ieee_quiet_nan)
    do j = 1, min(n, size(fields))
      if (verify(trim(fields(j)), "0123456789.-+eE") /= 0 &
          .or. len_trim(fields(j)) == 0) cycle
      read (fields(j), *, iostat=ios) x(j)
      if (ios /= 0) x(j) = ieee_value(1.0_dp, ieee_quiet_nan)
    end do
  end function to_real
end module csv_io
