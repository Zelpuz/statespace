!> Local level model for the Nile data (DK ch. 2) using the matrix API.
!> Produces the same output as example/locallevel.f90.
!>
!> Run from the repo root:  fpm run --example nile_local_level
program nile_local_level
  use statespace
  implicit none

  real(dp), allocatable :: y(:, :)
  type(ssm_rep_t) :: rep
  type(filter_result_t) :: fres
  type(smoother_result_t) :: sres
  integer :: info, t

  y = read_nile("data/nile.csv")

  ! y_t = alpha_t + eps_t,  alpha_t+1 = alpha_t + eta_t
  rep = ssm_rep(y, m=1, r=1)
  rep%Z = 1.0_dp
  rep%T = 1.0_dp
  rep%H = 15099.0_dp
  rep%Q = 1469.1_dp
  call rep%initialize_known([0.0_dp], reshape([1.0e7_dp], [1, 1]))

  call kalman_filter(rep, fres, info)
  if (info /= SS_OK) error stop "kalman_filter failed"
  call state_smoother(rep, fres, sres, info)
  if (info /= SS_OK) error stop "state_smoother failed"

  print '(*(g0, :, ","))', "t", "y", "a", "P", "v", "F", "alphahat", "V"
  do t = 1, rep%nobs
    print '(*(g0, :, ","))', t, y(1, t), fres%a(1, t), fres%P(1, 1, t), fres%v(1, t), &
      fres%F(1, 1, t), sres%alphahat(1, t), sres%V(1, 1, t)
  end do
  print '(a, g0)', "# loglik = ", fres%llf

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
end program nile_local_level
