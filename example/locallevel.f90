!> Example of local level model written with explicit loops and the DK2e notation
!> Contrast this with example/nile_local_level.f90, which uses the statespace lib
!> and produces the same result.
!> Run from repo root:  fpm run --example locallevel

program locallevel
    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307)
    character(len=*), parameter :: filename = "data/nile.csv"
    real(dp), allocatable :: datetime(:), y(:), a(:), P(:), v(:), F(:), K(:)
    real(dp), allocatable :: a_given(:), P_given(:)
    real(dp), allocatable :: r(:), N(:), v_up(:), L(:), expected_alpha(:)
    character(len=200) :: line
    integer :: unit, ios, len, t
    real(dp) :: a_0, P_0, var_epsilon, var_eta


    open(newunit=unit, file=filename, status="old", action="read")

    ! determine file length loop
    read(unit, "(A)", iostat=ios) line ! skip header row
    len = 0
    do
        read(unit, "(A)", iostat=ios) line
        if (ios /= 0) exit
        len = len + 1
    end do

    rewind(unit) ! go back to head of file
    read(unit, "(A)", iostat=ios) line ! skip header again

    allocate(datetime(len), y(len), a(len), P(len), v(len), F(len), K(len))
    allocate(a_given(len), P_given(len), r(len), N(len), v_up(len), L(len))
    allocate(expected_alpha(len))

    ! file read loop
    do t = 1, len
        read(unit, *) datetime(t), y(t)
    end do

    var_epsilon = 15099.0
    var_eta = 1469.1
    a_0 = 0.0
    P_0 = 1e7

    a(1) = a_0
    P(1) = P_0

    ! kalman filter loop
    do t = 1, len
        v(t) = y(t) - a(t)
        F(t) = P(t) + var_epsilon
        K(t) = P(t) / F(t)
        a_given(t) = a(t) + (K(t) * v(t))
        P_given(t) = (1 - K(t)) * P(t)
        if (t < len) then
            a(t+1) = a_given(t)
            P(t+1) = P_given(t) + var_eta
        end if
    end do

    r(len) = 0.0
    N(len) = 0.0

    ! state smoother loop
    do t = len, 1, -1
        L(t) = 1 - K(t)
        if (t > 1) then
            r(t-1) = v(t) / F(t) + L(t) * r(t)
            expected_alpha(t) = a(t) + P(t) * r(t-1)
            N(t-1) = 1 / F(t) + N(t) * L(t) ** 2
            v_up(t) = P(t) - N(t-1) * P(t) ** 2
        else
            expected_alpha(t) = a(t) + P(t) * (v(t) / F(t) + L(t) * r(t))
            v_up(t) = P(t) - (1 / F(t) + N(t) * L(t) ** 2) * P(t) ** 2
        end if
    end do

    ! print results
    print '(*(g0, :, ","))', "t", "y(t)", "a(t)", &
                             "v(t)", "a_given(t)", &
                             "P(t)", "F(t)", "K(t)", &
                             "P_given(t)", "r(t)", &
                             "expected_alpha(t)", "N(t)", &
                             "v_up(t)"
    do t = 1, len
        print '(*(g0, :, ","))', t, y(t), a(t), v(t), &
                                 a_given(t), P(t), F(t), &
                                 K(t), P_given(t), r(t), &
                                 expected_alpha(t), N(t), &
                                 v_up(t)
    end do

end program locallevel
