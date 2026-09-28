!> Kind parameters, constants, and status codes shared by the library.
module statespace_kinds
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private

  integer, parameter, public :: dp = real64

  real(dp), parameter, public :: log2pi = 1.8378770664093454835606594728112_dp

  !> Status codes returned through `info` arguments.
  integer, parameter, public :: SS_OK = 0
  integer, parameter, public :: SS_ERR_DIM = 1        !< inconsistent array dimensions
  integer, parameter, public :: SS_ERR_NOT_PD = 2     !< matrix not positive definite
  integer, parameter, public :: SS_ERR_INIT = 3       !< missing or invalid
                                                      !! initialization
  integer, parameter, public :: SS_ERR_UNSUPPORTED = 4 !< feature not implemented yet
  integer, parameter, public :: SS_ERR_SINGULAR = 5   !< singular linear system
  integer, parameter, public :: SS_ERR_NOT_STATIONARY = 6 !< T has an eigenvalue
                                                          !! |lambda| >= 1
  integer, parameter, public :: SS_ERR_NOT_CONVERGED = 7  !< an iteration did not
                                                          !! converge
end module statespace_kinds
