program check
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite, new_testsuite, testsuite_type
  use test_filter, only: collect_filter
  use test_mle, only: collect_mle
  use test_simsmooth, only: collect_simsmooth
  use test_smoothing, only: collect_smoothing
  use test_score, only: collect_score
  use test_diagnostics, only: collect_diagnostics
  use test_recursions, only: collect_recursions
  use test_structural, only: collect_structural
  use test_arima, only: collect_arima
  implicit none
  type(testsuite_type), allocatable :: testsuites(:)
  integer :: stat, is

  stat = 0
  testsuites = [ &
               new_testsuite("filter", collect_filter), &
               new_testsuite("mle", collect_mle), &
               new_testsuite("simsmooth", collect_simsmooth), &
               new_testsuite("smoothing", collect_smoothing), &
               new_testsuite("score", collect_score), &
               new_testsuite("diagnostics", collect_diagnostics), &
               new_testsuite("recursions", collect_recursions), &
               new_testsuite("structural", collect_structural), &
               new_testsuite("arima", collect_arima) &
               ]
  do is = 1, size(testsuites)
    write (error_unit, '("# Testing:", 1x, a)') testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write (error_unit, '(i0, 1x, "test(s) failed!")') stat
    error stop
  end if
end program check
