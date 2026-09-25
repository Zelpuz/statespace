!> Simulation and the simulation smoother.
module test_simsmooth
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use statespace
  use fixture_io, only: fixture_t, load_fixture, rep_from_fixture
  implicit none
  private

  public :: collect_simsmooth

  real(dp), parameter :: rtol = 1.0e-9_dp

contains

  subroutine collect_simsmooth(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("known_init", test_known), &
                new_unittest("time_varying", test_timevarying), &
                new_unittest("missing", test_missing), &
                new_unittest("stationary_init", test_stationary), &
                new_unittest("exact_diffuse", test_diffuse), &
                new_unittest("exact_diffuse_seasonal_missing", test_diffuse_seasonal), &
                new_unittest("draw_moments_missing", test_moments_missing), &
                new_unittest("draw_moments_mixed_init", test_moments_mixed) &
                ]
  end subroutine collect_simsmooth

  pure real(dp) function relerr(actual, expected)
    real(dp), intent(in) :: actual(:), expected(:)

    relerr = maxval(abs(actual - expected)) / max(1.0_dp, maxval(abs(expected)))
  end function relerr

  subroutine check_close(error, actual, expected, label)
    type(error_type), allocatable, intent(out) :: error
    real(dp), intent(in) :: actual(:), expected(:)
    character(len=*), intent(in) :: label
    character(len=32) :: buf

    write (buf, '(es10.3)') relerr(actual, expected)
    call check(error, relerr(actual, expected) <= rtol, label//": relative error "//trim(buf))
  end subroutine check_close

  !> Same variates as statsmodels must give the same draws. eps is compared at
  !> observed elements only: for missing elements statsmodels returns the
  !> unconditional draw eps+. statsmodels leaves a1 out of its generated series
  !> when there are no missing data (harmless, it cancels), so those are only
  !> compared when y has missing values.
  subroutine check_sim_fixture(error, base)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), intent(in) :: base
    type(fixture_t) :: fx, sx
    type(ssm_rep_t) :: rep
    type(simsmooth_result_t) :: sim
    real(dp), allocatable :: y(:, :), alpha(:, :), eps(:, :), eta(:, :), gen_state(:, :)
    logical, allocatable :: obs(:, :)
    integer :: info, n

    fx = load_fixture("test/fixtures/"//base//".txt")
    sx = load_fixture("test/fixtures/sim_"//base//".txt")
    rep = rep_from_fixture(fx)
    n = rep%nobs
    allocate (y(rep%k_endog, n), alpha(rep%k_states, n), eps(rep%k_endog, n), &
              eta(rep%k_posdef, n))

    call simulate(rep, sx%get1('u_init'), sx%get2('u_eps'), sx%get2('u_eta'), y, alpha, eps, &
                  eta, info)
    call check(error, info, SS_OK, "simulate info")
    if (allocated(error)) return
    if (any(ieee_is_nan(rep%y))) then
      gen_state = sx%get2('generated_state')
      call check_close(error, pack(y, .true.), sx%get1('generated_obs'), "generated_obs")
      if (allocated(error)) return
      call check_close(error, pack(alpha, .true.), pack(gen_state(:, 1:n), .true.), &
                       "generated_state")
      if (allocated(error)) return
    end if

    call simulation_smoother(rep, sx%get1('u_init'), sx%get2('u_eps'), sx%get2('u_eta'), sim, &
                             info)
    call check(error, info, SS_OK, "simulation_smoother info")
    if (allocated(error)) return
    call check_close(error, pack(sim%state, .true.), sx%get1('state'), "state")
    if (allocated(error)) return
    call check_close(error, pack(sim%eta, .true.), sx%get1('eta'), "eta")
    if (allocated(error)) return
    obs = .not. ieee_is_nan(rep%y)
    call check_close(error, pack(sim%eps, obs), pack(sx%get2('eps'), obs), "eps")
  end subroutine check_sim_fixture

  subroutine test_known(error)
    type(error_type), allocatable, intent(out) :: error

    call check_sim_fixture(error, "mv_invariant")
  end subroutine test_known

  subroutine test_timevarying(error)
    type(error_type), allocatable, intent(out) :: error

    call check_sim_fixture(error, "mv_timevarying")
  end subroutine test_timevarying

  subroutine test_missing(error)
    type(error_type), allocatable, intent(out) :: error

    call check_sim_fixture(error, "mv_missing")
  end subroutine test_missing

  subroutine test_stationary(error)
    type(error_type), allocatable, intent(out) :: error

    call check_sim_fixture(error, "mv_stationary")
  end subroutine test_stationary

  subroutine test_diffuse(error)
    type(error_type), allocatable, intent(out) :: error

    call check_sim_fixture(error, "mv_diffuse")
  end subroutine test_diffuse

  subroutine test_diffuse_seasonal(error)
    type(error_type), allocatable, intent(out) :: error

    call check_sim_fixture(error, "uc_trend_seasonal_exact")
  end subroutine test_diffuse_seasonal

  !> Missing elements with correlated H: the eps draws there must have the
  !> conditional moments, where statsmodels returns unconditional draws.
  subroutine test_moments_missing(error)
    type(error_type), allocatable, intent(out) :: error

    call check_moments(error, "mv_missing")
  end subroutine test_moments_missing

  !> Diffuse level with a stationary AR(1), so P_star is singular. statsmodels
  !> scales the AR variate by its variance instead of its standard deviation
  !> here (unchecked potrf failure), so its draws are no reference.
  subroutine test_moments_mixed(error)
    type(error_type), allocatable, intent(out) :: error

    call check_moments(error, "uc_level_ar1_mixed")
  end subroutine test_moments_mixed

  !> Draws from random variates must have mean alphahat and variance V, and
  !> the eps draws the smoothed eps moments. Tolerances are 5 standard errors.
  subroutine check_moments(error, base)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), intent(in) :: base
    integer, parameter :: ndraw = 4000
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    type(simsmooth_result_t) :: sim
    real(dp), allocatable :: u_init(:), u_eps(:, :), u_eta(:, :)
    real(dp), allocatable :: sa(:, :), sa2(:, :), se(:, :), se2(:, :), sd(:, :)
    integer :: info, k, i, t, p, m, n, seed_size
    integer, allocatable :: seed(:)
    real(dp) :: z, zmax
    character(len=32) :: buf

    call random_seed(size=seed_size)
    seed = [(12345 + 7 * k, k=1, seed_size)]
    call random_seed(put=seed)

    fx = load_fixture("test/fixtures/"//base//".txt")
    rep = rep_from_fixture(fx)
    p = rep%k_endog; m = rep%k_states; n = rep%nobs
    call kalman_filter(rep, fres, info)
    call state_smoother(rep, fres, sres, info)

    allocate (u_init(m), u_eps(p, n), u_eta(rep%k_posdef, n))
    allocate (sa(m, n), sa2(m, n), se(p, n), se2(p, n), source=0.0_dp)
    do k = 1, ndraw
      call draw_standard_normal(u_init)
      do t = 1, n
        call draw_standard_normal(u_eps(:, t))
        call draw_standard_normal(u_eta(:, t))
      end do
      call simulation_smoother(rep, u_init, u_eps, u_eta, sim, info)
      sa = sa + sim%state
      sa2 = sa2 + sim%state**2
      se = se + sim%eps
      se2 = se2 + sim%eps**2
    end do
    sa = sa / ndraw
    sa2 = sa2 / ndraw - sa**2
    se = se / ndraw
    se2 = se2 / ndraw - se**2

    ! Standardized errors of the sample means, and of the sample variances
    ! (Var of a sample variance ~ 2 sigma^4 / N for normal draws).
    zmax = 0.0_dp
    allocate (sd(m, n))
    do t = 1, n
      do i = 1, m
        z = abs(sa(i, t) - sres%alphahat(i, t)) / sqrt(sres%V(i, i, t) / ndraw)
        zmax = max(zmax, z)
        z = abs(sa2(i, t) - sres%V(i, i, t)) / (sqrt(2.0_dp / ndraw) * sres%V(i, i, t))
        zmax = max(zmax, z)
      end do
      do i = 1, p
        z = abs(se(i, t) - sres%epshat(i, t)) / sqrt(sres%epsvar(i, i, t) / ndraw)
        zmax = max(zmax, z)
        z = abs(se2(i, t) - sres%epsvar(i, i, t)) / (sqrt(2.0_dp / ndraw) * sres%epsvar(i, i, t))
        zmax = max(zmax, z)
      end do
    end do
    ! A few hundred statistics; P(max |z| > 5) is well below 1e-3.
    write (buf, '(f8.2)') zmax
    call check(error, zmax < 5.0_dp, "largest standardized error "//trim(buf))
  end subroutine check_moments
end module test_simsmooth
