!> Analytic score (DK 7.3.3) against numerical derivatives of the log
!> likelihood, and the EM algorithm (DK 7.3.4).
module test_score
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use statespace
  use fixture_io, only: fixture_t, load_fixture, rep_from_fixture
  implicit none
  private

  public :: collect_score

  integer, parameter :: TO_H = 1, TO_Q = 2, SCALE_H = 3, TO_T = 4, TO_R = 5

  !> Each parameter sets H(row, col) = H(col, row), Q(row, col) = Q(col, row),
  !> T(row, col), or scales every slice of a base H (for time-varying H).
  !> With `cholesky`, the H and Q entries (which must cover lower triangles)
  !> are optimized through Cholesky factors, keeping H and Q positive definite.
  type, extends(ssm_model_t) :: varmodel_t
    integer, allocatable :: target(:), row(:), col(:)
    real(dp), allocatable :: H0(:, :, :)
    logical :: cholesky = .false.
  contains
    procedure :: update => var_update
    procedure :: start_params => var_start_params
    procedure :: transform_params => var_transform
    procedure :: untransform_params => var_untransform
    procedure :: map_cholesky
  end type varmodel_t

contains

  subroutine collect_score(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("nile_exact_diffuse", test_nile), &
                new_unittest("correlated_H_diffuse_missing", test_mv_diffuse), &
                new_unittest("correlated_H_missing", test_mv_missing), &
                new_unittest("time_varying_H", test_mv_timevarying), &
                new_unittest("trend_seasonal_missing", test_seasonal), &
                new_unittest("mixed_init_stationary_block", test_mixed), &
                new_unittest("concentrated_scale", test_concentrated), &
                new_unittest("parameter_in_R", test_param_in_R), &
                new_unittest("unsupported_detected", test_unsupported), &
                new_unittest("em_nile", test_em_nile), &
                new_unittest("em_full_H_Q_missing", test_em_mv_missing), &
                new_unittest("em_diagonal_seasonal_missing", test_em_seasonal), &
                new_unittest("em_rejects_stationary_init", test_em_stationary), &
                new_unittest("hybrid_score_mask", test_hybrid_mask), &
                new_unittest("hybrid_gradient_fit", test_hybrid_fit) &
                ]
  end subroutine collect_score

  subroutine var_update(self, params)
    class(varmodel_t), intent(inout) :: self
    real(dp), intent(in) :: params(:)
    integer :: i

    do i = 1, size(params)
      select case (self%target(i))
      case (TO_H)
        self%rep%H(self%row(i), self%col(i), 1) = params(i)
        self%rep%H(self%col(i), self%row(i), 1) = params(i)
      case (TO_Q)
        self%rep%Q(self%row(i), self%col(i), 1) = params(i)
        self%rep%Q(self%col(i), self%row(i), 1) = params(i)
      case (SCALE_H)
        self%rep%H = params(i) * self%H0
      case (TO_T)
        self%rep%T(self%row(i), self%col(i), 1) = params(i)
      case (TO_R)
        self%rep%R(self%row(i), self%col(i), 1) = params(i)
      end select
    end do
  end subroutine var_update

  function var_start_params(self) result(params)
    class(varmodel_t), intent(in) :: self
    real(dp), allocatable :: params(:)

    allocate (params(self%k_params), source=1.0_dp)
  end function var_start_params

  !> Unconstrained Cholesky entries -> H and Q entries (H = L L', Q = L L').
  function var_transform(self, unconstrained) result(constrained)
    class(varmodel_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)

    constrained = self%map_cholesky(unconstrained, .true.)
  end function var_transform

  function var_untransform(self, constrained) result(unconstrained)
    class(varmodel_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)

    unconstrained = self%map_cholesky(constrained, .false.)
  end function var_untransform

  !> forward: Cholesky entries -> matrix entries; otherwise the inverse. Only
  !> active with `cholesky`; the identity otherwise.
  function map_cholesky(self, x, forward) result(y)
    class(varmodel_t), intent(in) :: self
    real(dp), intent(in) :: x(:)
    logical, intent(in) :: forward
    real(dp), allocatable :: y(:), A(:, :), L(:, :)
    integer :: target, i, n

    y = x
    if (.not. self%cholesky) return
    do target = TO_H, TO_Q
      n = merge(self%rep%k_endog, self%rep%k_posdef, target == TO_H)
      allocate (A(n, n), source=0.0_dp)
      do i = 1, size(x)
        if (self%target(i) == target) A(self%row(i), self%col(i)) = x(i)
      end do
      if (forward) then
        L = matmul(A, transpose(A))
      else
        A = A + transpose(A) - diag_of(A)
        L = lower_cholesky(A)
      end if
      do i = 1, size(x)
        if (self%target(i) == target) y(i) = L(self%row(i), self%col(i))
      end do
      deallocate (A)
    end do
  end function map_cholesky

  pure function diag_of(A) result(D)
    real(dp), intent(in) :: A(:, :)
    real(dp) :: D(size(A, 1), size(A, 2))
    integer :: i

    D = 0.0_dp
    do i = 1, size(A, 1)
      D(i, i) = A(i, i)
    end do
  end function diag_of

  function lower_cholesky(A) result(L)
    real(dp), intent(in) :: A(:, :)
    real(dp), allocatable :: L(:, :)
    integer :: i, j, n

    n = size(A, 1)
    allocate (L(n, n), source=0.0_dp)
    do j = 1, n
      L(j, j) = sqrt(A(j, j) - sum(L(j, 1:j - 1)**2))
      do i = j + 1, n
        L(i, j) = (A(i, j) - sum(L(i, 1:j - 1) * L(j, 1:j - 1))) / L(j, j)
      end do
    end do
  end function lower_cholesky

  !> A varmodel on a fixture's model; `spec` rows are (target, row, col).
  function varmodel(path, spec) result(mod)
    character(len=*), intent(in) :: path
    integer, intent(in) :: spec(:, :)
    type(varmodel_t) :: mod
    type(fixture_t) :: fx

    fx = load_fixture(path)
    mod%rep = rep_from_fixture(fx)
    mod%H0 = mod%rep%H
    mod%k_params = size(spec, 1)
    mod%target = spec(:, 1)
    mod%row = spec(:, 2)
    mod%col = spec(:, 3)
  end function varmodel

  !> Current parameter values of a varmodel, read back from its matrices.
  function current_params(mod) result(params)
    type(varmodel_t), intent(in) :: mod
    real(dp), allocatable :: params(:)
    integer :: i

    allocate (params(mod%k_params))
    do i = 1, mod%k_params
      select case (mod%target(i))
      case (TO_H)
        params(i) = mod%rep%H(mod%row(i), mod%col(i), 1)
      case (TO_Q)
        params(i) = mod%rep%Q(mod%row(i), mod%col(i), 1)
      case (SCALE_H)
        params(i) = 1.0_dp
      case (TO_T)
        params(i) = mod%rep%T(mod%row(i), mod%col(i), 1)
      case (TO_R)
        params(i) = mod%rep%R(mod%row(i), mod%col(i), 1)
      end select
    end do
  end function current_params

  !> Analytic score against central differences of the log likelihood, at
  !> the fixture's parameters moved by `shift` (so the score is not ~0).
  subroutine check_score(error, mod, shift)
    type(error_type), allocatable, intent(out) :: error
    class(ssm_model_t), intent(inout) :: mod
    real(dp), intent(in) :: shift(:)
    real(dp), allocatable :: p0(:), p(:), score(:), fd(:)
    real(dp) :: llf, h, err
    integer :: info, i
    character(len=64) :: buf

    select type (mod)
    type is (varmodel_t)
      p0 = current_params(mod) * shift
    class default
      p0 = shift
    end select
    allocate (score(size(p0)), fd(size(p0)))
    call analytic_score(mod, p0, score, llf, info)
    call check(error, info, SS_OK, "analytic_score info")
    if (allocated(error)) return
    call check(error, abs(llf - mod%loglike(p0, info)) <= 1.0e-8_dp * abs(llf), "llf")
    if (allocated(error)) return
    do i = 1, size(p0)
      h = 1.0e-4_dp * max(abs(p0(i)), 1.0e-2_dp)
      p = p0
      p(i) = p0(i) + h
      fd(i) = mod%loglike(p, info)
      p(i) = p0(i) - h
      fd(i) = (fd(i) - mod%loglike(p, info)) / (2.0_dp * h)
    end do
    err = maxval(abs(score - fd) / (abs(fd) + 1.0e-6_dp * maxval(abs(fd))))
    write (buf, '(es10.3)') err
    call check(error, err < 1.0e-5_dp, "score vs finite differences: relative error "//trim(buf))
  end subroutine check_score

  subroutine test_nile(error)
    type(error_type), allocatable, intent(out) :: error
    type(varmodel_t) :: mod

    mod = varmodel("test/fixtures/nile_llevel_exact.txt", reshape([TO_H, TO_Q, 1, 1, 1, 1], [2, 3]))
    call check_score(error, mod, [1.2_dp, 0.7_dp])
  end subroutine test_nile

  subroutine test_mv_diffuse(error)
    type(error_type), allocatable, intent(out) :: error
    type(varmodel_t) :: mod

    mod = varmodel("test/fixtures/mv_diffuse.txt", reshape( &
                   [TO_H, TO_H, TO_H, TO_Q, TO_Q, TO_Q, &
                    1, 2, 2, 1, 2, 2, &
                    1, 1, 2, 1, 1, 2], [6, 3]))
    call check_score(error, mod, [1.2_dp, 0.8_dp, 1.1_dp, 0.9_dp, 1.3_dp, 1.2_dp])
  end subroutine test_mv_diffuse

  subroutine test_mv_missing(error)
    type(error_type), allocatable, intent(out) :: error
    type(varmodel_t) :: mod

    mod = varmodel("test/fixtures/mv_missing.txt", reshape( &
                   [TO_H, TO_H, TO_H, TO_H, TO_Q, TO_Q, &
                    1, 2, 3, 3, 1, 2, &
                    1, 2, 3, 1, 1, 2], [6, 3]))
    call check_score(error, mod, [1.2_dp, 0.8_dp, 1.1_dp, 0.7_dp, 1.3_dp, 0.9_dp])
  end subroutine test_mv_missing

  subroutine test_mv_timevarying(error)
    type(error_type), allocatable, intent(out) :: error
    type(varmodel_t) :: mod

    mod = varmodel("test/fixtures/mv_timevarying.txt", reshape( &
                   [SCALE_H, TO_Q, TO_Q, TO_Q, &
                    1, 1, 2, 2, &
                    1, 1, 1, 2], [4, 3]))
    call check_score(error, mod, [1.3_dp, 0.8_dp, 1.2_dp, 1.1_dp])
  end subroutine test_mv_timevarying

  subroutine test_seasonal(error)
    type(error_type), allocatable, intent(out) :: error
    type(varmodel_t) :: mod

    mod = varmodel("test/fixtures/uc_trend_seasonal_exact.txt", reshape( &
                   [TO_H, TO_Q, TO_Q, TO_Q, &
                    1, 1, 2, 3, &
                    1, 1, 2, 3], [4, 3]))
    call check_score(error, mod, [1.2_dp, 0.8_dp, 1.5_dp, 0.7_dp])
  end subroutine test_seasonal

  !> Diffuse level with a stationary AR(1): the AR variance moves P_star.
  subroutine test_mixed(error)
    type(error_type), allocatable, intent(out) :: error
    type(varmodel_t) :: mod

    mod = varmodel("test/fixtures/uc_level_ar1_mixed.txt", reshape( &
                   [TO_H, TO_Q, TO_Q, &
                    1, 1, 2, &
                    1, 1, 2], [3, 3]))
    call check_score(error, mod, [1.2_dp, 0.8_dp, 1.4_dp])
  end subroutine test_mixed

  !> DK (7.16): a parameter in R, here the loading of the level shock on the
  !> slope of a local linear trend (exact diffuse).
  subroutine test_param_in_R(error)
    type(error_type), allocatable, intent(out) :: error
    type(varmodel_t) :: mod

    mod = varmodel("test/fixtures/nile_lltrend_exact.txt", reshape( &
                   [TO_H, TO_Q, TO_Q, TO_R, &
                    1, 1, 2, 2, &
                    1, 1, 2, 1], [4, 3]))
    mod%rep%R(2, 1, 1) = 0.3_dp
    call check_score(error, mod, [1.2_dp, 0.8_dp, 1.5_dp, 0.7_dp])
  end subroutine test_param_in_R

  !> Nile local level with sigma2_eps concentrated out: H = 1, param q = Q.
  subroutine test_concentrated(error)
    type(error_type), allocatable, intent(out) :: error
    type(varmodel_t) :: mod

    mod = varmodel("test/fixtures/nile_llevel_exact.txt", reshape([TO_Q, 1, 1], [1, 3]))
    mod%rep%H = 1.0_dp
    mod%rep%Q = 0.2_dp
    mod%concentrate_scale = .true.
    call check_score(error, mod, [0.5_dp])
  end subroutine test_concentrated

  !> A parameter moving T (AR coefficient) and a burn-in are outside the formula.
  subroutine test_unsupported(error)
    type(error_type), allocatable, intent(out) :: error
    type(varmodel_t) :: mod
    real(dp) :: score(2), llf
    integer :: info

    mod = varmodel("test/fixtures/nile_llevel_exact.txt", reshape([TO_H, TO_Q, 1, 1, 1, 1], [2, 3]))
    mod%rep%loglikelihood_burn = 1
    call analytic_score(mod, [15000.0_dp, 1500.0_dp], score, llf, info)
    call check(error, info, SS_ERR_UNSUPPORTED, "burn-in")
    if (allocated(error)) return

    mod%rep%loglikelihood_burn = 0
    mod%target(2) = TO_T
    call analytic_score(mod, [15000.0_dp, 0.9_dp], score, llf, info)
    call check(error, info, SS_ERR_UNSUPPORTED, "parameter in T")
  end subroutine test_unsupported
  !> Nile local level, exact diffuse, EM from var(y)/2 for both variances:
  !> monotone, and converging to the MLE (DK: 15099 and 1469.1).
  subroutine test_em_nile(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx, mx
    type(ssm_rep_t) :: rep
    real(dp), allocatable :: path(:)
    real(dp) :: llf, v
    integer :: info, niter
    character(len=64) :: buf

    fx = load_fixture("test/fixtures/nile_llevel_exact.txt")
    mx = load_fixture("test/fixtures/nile_llevel_mle_exact.txt")
    rep = rep_from_fixture(fx)
    associate (y => rep%y(1, :))
      v = sum((y - sum(y) / size(y))**2) / size(y)
    end associate
    rep%H = v / 2
    rep%Q = v / 2
    call em_variances(rep, 20000, 1.0e-10_dp, llf, niter, info, llf_path=path)
    call check(error, info, SS_OK, "em info")
    if (allocated(error)) return
    call check(error, all(path(2:) >= path(:niter - 1) - 1.0e-9_dp), "log likelihood decreased")
    if (allocated(error)) return
    write (buf, '(2f12.3, i7)') rep%H(1, 1, 1), rep%Q(1, 1, 1), niter
    call check(error, abs(rep%H(1, 1, 1) / 15098.518_dp - 1) < 1.0e-3_dp .and. &
               abs(rep%Q(1, 1, 1) / 1469.176_dp - 1) < 1.0e-3_dp, "EM estimates "//trim(buf))
    if (allocated(error)) return
    call check(error, abs(llf - sum(mx%get1('llf_tight'))) < 1.0e-6_dp, "EM log likelihood")
  end subroutine test_em_nile

  !> EM is monotone, and the MLE is an EM fixed point: fit by maximum
  !> likelihood (analytic score), then one EM step must leave the estimates
  !> essentially unchanged. (EM itself converges too slowly on these flat
  !> likelihoods to reach the MLE in a test.)
  subroutine check_em_fixed_point(error, mod, diagonal)
    type(error_type), allocatable, intent(out) :: error
    type(varmodel_t), intent(inout) :: mod
    logical, intent(in) :: diagonal
    type(varmodel_t) :: em
    type(fit_result_t) :: res
    type(fit_options_t) :: opts
    real(dp), allocatable :: path(:), p0(:)
    real(dp) :: llf, change
    integer :: info, niter
    character(len=64) :: buf

    em = mod
    call em_variances(em%rep, 200, 0.0_dp, llf, niter, info, diagonal_H=diagonal, &
                      diagonal_Q=diagonal, llf_path=path)
    call check(error, info, SS_OK, "em info")
    if (allocated(error)) return
    call check(error, all(path(2:) >= path(:niter - 1) - 1.0e-9_dp), "log likelihood decreased")
    if (allocated(error)) return

    opts%factr = 10.0_dp
    opts%pgtol = 1.0e-10_dp
    opts%gradient = GRADIENT_ANALYTIC
    mod%cholesky = .true.
    call fit(mod, res, start_params=current_params(mod), options=opts, info=info)
    call check(error, res%converged, "MLE: "//trim(res%message))
    if (allocated(error)) return
    p0 = current_params(mod)
    call em_step(mod%rep, llf, info, diagonal_H=diagonal, diagonal_Q=diagonal)
    call check(error, info, SS_OK, "em_step info")
    if (allocated(error)) return
    change = maxval(abs(current_params(mod) - p0) / max(abs(p0), 1.0e-8_dp))
    write (buf, '(es10.2)') change
    call check(error, change < 1.0e-4_dp, "EM step from the MLE moved it by "//trim(buf))
  end subroutine check_em_fixed_point

  subroutine test_em_mv_missing(error)
    type(error_type), allocatable, intent(out) :: error
    type(varmodel_t) :: mod

    mod = varmodel("test/fixtures/mv_missing.txt", reshape( &
                   [TO_H, TO_H, TO_H, TO_H, TO_H, TO_H, TO_Q, TO_Q, TO_Q, &
                    1, 2, 3, 2, 3, 3, 1, 2, 2, &
                    1, 1, 1, 2, 2, 3, 1, 1, 2], [9, 3]))
    call check_em_fixed_point(error, mod, .false.)
  end subroutine test_em_mv_missing

  subroutine test_em_seasonal(error)
    type(error_type), allocatable, intent(out) :: error
    type(varmodel_t) :: mod

    mod = varmodel("test/fixtures/uc_trend_seasonal_exact.txt", reshape( &
                   [TO_H, TO_Q, TO_Q, TO_Q, &
                    1, 1, 2, 3, &
                    1, 1, 2, 3], [4, 3]))
    call check_em_fixed_point(error, mod, .true.)
  end subroutine test_em_seasonal

  subroutine test_em_stationary(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ssm_rep_t) :: rep
    real(dp) :: llf
    integer :: info

    fx = load_fixture("test/fixtures/uc_level_ar1_mixed.txt")
    rep = rep_from_fixture(fx)
    call em_step(rep, llf, info)
    call check(error, info, SS_ERR_UNSUPPORTED, "stationary block depends on Q")
  end subroutine test_em_stationary

  !> Irregular + level + damped cycle on simulated data.
  function cycle_model() result(mod)
    type(structural_model_t) :: mod
    type(component_holder_t) :: comps(3)
    real(dp) :: y(1, 200), e(200), u(200)
    integer :: t, info

    call draw_standard_normal(e)
    call draw_standard_normal(u)
    y(1, 1) = 0.0_dp
    do t = 2, 200
      y(1, t) = y(1, t - 1) + 0.2_dp * e(t)
    end do
    y(1, :) = y(1, :) + 2.0_dp * sin(2 * acos(-1.0_dp) * [(t, t=1, 200)] / 20.0_dp) + u
    allocate (irregular_t :: comps(1)%c)
    allocate (level_t :: comps(2)%c)
    comps(3)%c = cycle_t(period_min=5.0_dp, period_max=50.0_dp)
    mod = structural_model(y, comps, info)
  end function cycle_model

  !> With `analytic`, variance parameters are scored and the cycle's
  !> frequency and damping (in T) are marked as not covered.
  subroutine test_hybrid_mask(error)
    type(error_type), allocatable, intent(out) :: error
    type(structural_model_t) :: mod
    real(dp) :: p0(5), score(5), llf, h, fd, pp(5), pm(5)
    logical :: analytic(5)
    integer :: info, i

    mod = cycle_model()
    p0 = [1.0_dp, 0.04_dp, 0.1_dp, 2 * acos(-1.0_dp) / 20, 0.9_dp]
    call analytic_score(mod, p0, score, llf, info)
    call check(error, info, SS_ERR_UNSUPPORTED, "without the mask")
    if (allocated(error)) return
    call analytic_score(mod, p0, score, llf, info, analytic)
    call check(error, info, SS_OK, "with the mask")
    if (allocated(error)) return
    call check(error, all(analytic .eqv. [.true., .true., .true., .false., .false.]), "mask")
    if (allocated(error)) return
    do i = 1, 3
      h = 1.0e-5_dp * p0(i)
      pp = p0; pp(i) = p0(i) + h
      pm = p0; pm(i) = p0(i) - h
      fd = (mod%loglike(pp, info) - mod%loglike(pm, info)) / (2 * h)
      call check(error, abs(score(i) - fd) <= 1.0e-5_dp * max(1.0_dp, abs(fd)), "score vs FD")
      if (allocated(error)) return
    end do
  end subroutine test_hybrid_mask

  !> A fit with the hybrid gradient reaches the same optimum as a fully
  !> numerical one.
  subroutine test_hybrid_fit(error)
    type(error_type), allocatable, intent(out) :: error
    type(structural_model_t) :: m1, m2
    type(fit_result_t) :: r1, r2
    type(fit_options_t) :: opts
    integer :: info

    m1 = cycle_model()
    m2 = m1
    opts%factr = 10.0_dp
    opts%pgtol = 1.0e-9_dp
    opts%compute_cov = .false.
    call fit(m1, r1, options=opts, info=info)
    call check(error, info, SS_OK, "hybrid fit")
    if (allocated(error)) return
    call check(error, r1%analytic_gradient, "hybrid gradient used")
    if (allocated(error)) return
    opts%gradient = GRADIENT_NUMERICAL
    call fit(m2, r2, options=opts, info=info)
    call check(error, abs(r1%llf - r2%llf) < 1.0e-6_dp, "same log likelihood")
    if (allocated(error)) return
    call check(error, maxval(abs(r1%params - r2%params) / max(abs(r2%params), 1.0e-3_dp)) < 1.0e-3_dp, &
               "same estimates")
  end subroutine test_hybrid_fit
end module test_score
