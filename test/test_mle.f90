!> Maximum likelihood estimation against statsmodels fixtures.
module test_mle
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use statespace
  use fixture_io, only: fixture_t, load_fixture, rep_from_fixture
  implicit none
  private

  public :: collect_mle

  !> Local level model (DK ch. 2) with params = [sigma2_eps, sigma2_eta],
  !> set up like statsmodels' UnobservedComponents('llevel') default:
  !> approximate diffuse init with kappa = 1e6 and the first observation burned.
  type, extends(ssm_model_t) :: local_level_t
  contains
    procedure :: update => ll_update
    procedure :: start_params => ll_start_params
    procedure :: transform_params => ll_transform
    procedure :: untransform_params => ll_untransform
  end type local_level_t

  !> AR(2) with params = [phi_1, phi_2, sigma2], stationary initialization,
  !> set up like statsmodels' SARIMAX(order=(2, 0, 0)).
  type, extends(ssm_model_t) :: ar2_t
  contains
    procedure :: update => ar2_update
    procedure :: start_params => ar2_start_params
    procedure :: transform_params => ar2_transform
    procedure :: untransform_params => ar2_untransform
  end type ar2_t

  !> AR(2) with the innovation variance concentrated out: params = [phi_1, phi_2].
  type, extends(ssm_model_t) :: ar2c_t
  contains
    procedure :: update => ar2c_update
    procedure :: start_params => ar2c_start_params
    procedure :: transform_params => ar2c_transform
    procedure :: untransform_params => ar2c_untransform
  end type ar2c_t

contains

  subroutine collect_mle(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ new_unittest("transform_roundtrip", test_transform), &
                 new_unittest("loglike_at_statsmodels_params", test_loglike), &
                 new_unittest("fit_from_statsmodels_start", test_fit_sm_start), &
                 new_unittest("fit_from_default_start", test_fit_default_start), &
                 new_unittest("fit_tight_tolerances", test_fit_tight), &
                 new_unittest("fit_exact_diffuse", test_fit_exact_diffuse), &
                 new_unittest("fit_analytic_vs_numerical_gradient", &
                              test_fit_gradients), new_unittest("estimation_bias", &
        test_estimation_bias), new_unittest("concentrated_nile", &
                                            test_concentrated_nile), &
                 new_unittest("concentrated_multivariate_diffuse", &
                              test_concentrated_mv), &
                 new_unittest("concentrated_ar2_fit", test_concentrated_ar2_fit), &
                 new_unittest("stationary_transform", test_stationary_transform), &
                 new_unittest("ar2_loglike", test_ar2_loglike), &
                 new_unittest("ar2_fit", test_ar2_fit), &
                 new_unittest("fit_many_matches_fit", test_fit_many) ]
  end subroutine collect_mle

  function local_level(y) result(model)
    real(dp), intent(in) :: y(:, :)
    type(local_level_t) :: model

    model%k_params = 2
    model%rep = ssm_rep(y, 1, 1)
    model%rep%Z = 1.0_dp
    model%rep%T = 1.0_dp
    call model%rep%initialize_approximate_diffuse(1.0e6_dp)
    model%rep%loglikelihood_burn = 1
  end function local_level

  subroutine ll_update(self, params)
    class(local_level_t), intent(inout) :: self
    real(dp), intent(in) :: params(:)

    self%rep%H(1, 1, 1) = params(1)
    self%rep%Q(1, 1, 1) = params(2)
  end subroutine ll_update

  function ll_start_params(self) result(params)
    class(local_level_t), intent(in) :: self
    real(dp), allocatable :: params(:)
    real(dp) :: v

    associate (y => self%rep%y(1, :))
      v = sum((y - sum(y) / size(y))**2) / size(y)
    end associate
    params = [v / 2, v / 2]
  end function ll_start_params

  function ll_transform(self, unconstrained) result(constrained)
    class(local_level_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)

    constrained = constrain_positive(unconstrained)
  end function ll_transform

  function ll_untransform(self, constrained) result(unconstrained)
    class(local_level_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)

    unconstrained = unconstrain_positive(constrained)
  end function ll_untransform

  function ar2(y) result(model)
    real(dp), intent(in) :: y(:, :)
    type(ar2_t) :: model

    model%k_params = 3
    model%rep = ssm_rep(y, 2, 1)
    model%rep%Z(1, 1, 1) = 1.0_dp
    model%rep%T(2, 1, 1) = 1.0_dp
    call model%rep%initialize_stationary()
  end function ar2

  subroutine ar2_update(self, params)
    class(ar2_t), intent(inout) :: self
    real(dp), intent(in) :: params(:)

    self%rep%T(1, :, 1) = params(1:2)
    self%rep%Q(1, 1, 1) = params(3)
  end subroutine ar2_update

  function ar2_start_params(self) result(params)
    class(ar2_t), intent(in) :: self
    real(dp), allocatable :: params(:)

    associate (y => self%rep%y(1, :))
      params = [0.0_dp, 0.0_dp, sum(y**2) / size(y)]
    end associate
  end function ar2_start_params

  function ar2_transform(self, unconstrained) result(constrained)
    class(ar2_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)

    constrained = [constrain_stationary(unconstrained(1:2)), &
                   constrain_positive(unconstrained(3))]
  end function ar2_transform

  function ar2_untransform(self, constrained) result(unconstrained)
    class(ar2_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)

    unconstrained = [unconstrain_stationary(constrained(1:2)), &
                     unconstrain_positive(constrained(3))]
  end function ar2_untransform

  subroutine ar2c_update(self, params)
    class(ar2c_t), intent(inout) :: self
    real(dp), intent(in) :: params(:)

    self%rep%T(1, :, 1) = params(1:2)
    self%rep%Q(1, 1, 1) = 1.0_dp
  end subroutine ar2c_update

  function ar2c_start_params(self) result(params)
    class(ar2c_t), intent(in) :: self
    real(dp), allocatable :: params(:)

    params = [0.0_dp, 0.0_dp]
  end function ar2c_start_params

  function ar2c_transform(self, unconstrained) result(constrained)
    class(ar2c_t), intent(in) :: self
    real(dp), intent(in) :: unconstrained(:)
    real(dp), allocatable :: constrained(:)

    constrained = constrain_stationary(unconstrained)
  end function ar2c_transform

  function ar2c_untransform(self, constrained) result(unconstrained)
    class(ar2c_t), intent(in) :: self
    real(dp), intent(in) :: constrained(:)
    real(dp), allocatable :: unconstrained(:)

    unconstrained = unconstrain_stationary(constrained)
  end function ar2c_untransform

  function nile_model() result(model)
    type(local_level_t) :: model
    type(fixture_t) :: fx

    fx = load_fixture("test/fixtures/nile_llevel_known.txt")
    model = local_level(fx%get2('y'))
  end function nile_model

  pure real(dp) function relerr(actual, expected)
    real(dp), intent(in) :: actual(:), expected(:)

    relerr = maxval(abs(actual - expected) / max(1.0_dp, abs(expected)))
  end function relerr

  subroutine check_rel(error, actual, expected, tol, label)
    type(error_type), allocatable, intent(out) :: error
    real(dp), intent(in) :: actual(:), expected(:), tol
    character(len=*), intent(in) :: label
    character(len=32) :: buf

    write (buf, '(es10.3)') relerr(actual, expected)
    call check(error, relerr(actual, expected) <= tol, &
               label//": relative error "//trim(buf))
  end subroutine check_rel

  subroutine test_transform(error)
    type(error_type), allocatable, intent(out) :: error
    type(local_level_t) :: model
    real(dp), parameter :: p(2) = [15099.0_dp, 1469.1_dp]

    model = nile_model()
    call check_rel(error, model%transform_params(model%untransform_params(p)), p, &
                   1.0e-14_dp, "roundtrip")
  end subroutine test_transform

  subroutine test_loglike(error)
    type(error_type), allocatable, intent(out) :: error
    type(local_level_t) :: model
    type(fixture_t) :: fx
    integer :: info
    real(dp) :: llf

    model = nile_model()
    fx = load_fixture("test/fixtures/nile_llevel_mle_approx.txt")
    llf = model%loglike(fx%get1('params'), info)
    call check(error, info, SS_OK, "loglike info")
    if (allocated(error)) return
    call check_rel(error, [llf], fx%get1('llf'), 1.0e-12_dp, "llf")
  end subroutine test_loglike

  !> Compare a fit with statsmodels' tightly converged optimum. `tol` is the
  !> relative tolerance on the parameters. The likelihood is flat near the
  !> optimum, so fits with default tolerances (ours and statsmodels') stop
  !> up to ~1% away in sigma2_eta; only the tight fit is checked closely.
  subroutine check_fit(error, model, res, info, tol)
    type(error_type), allocatable, intent(out) :: error
    type(local_level_t), intent(in) :: model
    type(fit_result_t), intent(in) :: res
    integer, intent(in) :: info
    real(dp), intent(in) :: tol
    type(fixture_t) :: fx
    real(dp) :: llf_default, llf_tight
    character(len=64) :: buf

    fx = load_fixture("test/fixtures/nile_llevel_mle_approx.txt")
    llf_default = sum(fx%get1('llf'))
    llf_tight = sum(fx%get1('llf_tight'))

    call check(error, info, SS_OK, "fit info")
    if (allocated(error)) return
    call check(error, res%converged, "not converged: "//trim(res%message))
    if (allocated(error)) return
    ! At least as good as statsmodels at default tolerances, no better than
    ! the tight optimum.
    write (buf, '(2es12.4)') res%llf - llf_default, res%llf - llf_tight
    call check(error, res%llf >= llf_default - 1.0e-7_dp &
               .and. res%llf <= llf_tight + 1.0e-8_dp, "llf out of range: "//trim(buf))
    if (allocated(error)) return
    call check_rel(error, res%params, fx%get1('params_tight'), tol, "params")
    if (allocated(error)) return
    call check_rel(error, res%bse, fx%get1('bse_tight'), tol, "bse")
    if (allocated(error)) return
    call check_rel(error, [res%aic, res%bic], &
                   [fx%get1('aic_tight'), fx%get1('bic_tight')], 1.0e-6_dp, "aic/bic")
    if (allocated(error)) return
    call check_rel(error, [model%rep%H(1, 1, 1), model%rep%Q(1, 1, 1)], res%params, &
                   0.0_dp, "model left at estimates")
  end subroutine check_fit

  subroutine test_fit_sm_start(error)
    type(error_type), allocatable, intent(out) :: error
    type(local_level_t) :: model
    type(fit_result_t) :: res
    type(fixture_t) :: fx
    integer :: info

    model = nile_model()
    fx = load_fixture("test/fixtures/nile_llevel_mle_approx.txt")
    call fit(model, res, start_params=fx%get1('start_params'), info=info)
    call check_fit(error, model, res, info, 2.0e-2_dp)
  end subroutine test_fit_sm_start

  subroutine test_fit_default_start(error)
    type(error_type), allocatable, intent(out) :: error
    type(local_level_t) :: model
    type(fit_result_t) :: res
    integer :: info

    model = nile_model()
    call fit(model, res, info=info)
    call check_fit(error, model, res, info, 2.0e-2_dp)
  end subroutine test_fit_default_start

  subroutine test_fit_tight(error)
    type(error_type), allocatable, intent(out) :: error
    type(local_level_t) :: model
    type(fit_result_t) :: res
    type(fit_options_t) :: opts
    integer :: info

    model = nile_model()
    opts%factr = 10.0_dp
    opts%pgtol = 1.0e-9_dp
    call fit(model, res, options=opts, info=info)
    call check_fit(error, model, res, info, 1.0e-4_dp)
  end subroutine test_fit_tight
  !> Exact diffuse Nile local level (DK 2.10): the optimum must match
  !> statsmodels' tightly converged fit and DK's 15099 / 1469.1.
  subroutine test_fit_exact_diffuse(error)
    type(error_type), allocatable, intent(out) :: error
    type(local_level_t) :: model
    type(fit_result_t) :: res
    type(fit_options_t) :: opts
    type(fixture_t) :: fx
    integer :: info

    model = nile_model()
    call model%rep%initialize_diffuse()
    model%rep%loglikelihood_burn = 0
    fx = load_fixture("test/fixtures/nile_llevel_mle_exact.txt")
    opts%factr = 10.0_dp
    opts%pgtol = 1.0e-9_dp
    call fit(model, res, options=opts, info=info)
    call check(error, info, SS_OK, "fit info")
    if (allocated(error)) return
    call check(error, res%converged, "not converged: "//trim(res%message))
    if (allocated(error)) return
    call check(error, abs(res%llf - sum(fx%get1('llf_tight'))) < 1.0e-8_dp, "llf")
    if (allocated(error)) return
    call check_rel(error, res%params, fx%get1('params_tight'), 1.0e-4_dp, "params")
    if (allocated(error)) return
    call check_rel(error, res%params, [15099.0_dp, 1469.1_dp], 1.0e-3_dp, &
                   "DK estimates")
    if (allocated(error)) return
    call check_rel(error, res%bse, fx%get1('bse_tight'), 1.0e-4_dp, "bse")
    if (allocated(error)) return
    call check_rel(error, [res%aic, res%bic], &
                   [fx%get1('aic_tight'), fx%get1('bic_tight')], 1.0e-6_dp, "aic/bic")
  end subroutine test_fit_exact_diffuse

  !> Local level with sigma2_eps concentrated out (H = 1, Q = q), against
  !> statsmodels' filter_concentrated: exact diffuse and known initialization,
  !> with and without missing data. The concentrated log likelihood must also
  !> equal the ordinary one of the model scaled by sigma2_hat.
  subroutine test_concentrated_nile(error)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), parameter :: inits(2) = ["diffuse", "known  "]
    character(len=*), parameter :: tags(2) = ["full   ", "missing"]
    type(fixture_t) :: fx, cx
    type(ssm_rep_t) :: rep
    real(dp) :: llf, scale
    integer :: info, i, j

    fx = load_fixture("test/fixtures/nile_llevel_known.txt")
    cx = load_fixture("test/fixtures/concentrated.txt")
    do i = 1, 2
      do j = 1, 2
        if (j == 1) then
          rep = ssm_rep(fx%get2('y'), 1, 1)
        else
          rep = ssm_rep(cx%get2('y_missing'), 1, 1)
        end if
        rep%Z = 1.0_dp
        rep%T = 1.0_dp
        rep%H = 1.0_dp
        rep%Q = sum(cx%get1('q'))
        if (i == 1) then
          call rep%initialize_diffuse()
        else
          call rep%initialize_known([1000.0_dp], reshape([100.0_dp], [1, 1]))
        end if
        llf = loglike_concentrated(rep, scale, info)
        call check(error, info, SS_OK, "loglike_concentrated info")
        if (allocated(error)) return
        call check_rel(error, [llf, scale], &
                       [sum(cx%get1(trim(inits(i))//'_'//trim(tags(j))//'_llf')), &
                        sum(cx%get1(trim(inits(i))//'_'//trim(tags(j))//'_scale'))], &
                       1.0e-9_dp, trim(inits(i))//' '//trim(tags(j)))
        if (allocated(error)) return
        call rep%scale_by(scale)
        call check_rel(error, [loglike(rep, info)], [llf], 1.0e-10_dp, &
                       "scaled model llf")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_concentrated_nile

  !> p = 2 with correlated H, exact diffuse, and a missing element in the
  !> diffuse period, for the univariate and multivariate diffuse filters.
  subroutine test_concentrated_mv(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx, cx
    type(ssm_rep_t) :: rep
    real(dp) :: llf, scale
    integer :: info, k

    fx = load_fixture("test/fixtures/mv_diffuse.txt")
    cx = load_fixture("test/fixtures/concentrated.txt")
    do k = 1, 2
      rep = rep_from_fixture(fx)
      if (k == 2) rep%diffuse_method = DIFFUSE_MULTIVARIATE
      llf = loglike_concentrated(rep, scale, info)
      call check(error, info, SS_OK, "loglike_concentrated info")
      if (allocated(error)) return
      call check_rel(error, [llf, scale], [sum(cx%get1('mv_diffuse_llf')), &
                                           sum(cx%get1('mv_diffuse_scale'))], &
                     1.0e-9_dp, "mv_diffuse")
      if (allocated(error)) return
    end do
  end subroutine test_concentrated_mv

  subroutine test_concentrated_ar2_fit(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx, cx
    type(ar2c_t) :: model
    type(fit_result_t) :: res
    type(fit_options_t) :: opts
    integer :: info

    fx = load_fixture("test/fixtures/ar2.txt")
    cx = load_fixture("test/fixtures/concentrated.txt")
    model%k_params = 2
    model%concentrate_scale = .true.
    model%rep = ssm_rep(fx%get2('y'), 2, 1)
    model%rep%Z(1, 1, 1) = 1.0_dp
    model%rep%T(2, 1, 1) = 1.0_dp
    call model%rep%initialize_stationary()
    opts%factr = 10.0_dp
    opts%pgtol = 1.0e-9_dp
    call fit(model, res, options=opts, info=info)
    call check(error, info, SS_OK, "fit info")
    if (allocated(error)) return
    call check(error, abs(res%llf - sum(cx%get1('ar2_llf'))) < 1.0e-8_dp, "llf")
    if (allocated(error)) return
    call check_rel(error, res%params, cx%get1('ar2_params'), 1.0e-4_dp, "params")
    if (allocated(error)) return
    call check_rel(error, [res%scale], cx%get1('ar2_scale'), 1.0e-4_dp, "scale")
    if (allocated(error)) return
    call check_rel(error, [res%aic, res%bic], &
                   [sum(cx%get1('ar2_aic')), sum(cx%get1('ar2_bic'))], 1.0e-6_dp, &
                   "aic/bic")
  end subroutine test_concentrated_ar2_fit

  !> The analytic score (DK 7.3.3) must reach the same optimum as numerical
  !> gradients, with far fewer likelihood evaluations.
  subroutine test_fit_gradients(error)
    type(error_type), allocatable, intent(out) :: error
    type(local_level_t) :: model
    type(fit_result_t) :: ra, rn
    type(fit_options_t) :: opts
    integer :: info
    character(len=64) :: buf

    model = nile_model()
    call model%rep%initialize_diffuse()
    model%rep%loglikelihood_burn = 0
    opts%factr = 10.0_dp
    opts%pgtol = 1.0e-9_dp
    call fit(model, ra, options=opts, info=info)
    call check(error, info, SS_OK, "analytic fit info")
    if (allocated(error)) return
    call check(error, ra%analytic_gradient, "analytic score used")
    if (allocated(error)) return
    opts%gradient = GRADIENT_NUMERICAL
    call fit(model, rn, options=opts, info=info)
    call check(error, .not. rn%analytic_gradient, "numerical gradient used")
    if (allocated(error)) return
    call check_rel(error, ra%params, rn%params, 1.0e-5_dp, "same optimum")
    if (allocated(error)) return
    write (buf, '(i0, " vs ", i0)') ra%nfev, rn%nfev
    call check(error, ra%nfev < rn%nfev, "likelihood evaluations "//trim(buf))
  end subroutine test_fit_gradients

  !> DK 7.3.7 bias estimate (7.20): zero with zero parameter covariance,
  !> reproducible, no failed draws on the unconstrained scale, and simulation
  !> noise (antithetic draws) small next to the estimated bias.
  subroutine test_estimation_bias(error)
    type(error_type), allocatable, intent(out) :: error
    type(local_level_t) :: model
    type(fit_result_t) :: res, fixed
    type(fit_options_t) :: opts
    type(filter_result_t) :: fres
    type(smoother_result_t) :: sres
    real(dp), allocatable :: bias(:, :), bias2(:, :), biasV(:, :, :), rel(:)
    integer, allocatable :: seed(:)
    integer :: info, nseed, failed, k
    character(len=64) :: buf

    model = nile_model()
    call model%rep%initialize_diffuse()
    model%rep%loglikelihood_burn = 0
    opts%factr = 10.0_dp
    opts%pgtol = 1.0e-9_dp
    call fit(model, res, options=opts, info=info)
    call model%smooth(res%params, fres, sres, info)
    allocate (bias(1, model%rep%nobs), bias2(1, model%rep%nobs), &
              biasV(1, 1, model%rep%nobs))

    fixed = res
    fixed%cov_params = 0.0_dp * res%cov_params
    call estimation_bias(model, fixed, 4, bias, info, biasV)
    call check(error, info, SS_OK, "info")
    if (allocated(error)) return
    call check(error, maxval(abs(bias)) < 1.0e-9_dp * maxval(abs(sres%alphahat)) &
               .and. maxval(abs(biasV)) < 1.0e-9_dp * maxval(sres%V), &
               "zero covariance gives zero bias")
    if (allocated(error)) return
    call estimation_bias(model, res, 5, bias, info)
    call check(error, info, SS_ERR_DIM, "antithetic draws need an even N")
    if (allocated(error)) return

    call random_seed(size=nseed)
    seed = [(97 + 13 * k, k=1, nseed)]
    call random_seed(put=seed)
    call estimation_bias(model, res, 500, bias, info, biasV, failed=failed)
    call check(error, info, SS_OK, "info")
    if (allocated(error)) return
    call random_seed(put=seed)
    call estimation_bias(model, res, 500, bias2, info)
    call check_rel(error, pack(bias2, .true.), pack(bias, .true.), 0.0_dp, &
                   "reproducible")
    if (allocated(error)) return
    call check(error, failed == 0, "no failed draws on the unconstrained scale")
    if (allocated(error)) return

    ! Simulation noise, from a second seed, must be small next to the bias.
    seed = [(31 + 7 * k, k=1, nseed)]
    call random_seed(put=seed)
    call estimation_bias(model, res, 500, bias2, info)
    rel = abs(bias(1, :)) / sqrt(sres%V(1, 1, :))
    write (buf, '(3f9.5)') maxval(rel), sum(rel) / size(rel), &
      maxval(abs(bias2(1, :) - bias(1, :)) / sqrt(sres%V(1, 1, :)))
    call check(error, maxval(abs(bias2(1, :) - bias(1, :)) / sqrt(sres%V(1, 1, :))) &
               < 0.2_dp * maxval(rel), "simulation noise vs bias: "//trim(buf))
  end subroutine test_estimation_bias

  subroutine test_stationary_transform(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    real(dp), allocatable :: x(:)

    fx = load_fixture("test/fixtures/ar2.txt")
    x = fx%get1('transform_in')
    call check_rel(error, constrain_stationary(x), fx%get1('transform_out'), &
                   1.0e-14_dp, "matches statsmodels")
    if (allocated(error)) return
    call check_rel(error, unconstrain_stationary(constrain_stationary(x)), x, &
                   1.0e-12_dp, "roundtrip")
  end subroutine test_stationary_transform

  subroutine test_ar2_loglike(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ar2_t) :: model
    type(filter_result_t) :: fres
    integer :: info

    fx = load_fixture("test/fixtures/ar2.txt")
    model = ar2(fx%get2('y'))
    call model%filter(fx%get1('params'), fres, info)
    call check(error, info, SS_OK, "filter info")
    if (allocated(error)) return
    call check_rel(error, fres%llf_obs, fx%get1('llf_obs'), 1.0e-10_dp, "llf_obs")
    if (allocated(error)) return
    call check_rel(error, [model%loglike(fx%get1('params'), info)], fx%get1('llf'), &
                   1.0e-12_dp, "llf")
  end subroutine test_ar2_loglike

  subroutine test_ar2_fit(error)
    type(error_type), allocatable, intent(out) :: error
    type(fixture_t) :: fx
    type(ar2_t) :: model
    type(fit_result_t) :: res
    type(fit_options_t) :: opts
    integer :: info

    fx = load_fixture("test/fixtures/ar2.txt")
    model = ar2(fx%get2('y'))
    opts%factr = 10.0_dp
    opts%pgtol = 1.0e-9_dp
    call fit(model, res, options=opts, info=info)
    call check(error, info, SS_OK, "fit info")
    if (allocated(error)) return
    call check(error, res%analytic_gradient, &
               "hybrid gradient: sigma2 analytic, AR (in T) numerical")
    if (allocated(error)) return
    call check(error, res%converged, "not converged: "//trim(res%message))
    if (allocated(error)) return
    call check(error, abs(res%llf - sum(fx%get1('llf_tight'))) < 1.0e-8_dp, "llf")
    if (allocated(error)) return
    call check_rel(error, res%params, fx%get1('params_tight'), 1.0e-4_dp, "params")
    if (allocated(error)) return
    call check_rel(error, res%bse, fx%get1('bse_tight'), 1.0e-3_dp, "bse")
  end subroutine test_ar2_fit

  !> fit_many (parallel with OpenMP) gives the same results as fit.
  subroutine test_fit_many(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: ns = 8, n = 80
    type(component_holder_t) :: comps(2)
    type(structural_model_t) :: models(ns), one
    type(fit_result_t) :: res(ns), r1
    real(dp) :: y(1, n, ns), e(n), u(n)
    integer :: info(ns), i, t, stat

    allocate (irregular_t :: comps(1)%c)
    allocate (level_t :: comps(2)%c)
    do i = 1, ns
      call draw_standard_normal(e)
      call draw_standard_normal(u)
      y(1, 1, i) = e(1)
      do t = 2, n
        y(1, t, i) = y(1, t - 1, i) + 0.3_dp * e(t)
      end do
      y(1, :, i) = y(1, :, i) + u
      models(i) = structural_model(y(:, :, i), comps, stat)
    end do
    call fit_many(models, res, info)
    do i = 1, ns
      one = structural_model(y(:, :, i), comps, stat)
      call fit(one, r1, info=stat)
      call check(error, info(i), stat, "info")
      if (allocated(error)) return
      call check(error, all(res(i)%params == r1%params) .and. res(i)%llf == r1%llf, &
                 "identical to fit")
      if (allocated(error)) return
    end do
  end subroutine test_fit_many
end module test_mle
