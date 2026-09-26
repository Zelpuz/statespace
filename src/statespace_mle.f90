!> Maximum likelihood estimation (DK 7.3) with L-BFGS-B.
!>
!> The optimizer minimizes -loglike / nobs over the unconstrained parameters,
!> as statsmodels does. The gradient is the analytic score (DK 7.3.3, one
!> smoother pass; see `statespace_score`) for the parameters it covers, with
!> central differences for the others (parameters in Z or T, as DK
!> recommend), or entirely by central differences when it covers none. Standard errors come from a numerical
!> Hessian of the log likelihood in the unconstrained parameters, mapped to
!> the constrained ones by the delta method (DK 7.3.6).
!>
!> `fit_many` fits independent models in parallel when the library is built
!> with OpenMP (`--flag -fopenmp`); otherwise it runs them in turn.
module statespace_mle
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use lbfgsb_module, only: setulb
  use statespace_kinds, only: dp, SS_OK, SS_ERR_NOT_PD
  use statespace_linalg, only: chol_inv
  use statespace_kinds, only: SS_ERR_UNSUPPORTED, SS_ERR_DIM
  use statespace_model, only: ssm_model_t
  use statespace_score, only: analytic_score
  use statespace_filter, only: filter_result_t
  use statespace_smoother, only: smoother_result_t
  use statespace_simsmooth, only: draw_standard_normal, psd_sqrt
  implicit none
  private

  public :: fit_options_t, fit_result_t, fit, fit_many, numerical_hessian, estimation_bias

  !> Gradient methods for `fit`.
  integer, parameter, public :: GRADIENT_AUTO = 0       !< analytic if supported, else numerical
  integer, parameter, public :: GRADIENT_NUMERICAL = 1
  integer, parameter, public :: GRADIENT_ANALYTIC = 2

  type :: fit_options_t
    integer :: maxiter = 500
    integer :: m = 10                 !< L-BFGS memory (number of corrections)
    real(dp) :: factr = 1.0e7_dp      !< L-BFGS-B relative reduction tolerance / machine eps
    real(dp) :: pgtol = 1.0e-5_dp     !< projected gradient tolerance
    logical :: compute_cov = .true.
    integer :: iprint = -1            !< L-BFGS-B output level; < 0 is silent
    integer :: gradient = GRADIENT_AUTO
  end type fit_options_t

  type :: fit_result_t
    real(dp), allocatable :: params(:)        !< constrained estimates
    real(dp), allocatable :: cov_params(:, :) !< delta method from the unconstrained Hessian
    real(dp), allocatable :: bse(:)           !< standard errors
    real(dp) :: llf = 0.0_dp
    real(dp) :: scale = 1.0_dp                !< estimated scale if concentrated out
    real(dp) :: aic = 0.0_dp, bic = 0.0_dp
    integer :: niter = 0, nfev = 0
    logical :: converged = .false.
    logical :: analytic_gradient = .false.    !< whether the analytic score was used (for at least
                                              !< one parameter; see analytic_gradient)
    character(len=60) :: message = ''
  end type fit_result_t

  !> Objective value used where the likelihood cannot be evaluated, e.g. a
  !> non-positive-definite F_t. It makes L-BFGS-B's line search back off.
  real(dp), parameter :: bad_objective = 1.0e20_dp

contains

  !> Fit `model` by maximum likelihood. On return the model holds the estimates.
  !> If the estimates are available but the Hessian is not negative definite,
  !> `info` is SS_ERR_NOT_PD and `res%cov_params`/`res%bse` are not set.
  subroutine fit(model, res, start_params, options, info)
    class(ssm_model_t), intent(inout) :: model
    type(fit_result_t), intent(out) :: res
    real(dp), intent(in), optional :: start_params(:)
    type(fit_options_t), intent(in), optional :: options
    integer, intent(out) :: info
    type(fit_options_t) :: opts
    character(len=60) :: task, csave
    logical :: lsave(4)
    integer :: isave(44)
    real(dp) :: dsave(29), f, nobs, logdet
    real(dp), allocatable :: x(:), g(:), l(:), u(:), wa(:), J(:, :)
    integer, allocatable :: nbd(:), iwa(:)
    integer :: i, k, m, n_eff, k_eff, nfd
    logical :: use_analytic

    if (present(options)) opts = options
    k = model%k_params
    m = opts%m
    nobs = real(model%rep%nobs, dp)

    if (present(start_params)) then
      x = model%untransform_params(start_params)
    else
      x = model%untransform_params(model%start_params())
    end if
    allocate (g(k), l(k), u(k), nbd(k), iwa(3 * k), wa(2 * m * k + 5 * k + 11 * m * m + 8 * m))
    l = 0.0_dp
    u = 0.0_dp
    nbd = 0

    use_analytic = opts%gradient /= GRADIENT_NUMERICAL
    task = 'START'
    do
      call setulb(k, m, x, l, u, nbd, f, g, opts%factr, opts%pgtol, wa, iwa, task, &
                  opts%iprint, csave, lsave, isave, dsave)
      if (task(1:2) == 'FG') then
        if (use_analytic) then
          call analytic_gradient(model, x, nobs, f, g, info, nfd)
          if (info == SS_ERR_UNSUPPORTED .and. opts%gradient == GRADIENT_AUTO) then
            use_analytic = .false.
          else if (info == SS_ERR_UNSUPPORTED) then
            return
          else
            res%nfev = res%nfev + 1 + nfd
          end if
        end if
        if (.not. use_analytic) then
          f = objective(model, x, nobs)
          call gradient(model, x, nobs, g)
          res%nfev = res%nfev + 1 + 2 * k
        end if
      else if (task(1:5) == 'NEW_X') then
        if (isave(30) >= opts%maxiter) then
          task = 'STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT'
          exit
        end if
      else
        exit
      end if
    end do

    res%analytic_gradient = use_analytic
    res%niter = isave(30)
    res%message = task
    res%converged = task(1:4) == 'CONV'
    res%params = model%transform_params(x)
    res%llf = model%loglike(res%params, info)
    if (info /= SS_OK) return

    n_eff = model%rep%nobs - model%rep%loglikelihood_burn
    res%scale = model%scale
    ! Diffuse initial states and a concentrated scale count as estimated
    ! parameters, as in statsmodels.
    k_eff = k + model%rep%k_diffuse()
    if (model%concentrate_scale) k_eff = k_eff + 1
    res%aic = -2.0_dp * res%llf + 2.0_dp * k_eff
    res%bic = -2.0_dp * res%llf + k_eff * log(real(n_eff, dp))

    if (opts%compute_cov) then
      ! Hessian in the unconstrained psi (DK 7.3.6), where every step is a
      ! valid parameter even for a variance near zero, then the delta method:
      ! cov(theta) = J cov(psi) J', J = d theta / d psi.
      x = model%untransform_params(res%params)
      res%cov_params = -numerical_hessian(model, x, info, transformed=.false.)
      if (info == SS_OK) call chol_inv(res%cov_params, logdet, info)
      if (info == SS_OK) then
        J = transform_jacobian(model, x)
        res%cov_params = matmul(J, matmul(res%cov_params, transpose(J)))
      end if
      if (info /= SS_OK) then
        info = SS_ERR_NOT_PD
        deallocate (res%cov_params)
      else
        allocate (res%bse(k))
        do i = 1, k
          res%bse(i) = sqrt(res%cov_params(i, i))
        end do
      end if
    end if

    ! Leave the model at the estimates (the Hessian evaluations moved it).
    call model%update(res%params)
    model%scale = res%scale
  end subroutine fit

  !> Fit each of `models` (e.g. one model per series) as `fit` does, in
  !> parallel over the models with OpenMP. Each model is only touched by its
  !> own thread. For large models, set OPENBLAS_NUM_THREADS=1 (or the
  !> equivalent) to avoid nested threading in BLAS.
  subroutine fit_many(models, res, info, options)
    class(ssm_model_t), intent(inout) :: models(:)
    type(fit_result_t), intent(out) :: res(:)
    integer, intent(out) :: info(:)
    type(fit_options_t), intent(in), optional :: options
    type(fit_options_t) :: opts
    integer :: i

    if (present(options)) opts = options
    !$omp parallel do schedule(dynamic)
    do i = 1, size(models)
      call fit(models(i), res(i), options=opts, info=info(i))
    end do
    !$omp end parallel do
  end subroutine fit_many

  !> -loglike / nobs at unconstrained parameters x.
  real(dp) function objective(model, x, nobs) result(f)
    class(ssm_model_t), intent(inout) :: model
    real(dp), intent(in) :: x(:), nobs
    integer :: info
    real(dp) :: llf

    llf = model%loglike(x, info, transformed=.false.)
    if (info /= SS_OK .or. .not. ieee_is_finite(llf)) then
      f = bad_objective
    else
      f = -llf / nobs
    end if
  end function objective

  !> Bias from estimating the parameters (DK 7.3.7, eq. 7.20). Treating
  !> psi_hat as the true value, psi_hat is approximately N(psi, Omega)
  !> (DK 7.3.6), so draw psi^(i) from N(psi_hat, Omega), i = 1..N, and
  !> estimate the bias of the smoothed state by
  !>
  !>     B = (1/N) sum_i alphahat(psi^(i)) - alphahat(psi_hat)
  !>
  !> (and likewise for V with `bias_V`). As in DK, psi is the unconstrained
  !> parameter vector (e.g. log sigma^2 / 2), with Omega the inverse Hessian
  !> there; it is obtained from `fres%cov_params` by the delta method, which
  !> is exact at the maximum. With `antithetic` (default), N/2 draws are made
  !> and each is paired with psi^(N-i+1) = 2 psi_hat - psi^(i), balancing the
  !> sample for location (DK 7.3.7); `ndraw` must then be even. Draws where
  !> the model cannot be evaluated are skipped with their pair and counted in
  !> `failed`. Uses the random number generator (`random_seed` for
  !> reproducibility). The model is left at psi_hat.
  subroutine estimation_bias(model, fres, ndraw, bias_alpha, info, bias_V, antithetic, failed)
    class(ssm_model_t), intent(inout) :: model
    type(fit_result_t), intent(in) :: fres
    integer, intent(in) :: ndraw
    real(dp), intent(out) :: bias_alpha(:, :)          !< (m, n)
    integer, intent(out) :: info
    real(dp), intent(out), optional :: bias_V(:, :, :)   !< (m, m, n)
    logical, intent(in), optional :: antithetic
    integer, intent(out), optional :: failed
    type(filter_result_t) :: filt
    type(smoother_result_t) :: s0, sp, sm
    real(dp), allocatable :: C(:, :), u(:), x0(:), J(:, :), Jinv(:, :), sumV(:, :, :)
    real(dp) :: logdet
    logical :: anti
    integer :: k, nok, nfail, stat, stat2

    anti = .true.
    if (present(antithetic)) anti = antithetic
    info = SS_ERR_DIM
    if (ndraw < 1 .or. (anti .and. mod(ndraw, 2) /= 0)) return
    k = size(fres%params)
    call model%smooth(fres%params, filt, s0, info)
    if (info /= SS_OK) return

    ! Omega on the unconstrained scale: J^-1 cov_params J^-T, J = d psi_c / d x.
    x0 = model%untransform_params(fres%params)
    J = transform_jacobian(model, x0)
    Jinv = matmul(transpose(J), J)
    call chol_inv(Jinv, logdet, info)
    if (info /= SS_OK) return
    Jinv = matmul(Jinv, transpose(J))
    C = psd_sqrt(matmul(Jinv, matmul(fres%cov_params, transpose(Jinv))))

    allocate (u(k))
    bias_alpha = 0.0_dp
    allocate (sumV(size(s0%V, 1), size(s0%V, 2), size(s0%V, 3)), source=0.0_dp)
    nok = 0
    nfail = 0
    do while (nok < ndraw)
      if (nfail > 10 * ndraw) then
        info = SS_ERR_NOT_PD        ! nearly every draw invalid
        exit
      end if
      call draw_standard_normal(u)
      call model%smooth(model%transform_params(x0 + matmul(C, u)), filt, sp, stat, &
                        transformed=.true.)
      if (anti) then
        call model%smooth(model%transform_params(x0 - matmul(C, u)), filt, sm, stat2, &
                          transformed=.true.)
      else
        stat2 = SS_OK
      end if
      if (stat /= SS_OK .or. stat2 /= SS_OK) then
        nfail = nfail + merge(2, 1, anti)
        cycle
      end if
      bias_alpha = bias_alpha + sp%alphahat
      sumV = sumV + sp%V
      nok = nok + 1
      if (anti) then
        bias_alpha = bias_alpha + sm%alphahat
        sumV = sumV + sm%V
        nok = nok + 1
      end if
    end do
    if (nok == ndraw) info = SS_OK
    if (nok > 0) then
      bias_alpha = bias_alpha / nok - s0%alphahat
      if (present(bias_V)) bias_V = sumV / nok - s0%V
    end if
    if (present(failed)) failed = nfail
    call model%smooth(fres%params, filt, s0, stat)   ! leave the model at psi_hat
  end subroutine estimation_bias

  !> Objective and gradient with the analytic score (DK 7.3.3) where it
  !> applies: g = -J' score / nobs, J = d psi / d x. Each unconstrained x_i
  !> that moves a parameter outside the score formula's scope (e.g. one in
  !> Z or T) gets a central difference of the objective instead. Returns
  !> SS_ERR_UNSUPPORTED when no parameter is covered; `nfd` counts the extra
  !> objective evaluations.
  subroutine analytic_gradient(model, x, nobs, f, g, info, nfd)
    class(ssm_model_t), intent(inout) :: model
    real(dp), intent(in), contiguous :: x(:)
    real(dp), intent(in) :: nobs
    real(dp), intent(out) :: f
    real(dp), intent(out), contiguous :: g(:)
    integer, intent(out) :: info
    integer, intent(out) :: nfd
    real(dp), allocatable :: psi(:), score(:), J(:, :)
    logical, allocatable :: analytic(:)
    real(dp) :: llf
    integer :: i

    nfd = 0
    psi = model%transform_params(x)
    allocate (score(size(psi)), analytic(size(psi)))
    call analytic_score(model, psi, score, llf, info, analytic)
    if (info == SS_ERR_UNSUPPORTED) return
    if (info /= SS_OK .or. .not. ieee_is_finite(llf)) then
      f = bad_objective
      g = 0.0_dp
      info = SS_OK
      return
    end if
    if (.not. any(analytic)) then
      info = SS_ERR_UNSUPPORTED
      return
    end if
    f = -llf / nobs
    J = transform_jacobian(model, x)
    g = -matmul(transpose(J), score) / nobs
    do i = 1, size(x)
      if (any(J(:, i) /= 0.0_dp .and. .not. analytic)) then
        g(i) = gradient_component(model, x, nobs, i)
        nfd = nfd + 2
      end if
    end do
  end subroutine analytic_gradient

  !> Central-difference gradient of `objective`.
  subroutine gradient(model, x, nobs, g)
    class(ssm_model_t), intent(inout) :: model
    real(dp), intent(in), contiguous :: x(:)
    real(dp), intent(in) :: nobs
    real(dp), intent(out), contiguous :: g(:)
    integer :: i

    do i = 1, size(x)
      g(i) = gradient_component(model, x, nobs, i)
    end do
  end subroutine gradient

  !> Central difference of `objective` in x_i; 0 if either side fails.
  real(dp) function gradient_component(model, x, nobs, i) result(gi)
    class(ssm_model_t), intent(inout) :: model
    real(dp), intent(in), contiguous :: x(:)
    real(dp), intent(in) :: nobs
    integer, intent(in) :: i
    real(dp), allocatable :: xh(:)
    real(dp) :: h, fp, fm

    xh = x
    h = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * max(abs(x(i)), 1.0_dp)
    xh(i) = x(i) + h
    fp = objective(model, xh, nobs)
    xh(i) = x(i) - h
    fm = objective(model, xh, nobs)
    if (fp >= bad_objective .or. fm >= bad_objective) then
      gi = 0.0_dp
    else
      gi = (fp - fm) / (2.0_dp * h)
    end if
  end function gradient_component

  !> Central-difference Jacobian d theta / d psi of `transform_params` at
  !> unconstrained x.
  function transform_jacobian(model, x) result(J)
    class(ssm_model_t), intent(in) :: model
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: J(:, :)
    real(dp), allocatable :: xh(:)
    real(dp) :: h
    integer :: i

    allocate (J(size(x), size(x)))
    xh = x
    do i = 1, size(x)
      h = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * max(abs(x(i)), 1.0_dp)
      xh(i) = x(i) + h
      J(:, i) = model%transform_params(xh)
      xh(i) = x(i) - h
      J(:, i) = (J(:, i) - model%transform_params(xh)) / (2.0_dp * h)
      xh(i) = x(i)
    end do
  end function transform_jacobian

  !> Central-difference Hessian of the log likelihood at `params`, which are
  !> constrained unless `transformed` is false. If any evaluation fails,
  !> `info` is set and the result is not valid.
  function numerical_hessian(model, params, info, transformed) result(hess)
    class(ssm_model_t), intent(inout) :: model
    real(dp), intent(in) :: params(:)
    integer, intent(out) :: info
    logical, intent(in), optional :: transformed
    real(dp), allocatable :: hess(:, :)
    real(dp), allocatable :: h(:)
    real(dp) :: f0
    integer :: i, j, k
    logical :: tr

    tr = .true.
    if (present(transformed)) tr = transformed

    k = size(params)
    allocate (hess(k, k))
    h = epsilon(1.0_dp)**0.25_dp * max(abs(params), 1.0_dp)
    f0 = model%loglike(params, info, transformed=tr)
    if (info /= SS_OK) return

    do i = 1, k
      hess(i, i) = (llf_at([i], [h(i)]) - 2.0_dp * f0 + llf_at([i], [-h(i)])) / h(i)**2
      if (info /= SS_OK) return
      do j = 1, i - 1
        hess(i, j) = (llf_at([i, j], [h(i), h(j)]) - llf_at([i, j], [h(i), -h(j)]) &
                      - llf_at([i, j], [-h(i), h(j)]) + llf_at([i, j], [-h(i), -h(j)])) &
                     / (4.0_dp * h(i) * h(j))
        if (info /= SS_OK) return
        hess(j, i) = hess(i, j)
      end do
    end do

  contains

    !> Log likelihood with params(idx) shifted by step.
    real(dp) function llf_at(idx, step)
      integer, intent(in) :: idx(:)
      real(dp), intent(in) :: step(:)
      real(dp), allocatable :: p(:)
      integer :: stat

      p = params
      p(idx) = p(idx) + step
      llf_at = model%loglike(p, stat, transformed=tr)
      if (stat /= SS_OK) info = stat
    end function llf_at
  end function numerical_hessian
end module statespace_mle
