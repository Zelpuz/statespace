!> Goodness of fit (DK 7.4) and diagnostic checking (DK 7.5).
!>
!> Residuals:
!>   * `standardized_residuals`: one-step errors e_t = L_t^-1 v_t over the
!>     observed elements, L_t the lower Cholesky factor of F_t (as in
!>     statsmodels); NaN for missing elements and diffuse periods.
!>   * `auxiliary_residuals`: smoothed disturbances standardized by their own
!>     variances, Var(epshat_t) = H_t - Var(eps_t|Y), Var(etahat_t) =
!>     Q_t - Var(eta_t|Y) (DK 2.12, 7.5), for detecting outliers and breaks.
!> Tests on a residual series (NaNs skipped), as in statsmodels:
!>   * `ljung_box`: serial correlation, Q(k) = n(n+2) sum_j rho_j^2 / (n - j).
!>   * `jarque_bera`: normality from skewness and kurtosis.
!>   * `breakvar_test`: heteroskedasticity H(h) = sum of the last h squared
!>     residuals over the sum of the first h.
!> `diagnostic_start` gives the first period to use, after the burn-in and
!> the diffuse period.
module statespace_diagnostics
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use statespace_kinds, only: dp, SS_OK, SS_ERR_DIM
  use statespace_rep, only: ssm_rep_t, tidx
  use statespace_filter, only: filter_result_t, kalman_filter, steady_state
  use statespace_smoothing, only: conventional_gain
  use statespace_linalg, only: ldl_psd
  use statespace_smoother, only: smoother_result_t
  use statespace_simsmooth, only: psd_sqrt
  use statespace_special, only: chi2_sf, f_cdf
  implicit none
  private

  public :: standardized_residuals, auxiliary_residuals, diagnostic_start
  public :: ljung_box, jarque_bera, breakvar_test, r2_diffuse, prediction_error_variance
  public :: auxiliary_residuals_vector, de_jong_penzer, least_squares_residuals

contains

  !> First period of the residuals to use in tests: after the log likelihood
  !> burn-in and the diffuse period.
  integer function diagnostic_start(rep, fres)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres

    diagnostic_start = max(rep%loglikelihood_burn, fres%nobs_diffuse) + 1
  end function diagnostic_start

  !> e_t = L_t^-1 v_t over the observed elements (F_oo = L L').
  subroutine standardized_residuals(fres, e)
    type(filter_result_t), intent(in) :: fres
    real(dp), intent(out), contiguous :: e(:, :)        !< (p, n)
    real(dp), allocatable :: L(:, :)
    integer, allocatable :: idx(:)
    real(dp) :: x
    integer :: t, i, j, n_o

    e = ieee_value(1.0_dp, ieee_quiet_nan)
    do t = fres%nobs_diffuse + 1, fres%nobs
      idx = pack([(i, i=1, fres%k_endog)], .not. ieee_is_nan(fres%v(:, t)))
      n_o = size(idx)
      if (n_o == 0) cycle
      L = psd_sqrt(fres%F(idx, idx, t))
      ! Forward substitution L e = v
      do i = 1, n_o
        x = fres%v(idx(i), t)
        do j = 1, i - 1
          x = x - L(i, j) * e(idx(j), t)
        end do
        e(idx(i), t) = x / L(i, i)
      end do
    end do
  end subroutine standardized_residuals

  !> Auxiliary residuals: epshat_t,i / sqrt(Var(epshat_t)_ii) and
  !> etahat_t,i / sqrt(Var(etahat_t)_ii); NaN where that variance is zero.
  subroutine auxiliary_residuals(rep, sres, eps_std, eta_std)
    type(ssm_rep_t), intent(in) :: rep
    type(smoother_result_t), intent(in) :: sres
    real(dp), intent(out), contiguous :: eps_std(:, :)   !< (p, n)
    real(dp), intent(out), contiguous :: eta_std(:, :)   !< (r, n)
    real(dp) :: v, nan
    integer :: t, i, ih, iq

    nan = ieee_value(1.0_dp, ieee_quiet_nan)
    do t = 1, rep%nobs
      ih = tidx(size(rep%H, 3), t)
      iq = tidx(size(rep%Q, 3), t)
      do i = 1, rep%k_endog
        v = rep%H(i, i, ih) - sres%epsvar(i, i, t)
        eps_std(i, t) = merge(sres%epshat(i, t) / sqrt(max(v, tiny(1.0_dp))), nan, v > 0.0_dp)
      end do
      do i = 1, rep%k_posdef
        v = rep%Q(i, i, iq) - sres%etavar(i, i, t)
        eta_std(i, t) = merge(sres%etahat(i, t) / sqrt(max(v, tiny(1.0_dp))), nan, v > 0.0_dp)
      end do
    end do
  end subroutine auxiliary_residuals

  !> Least squares residuals (DK 6.2.4): with regression coefficients in
  !> states first..last (diffuse, constant), v+_t = y_t - Z_t a_t - X_t beta_hat
  !> with beta_hat the estimate from the whole sample, computed as the
  !> innovations of the model with those states fixed at beta_hat (the
  !> constructed measurement equation y_t - X_t beta_hat = Z_t alpha_t + eps_t).
  !> The recursive residuals of DK 6.2.4 are the ordinary innovations
  !> `fres%v` of the model itself after the diffuse period.
  subroutine least_squares_residuals(rep, first, last, vplus, info)
    type(ssm_rep_t), intent(in) :: rep
    integer, intent(in) :: first, last
    real(dp), intent(out), contiguous :: vplus(:, :)     !< (p, n)
    integer, intent(out) :: info
    type(ssm_rep_t) :: fixed
    type(filter_result_t) :: fres
    real(dp), allocatable :: a1(:), Pstar(:, :), Pinf(:, :)
    integer :: m

    call kalman_filter(rep, fres, info)
    if (info /= SS_OK) return
    m = rep%k_states
    allocate (a1(m), Pstar(m, m), Pinf(m, m))
    call rep%initial_state(a1, Pstar, Pinf, info)
    if (info /= SS_OK) return
    ! beta_hat: the filtered estimate at n (constant coefficients).
    a1(first:last) = fres%att(first:last, rep%nobs)
    Pstar(first:last, :) = 0.0_dp; Pstar(:, first:last) = 0.0_dp
    Pinf(first:last, :) = 0.0_dp; Pinf(:, first:last) = 0.0_dp
    fixed = rep
    call fixed%initialize_general(a1, Pstar, Pinf)
    call kalman_filter(fixed, fres, info)
    if (info /= SS_OK) return
    vplus = fres%v
  end subroutine least_squares_residuals

  !> Auxiliary residuals standardized as vectors (DK 7.5):
  !> epshat^s_t = B_t epshat_t with [Var(epshat_t)]^-1 = B_t' B_t, and likewise
  !> for eta, using Var(epshat_t) = H_t - Var(eps_t|Y) and Var(etahat_t) =
  !> Q_t - Var(eta_t|Y). B_t = L^-1 from Var = L L' (Cholesky); elements
  !> with zero variance (e.g. no information) are NaN.
  subroutine auxiliary_residuals_vector(rep, sres, eps_std, eta_std)
    type(ssm_rep_t), intent(in) :: rep
    type(smoother_result_t), intent(in) :: sres
    real(dp), intent(out), contiguous :: eps_std(:, :)   !< (p, n)
    real(dp), intent(out), contiguous :: eta_std(:, :)   !< (r, n)
    integer :: t, ih, iq

    do t = 1, rep%nobs
      ih = tidx(size(rep%H, 3), t)
      iq = tidx(size(rep%Q, 3), t)
      eps_std(:, t) = whiten(rep%H(:, :, ih) - sres%epsvar(:, :, t), sres%epshat(:, t))
      eta_std(:, t) = whiten(rep%Q(:, :, iq) - sres%etavar(:, :, t), sres%etahat(:, t))
    end do
  end subroutine auxiliary_residuals_vector

  !> de Jong and Penzer (1998) statistics (DK 7.5): r_i,t / sqrt(N_ii,t) for
  !> the state equation, i = 1..m, t = 1..n, and e_i,t / sqrt(D_ii,t) for the
  !> observation equation with e_t = F_t^-1 v_t - K_t' r_t and
  !> D_t = F_t^-1 + K_t' N_t K_t (DK 4.5.3). NaN where undefined: zero
  !> variance, missing elements, and the diffuse period for e.
  subroutine de_jong_penzer(rep, fres, sres, r_stat, e_stat)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    type(smoother_result_t), intent(in) :: sres
    real(dp), intent(out), contiguous :: r_stat(:, :)    !< (m, n)
    real(dp), intent(out), contiguous :: e_stat(:, :)    !< (p, n)
    real(dp), allocatable :: K(:, :), Finv(:, :), e(:), D(:, :), vz(:)
    real(dp) :: nan
    integer :: t, i, info

    nan = ieee_value(1.0_dp, ieee_quiet_nan)
    allocate (K(rep%k_states, rep%k_endog), Finv(rep%k_endog, rep%k_endog))
    e_stat = nan
    do t = 1, rep%nobs
      do i = 1, rep%k_states
        if (sres%N(i, i, t) > 0.0_dp) then
          r_stat(i, t) = sres%r(i, t) / sqrt(sres%N(i, i, t))
        else
          r_stat(i, t) = nan
        end if
      end do
      call conventional_gain(rep, fres, t, K, Finv, info)
      if (info /= SS_OK) cycle
      vz = merge(0.0_dp, fres%v(:, t), ieee_is_nan(fres%v(:, t)))
      e = matmul(Finv, vz) - matmul(transpose(K), sres%r(:, t))
      D = Finv + matmul(transpose(K), matmul(sres%N(:, :, t), K))
      do i = 1, rep%k_endog
        if (.not. ieee_is_nan(rep%y(i, t)) .and. D(i, i) > 0.0_dp) e_stat(i, t) = e(i) / sqrt(D(i, i))
      end do
    end do
  end subroutine de_jong_penzer

  !> L^-1 x for V = L D L' factored as L D^(1/2) (the Cholesky factor when V is
  !> positive definite); NaN for elements whose pivot is zero.
  function whiten(V, x) result(z)
    real(dp), intent(in), contiguous :: V(:, :), x(:)
    real(dp) :: z(size(x))
    real(dp), allocatable :: L(:, :), D(:)
    integer :: i, j

    allocate (L(size(x), size(x)), D(size(x)))
    call ldl_psd(V, L, D)
    ! Forward substitution with the unit lower L, then scale by D^(-1/2).
    do i = 1, size(x)
      z(i) = x(i)
      do j = 1, i - 1
        z(i) = z(i) - L(i, j) * z(j)
      end do
    end do
    do i = 1, size(x)
      if (D(i) > 0.0_dp) then
        z(i) = z(i) / sqrt(D(i))
      else
        z(i) = ieee_value(1.0_dp, ieee_quiet_nan)
      end if
    end do
  end function whiten

  !> Ljung-Box statistics Q(k) and p-values for lags k = 1..size(stat), from
  !> the sample autocorrelations of x (demeaned, divided by the full-sample
  !> sum of squares). Pairs with a NaN are skipped. `model_df` is subtracted
  !> from the degrees of freedom.
  subroutine ljung_box(x, stat, pvalue, model_df)
    real(dp), intent(in), contiguous :: x(:)
    real(dp), intent(out), contiguous :: stat(:), pvalue(:)
    integer, intent(in), optional :: model_df
    logical, allocatable :: ok(:)
    real(dp), allocatable :: z(:)
    real(dp) :: mean, denom, rho, total
    integer :: n, k, df0

    df0 = 0
    if (present(model_df)) df0 = model_df
    ok = .not. ieee_is_nan(x)
    n = count(ok)
    mean = sum(x, mask=ok) / n
    z = merge(x - mean, 0.0_dp, ok)
    denom = sum(z**2)
    total = 0.0_dp
    do k = 1, size(stat)
      rho = sum(z(1:size(x) - k) * z(k + 1:)) / denom
      total = total + rho**2 / (n - k)
      stat(k) = n * (n + 2.0_dp) * total
      if (k > df0) then
        pvalue(k) = chi2_sf(stat(k), real(k - df0, dp))
      else
        pvalue(k) = ieee_value(1.0_dp, ieee_quiet_nan)
      end if
    end do
  end subroutine ljung_box

  !> Jarque-Bera normality test on the non-NaN elements of x:
  !> JB = n/6 (S^2 + (K - 3)^2 / 4), with the (biased) sample skewness S and
  !> kurtosis K, and p-value P(chi2(2) > JB) = exp(-JB/2).
  subroutine jarque_bera(x, jb, pvalue, skew, kurtosis)
    real(dp), intent(in), contiguous :: x(:)
    real(dp), intent(out) :: jb, pvalue, skew, kurtosis
    real(dp), allocatable :: z(:)
    real(dp) :: m2
    integer :: n

    z = pack(x, .not. ieee_is_nan(x))
    n = size(z)
    z = z - sum(z) / n
    m2 = sum(z**2) / n
    skew = (sum(z**3) / n) / m2**1.5_dp
    kurtosis = (sum(z**4) / n) / m2**2
    jb = n / 6.0_dp * (skew**2 + (kurtosis - 3.0_dp)**2 / 4.0_dp)
    pvalue = exp(-0.5_dp * jb)
  end subroutine jarque_bera

  !> Heteroskedasticity test H(h) (DK 2.12): the sum of the last h squared
  !> elements of x over the sum of the first h (NaNs skipped), with a
  !> two-sided p-value from F(dfn, dfd) applied to H dfd / dfn (dfn, dfd the
  !> non-missing counts). Default h = round(n / 3), rounding half to even as
  !> statsmodels does.
  subroutine breakvar_test(x, stat, pvalue, h)
    real(dp), intent(in), contiguous :: x(:)
    real(dp), intent(out) :: stat, pvalue
    integer, intent(in), optional :: h
    real(dp), allocatable :: last(:), first(:)
    real(dp) :: lower
    integer :: hh, n, dfn, dfd

    n = size(x)
    if (present(h)) then
      hh = h
    else
      hh = nint_even(n / 3.0_dp)
    end if
    last = x(n - hh + 1:)
    first = x(1:hh)
    dfn = count(.not. ieee_is_nan(last))
    dfd = count(.not. ieee_is_nan(first))
    stat = sum(last**2, mask=.not. ieee_is_nan(last)) / sum(first**2, mask=.not. ieee_is_nan(first))
    lower = f_cdf(stat * dfd / dfn, real(dfn, dp), real(dfd, dp))
    pvalue = 2.0_dp * min(lower, 1.0_dp - lower)
  end subroutine breakvar_test

  !> Coefficient of determination against a random walk with drift
  !> (Harvey 1989): R2_D = 1 - SSE / sum_t (Delta y_t - mean Delta y)^2, with
  !> SSE the sum of squared one-step errors after the diffuse period and
  !> burn-in. For a univariate series (element `i`, default 1).
  real(dp) function r2_diffuse(rep, fres, i) result(r2)
    type(ssm_rep_t), intent(in) :: rep
    type(filter_result_t), intent(in) :: fres
    integer, intent(in), optional :: i
    real(dp), allocatable :: dy(:), v(:)
    integer :: k, d

    k = 1
    if (present(i)) k = i
    d = diagnostic_start(rep, fres)
    v = fres%v(k, d:)
    dy = rep%y(k, 2:) - rep%y(k, :rep%nobs - 1)
    dy = pack(dy, .not. ieee_is_nan(dy))
    r2 = 1.0_dp - sum(v**2, mask=.not. ieee_is_nan(v)) / sum((dy - sum(dy) / size(dy))**2)
  end function r2_diffuse

  !> Prediction error variance (DK 7.4): the steady-state F of a
  !> time-invariant model (DK 2.11, 4.3.4; see `steady_state`). Time-varying
  !> models give SS_ERR_UNSUPPORTED.
  subroutine prediction_error_variance(rep, F, info)
    type(ssm_rep_t), intent(in) :: rep
    real(dp), intent(out), contiguous :: F(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: P(:, :)

    info = SS_ERR_DIM
    if (any(shape(F) /= [rep%k_endog, rep%k_endog])) return
    allocate (P(rep%k_states, rep%k_states))
    call steady_state(rep, P, F, info)
  end subroutine prediction_error_variance

  !> Round to the nearest integer, halves to even (numpy's round).
  integer function nint_even(x)
    real(dp), intent(in) :: x
    real(dp) :: f

    nint_even = floor(x)
    f = x - nint_even
    if (f > 0.5_dp .or. (f == 0.5_dp .and. mod(nint_even, 2) /= 0)) nint_even = nint_even + 1
  end function nint_even
end module statespace_diagnostics
