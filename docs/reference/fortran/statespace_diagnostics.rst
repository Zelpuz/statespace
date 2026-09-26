``statespace_diagnostics``
==========================

Goodness of fit (DK §7.4) and diagnostic checking (DK §2.12, §7.5).

Residuals
---------

.. code-block:: fortran

   subroutine standardized_residuals(fres, e)
     real(dp), intent(out) :: e(:, :)                   ! (p, n)
   subroutine auxiliary_residuals(rep, sres, eps_std, eta_std)
     real(dp), intent(out) :: eps_std(:, :), eta_std(:, :)   ! (p, n), (r, n)
   subroutine auxiliary_residuals_vector(rep, sres, eps_std, eta_std)
   subroutine de_jong_penzer(rep, fres, sres, r_stat, e_stat)
     real(dp), intent(out) :: r_stat(:, :), e_stat(:, :)     ! (m, n), (p, n)
   subroutine least_squares_residuals(rep, first, last, vplus, info)
     real(dp), intent(out) :: vplus(:, :)               ! (p, n)
   integer function diagnostic_start(rep, fres)

``standardized_residuals``
   :math:`e_t = L_t^{-1} v_t` over the observed elements, with
   :math:`F_{t,oo} = L_t L_t'` (Cholesky, as in statsmodels); NaN for
   missing elements and in the diffuse periods (statsmodels reports 0).
``auxiliary_residuals``
   The smoothed disturbances standardized element by element by
   :math:`Var(\hat\varepsilon_t) = H_t - Var(\varepsilon_t | Y)` and
   :math:`Var(\hat\eta_t) = Q_t - Var(\eta_t | Y)`, for detecting outliers
   and structural breaks (DK §2.12.2, §7.5). NaN where the variance is
   zero.
``auxiliary_residuals_vector``
   The same standardized as vectors, :math:`B_t \hat\varepsilon_t` with
   :math:`B_t' B_t = Var(\hat\varepsilon_t)^{-1}` (DK §7.5).
``de_jong_penzer``
   The shock statistics of de Jong and Penzer (1998) (DK §7.5):
   :math:`r_{i,t} / \sqrt{N_{ii,t}}` for the state equation and
   :math:`e_{i,t} / \sqrt{D_{ii,t}}` for the observation equation, with
   :math:`e_t = F_t^{-1} v_t - K_t' r_t` and :math:`D_t = F_t^{-1} + K_t'
   N_t K_t` (DK §4.5.3). NaN where undefined.
``least_squares_residuals``
   DK §6.2.4: with regression coefficients in states ``first..last``
   (diffuse, constant), the residuals at their full-sample estimate
   :math:`\hat\beta`, computed as the innovations of the model with those
   states fixed at :math:`\hat\beta`. The innovations of the original model
   are the recursive residuals.
``diagnostic_start``
   First period of the residuals to use in tests: after the burn-in and the
   diffuse periods.

Tests
-----

On a residual series, NaNs skipped, as in statsmodels (DK §2.12.1):

.. code-block:: fortran

   subroutine ljung_box(x, stat, pvalue, model_df)
     real(dp), intent(in) :: x(:)
     real(dp), intent(out) :: stat(:), pvalue(:)     ! lags 1..size(stat)
     integer, intent(in), optional :: model_df
   subroutine jarque_bera(x, jb, pvalue, skew, kurtosis)
   subroutine breakvar_test(x, stat, pvalue, h)
     integer, intent(in), optional :: h              ! default n / 3

``ljung_box``
   :math:`Q(k) = n(n+2) \sum_{j=1}^k \rho_j^2 / (n - j)` against
   :math:`\chi^2_{k - model\_df}`. Pairs with a NaN are skipped.
``jarque_bera``
   Normality from skewness and kurtosis, against :math:`\chi^2_2`.
``breakvar_test``
   Heteroskedasticity: H(h), the sum of the last h squared residuals over
   the sum of the first h, against F(h, h); two-sided p-value.

Goodness of fit
---------------

.. code-block:: fortran

   real(dp) function r2_diffuse(rep, fres, i) result(r2)
   subroutine prediction_error_variance(rep, F, info)
     real(dp), intent(out) :: F(:, :)

``r2_diffuse``
   :math:`R^2_D = 1 - SSE / \sum_t (\Delta y_t - \overline{\Delta y})^2`,
   against a random walk with drift (Harvey 1989), for series i.
``prediction_error_variance``
   The steady-state :math:`\bar F` of a time-invariant model (DK §4.3.4,
   §7.4; see ``steady_state``). DK's printed values in §8.2 differ by the
   factor (n - d)/n; see :doc:`../../examples/dk_8_2`.
