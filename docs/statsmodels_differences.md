# Differences from statsmodels

This library follows Durbin and Koopman, *Time Series Analysis by State Space Methods*
(2nd ed., 2012). Where statsmodels' `tsa.statespace` makes a different choice,
we follow DK. This page lists those choices, the features only one of the two libraries
has, and the statsmodels bugs we found while building the test fixtures.

The comparison is against **statsmodels 0.15.0** (the version in `.venv`). Section
numbers refer to DK.

## Methods in both libraries

### Filter and smoother output

| | This library | statsmodels |
|---|---|---|
| v, F, F∞ in diffuse and univariate periods | In the original coordinates, the same as in conventional periods. The per-element values of the univariate recursion are also stored, in the `uv_*` arrays | In the transformed (element-by-element) coordinates |
| K and F⁻¹ in univariate periods | NaN, because the univariate recursion computes neither | The values from the transformed coordinates |
| F∞ for a missing element | Z P∞ Z' like any other element | 0 |
| Smoothed ε in diffuse and univariate periods | Full matrices, from the exact identity ε_o = y_o − d_o − Z_o α_t | Diagonal variances only, in the transformed coordinates |
| Smoothed ε for a missing element | The true conditional moments: mean H u, variance H − H D H (DK 4.5). EM needs these | Mean 0 and variance H, the unconditional moments |
| Standardized residuals where undefined (missing, or F = 0 in the diffuse period) | NaN | 0 |
| Smoothing weights (4.8.3) in diffuse periods | Computed, including the weights on c and a₁ | NaN |

The two libraries agree everywhere else. Across the fixtures, the largest relative
difference is about 1e-10.

### Filter variants

- **Exact diffuse filter (5.2).** statsmodels always uses the univariate form.
  Ours also uses it by default. `DIFFUSE_MULTIVARIATE` switches to DK's multivariate
  form, which falls back to the univariate form in periods where F∞ is singular but
  nonzero. `fres%method(t)` records which recursion each period used.
- **Univariate treatment with correlated H (6.4.3).** Ours works with H_t that is
  non-diagonal, time-varying and partly missing: each period's observed block gets its
  own LDL' factorization.
- **Diffuse tolerances.** These match statsmodels: 1e-10 on F∞ and on ‖P∞‖²_F. Elements
  with F\* ≤ 1e-10 are skipped in the diffuse period, as in statsmodels.
- **Steady state.** statsmodels' conventional filter switches to a steady-state
  shortcut once P_t converges. Ours doesn't yet (milestone 8), which explains the
  ~1e-10 differences in some time-invariant multivariate fixtures.

### Default initialization of the built-in models

| | This library | statsmodels UC and SARIMAX |
|---|---|---|
| Nonstationary states | Exact diffuse (5.2) | Approximate diffuse, with variance 1e6 (1e10 for SARIMAX with state regression). `use_exact_diffuse=True` selects exact diffuse |
| Damped cycle | Stationary, as in DK 3.2.4 | Diffuse, like the other UC states |
| Log likelihood | Every observation contributes | With the approximate diffuse default, `loglikelihood_burn` equals the number of diffuse states and those observations are left out of `llf`. `llf_obs` still has every period, so the tests compare `llf_obs` |

### Parameter transforms

| Parameter | This library | statsmodels |
|---|---|---|
| Variances | σ² = exp(2ψ) (DK 7.3.2) | σ² = x² |
| Full covariance matrices | Σ = LL', with L's diagonal exp(ψ) | UC has no full covariances |
| Cycle frequency λ in (λ_lo, λ_hi) | mid + half · ψ/√(1+ψ²) | A logistic between the bounds |
| Cycle damping ρ in (0, 1) | ψ/√(1+ψ²), mapped onto (0, 1) | A logistic |
| AR and MA coefficients | Monahan's transform; MA = −(the AR transform) | The same, so the values agree to 1e-14 |

The constrained parameters mean the same thing in both libraries, so fitted values and
log likelihoods agree. The unconstrained values, and the optimizer path, differ.

The cycle period bounds default to 2 up to the number of observations. statsmodels uses
1.5 to 12 years when the data's frequency is known, and (2, ∞) otherwise.

### Estimation

- **Optimizer.** Both use L-BFGS-B: ours is the Fortran `jacobwilliams/lbfgsb`, and
  statsmodels uses SciPy's. The defaults are statsmodels': factr = 1e7, pgtol = 1e-5,
  m = 10, minimizing −llf/n.
- **Gradient.**
  - Ours defaults to DK's analytic score (7.14)/(7.16), which comes from one smoother
    pass through the Fisher identity. It covers parameters in H, R, Q and a
    stationary P\*. Parameters that move Z or T fall back to finite differences, as DK
    recommend.
  - statsmodels defaults to complex-step differentiation. Its analytic score is
    `_score_harvey` (Harvey 1989), which differentiates the filter recursions.
- **Standard errors.** Ours come from a finite-difference Hessian of the log
  likelihood (DK 7.3.6). statsmodels defaults to the outer product of gradients
  (`cov_type="opg"`). Its `cov_type="approx"` is closer to ours: a numerical
  Hessian, by complex step.
- **AIC/BIC.** We follow statsmodels: k = k_params + k_diffuse (+1 with a concentrated
  scale), and n = nobs − burn. DK (7.4) also divide by n; neither library does.

### Built-in models

**Unobserved components** (`structural_model_t` vs `UnobservedComponents`)

- Ours is multivariate: SUTSE (DK 3.3) with diagonal, full or no disturbance
  covariance per component, and signal loadings for common levels, latent risk and
  dynamic factors. statsmodels' UC is univariate.
- Ours has the Harrison–Stevens seasonal (3.2.2) alongside the dummy and trigonometric
  forms. statsmodels' `freq_seasonal` can keep fewer than s/2 harmonics; ours always
  uses all of them.
- Interventions (step, pulse, slope) and regression effects are components in ours.
  Regression coefficients are diffuse states, either fixed or random walks. statsmodels
  handles these through `exog`, estimated either by MLE or as states.
- statsmodels' `autoregressive` component corresponds to adding an `arima_t` component
  in ours.

**ARIMA** (`arima_t` vs `SARIMAX`)

- The state layout is the same, DK 3.4 with the differences held as states, so the
  states agree one by one.
- Ours allows seasonal differencing up to D = 1; statsmodels allows any D.
- statsmodels has several options ours lacks: `trend` polynomials,
  `simple_differencing`, `hamilton_representation` and `measurement_error`. In ours, a
  trend comes from a `regression_t` component and measurement error from an
  `irregular_t` component.
- SARIMAX with `time_varying_regression` corresponds to `regression_t` with
  `random_walk`.

### Forecasting

Our `forecast` supports only time-invariant models. For a time-varying model, append
NaNs to y and filter. statsmodels' `get_forecast` handles time-varying models, given
future exog.

## Only in this library

| Feature | DK |
|---|---|
| Augmented Kalman filter and smoother, with δ̂ and Var(δ̂). statsmodels defines `FILTER_AUGMENTED` but doesn't implement it | 5.7 |
| Square-root filter and smoother, by Householder QR. statsmodels defines `FILTER_SQUARE_ROOT` but doesn't implement it | 6.3 |
| Multivariate exact initial filter and smoother | 5.2–5.3 |
| General initialization from any (a₁, P\*, P∞), as well as block-wise mixtures | 5.1 |
| Fast, two-filter and Whittle smoothers | 4.6.2–4.6.4 |
| Updating smoothed estimates, and fixed-point and fixed-lag smoothers | 4.4.5–4.4.6 |
| Filtering weights | 4.8.2 |
| de Jong–Shephard disturbance simulation smoother, which also draws α₁ | 4.9.3 |
| Dense matrix form of the likelihood and smoother, used as a test oracle | 4.13 |
| Linear restrictions on the states | 6.6 |
| Likelihood with elements of α₁ fixed but unknown | 7.2.4 |
| Marginal likelihood for regression effects (Francke, Koopman and de Vos 2010) | 7.2.6 |
| EM for H and Q in any model. statsmodels has EM only in `DynamicFactorMQ` | 7.3.4 |
| Effect of parameter estimation errors on smoothed states (bias estimate, antithetic draws) | 7.3.7 |
| Auxiliary residuals, element-wise and vector-standardized; de Jong–Penzer statistics | 7.5 |
| Least squares residuals | 6.2.4 |
| R²_D and prediction error variance | 7.4 |
| Harrison–Stevens seasonal; multivariate structural models (SUTSE, common levels, latent risk) | 3.2–3.3 |
| Continuous-time local level and smooth trend, and weighted irregulars for unequal spacing | 3.8 |
| Discrete and continuous smoothing splines | 3.9 |

## Only in statsmodels

| Feature | Notes |
|---|---|
| Chandrasekhar recursions (`FILTER_CHANDRASEKHAR`) | Not in DK |
| Steady-state switch in the conventional filter | Planned for milestone 8 |
| Memory-conservation options, filter timing, choice of matrix inversion method, forced symmetry | Ours always uses Cholesky for F and stores every output |
| Chan–Jeliazkov (CFA) simulation smoother | Not in DK |
| News and revisions (`news`), smoothed-state decomposition and gains, impulse responses | |
| `append`/`extend`/`apply` for new data, `fix_params`, `fit_constrained` | |
| Forecasting for time-varying models | See [Forecasting](#forecasting) |
| VARMAX, DynamicFactorMQ (mixed frequency), ETS `ExponentialSmoothing`, `RecursiveLS` with CUSUM tests | Our DK 3.5 exponential smoothing is only tested as equivalent to the local level filter. Recursive residuals are the innovations of a regression model |
| OPG (the default) and robust covariance types | Ours has only the numerical Hessian |
| Results objects: summaries, plots, pandas indexes | Planned for the Python layer (milestone 9) |

## Known statsmodels bugs

These were found while comparing against statsmodels and confirmed against DK or
against an independent route to the same result. We have not reported them upstream.

1. **The exact diffuse smoother uses the wrong transition matrix when T varies.**
   - **Where:** `_smoothers/_univariate_diffuse.pyx.in`, around line 416. The source has
     a `TODO` asking whether this is the right transition matrix when the matrices vary
     over time.
   - **Effect:** r and r⁽¹⁾ are propagated with `model._transition`, the current
     period's T, instead of the previous period's. When T changes inside the diffuse
     period, the smoothed states are wrong.
   - **Evidence:** in the `mv_diffuse_timevarying` fixture, statsmodels' α̂ stays 1.46
     away from the large-κ approximate-diffuse limit for every κ, while ours converges
     at the rate 1/κ (test `exact_diffuse_vs_approximate`). Our tests use that fixture
     only for filter output, and for smoother output after the diffuse period.
2. **The simulation smoother mishandles a singular P\*.**
   - **Where:** `_simulation_smoother.pyx.in`, around line 640. `cholesky()` calls
     `potrf` without checking `info`.
   - **Effect:** with a mixed diffuse and stationary initialization, P\* is singular,
     `potrf` fails, and the initial-state variate ends up scaled by the variance
     instead of the standard deviation. For the AR state in our mixed fixture that
     gives 398.47 where the standard deviation is expected.
   - **Evidence:** our test checks this case through the moments of 4000 draws
     instead of statsmodels' draws.
3. **The simulation smoother draws unconditional ε for missing elements.** For a
   missing element of y_t, the drawn ε comes from N(0, H) and ignores the observed
   elements it is correlated with. This follows from the smoother convention for
   missing elements (mean 0, variance H) described above. Ours draws from the
   conditional distribution.

Two more differences can look like bugs but are conventions:

- **`loglikelihood_burn`.** UC and SARIMAX leave the first observations out of `llf`
  by default, and a fitted SARIMAX's `llf` leaves out the first d. Use `llf_obs` to
  compare.
- **Diffuse damped cycle.** statsmodels starts a damped (stationary) cycle as diffuse.
  DK 3.2.4 give its stationary distribution, which changes the log likelihood; in one
  fixture, −65.72 against statsmodels' −63.60.
