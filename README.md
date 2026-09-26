# statespace

Linear Gaussian state space models, following Part I of Durbin and Koopman,
*Time Series Analysis by State Space Methods* (2nd ed., 2012). The core is a Fortran
library; `ssfortran` is its Python package.

- **Filtering and smoothing:**
  - the Kalman filter with a steady-state shortcut
  - exact diffuse initialization, in both univariate and multivariate forms
  - univariate treatment, including correlated and time-varying H
  - state, disturbance, fast, classical, two-filter, fixed-point and fixed-lag smoothers
  - augmented and square-root filters
  - collapsing large observation vectors, and linear restrictions
- **Likelihood and estimation:**
  - exact, diffuse, concentrated and marginal log likelihoods
  - maximum likelihood with L-BFGS-B, using DK's analytic score where it applies and
    numerical derivatives elsewhere
  - EM for variance parameters
  - standard errors, and the effect of parameter estimation on the smoothed states
- **Simulation, forecasting and diagnostics:** simulation smoothers, forecasts,
  standardized and auxiliary residuals, and residual tests.
- **Built-in models (DK ch. 3):** irregular, level, trend, seasonal (dummy,
  trigonometric, Harrison–Stevens), cycle, regression and intervention effects, ARIMA,
  continuous-time components and splines. Components apply to one series or several
  (SUTSE), and can load on common signals.

The documentation (user guide, examples, design notes, and Python, Fortran and C
references) builds with `make -C docs html`; see [docs/install.rst](docs/install.rst).
The illustrations of DK chapter 8 are reproduced in `example/`. Differences from
statsmodels' `tsa.statespace` are listed in
[docs/statsmodels_differences.md](docs/statsmodels_differences.md).

## Fortran

The library builds with [fpm](https://fpm.fortran-lang.org) and needs gfortran,
LAPACK and BLAS.

```sh
fpm test --profile release                              # the test suite
fpm run --profile release --example nile_mle            # a model defined in Fortran
fpm run --profile release --example dk_8_2_seatbelt     # needs data/fetch_dk_data.py first
```

Models are defined by extending `ssm_model_t` (see `example/nile_mle.f90`), or
assembled from components with `structural_model` (see `example/dk_8_2_seatbelt.f90`).
`fit_many` fits independent models in parallel when built with
`--flag -fopenmp --link-flag -fopenmp`.

## Python

`pip` builds the shared library with CMake through scikit-build-core. It needs
gfortran, LAPACK and BLAS installed.

```sh
pip install .              # or: pip install .[test] && pytest
```

```python
import numpy as np
import ssfortran as ss

y = np.loadtxt("data/nile.csv", delimiter=",", skiprows=1)[:, 1]

# Built-in components
mod = ss.StructuralModel(y, [ss.Irregular(), ss.Level()])
res = mod.fit()
print(res.summary())
smoothed = res.smooth().smoothed_state

# Matrix-level: fill the system matrices yourself
rep = ss.Representation(y, k_states=1)
rep["design"] = [[1.0]]; rep["transition"] = [[1.0]]; rep["selection"] = [[1.0]]
rep["obs_cov"] = [[15099.0]]; rep["state_cov"] = [[1469.1]]
rep.initialize_diffuse()
print(rep.loglike())
```

Fitting results have `summary()` (estimates with z-tests and intervals, fit
statistics, residual tests), `fittedvalues`, `resid`, `get_forecast(steps)` with
`conf_int()`, `components()` and `estimation_bias()`. A `Representation` also offers:
- the other smoothers of DK ch. 4: fast, classical, two-filter, Whittle, fixed-point,
  fixed-lag and updating
- smoothed covariances between periods, and filtering and smoothing weights
- the square-root and augmented filters (with the marginal likelihood)
- EM, collapsing, and linear state restrictions
- auxiliary residuals, de Jong–Penzer statistics, least squares residuals and R²_D
- the mean-correction and de Jong–Shephard simulation smoothers

A model with parameters can be defined in three ways:

| | How | Speed |
|---|---|---|
| `StructuralModel` | built-in components | all in Fortran; works with `fit_many` |
| `MappedModel` | declare which matrix entries each parameter sets | all in Fortran; works with `fit_many` |
| `MLEModel` | subclass and write `update(params)`, as in statsmodels | calls Python for each likelihood evaluation |

For the timing comparison with statsmodels from Python, run
`bench/bench_python.py`. It shows 2–8× speedups for single models, and 36× for 1000
series with `fit_many`. `example/bench.f90` and `bench/bench_statsmodels.py` compare
the Fortran core directly.

To develop without installing:

```sh
cmake -S . -B build/cmake -G Ninja && cmake --build build/cmake
PYTHONPATH=python pytest python/tests
```

## License

MIT (see [LICENSE](LICENSE)); third-party notices are in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md), and citations in
[docs/references.md](docs/references.md).
