# statespace

Linear Gaussian state space models, following Part I of Durbin and Koopman, *Time Series Analysis by State Space Methods* (2nd ed., 2012). Computation is done in Fortran. Users can write models in Fortran or Python.

- **Filtering and smoothing:**
    - the Kalman filter with a steady-state shortcut
    - exact diffuse initialization, in both univariate and multivariate forms
    - univariate treatment, including correlated and time-varying H
    - state, disturbance, fast, classical, two-filter, fixed-point and fixed-lag smoothers
    - augmented and square-root filters
    - collapsing large observation vectors, and linear restrictions
- **Likelihood and estimation:**
    - exact, diffuse, concentrated and marginal log likelihoods
    - maximum likelihood with L-BFGS-B, using DK's analytic score where it applies and numerical derivatives elsewhere
    - EM for variance parameters
    - standard errors, and the effect of parameter estimation on the smoothed states
- **Simulation, forecasting and diagnostics:** simulation smoothers, forecasts, standardized and auxiliary residuals, and residual tests.
- **Built-in models (DK ch. 3):** irregular, level, trend, seasonal (dummy, trigonometric, Harrison–Stevens), cycle, regression and intervention effects, ARIMA, continuous-time components and splines. Components apply to one series or several (SUTSE), and can load on common signals.

**Documentation:** <https://zelpuz.github.io/statespace/>, with a user guide, the examples of DK chapter 8, design notes, the Python and Fortran API references, and the [differences from statsmodels](https://zelpuz.github.io/statespace/statsmodels_differences.html).

## Fortran

The library builds with [fpm](https://fpm.fortran-lang.org) and needs gfortran, LAPACK and BLAS.

```sh
fpm test --profile release                              # the test suite
fpm run --profile release --example nile_mle            # a model defined in Fortran
fpm run --profile release --example dk_8_2_seatbelt     # needs data/fetch_dk_data.py first
```

Models are defined by extending `ssm_model_t` (see `example/nile_mle.f90`), or assembled from components with `structural_model` (see `example/dk_8_2_seatbelt.f90`). `fit_many` fits independent models in parallel when built with `--flag -fopenmp --link-flag -fopenmp`.

## Python

```sh
pip install ssfortran
# or
uv add ssfortran           # in a uv project; or: uv pip install ssfortran
```

Wheels for Linux (x86_64 and aarch64, glibc 2.28 or later) include the compiled library, gfortran's runtime and OpenBLAS, so they need no compiler. On other platforms pip and uv build from the source distribution, which needs gfortran, LAPACK and BLAS. The package needs Python 3.10 or later and numpy.

```python
import numpy as np
import ssfortran as ss

y = np.loadtxt("data/nile.csv", delimiter=",", skiprows=1)[:, 1]  # from this repository

# Build a model with built-in components and defaults...
mod = ss.StructuralModel(y, [ss.Irregular(), ss.Level()])
res = mod.fit()
print(res.summary())
smoothed = res.smooth().smoothed_state

# ... or define the matrices and starting parameters yourself.
mod = ss.MappedModel(
    y,
    k_states=1,
    k_params=2,
    param_names=["sigma2.irregular", "sigma2.level"],
    start_params=[np.var(y) / 2, np.var(y) / 2],
)
mod["design"] = mod["transition"] = mod["selection"] = [[1.0]]
mod.initialize_diffuse()
mod.map(0, "obs_cov", 0, 0).map(1, "state_cov", 0, 0)
mod.constrain([0, 1], "positive")
res = mod.fit()
```

Fitting results have `summary()` (estimates with z-tests and intervals, fit statistics, residual tests), `fittedvalues`, `resid`, `get_forecast(steps)` with `conf_int()`, `components()` and `estimation_bias()`. Every model builds a `Representation`, the system at given parameter values; `mod.representation(res.params)` returns it, and it offers:
 - the other smoothers of DK ch. 4: fast, classical, two-filter, Whittle, fixed-point, fixed-lag and updating
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

For the timing comparison with statsmodels from Python, run `bench/bench_python.py`. It shows 2–8× speedups for single models, and 36× for 1000 series with `fit_many`. `example/bench.f90` and `bench/bench_statsmodels.py` compare the Fortran core directly.

## Contributing

Contribute via issues or pull requests. LLM coding assistance is permitted but submissions should have a responsible human user attached.

To develop without installing:

```sh
cmake -S . -B build/cmake -G Ninja && cmake --build build/cmake
pytest            # uses python/src and build/cmake (pyproject.toml)
```

To build the documentation: `pip install --group docs`, then `make -C docs html`.

Code is formatted to 88 columns. `ruff format` and `ruff check` cover the Python (settings in `pyproject.toml`); [Fortitude](https://fortitude.readthedocs.io)'s `fortitude check` covers the Fortran (settings in `fpm.toml`), whose long lines are wrapped by hand.

## License

MIT (see [LICENSE](https://github.com/Zelpuz/statespace/blob/master/LICENSE)); third-party notices are in [THIRD_PARTY_NOTICES.md](https://github.com/Zelpuz/statespace/blob/master/THIRD_PARTY_NOTICES.md), and citations in [docs/references.md](https://github.com/Zelpuz/statespace/blob/master/docs/references.md).

## LLM disclosure

Initial (and possibly subsequent) writes of this package relied on LLM coding assistance. The work is provided as-is, without warranty, as outlined in the `LICENSE` file.
