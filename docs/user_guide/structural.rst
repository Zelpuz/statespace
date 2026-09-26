Structural and ARIMA models
===========================

:class:`~ssfortran.StructuralModel` builds a model from components
(DK ch. 3). The state vector stacks the components' states in the order
given; the system matrices are block diagonal (DK eq. 3.10); each block is
initialized as its component needs (DK §5.6).

Components
----------

=====================================  ==========================================  ======
component                              model                                       DK
=====================================  ==========================================  ======
:class:`~ssfortran.Irregular`          :math:`\varepsilon_t`                       §3.2
:class:`~ssfortran.Level`              :math:`\mu_{t+1} = \mu_t + \xi_t`           ch. 2
:class:`~ssfortran.Trend`              local linear or smooth trend                §3.2.1
:class:`~ssfortran.Seasonal`           dummy, trigonometric or Harrison-Stevens    §3.2.2
:class:`~ssfortran.Cycle`              damped or undamped stochastic cycle         §3.2.4
:class:`~ssfortran.Regression`         fixed or random-walk coefficients           §3.2.5
:class:`~ssfortran.ARIMA`              ARIMA(p, d, q)(P, D, Q)_s                   §3.4
:class:`~ssfortran.ContinuousLevel`    level at arbitrary times                    §3.8.1
:class:`~ssfortran.ContinuousTrend`    smooth trend at arbitrary times; splines    §3.8.2
=====================================  ==========================================  ======

The basic structural model, level plus seasonal plus irregular, for
quarterly data:

>>> rng = np.random.default_rng(7)
>>> n = 120
>>> y = (np.cumsum(0.3 * rng.standard_normal(n)) + np.tile([2.0, -1.0, 0.5, -1.5], n // 4)
...      + 0.5 * rng.standard_normal(n))
>>> mod = ss.StructuralModel(y, [ss.Irregular(), ss.Level(), ss.Seasonal(4)])
>>> mod.param_names
['sigma2.irregular', 'sigma2.level', 'sigma2.seasonal']
>>> res = mod.fit()
>>> comps = res.components()
>>> sorted(comps)
['irregular', 'level', 'seasonal']
>>> np.allclose(comps["irregular"] + comps["level"] + comps["seasonal"], y)
True

The components add up to the data: the irregular is the smoothed
:math:`\hat\varepsilon_t`, and the others are :math:`Z \hat\alpha` over their
states.

A damped cycle starts from its unconditional distribution, as DK §3.2.4
give; statsmodels starts it diffuse (see
:doc:`../design/cycle_initialization`). Its period is kept within
``period_bounds``.

Regression and interventions
----------------------------

Regressors enter through :class:`~ssfortran.Regression`, with the
coefficients as diffuse states, fixed or random walks (DK §3.6).
Interventions (DK §3.2.5) are regressors: a step
:math:`w_t = 1\{t \ge \tau\}` shifts the level, a pulse is an outlier, a
ramp changes the slope.

>>> step = (np.arange(n) >= 60).astype(float)
>>> mod = ss.StructuralModel(y + 3 * step, [ss.Irregular(), ss.Level(),
...                                         ss.Seasonal(4), ss.Regression(step)])
>>> res = mod.fit()
>>> beta = res.smooth().smoothed_state[-1, -1]      # the last state is the coefficient
>>> bool(abs(beta - 3.0) < 1.0)
True

ARIMA
-----

:class:`~ssfortran.ARIMA` puts an ARIMA model in state space form
(DK §3.4), with the differences as diffuse states and the ARMA part
stationary. The state layout equals statsmodels' SARIMAX. A constant or
regressors make a regression with ARMA errors (DK §3.6.2):

>>> z = rng.standard_normal(n)
>>> mod = ss.StructuralModel(y, [ss.Regression(np.column_stack([np.ones(n), z])),
...                              ss.ARIMA(order=(1, 0, 1))])
>>> mod.param_names
['ar.L1', 'ma.L1', 'sigma2']

Several series
--------------

With p > 1 each component applies to every series: the series have their
own states and correlated disturbances (seemingly unrelated time series
equations, DK §3.3). ``cov="full"`` estimates the full disturbance variance
through its Cholesky factor. :class:`~ssfortran.Regression` and
:class:`~ssfortran.ARIMA` apply to one series, chosen by ``series``.

Signal loadings, :math:`Z = \Lambda Z_{sig}` (``loading``,
``loading_free``), give common levels, latent risk and dynamic factor
models (DK §3.3.2, §3.3.3, §3.7). Components then describe the signals,
except those with ``at_observations=True``, such as series-specific
intercepts. The example of DK §8.3 uses a common level:
:doc:`../examples/dk_8_3`.

Continuous time and splines
---------------------------

:class:`~ssfortran.ContinuousLevel` and :class:`~ssfortran.ContinuousTrend`
take the observation times, which may be irregular or repeated
(DK §3.8). A continuous trend with an irregular gives the cubic smoothing
spline, with smoothing parameter
:math:`\sigma^2_\varepsilon / \sigma^2` estimated by maximum likelihood
(DK §3.9.2); see :doc:`../examples/dk_8_5`.
