Forecasting and simulation
==========================

Forecasting
-----------

Forecasts are predictions with missing observations after the sample
(DK §4.11). :meth:`~ssfortran.FitResults.get_forecast` returns the means and
variances for time-invariant models:

>>> y = np.loadtxt("data/nile.csv", delimiter=",", skiprows=1)[:, 1]
>>> res = ss.StructuralModel(y, [ss.Irregular(), ss.Level()]).fit()
>>> fc = res.get_forecast(3)
>>> fc.predicted_mean.shape, fc.var_pred_mean.shape
((3,), (3, 1, 1))
>>> lower, upper = fc.conf_int(alpha=0.1)

For models with time-varying matrices, append NaN observations to the data,
with the matrices for the forecast periods, and read ``forecasts`` and
``forecasts_error_cov`` from the filter.

Simulation
----------

:meth:`~ssfortran.Representation.simulate` draws observations, states and
disturbances from the model. The random input is standard normal variates,
from a numpy generator or given directly, so draws are reproducible and can
match statsmodels' for the same variates.

The simulation smoothers draw states and disturbances from their
distribution given the data:

* :meth:`~ssfortran.Representation.simulation_smoother`: the
  mean-correction method of Durbin and Koopman (2002) (DK §4.9.2), also
  with diffuse states (DK §5.5);
* :meth:`~ssfortran.Representation.djs_simulation_smoother`: the method of
  de Jong and Shephard (1995) (DK §4.9.3).

>>> rep = res.model.representation(res.params)
>>> draws = np.array([rep.simulation_smoother(rng=i)[0] for i in range(200)])
>>> draws.shape
(200, 1, 100)

The mean of the draws approaches the smoothed state as their number grows.
For a missing element the drawn observation disturbance is conditional on
the observed elements; statsmodels draws it from its unconditional
distribution.

The effect of estimating the parameters
---------------------------------------

The smoothed state depends on the estimated parameters.
:meth:`~ssfortran.FitResults.estimation_bias` estimates the bias from
treating them as known (DK §7.3.7):

>>> bias = res.estimation_bias(ndraw=200, seed=1)
>>> bias.shape
(1, 100)
