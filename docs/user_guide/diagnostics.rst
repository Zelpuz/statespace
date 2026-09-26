Diagnostic checking
===================

Residuals
---------

The standardized one-step prediction errors
:math:`e_t = L_t^{-1} v_t`, :math:`F_t = L_t L_t'`, are independent standard
normal when the model is correct (DK §2.12, §7.5). They are NaN in the
diffuse periods and where observations are missing.

>>> y = np.loadtxt("data/nile.csv", delimiter=",", skiprows=1)[:, 1]
>>> res = ss.StructuralModel(y, [ss.Irregular(), ss.Level()]).fit()
>>> e = res.standardized_residuals
>>> bool(np.isnan(e[0]))
True

Tests
-----

:meth:`~ssfortran.FitResults.diagnostics` applies three tests to the
residuals after the diffuse periods (DK §2.12.1): Ljung-Box for serial
correlation, Jarque-Bera for normality and H for heteroskedasticity. They
are also available on any series in :mod:`ssfortran.diagnostics`, and
:meth:`~ssfortran.FitResults.summary` reports them:

>>> print(res.summary())                                 # doctest: +ELLIPSIS
============...
Model:          StructuralModel         Log likelihood:       -633.4646
Observations:   100                     AIC:                  1272.9291
Diffuse periods:1                       BIC:                  1280.7446
============...

Auxiliary residuals
-------------------

The smoothed disturbances, standardized by their variances
(:meth:`~ssfortran.Representation.auxiliary_residuals`), point to outliers
(observation disturbances) and structural breaks (state disturbances)
(DK §2.12.2, §7.5). The statistics of de Jong and Penzer (1998)
(:meth:`~ssfortran.Representation.de_jong_penzer`) do the same for every
state. In the Nile series the auxiliary level residual is largest in 1899,
when the Aswan dam was built:

>>> eps, eta = res.model.representation(res.params).auxiliary_residuals()
>>> int(1871 + np.nanargmax(np.abs(eta[0])))
1898

The state disturbance of 1898 moves the level of 1899.

Goodness of fit
---------------

:meth:`~ssfortran.Representation.r2_diffuse` compares the one-step
prediction errors with those of a random walk with drift (Harvey 1989);
:meth:`~ssfortran.Representation.steady_state` gives the prediction error
variance of DK §7.4.
