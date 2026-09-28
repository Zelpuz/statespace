Quickstart
==========

The local level model for the flow of the Nile (DK ch. 2),

.. math::

   y_t = \mu_t + \varepsilon_t, \quad \varepsilon_t \sim N(0, \sigma^2_\varepsilon),
   \qquad \mu_{t+1} = \mu_t + \eta_t, \quad \eta_t \sim N(0, \sigma^2_\eta),

estimated by maximum likelihood with an exact diffuse initial level. DK
report :math:`\hat\sigma^2_\varepsilon = 15099` and
:math:`\hat\sigma^2_\eta = 1469.1`; the likelihood is flat near its
maximum, so optimizers agree to about four digits.

Python
------

A built-in model:

>>> import numpy as np
>>> import ssfortran as ss
>>> y = np.loadtxt("data/nile.csv", delimiter=",", skiprows=1)[:, 1]
>>> mod = ss.StructuralModel(y, [ss.Irregular(), ss.Level()])
>>> res = mod.fit(factr=10, pgtol=1e-9)
>>> np.round(res.params)
array([15099.,  1469.])
>>> round(res.llf, 4)
-633.4646

The smoothed level and its variance:

>>> sm = res.smooth()
>>> level = sm.smoothed_state[0]
>>> np.round(level[:3])
array([1112., 1111., 1105.])

The same model declared from its matrices, with the two variances as
parameters (:doc:`user_guide/custom_models`):

>>> mod = ss.MappedModel(y, k_states=1, k_params=2,
...                      param_names=["sigma2.irregular", "sigma2.level"],
...                      start_params=[np.var(y) / 2, np.var(y) / 2])
>>> mod["design"] = mod["transition"] = mod["selection"] = [[1.0]]
>>> mod.initialize_diffuse()
>>> _ = mod.map(0, "obs_cov", 0, 0).map(1, "state_cov", 0, 0)
>>> _ = mod.constrain([0, 1], "positive")
>>> res = mod.fit(factr=10, pgtol=1e-9)
>>> np.round(res.params)
array([15099.,  1469.])

Every model builds its :class:`~ssfortran.Representation`, the system at
given parameter values. ``mod.representation(res.params)`` returns it, to
run the other algorithms of DK Part I (:doc:`user_guide/filtering`).

Fortran
-------

The same model defined by extending ``ssm_model_t``, from
``example/nile_mle.f90``:

.. literalinclude:: ../example/nile_mle.f90
   :language: fortran
   :lines: 79-

Output::

   $ fpm run --profile release --example nile_mle
   CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL
   iterations: 10  likelihood evaluations: 12
   log likelihood:    -633.464564
   AIC:  1272.9291  BIC:  1280.7446
                             estimate       std err
       sigma2.irregular    15098.5183     3145.5481
           sigma2.level     1469.1764     1280.3752

Next
----

* :doc:`user_guide/index` explains the model, the algorithms and the
  choices behind them.
* :doc:`examples/index` reproduces the illustrations of DK ch. 8.
* :doc:`reference/index` lists every class, routine and argument.
