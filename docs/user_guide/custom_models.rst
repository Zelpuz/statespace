Defining models
===============

Models outside the built-in components are defined in one of three ways.

.. list-table::
   :header-rows: 1

   * - kind
     - how
     - estimation
   * - :class:`~ssfortran.MappedModel`
     - declare which entries each parameter sets
     - in Fortran; works with ``fit_many``
   * - :class:`~ssfortran.MLEModel`
     - subclass and write ``update(params)``
     - calls Python at every evaluation
   * - Fortran ``ssm_model_t``
     - extend the type
     - in Fortran

Declared models
---------------

A :class:`~ssfortran.MappedModel` starts from a
:class:`~ssfortran.Representation` holding the fixed entries and the
initialization. :meth:`~ssfortran.MappedModel.map` lets a parameter set an
entry, :meth:`~ssfortran.MappedModel.cov` a covariance block, and
:meth:`~ssfortran.MappedModel.constrain` sets the transforms. An ARMA(1, 1)
in DK's form (DK §3.4):

>>> y = np.random.default_rng(8).standard_normal(200)
>>> rep = ss.Representation(y, k_states=2, k_posdef=1)
>>> rep["design"] = [[1.0, 0.0]]
>>> rep["transition"] = [[0.0, 1.0], [0.0, 0.0]]
>>> rep["selection"] = [[1.0], [0.0]]
>>> rep.initialize_stationary()
>>> arma = ss.MappedModel(rep, 3, ["ar.L1", "ma.L1", "sigma2"], start_params=[0.0, 0.0, 1.0])
>>> _ = arma.map(0, "transition", 0, 0)          # T[0, 0] = phi
>>> _ = arma.map(1, "selection", 1, 0)           # R[1, 0] = theta
>>> _ = arma.map(2, "state_cov", 0, 0)           # Q = sigma2
>>> _ = arma.constrain(0, "stationary").constrain(1, "invertible").constrain(2, "positive")
>>> res = arma.fit()
>>> res.param_names
['ar.L1', 'ma.L1', 'sigma2']

It gives the same likelihood as ``ss.ARIMA(order=(1, 0, 1))``. Declared
models run entirely in Fortran, so they are as fast as the built-in ones
and :func:`~ssfortran.fit_many` can fit them in parallel.

Starting from components
------------------------

The built-in components can supply the fixed part of a declared model.
:meth:`~ssfortran.Model.representation` returns a structural model's
representation at given parameters: the components' system matrices and
initialization, with their states in the order the components were given.
A :class:`~ssfortran.MappedModel` copies it as its template; mapped
entries overwrite it at each evaluation. Mapping the structural model's
own parameters reproduces it:

>>> y = np.loadtxt("data/nile.csv", delimiter=",", skiprows=1)[:, 1]
>>> st = ss.StructuralModel(y, [ss.Irregular(), ss.Trend()])
>>> st.param_names
['sigma2.irregular', 'sigma2.level', 'sigma2.slope']
>>> rep = st.representation(st.start_params)
>>> mod = ss.MappedModel(rep, 3, st.param_names, start_params=st.start_params)
>>> _ = mod.map(0, "obs_cov", 0, 0)       # H = sigma2.irregular
>>> _ = mod.map(1, "state_cov", 0, 0)     # Q[0, 0] = sigma2.level
>>> _ = mod.map(2, "state_cov", 1, 1)     # Q[1, 1] = sigma2.slope
>>> _ = mod.constrain([0, 1, 2], "positive")
>>> bool(np.isclose(mod.fit().llf, st.fit().llf))
True

From there, further ``map`` calls add parameters of your own; give the
model more parameters and names to match. ``rep["transition"]`` and the
other matrices show the layout to map into.

A map sets an entry to a multiple of one parameter, which covers the
irregular, level, trend, seasonal and regression components and
non-seasonal ARIMA. The cycle (:math:`\rho \cos\lambda` in :math:`T`) and
seasonal ARIMA (products of coefficients) are not of that form. For
those, use an :class:`~ssfortran.MLEModel` whose ``update`` copies the
matrices of ``st.representation(params)`` and changes what it needs;
building the representation at each evaluation costs some speed.

Python models
-------------

:class:`~ssfortran.MLEModel` follows statsmodels: the subclass sets the
fixed matrices in ``__init__`` and overrides ``update``,
``start_params`` and, optionally, the transforms and names. See the
example in :class:`~ssfortran.MLEModel`.

The library calls ``update`` at every likelihood evaluation, and the score
calls it 2k more times per gradient for its matrix derivatives. For the
local level model this is still about twice as fast as statsmodels (see
:doc:`performance`). An exception raised in ``update`` propagates to the
caller, and the model remains usable. Python models cannot be fitted with
``fit_many``.

Fortran models
--------------

In Fortran, extend ``ssm_model_t``, set up ``rep`` in a constructor, and
implement ``update`` and ``start_params``; see
:doc:`../reference/fortran/statespace_model`. New components for
structural models extend ``component_t``
(:doc:`../reference/fortran/statespace_components`).
