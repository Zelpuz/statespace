Estimation
==========

Maximum likelihood
------------------

A model maps a parameter vector :math:`\psi` to the system matrices;
:meth:`~ssfortran.Model.fit` maximizes the log likelihood over it with
L-BFGS-B (DK §7.3.2). The optimizer works on unconstrained values, which
the model transforms: variances by :math:`\sigma^2 = \exp(2\psi)`
(DK §7.3.2), bounded parameters such as a cycle's damping by
:math:`m + h\psi / \sqrt{1 + \psi^2}`, AR and MA coefficients so that the
polynomials stay stationary and invertible (Monahan 1984).

>>> y = np.loadtxt("data/nile.csv", delimiter=",", skiprows=1)[:, 1]
>>> mod = ss.StructuralModel(y, [ss.Irregular(), ss.Level()])
>>> res = mod.fit()
>>> res.converged, res.analytic_gradient
(True, True)

The defaults (``factr = 1e7``, ``pgtol = 1e-5``) are statsmodels'. The
likelihood of small models is often flat near the maximum, so these stop
up to about 1% from the optimum in some parameters; ``factr = 10`` and
``pgtol = 1e-9`` give a tight optimum.

The gradient
------------

The score (DK §7.3.3) comes from one smoother pass. By the Fisher identity,
for parameters in :math:`H_t`, :math:`R_t Q_t R_t'` and the stationary part
of the initial variance (DK eq. 7.16),

.. math::

   \frac{\partial \log L}{\partial \psi_i}
   = \tfrac12 \sum_t \mathrm{tr}\Big(\frac{\partial H_t}{\partial \psi_i}
     H_t^{-1}\big(\hat\varepsilon_t \hat\varepsilon_t' + Var(\varepsilon_t | Y)
     - H_t\big) H_t^{-1}\Big)
   + \tfrac12 \sum_t \mathrm{tr}\Big(\frac{\partial (R_t Q_t R_t')}{\partial \psi_i}
     (r_t r_t' - N_t)\Big) + \dots

Parameters that move Z or T (AR coefficients, a cycle's frequency, factor
loadings) get central differences, as DK recommend. The choice is made per
parameter, so a model with both kinds uses the score for its variances.
``res.analytic_gradient`` reports whether the score was used. See
:doc:`../design/estimation`.

Standard errors
---------------

``res.cov_params`` is the inverse of the negative Hessian of the log
likelihood in the unconstrained parameters, mapped to the constrained ones
by the delta method (DK §7.3.6). The Hessian is taken on the unconstrained
scale because steps there are always valid, also for a variance near zero.

Concentrated scale
------------------

With ``concentrate_scale = True``, H, Q and the proper part of the initial
variance are relative to a scale :math:`\sigma^2`, estimated in closed form
at each evaluation (DK §2.10.2); the model then has one parameter fewer.

>>> rep = ss.Representation(y, k_states=1)
>>> rep["design"] = rep["transition"] = rep["selection"] = [[1.0]]
>>> rep["obs_cov"] = [[1.0]]
>>> rep.initialize_diffuse()
>>> q = ss.MappedModel(rep, 1, ["q"], start_params=[0.1])
>>> _ = q.map(0, "state_cov", 0, 0).constrain(0, "positive")
>>> q.concentrate_scale = True
>>> rq = q.fit(factr=10, pgtol=1e-9)
>>> round(float(rq.params[0]), 3), round(rq.scale, -1)
(0.097, 15100.0)

This is the signal-to-noise ratio :math:`q = \sigma^2_\eta /
\sigma^2_\varepsilon` of DK §2.10.2.

EM
--

:meth:`~ssfortran.Representation.em` estimates H and Q by the EM algorithm
(DK §7.3.4). Each iteration raises the likelihood; convergence is slow near
the optimum, so EM serves best to find starting values.

Information criteria
--------------------

``res.aic`` and ``res.bic`` count the diffuse states and a concentrated
scale as parameters, as statsmodels does (DK §7.4 also divide by n).

Uncertainty from estimation
---------------------------

Smoothed states treat :math:`\hat\psi` as known.
:meth:`~ssfortran.FitResults.estimation_bias` estimates the bias this
causes (DK §7.3.7, eq. 7.20) by drawing parameters from their estimated
distribution.

Many series
-----------

:func:`~ssfortran.fit_many` fits independent built-in or mapped models in
parallel, releasing the GIL. Models defined by Python callbacks
(:class:`~ssfortran.MLEModel`) cannot run on the library's threads.
