Estimation
==========

Variance transform
------------------

**Decision.** Variances are estimated through :math:`\sigma^2 =
\exp(2\psi)`, DK's transform (§7.3.2). Parameters within bounds use DK's
bounded transform, shifted to the interval: :math:`m + h \psi / \sqrt{1 +
\psi^2}`. AR and MA coefficients use Monahan's (1984) transform, as in
statsmodels.

**Why.** DK's transforms; :math:`\exp(2\psi)` also keeps the optimizer away
from exactly zero variances, where :math:`x^2` has a flat point.

**Alternatives.** statsmodels uses :math:`\sigma^2 = x^2` and a logistic for
bounded parameters. The constrained estimates agree; the paths of the
optimizer differ.

The gradient: a hybrid score
----------------------------

**Decision.** The gradient uses DK's analytic score (§7.3.3, eq. 7.16;
Koopman and Shephard 1992) for each parameter that moves only H, R, Q or a
stationary part of the initial variance, and central differences for each
parameter that moves Z or T. The matrix derivatives the score needs are
central differences of what the model's ``update`` sets, so models need no
derivative code.

**Why.** The score costs one filter and smoother pass for all covered
parameters, against two likelihood evaluations per parameter for
differences; for the Nile model the fit drops from 70 to 14 likelihood
evaluations. DK recommend numerical derivatives for parameters in Z and T,
where the score formula would need more terms. Deciding per parameter lets a
model with a cycle or ARMA part still use the score for its variances;
before, one such parameter made the whole gradient numerical.

**Alternatives.** statsmodels differentiates by complex steps by default and
offers Harvey's (1989) score, which differentiates the filter recursions.

**Where.** ``analytic_score`` (with its ``analytic`` mask) in
``statespace_score``; ``analytic_gradient`` in ``statespace_mle``.

Standard errors
---------------

**Decision.** The covariance of the estimates is the inverse of the
negative Hessian of the log likelihood in the unconstrained parameters,
mapped to the constrained ones by the delta method (DK §7.3.6).

**Why.** The Hessian in the constrained parameters failed for the DK §8.2
model: its seasonal variance is 5e-7, and a finite-difference step of 1e-4
made it negative. Steps on the unconstrained scale are always valid. At the
maximum the two give the same covariance.

**Alternatives.** statsmodels defaults to the outer product of gradients.

Estimates near zero
-------------------

Variance estimates can be near zero (5e-7 above) and the likelihood flat in
them. The default tolerances (statsmodels') stop up to about 1% from the
optimum in such directions; the examples and tests use ``factr = 10``,
``pgtol = 1e-9``.

Information criteria
--------------------

**Decision.** AIC and BIC count the diffuse states and a concentrated scale
as parameters and use n minus the burn-in, as statsmodels does; they are not
divided by n.

**Why.** Comparability with statsmodels. DK §7.4 divide by n; multiply
their values by n to compare.

EM
--

**Decision.** EM (DK §7.3.4) estimates time-invariant H and Q, optionally
diagonal, and requires an initialization that does not depend on them and
no burn-in.

**Why.** Under these conditions the M-step is closed form. With a
stationary initialization :math:`P_1` depends on Q and the M-step is not.

The effect of parameter estimation
----------------------------------

**Decision.** ``estimation_bias`` implements DK's bias estimate (§7.3.7,
eq. 7.20): draw parameters from their estimated distribution on the
unconstrained scale, in antithetic pairs, and average the change in the
smoothed state.

**Why.** It follows DK. An earlier version followed Hamilton's mean squared
error decomposition; drawing in the constrained space then put about 11% of
the draws for the Nile model at negative variances.
