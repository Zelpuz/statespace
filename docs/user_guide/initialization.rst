Initialization
==============

The initial state is (DK eq. 5.2)

.. math::

   \alpha_1 = a + A\delta + R_0 \eta_0, \qquad \eta_0 \sim N(0, Q_0),

with :math:`\delta` a vector of diffuse elements: its variance is
:math:`\kappa I` with :math:`\kappa \to \infty`. Then
:math:`\alpha_1 \sim N(a, P_* + \kappa P_\infty)` with
:math:`P_* = R_0 Q_0 R_0'` and :math:`P_\infty = A A'`.

Kinds
-----

:meth:`~ssfortran.Representation.initialize_known`
   :math:`N(a_1, P_1)`.
:meth:`~ssfortran.Representation.initialize_diffuse`
   Exact diffuse, :math:`P_\infty = I` (DK §5.2).
:meth:`~ssfortran.Representation.initialize_approximate_diffuse`
   :math:`N(0, \kappa I)` with a large finite :math:`\kappa`, default 1e6.
:meth:`~ssfortran.Representation.initialize_stationary`
   The unconditional distribution: mean :math:`(I - T)^{-1} c`, variance
   solving :math:`P = T P T' + R Q R'` (DK §5.6.2).
:meth:`~ssfortran.Representation.initialize`
   The general :math:`N(a_1, P_* + \kappa P_\infty)`.
:meth:`~ssfortran.Representation.initialize_block`
   One of the above for a block of states.

A stationary initialization is recomputed from the current matrices at
every filter run, so it follows the parameters during estimation.

Blocks
------

Models usually mix kinds, for example a diffuse trend with a stationary
ARMA disturbance (DK §5.6):

>>> rep = ss.Representation(np.zeros(20), k_states=3)
>>> rep["transition"] = [[1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.7]]
>>> rep.initialize_block(0, 2, ss.INIT_DIFFUSE)        # level and slope
>>> rep.initialize_block(2, 3, ss.INIT_STATIONARY)     # AR(1) disturbance

A stationary block must not depend on states outside it through T. The
built-in models choose the kinds themselves.

Exact and approximate diffuse
-----------------------------

The exact initial Kalman filter (DK §5.2) handles :math:`\kappa \to \infty`
exactly: it carries :math:`P_\infty` until it vanishes, after d periods, and
the log likelihood is then the diffuse log likelihood (DK §7.2.2). The
approximate method runs the ordinary filter with a large :math:`\kappa`,
which loses precision as :math:`\kappa` grows and gives a likelihood that
depends on :math:`\kappa`. The built-in models use the exact method; see
:doc:`../design/diffuse_initialization`.

The diffuse log likelihood accounts for the diffuse states, so no
observations need to be left out. statsmodels' structural and ARIMA models
default to the approximate method and leave out the first d observations
(``loglikelihood_burn``); compare per-period contributions, ``llf_obs``,
with them.

Diffuse elements as regression coefficients
--------------------------------------------

The augmented filter (DK §5.7,
:meth:`~ssfortran.Representation.augmented`) treats :math:`\delta` as fixed
unknown coefficients and estimates them by generalized least squares. It
gives the same diffuse log likelihood, the estimate :math:`\hat\delta` with
its variance, and two further likelihoods: with :math:`\delta` fixed but
unknown (DK §7.2.4) and the marginal likelihood (DK §7.2.6), which
``marginal_likelihood = True`` makes the default:

>>> y = np.loadtxt("data/nile.csv", delimiter=",", skiprows=1)[:, 1]
>>> rep = ss.Representation(y, k_states=1)
>>> rep["design"] = rep["transition"] = rep["selection"] = [[1.0]]
>>> rep["obs_cov"], rep["state_cov"] = [[15099.0]], [[1469.1]]
>>> rep.initialize_diffuse()
>>> aug = rep.augmented()
>>> round(aug.llf, 4) == round(rep.loglike(), 4)
True
>>> aug.delta.shape
(1,)
