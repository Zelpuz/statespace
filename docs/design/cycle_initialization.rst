Initialization of a damped cycle
================================

**Decision.** A damped cycle (:math:`\rho < 1`) starts from its
unconditional distribution, :math:`N(0, \sigma^2_\kappa / (1 - \rho^2) I)`.

**Why.** A damped cycle is stationary, and DK §3.2.4 give this distribution
for it; DK §5.6 initialize stationary components from their unconditional
distribution.

**Alternatives.** statsmodels' ``UnobservedComponents`` starts every state,
the damped cycle included, as diffuse. Both are valid models but give
different log likelihoods: in one of our test fixtures, -65.72 with the
stationary start against -63.60 with the diffuse one. To compare with
statsmodels, the fixture re-initializes its cycle as stationary.

**Where.** ``cycle_init_blocks`` in ``statespace_structural``.
