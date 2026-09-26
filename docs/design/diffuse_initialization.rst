Diffuse initialization
======================

Exact diffuse by default
------------------------

**Decision.** The built-in models initialize nonstationary states as exact
diffuse (DK §5.2) and report the diffuse log likelihood (DK §7.2.2), with no
observations left out.

**Why.** The exact treatment is the limit the approximate one only
approaches: with a large finite :math:`\kappa` the filter loses precision
as :math:`\kappa` grows and the log likelihood depends on :math:`\kappa`. DK
recommend the exact treatment.

**Alternatives.** statsmodels' structural and ARIMA models default to an
approximate diffuse variance (1e6, or 1e10 with regression states) and
leave the first d observations out of the log likelihood
(``loglikelihood_burn``). The approximate method remains available through
``initialize_approximate_diffuse``. Tests compare per-period contributions,
which do not depend on the burn-in.

**Where.** ``init_blocks`` of each component; ``statespace_filter``.

Univariate treatment of the diffuse periods
-------------------------------------------

**Decision.** The diffuse periods are filtered element by element
(DK §5.2.5 with §6.4) by default; the multivariate exact initial filter
(DK §5.2) is an option, with the univariate step where :math:`F_\infty` is
singular but nonzero.

**Why.** The univariate form needs no inverse of :math:`F_\infty`, which is
often singular, and works for any pattern of missing values. The
multivariate form is DK's main presentation and is kept for completeness
and as a check.

**Tolerances.** :math:`F_\infty` and :math:`\|P_\infty\|_F^2` below 1e-10
count as zero, as in statsmodels, and elements with :math:`F_*` below 1e-10
are skipped in diffuse periods.

The marginal likelihood
-----------------------

**Decision.** ``marginal_likelihood`` reports the marginal likelihood of
Francke, Koopman and de Vos (2010) (DK §7.2.6):
:math:`\log L_M = \log L_d + \tfrac{k}{2}\log 2\pi + \tfrac12 \log|S_*|`,
with :math:`S_*` from their recursion (21).

**Why.** The diffuse likelihood is not invariant to how the diffuse
directions are parameterized when parameters enter Z or T; the marginal
likelihood is. The correction depends on Z, T and the diffuse directions
only, so for variance parameters the estimates are the same.

**Where.** ``marginal_correction`` in ``statespace_filter``; ``llf_marginal``
of the augmented filter.
