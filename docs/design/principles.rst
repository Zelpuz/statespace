Principles and testing
======================

Durbin and Koopman as the source of truth
-----------------------------------------

**Decision.** The algorithms, their names and their defaults follow Durbin
and Koopman (2012), Part I. Where statsmodels differs, we follow DK.

**Why.** One reference keeps the library consistent: the recursions, the
initialization, the transforms and the diagnostics come from the same
derivations, and a disagreement has one arbiter. DK are also the authors of
most of the methods.

**Alternatives.** Matching statsmodels in every detail would ease
comparisons, but would carry over its choices where DK differ (the variance
transform, the damped cycle) and its known errors (see
:doc:`../statsmodels_differences`).

**Where.** Comments in the source cite DK sections and equations, as do
these documents.

statsmodels as a reference, not an authority
--------------------------------------------

**Decision.** statsmodels' ``tsa.statespace`` provides test data: the
fixtures in ``test/fixtures`` are its output for the same models, written by
``test/fixtures/make_fixtures.py``. Where it follows DK the results must
agree, typically to 1e-9 in relative terms; where it does not, the tests
check our results another way.

**Why.** An independent implementation catches errors that internal checks
miss. Using it as the only check would import its errors: three were found
this way (see :doc:`../statsmodels_differences`).

Checks where no reference exists
--------------------------------

Several algorithms have no counterpart in statsmodels. Each is checked
against another route to the same answer:

* the dense matrix form of the model (DK §4.13) against the recursions;
* the exact diffuse filter against the approximate one as
  :math:`\kappa \to \infty`, and against the augmented filter (DK §5.7);
* the univariate against the conventional filter;
* the square root, fast, classical and two-filter smoothers against the
  state smoother;
* collapsing against the full model;
* the marginal likelihood against the invariance and the examples of
  Francke, Koopman and de Vos (2010);
* the analytic score against finite differences;
* EM at the maximum likelihood estimate as a fixed point;
* the worked examples of DK §5.6 in closed form;
* the smoothing spline against SciPy's, and least squares residuals against
  ordinary least squares;
* the illustrations of DK ch. 8 against DK's printed values.

Completeness over convenience
-----------------------------

**Decision.** Every algorithm of DK Part I is implemented, including those
DK present for insight rather than use: the Whittle recursion (DK §4.6.3)
loses accuracy on long series and the dense form (DK §4.13) costs
:math:`O((nm)^3)`. They are documented as such.

**Why.** The library serves as a companion to the book; the extra
algorithms also check the main ones.
