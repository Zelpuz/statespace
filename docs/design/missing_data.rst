Missing observations
====================

**Decision.** Any subset of a period's observations may be missing (NaN).
The filter uses the observed elements (DK §4.10). For a missing element of
:math:`y_t`, the smoothed disturbance and its variance are the moments of
that element of :math:`\varepsilon_t` given the observed data,

.. math::

   E(\varepsilon_t | Y_n) = H_t u_t, \qquad
   Var(\varepsilon_t | Y_n) = H_t - H_t D_t H_t,

which are not zero and :math:`H_t` when :math:`H_t` has correlations. The
simulation smoothers draw missing elements from the same conditional
distribution. :math:`F_t` covers all elements, missing or not;
:math:`F_t^{-1}` covers the observed ones, zero-padded.

**Why.** EM (DK §7.3.4) and the analytic score (DK §7.3.3) take expectations
over the disturbances given the data; with correlated H, the observed
elements carry information about the missing ones, and the correct moments
are needed for both to be right. Zero-padding :math:`F_t^{-1}` keeps the
full-size smoother formulas exact without selecting rows of Z at each step.

**Alternatives.** statsmodels reports 0 and :math:`H_t` for missing
elements and draws them unconditionally in its simulation smoother. The
two agree on observed elements, which the tests compare.

**Standardized residuals** are NaN for missing elements and in the diffuse
periods, where statsmodels reports 0, so that tests on the residuals skip
them rather than count zeros.

**Where.** ``observed_inverse`` in ``statespace_filter``; the smoother's
disturbance step.
