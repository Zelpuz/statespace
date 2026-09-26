Output coordinates of the univariate treatment
==============================================

**Decision.** In periods processed element by element (DK §6.4), the
filter reports :math:`v_t`, :math:`F_t` and :math:`F_{\infty,t}` in the
original coordinates of :math:`y_t`, as in conventional periods. The
per-element quantities of the univariate recursion, in the coordinates
transformed by :math:`L^{-1}` from :math:`H_{oo} = L D L'`, are kept
separately (the ``uv_*`` arrays of ``filter_result_t``) for the smoother.
Smoothed observation disturbances are also in original coordinates, as
full matrices: from the identity
:math:`\varepsilon_{t,o} = y_{t,o} - d_{t,o} - Z_{t,o}\alpha_t`,
:math:`\hat\varepsilon_{t,o} = y_{t,o} - d_{t,o} - Z_{t,o}\hat\alpha_t` with
variance :math:`Z_{t,o} V_t Z_{t,o}'`.

**Why.** A user reading :math:`v_t` and :math:`F_t` should not have to know
which periods the filter treated element by element, which depends on the
initialization and on the filter options. Original coordinates make the
output of every period comparable, and make standardized residuals and
diagnostics correct across the diffuse periods.

**Alternatives.** statsmodels reports these quantities in the transformed
coordinates during diffuse and univariate periods, with diagonal
disturbance variances only. For p = 1 the two coincide; the tests compare
the transformed quantities (our ``uv_*`` arrays) with statsmodels' in those
periods and the original ones elsewhere.

**Cost.** Computing :math:`F_t` in original coordinates adds a matrix
product per univariate period; this is why the stored univariate filter
output is slower than statsmodels' (see :doc:`performance`), while the
log likelihood, which does not store it, is faster.

**Where.** ``observation_moments`` and ``store_elements`` in
``statespace_filter``; ``measurement_disturbance`` in
``statespace_smoother``.
