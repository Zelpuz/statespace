The steady state
================

A relative convergence test
---------------------------

**Decision.** In a time-invariant model the conventional filter stops
updating :math:`P_t`, :math:`F_t` and :math:`K_t` once
:math:`\|P_{t+1} - P_t\|_F \le \tau \|P_{t+1}\|_F`, with
:math:`\tau` = ``tol_steady`` = 1e-15 by default (DK §4.3.4). It resumes the
full recursion at a missing observation and never starts during the
diffuse periods.

**Why.** At :math:`\tau` = 1e-15, a few units in the last place, the
shortcut changes results by less than rounding; the tests compare it with
the full recursion, also with missing data, time-varying intercepts and a
concentrated scale. It reduces a step of the log likelihood to a few
matrix-vector products: 4-14 times faster in the benchmarks.

**Alternatives.** statsmodels switches when
:math:`\|P_{t+1} - P_t\|_F^2 < 10^{-19}`, an absolute threshold: for data
on a large scale it is never reached, and on a small scale it can be
reached early. This causes the differences near 1e-10 between the two
libraries in some time-invariant fixtures.

**Where.** ``check_steady``, ``use_steady`` and ``steady_step`` in
``statespace_filter``; ``tol_steady`` and ``t_steady``.

Computing the steady state
--------------------------

**Decision.** ``steady_state`` solves for :math:`\bar P` by the
structure-preserving doubling algorithm (Chu, Fan, Lin and Wang 2004), whose
k-th step equals step :math:`2^k` of the Riccati recursion from
:math:`P = 0`; with singular H it iterates the recursion.

**Why.** The recursion converges slowly exactly where the steady state is
of interest: a fixed regression coefficient is learned like 1/t, and a
nearly fixed seasonal (the DK §8.2 model has a seasonal variance of 5e-7)
over thousands of periods. Doubling reaches those limits in some 40 steps.

The prediction error variance
-----------------------------

**Decision.** ``prediction_error_variance`` returns the steady-state
:math:`\bar F`, the goodness-of-fit measure of DK §7.4.

**Note.** DK's printed value in §8.2, 0.00586717, is :math:`\bar F` times
(n - d)/n; see :doc:`../examples/dk_8_2`.
