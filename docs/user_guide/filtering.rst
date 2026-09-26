Filtering and smoothing
=======================

The Kalman filter
-----------------

The filter (DK §4.3) computes, for t = 1, ..., n, the one-step prediction
error and its variance,

.. math::

   v_t = y_t - d_t - Z_t a_t, \qquad F_t = Z_t P_t Z_t' + H_t,

the filtered state (DK eq. 4.24)

.. math::

   a_{t|t} = a_t + P_t Z_t' F_t^{-1} v_t, \qquad
   P_{t|t} = P_t - P_t Z_t' F_t^{-1} Z_t P_t,

and the prediction :math:`a_{t+1} = c_t + T_t a_{t|t}`,
:math:`P_{t+1} = T_t P_{t|t} T_t' + R_t Q_t R_t'`. The log likelihood is
(DK eq. 7.2)

.. math::

   \log L = -\frac{np}{2} \log 2\pi
            - \frac12 \sum_t \big(\log|F_t| + v_t' F_t^{-1} v_t\big).

>>> y = np.loadtxt("data/nile.csv", delimiter=",", skiprows=1)[:, 1]
>>> rep = ss.Representation(y, k_states=1)
>>> rep["design"] = rep["transition"] = rep["selection"] = [[1.0]]
>>> rep["obs_cov"], rep["state_cov"] = [[15099.0]], [[1469.1]]
>>> rep.initialize_diffuse()
>>> f = rep.filter()
>>> f.nobs_diffuse, round(f.llf, 4)
(1, -633.4646)

:meth:`~ssfortran.Representation.loglike` returns the same number without
storing the filter output; estimation uses it.

Univariate treatment
--------------------

With ``filter_method = FILTER_UNIVARIATE`` the filter processes the elements
of :math:`y_t` one at a time (DK §6.4). It needs no matrix inverse and
handles singular :math:`F_t`. When :math:`H_t` is not diagonal, each
period's observed block is first factorized, :math:`H_{oo} = L D L'`, and
the observation equation transformed by :math:`L^{-1}` (DK §6.4.3); this
works with time-varying H and any pattern of missing values.

The reported :math:`v_t` and :math:`F_t` are in the original coordinates in
every period; see :doc:`../design/output_coordinates`.

Diffuse periods
---------------

With diffuse states the exact initial filter (DK §5.2) runs until
:math:`P_\infty` vanishes, after ``nobs_diffuse`` periods. By default it
processes those periods element by element (DK §5.2.5).
``diffuse_method = DIFFUSE_MULTIVARIATE`` uses the multivariate form where
:math:`F_\infty` is zero or nonsingular. The per-period choice is in the
filter output's method array (Fortran).

Steady state
------------

In a time-invariant model :math:`P_t` converges. Once the relative change
falls below ``tol_steady`` (default 1e-15), the filter holds :math:`P_t`,
:math:`F_t` and :math:`K_t` fixed and updates only the state means
(DK §4.3.4), until a missing observation. ``t_steady`` reports where this
began. :meth:`~ssfortran.Representation.steady_state` computes the limit
directly:

>>> P, F = rep.steady_state()
>>> round(float(F[0, 0]), 2)
20600.26

See :doc:`../design/steady_state`.

State and disturbance smoothing
-------------------------------

The smoother (DK §4.4-4.5) runs backwards through

.. math::

   r_{t-1} = Z_t' F_t^{-1} v_t + L_t' r_t, \qquad
   N_{t-1} = Z_t' F_t^{-1} Z_t + L_t' N_t L_t,

with :math:`L_t = T_t - K_t Z_t`, and gives the smoothed state
:math:`\hat\alpha_t = a_t + P_t r_{t-1}` (DK eq. 4.39), its variance
:math:`V_t = P_t - P_t N_{t-1} P_t` (DK eq. 4.43), and the smoothed
disturbances with their variances.

>>> sm = rep.smooth()
>>> sm.smoothed_state.shape, sm.smoothed_state_cov.shape
((1, 100), (1, 1, 100))

Other smoothers
---------------

DK §4.4-4.8 give further algorithms; each is a method of
:class:`~ssfortran.Representation`:

* :meth:`~ssfortran.Representation.fast_smoother` (DK §4.6.2): means only;
* :meth:`~ssfortran.Representation.classical_smoother` (DK §4.6.1),
  :meth:`~ssfortran.Representation.two_filter_smoother` (DK §4.6.4),
  :meth:`~ssfortran.Representation.whittle_smoother` (DK §4.6.3);
* :meth:`~ssfortran.Representation.fixed_point_smoother` and
  :meth:`~ssfortran.Representation.fixed_lag_smoother` (DK §4.4.6), and
  :meth:`~ssfortran.Representation.update_smoothed` (DK §4.4.5) for
  observations that arrive later;
* :meth:`~ssfortran.Representation.smoothed_state_cov_between` and
  :meth:`~ssfortran.Representation.smoothed_state_autocov` (DK §4.7);
* :meth:`~ssfortran.Representation.filtered_state_weights` and
  :meth:`~ssfortran.Representation.smoothed_state_weights` (DK §4.8);
* ``smooth(method="sqrt")``, the square root filter and smoother (DK §6.3).

They agree with the standard smoother, which the tests check; the Whittle
recursion loses accuracy on long series, as DK note.

>>> np.allclose(rep.fast_smoother(), sm.smoothed_state)
True
