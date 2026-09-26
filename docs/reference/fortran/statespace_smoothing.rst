``statespace_smoothing``
========================

Further smoothing results (DK §4.4.5-4.8): updating smoothed estimates,
fixed-point and fixed-lag smoothing, other state smoothers, covariances of
the smoothed state across time, and filtering and smoothing weights.

The recursions are written in terms of the conventional :math:`K_t` and
:math:`F_t^{-1}`; ``conventional_gain`` forms them in any period after the
diffuse periods, whichever filter method was used. Unless stated, the
routines are not defined with diffuse states (``info`` is
``SS_ERR_UNSUPPORTED``).

Other state smoothers
---------------------

All give the same :math:`\hat\alpha_t` as ``state_smoother``; they illustrate
the algorithms of DK §4.6 and serve as checks.

.. code-block:: fortran

   subroutine fast_state_smoother(rep, fres, alphahat, info)
   subroutine classical_state_smoother(rep, fres, alphahat, V, info)
   subroutine two_filter_smoother(rep, fres, alphahat, V, info)
   subroutine whittle_smoother(rep, fres, alphahat, info)

``fast_state_smoother``
   DK §4.6.2: a backward pass for :math:`r_t` only, then
   :math:`\hat\alpha_1 = a_1 + P_1 r_0` (plus
   :math:`P_{\infty,1} r_0^{(1)}` with diffuse states) and
   :math:`\hat\alpha_{t+1} = c_t + T_t \hat\alpha_t + R_t Q_t R_t' r_t`.
   Defined with diffuse states.
``classical_state_smoother``
   The Rauch-Tung-Striebel form (DK §4.6.1):
   :math:`J_t = P_{t|t} T_t' P_{t+1}^{-1}`,
   :math:`\hat\alpha_t = a_{t|t} + J_t (\hat\alpha_{t+1} - a_{t+1})`,
   :math:`V_t = P_{t|t} + J_t (V_{t+1} - P_{t+1}) J_t'`, with
   :math:`P_{t+1}^{-1}` applied as a pseudo-inverse.
``two_filter_smoother``
   DK §4.6.4: a backward information filter for
   :math:`y_t, \dots, y_n` combined with the forward filter.
``whittle_smoother``
   DK §4.6.3: the backward recursion from the first-order conditions of the
   joint density of :math:`\alpha` and Y, started from
   :math:`\hat\alpha_n = a_{n|n}`. It needs :math:`T_t`, :math:`R_t Q_t R_t'`
   and :math:`H_t` nonsingular and, as DK note, loses accuracy
   geometrically going back; usable for short series only.

Updating, fixed-point and fixed-lag smoothing
---------------------------------------------

.. code-block:: fortran

   subroutine update_smoothed(rep, fres, n0, alphahat, V, info)
     real(dp), intent(inout) :: alphahat(:, :), V(:, :, :)    ! (m, n), (m, m, n)
   subroutine fixed_point_smoother(rep, fres, t, path_a, path_V, info)
     real(dp), intent(out) :: path_a(:, :), path_V(:, :, :)   ! (m, n-t+1), (m, m, n-t+1)
   subroutine fixed_lag_smoother(rep, fres, j, alphahat, V, info)
     real(dp), intent(out) :: alphahat(:, :), V(:, :, :)      ! (m, n), (m, m, n)

``update_smoothed``
   DK §4.4.5: columns 1..n0 hold estimates given :math:`y_1, \dots,
   y_{n_0}` on entry and given all n observations on exit, by the
   recursions (DK eqs. 4.49-4.50) for each new observation.
``fixed_point_smoother``
   DK §4.4.6: :math:`\hat\alpha_{t|k}` and :math:`V_{t|k}` for fixed t and
   k = t, ..., n, in column k - t + 1.
``fixed_lag_smoother``
   DK §4.4.6: :math:`\hat\alpha_{s|s+j}` and :math:`V_{s|s+j}` (DK eq. 4.51)
   for fixed j, in column s = 1, ..., n - j; later columns are NaN.

Covariances across time
-----------------------

.. code-block:: fortran

   subroutine smoothed_state_cov_between(rep, fres, sres, t, j, C, info)
     real(dp), intent(out) :: C(:, :)                 ! (m, m)
   subroutine smoothed_state_autocov(rep, fres, sres, acov, info)
     real(dp), intent(out) :: acov(:, :, :)           ! (m, m, n-1)

``smoothed_state_cov_between``
   :math:`Cov(\alpha_t, \alpha_j | Y_n)` (DK §4.7); for j > t,
   :math:`P_t L_t' \cdots L_{j-1}' (I - N_{j-1} P_j)`. Periods t, ..., j-1
   must be past the diffuse periods.
``smoothed_state_autocov``
   ``acov(:, :, t)`` is :math:`Cov(\alpha_{t+1}, \alpha_t | Y_n)`; NaN where
   the diffuse periods are involved.

Weights
-------

.. code-block:: fortran

   subroutine filtered_state_weights(rep, fres, Wa, Watt, info)
     real(dp), intent(out) :: Wa(:, :, :, :), Watt(:, :, :, :)   ! (m, p, n, n)
   subroutine smoothed_state_weights(rep, W, C, A, info, j_list)
     real(dp), intent(out) :: W(:, :, :, :)     ! (m, p, n, n)
     real(dp), intent(out) :: C(:, :, :, :)     ! (m, m, n, n)
     real(dp), intent(out) :: A(:, :, :)        ! (m, m, n)
     integer, intent(in), optional :: j_list(:)

``filtered_state_weights``
   DK §4.8.2 (Table 4.5): ``Wa(:, :, t, j)`` is the weight of :math:`y_j`
   in :math:`a_t`, :math:`L_{t-1} \cdots L_{j+1} K_j` for j < t;
   ``Watt`` the same for :math:`a_{t|t}`. Intercepts and the initial mean
   aside.
``smoothed_state_weights``
   DK §4.8.3, extended to the intercepts and the prior mean:
   :math:`\hat\alpha_t = \sum_j W_{tj} (y_j - d_j) + \sum_j C_{tj} c_j +
   A_t a_1`. Defined in diffuse periods. Costs one smoother run per
   observation; ``j_list`` restricts the computation to some j.

Building blocks
---------------

.. code-block:: fortran

   subroutine innovation_transition(rep, fres, t, L, info)
   subroutine conventional_gain(rep, fres, t, K, Finv, info)

``innovation_transition``
   :math:`L_t` with :math:`a_{t+1} - \alpha_{t+1} = L_t (a_t - \alpha_t) +
   \dots`: :math:`T_t - K_t Z_t` in conventional periods and
   :math:`T_t L_{t,p} \cdots L_{t,1}`, :math:`L_{t,i} = I - K_{t,i} Z^*_i`,
   in univariate ones (DK §4.3, §6.4).
``conventional_gain``
   :math:`K_t = T_t P_t Z_t' F_t^{-1}` and :math:`F_t^{-1}` over the
   observed elements for a period after the diffuse periods, formed from the
   stored :math:`P_t` and :math:`F_t` where the filter did not store them.
