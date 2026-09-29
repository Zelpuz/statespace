``statespace_arima``
====================

ARMA and ARIMA components (DK §3.4, §3.6.2, §5.6.2-5.6.4).

The stationary part :math:`y^*_t = \Delta^d \Delta_s^D y_t` follows the
ARMA model :math:`\phi(L)\Phi(L^s) y^*_t = \theta(L)\Theta(L^s) \zeta_t`,
put in DK's form (DK eqs. 3.19-3.20): with
:math:`r = \max(p + sP, q + sQ + 1)`,

.. math::

   T = \begin{bmatrix} \phi^* & I \\ & 0 \end{bmatrix}, \quad
   R = (1, \theta^*_1, \dots, \theta^*_{r-1})', \quad Z = e_1',

where :math:`\phi^*, \theta^*` are the coefficients of the multiplied
polynomials. For d > 0 the state is extended by the cumulative differences
(DK §3.4), and for D = 1 by the s lags of :math:`\Delta^d y`. The
differenced states are diffuse and the ARMA block stationary (DK §5.6.3).
The layout equals statsmodels' SARIMAX, so the states can be compared one
by one.

A constant or regressors, as in regression with ARMA errors (DK §3.6.2,
§5.6.4), are added as a ``regression_t`` component.

``arima_t``
-----------

.. code-block:: fortran

   type, extends(component_t) :: arima_t
     integer :: ar = 0, d = 0, ma = 0            ! p, d, q
     integer :: sar = 0, sd = 0, sma = 0         ! P, D, Q (D at most 1)
     integer :: s = 0                            ! seasonal period
     integer :: series = 1
     logical :: enforce_stationarity = .true.
     logical :: enforce_invertibility = .true.
   end type

Parameters: AR (p), seasonal AR (P), MA (q), seasonal MA (Q), then the
variance of :math:`\zeta_t`, named ``ar.L1, ..., ar.S.L12, ..., ma.L1, ...,
sigma2``. AR polynomials are kept stationary and MA polynomials invertible
by Monahan's (1984) transform, as in statsmodels; the variance uses
:math:`\exp(2\psi)`.

Example
-------

An ARMA(1, 1) with a constant (DK §5.6.4):

.. code-block:: fortran

   type(component_holder_t) :: comps(2)
   comps(1)%c = regression_t(x=reshape(spread(1.0_dp, 1, n), [n, 1]))
   comps(2)%c = arima_t(ar=1, ma=1)
   model = structural_model(y, comps, info)
