``statespace_special``
======================

Special functions for the p-values of the diagnostic tests.

.. code-block:: fortran

   real(dp) function gammp(a, x)       ! regularized lower incomplete gamma P(a, x)
   real(dp) function gammq(a, x)       ! Q(a, x) = 1 - P(a, x), accurate in the upper tail
   real(dp) function betai(a, b, x)    ! regularized incomplete beta I_x(a, b)
   real(dp) function chi2_sf(x, df)    ! P(chi2(df) > x)
   real(dp) function f_cdf(x, d1, d2)  ! P(F(d1, d2) <= x)

The incomplete gamma function uses its power series below
:math:`x = a + 1` and its continued fraction above, the incomplete beta
function its continued fraction; the continued fractions are evaluated by
the modified Lentz method (Thompson and Barnett 1986). The results match
SciPy to 1e-12 in the tests.
