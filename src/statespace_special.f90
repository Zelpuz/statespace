!> Special functions for test p-values: the regularized incomplete gamma and
!> beta functions, giving chi-square and F tail probabilities.
module statespace_special
  use statespace_kinds, only: dp
  implicit none
  private

  public :: gammp, gammq, betai, chi2_sf, f_cdf

  integer, parameter :: max_iter = 10000
  real(dp), parameter :: eps = 1.0e-15_dp
  real(dp), parameter :: fpmin = tiny(1.0_dp) / eps

contains

  !> P(chi2(df) > x).
  real(dp) function chi2_sf(x, df)
    real(dp), intent(in) :: x, df

    chi2_sf = gammq(0.5_dp * df, 0.5_dp * x)
  end function chi2_sf

  !> P(F(d1, d2) <= x).
  real(dp) function f_cdf(x, d1, d2)
    real(dp), intent(in) :: x, d1, d2

    if (x <= 0.0_dp) then
      f_cdf = 0.0_dp
    else
      f_cdf = betai(0.5_dp * d1, 0.5_dp * d2, d1 * x / (d1 * x + d2))
    end if
  end function f_cdf

  !> Regularized lower incomplete gamma function P(a, x).
  real(dp) function gammp(a, x)
    real(dp), intent(in) :: a, x

    if (x <= 0.0_dp) then
      gammp = 0.0_dp
    else if (x < a + 1.0_dp) then
      gammp = gamma_series(a, x)
    else
      gammp = 1.0_dp - gamma_cfrac(a, x)
    end if
  end function gammp

  !> Regularized upper incomplete gamma function Q(a, x) = 1 - P(a, x),
  !> computed directly in the upper tail for accuracy.
  real(dp) function gammq(a, x)
    real(dp), intent(in) :: a, x

    if (x <= 0.0_dp) then
      gammq = 1.0_dp
    else if (x < a + 1.0_dp) then
      gammq = 1.0_dp - gamma_series(a, x)
    else
      gammq = gamma_cfrac(a, x)
    end if
  end function gammq

  !> P(a, x) by its series.
  real(dp) function gamma_series(a, x) result(p)
    real(dp), intent(in) :: a, x
    real(dp) :: ap, del, total
    integer :: n

    ap = a
    total = 1.0_dp / a
    del = total
    do n = 1, max_iter
      ap = ap + 1.0_dp
      del = del * x / ap
      total = total + del
      if (abs(del) < abs(total) * eps) exit
    end do
    p = total * exp(-x + a * log(x) - log_gamma(a))
  end function gamma_series

  !> Q(a, x) by its continued fraction (modified Lentz).
  real(dp) function gamma_cfrac(a, x) result(q)
    real(dp), intent(in) :: a, x
    real(dp) :: b, c, d, h, an, del
    integer :: i

    b = x + 1.0_dp - a
    c = 1.0_dp / fpmin
    d = 1.0_dp / b
    h = d
    do i = 1, max_iter
      an = -i * (i - a)
      b = b + 2.0_dp
      d = an * d + b
      if (abs(d) < fpmin) d = fpmin
      c = b + an / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      del = d * c
      h = h * del
      if (abs(del - 1.0_dp) < eps) exit
    end do
    q = exp(-x + a * log(x) - log_gamma(a)) * h
  end function gamma_cfrac

  !> Regularized incomplete beta function I_x(a, b).
  real(dp) function betai(a, b, x)
    real(dp), intent(in) :: a, b, x
    real(dp) :: bt

    if (x <= 0.0_dp) then
      betai = 0.0_dp
      return
    else if (x >= 1.0_dp) then
      betai = 1.0_dp
      return
    end if
    bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * log(x) &
             + b * log(1.0_dp - x))
    if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
      betai = bt * beta_cfrac(a, b, x) / a
    else
      betai = 1.0_dp - bt * beta_cfrac(b, a, 1.0_dp - x) / b
    end if
  end function betai

  !> Continued fraction for the incomplete beta function (modified Lentz).
  real(dp) function beta_cfrac(a, b, x) result(h)
    real(dp), intent(in) :: a, b, x
    real(dp) :: qab, qap, qam, c, d, aa, del
    integer :: m, m2

    qab = a + b
    qap = a + 1.0_dp
    qam = a - 1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab * x / qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp / d
    h = d
    do m = 1, max_iter
      m2 = 2 * m
      aa = m * (b - m) * x / ((qam + m2) * (a + m2))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      h = h * d * c
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      del = d * c
      h = h * del
      if (abs(del - 1.0_dp) < eps) exit
    end do
  end function beta_cfrac
end module statespace_special
