!> Thin wrappers over BLAS/LAPACK.
!>
!> All wrappers take assumed-shape arrays and pass them straight to BLAS, so
!> callers must pass contiguous arrays (whole arrays or `x(:,:,t)` slices).
!> `gemm` and `gemv` are declared `contiguous`, and for small operands (the
!> common case in state space models) use inline loops instead of BLAS,
!> whose per-call overhead dominates at those sizes.
module statespace_linalg
  use statespace_kinds, only: dp, SS_OK, SS_ERR_NOT_PD, SS_ERR_SINGULAR, &
                              SS_ERR_NOT_STATIONARY
  implicit none
  private

  public :: gemm, gemv, chol_inv, symmetrize, eye, solve, solve_lyapunov
  public :: ldl_psd, solve_unit_lower, is_diagonal, psd_solve, tria

  !> Operands with m n k (gemm) or m n (gemv) up to this size use inline loops.
  integer, parameter :: small_gemm = 1000, small_gemv = 400
  !> chol_inv uses inline code up to this order.
  integer, parameter :: small_chol = 16

  interface
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      import :: dp
      character, intent(in) :: transa, transb
      integer, intent(in) :: m, n, k, lda, ldb, ldc
      real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
      real(dp), intent(inout) :: c(ldc, *)
    end subroutine dgemm

    subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import :: dp
      character, intent(in) :: trans
      integer, intent(in) :: m, n, lda, incx, incy
      real(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
      real(dp), intent(inout) :: y(*)
    end subroutine dgemv

    subroutine dpotrf(uplo, n, a, lda, info)
      import :: dp
      character, intent(in) :: uplo
      integer, intent(in) :: n, lda
      real(dp), intent(inout) :: a(lda, *)
      integer, intent(out) :: info
    end subroutine dpotrf

    subroutine dpotri(uplo, n, a, lda, info)
      import :: dp
      character, intent(in) :: uplo
      integer, intent(in) :: n, lda
      real(dp), intent(inout) :: a(lda, *)
      integer, intent(out) :: info
    end subroutine dpotri

    subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      import :: dp
      character, intent(in) :: side, uplo, transa, diag
      integer, intent(in) :: m, n, lda, ldb
      real(dp), intent(in) :: alpha, a(lda, *)
      real(dp), intent(inout) :: b(ldb, *)
    end subroutine dtrsm

    subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
      import :: dp
      integer, intent(in) :: m, n, lda, lwork
      real(dp), intent(inout) :: a(lda, *)
      real(dp), intent(out) :: tau(*), work(*)
      integer, intent(out) :: info
    end subroutine dgeqrf

    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import :: dp
      integer, intent(in) :: n, nrhs, lda, ldb
      real(dp), intent(inout) :: a(lda, *), b(ldb, *)
      integer, intent(out) :: ipiv(*), info
    end subroutine dgesv
  end interface

contains

  !> C <- alpha * op(A) * op(B) + beta * C, with op(X) = X or X'. As in
  !> BLAS, C is not read when beta = 0.
  subroutine gemm(transa, transb, alpha, A, B, beta, C)
    character, intent(in) :: transa, transb
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(in), contiguous :: A(:, :), B(:, :)
    real(dp), intent(inout), contiguous :: C(:, :)
    integer :: k, i, j, l

    if (transa == 'N') then
      k = size(A, 2)
    else
      k = size(A, 1)
    end if
    if (size(C, 1) * size(C, 2) * k > small_gemm) then
      call dgemm(transa, transb, size(C, 1), size(C, 2), k, alpha, A, &
                 max(1, size(A, 1)), B, max(1, size(B, 1)), beta, C, max(1, size(C, 1)))
      return
    end if

    if (beta == 0.0_dp) then
      C = 0.0_dp
    else if (beta /= 1.0_dp) then
      C = beta * C
    end if
    if (transa == 'N') then
      if (transb == 'N') then
        do j = 1, size(C, 2)
          do l = 1, k
            C(:, j) = C(:, j) + (alpha * B(l, j)) * A(:, l)
          end do
        end do
      else
        do j = 1, size(C, 2)
          do l = 1, k
            C(:, j) = C(:, j) + (alpha * B(j, l)) * A(:, l)
          end do
        end do
      end if
    else
      if (transb == 'N') then
        do j = 1, size(C, 2)
          do i = 1, size(C, 1)
            C(i, j) = C(i, j) + alpha * dot_product(A(:, i), B(:, j))
          end do
        end do
      else
        do j = 1, size(C, 2)
          do i = 1, size(C, 1)
            C(i, j) = C(i, j) + alpha * dot_product(A(:, i), B(j, :))
          end do
        end do
      end if
    end if
  end subroutine gemm

  !> y <- alpha * op(A) * x + beta * y. As in BLAS, y is not read when
  !> beta = 0.
  subroutine gemv(trans, alpha, A, x, beta, y)
    character, intent(in) :: trans
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(in), contiguous :: A(:, :), x(:)
    real(dp), intent(inout), contiguous :: y(:)
    integer :: i, j

    if (size(A, 1) * size(A, 2) > small_gemv) then
      call dgemv(trans, size(A, 1), size(A, 2), alpha, A, max(1, size(A, 1)), x, 1, &
                 beta, y, 1)
      return
    end if

    if (beta == 0.0_dp) then
      y = 0.0_dp
    else if (beta /= 1.0_dp) then
      y = beta * y
    end if
    if (trans == 'N') then
      do j = 1, size(A, 2)
        y = y + (alpha * x(j)) * A(:, j)
      end do
    else
      do i = 1, size(A, 2)
        y(i) = y(i) + alpha * dot_product(A(:, i), x)
      end do
    end if
  end subroutine gemv

  !> Invert a symmetric positive definite matrix in place via Cholesky and
  !> return log|A|. Both triangles of the result are filled.
  subroutine chol_inv(A, logdet, info)
    real(dp), intent(inout), contiguous :: A(:, :)
    real(dp), intent(out) :: logdet
    integer, intent(out) :: info
    integer :: i, j, n

    n = size(A, 1)
    info = SS_OK
    logdet = 0.0_dp
    if (n == 0) return
    if (n == 1) then
      if (.not. (A(1, 1) > 0.0_dp)) then
        info = SS_ERR_NOT_PD
        return
      end if
      logdet = log(A(1, 1))
      A(1, 1) = 1.0_dp / A(1, 1)
      return
    end if

    if (n <= small_chol) then
      call chol_inv_small(A, logdet, info)
      return
    end if

    call dpotrf('L', n, A, n, i)
    if (i /= 0) then
      info = SS_ERR_NOT_PD
      return
    end if
    do i = 1, n
      logdet = logdet + log(A(i, i))
    end do
    logdet = 2.0_dp * logdet
    call dpotri('L', n, A, n, i)
    if (i /= 0) then
      info = SS_ERR_NOT_PD
      return
    end if
    do j = 2, n
      do i = 1, j - 1
        A(i, j) = A(j, i)
      end do
    end do
  end subroutine chol_inv

  !> chol_inv for small matrices, without LAPACK's per-call overhead:
  !> A = L L', W = L^-1, A^-1 = W' W.
  subroutine chol_inv_small(A, logdet, info)
    real(dp), intent(inout), contiguous :: A(:, :)
    real(dp), intent(out) :: logdet
    integer, intent(out) :: info
    real(dp) :: L(size(A, 1), size(A, 1)), W(size(A, 1), size(A, 1)), d
    integer :: i, j, n

    n = size(A, 1)
    info = SS_OK
    logdet = 0.0_dp
    L = 0.0_dp
    do j = 1, n
      d = A(j, j) - sum(L(j, 1:j - 1)**2)
      if (.not. (d > 0.0_dp)) then
        info = SS_ERR_NOT_PD
        return
      end if
      L(j, j) = sqrt(d)
      logdet = logdet + log(d)
      do i = j + 1, n
        L(i, j) = (A(i, j) - sum(L(i, 1:j - 1) * L(j, 1:j - 1))) / L(j, j)
      end do
    end do
    W = 0.0_dp
    do j = 1, n
      W(j, j) = 1.0_dp / L(j, j)
      do i = j + 1, n
        W(i, j) = -sum(L(i, j:i - 1) * W(j:i - 1, j)) / L(i, i)
      end do
    end do
    do j = 1, n
      do i = j, n
        A(i, j) = sum(W(i:n, i) * W(i:n, j))
        A(j, i) = A(i, j)
      end do
    end do
  end subroutine chol_inv_small

  !> Solve A X = B for square A. B is overwritten with X; A is not modified.
  subroutine solve(A, B, info)
    real(dp), intent(in), contiguous :: A(:, :)
    real(dp), intent(inout), contiguous :: B(:, :)
    integer, intent(out) :: info
    real(dp), allocatable :: LU(:, :)
    integer, allocatable :: ipiv(:)
    integer :: n, stat

    n = size(A, 1)
    LU = A
    allocate (ipiv(n))
    call dgesv(n, size(B, 2), LU, max(1, n), ipiv, B, max(1, n), stat)
    info = SS_OK
    if (stat /= 0) info = SS_ERR_SINGULAR
  end subroutine solve

  !> Solve the discrete Lyapunov equation P = T P T' + V by doubling:
  !>
  !>     P_0 = V, A_0 = T;  P_k+1 = P_k + A_k P_k A_k',  A_k+1 = A_k A_k
  !>
  !> so that P_k = sum_{j < 2^k} T^j V T'^j. Fails with SS_ERR_NOT_STATIONARY
  !> if T has an eigenvalue on or outside the unit circle.
  subroutine solve_lyapunov(T, V, P, info)
    real(dp), intent(in), contiguous :: T(:, :), V(:, :)
    real(dp), intent(out), contiguous :: P(:, :)
    integer, intent(out) :: info
    integer, parameter :: max_doublings = 100
    real(dp), allocatable :: A(:, :), AP(:, :), Anext(:, :)
    integer :: k, m

    m = size(T, 1)
    allocate (AP(m, m), Anext(m, m))
    A = T
    P = V
    info = SS_OK
    do k = 1, max_doublings
      ! P <- P + A P A'
      call gemm('N', 'N', 1.0_dp, A, P, 0.0_dp, AP)
      call gemm('N', 'T', 1.0_dp, AP, A, 1.0_dp, P)
      call symmetrize(P)
      ! A <- A A. Once A = T^(2^k) is negligible, the remaining terms of the
      ! sum are O(eps^2) relative to P.
      call gemm('N', 'N', 1.0_dp, A, A, 0.0_dp, Anext)
      A = Anext
      if (maxval(abs(A)) <= epsilon(1.0_dp)) return
      if (.not. (maxval(abs(A)) < huge(1.0_dp)**0.25_dp)) exit
    end do
    info = SS_ERR_NOT_STATIONARY
  end subroutine solve_lyapunov

  !> A = L D L' for symmetric positive semi-definite A, with L unit lower
  !> triangular and D >= 0 diagonal. Pivots below a relative tolerance are set
  !> to zero together with the matching column of L, so singular A (e.g. some
  !> exact observations) is allowed.
  subroutine ldl_psd(A, L, D)
    real(dp), intent(in), contiguous :: A(:, :)
    real(dp), intent(out), contiguous :: L(:, :), D(:)
    real(dp) :: tol
    integer :: i, j, n

    n = size(A, 1)
    L = 0.0_dp
    tol = n * epsilon(1.0_dp) * maxval([(abs(A(i, i)), i=1, n), 0.0_dp])
    do j = 1, n
      L(j, j) = 1.0_dp
      D(j) = A(j, j) - sum(L(j, 1:j - 1)**2 * D(1:j - 1))
      if (D(j) <= tol) then
        D(j) = 0.0_dp
        cycle
      end if
      do i = j + 1, n
        L(i, j) = (A(i, j) - sum(L(i, 1:j - 1) * L(j, 1:j - 1) * D(1:j - 1))) / D(j)
      end do
    end do
  end subroutine ldl_psd

  !> B <- L^-1 B (or L'^-1 B if `transpose`) for unit lower triangular L.
  subroutine solve_unit_lower(L, B, transpose)
    real(dp), intent(in), contiguous :: L(:, :)
    real(dp), intent(inout), contiguous :: B(:, :)
    logical, intent(in), optional :: transpose
    character :: trans

    if (size(B) == 0) return
    trans = 'N'
    if (present(transpose)) then
      if (transpose) trans = 'T'
    end if
    call dtrsm('L', 'L', trans, 'U', size(B, 1), size(B, 2), 1.0_dp, L, &
               max(1, size(L, 1)), B, max(1, size(B, 1)))
  end subroutine solve_unit_lower

  !> X = A^+ B for symmetric positive semi-definite A, via A = L D L' with
  !> zero pivots treated as a pseudo-inverse.
  function psd_solve(A, B) result(X)
    real(dp), intent(in), contiguous :: A(:, :), B(:, :)
    real(dp), allocatable :: X(:, :)
    real(dp), allocatable :: L(:, :), D(:)
    integer :: i

    allocate (L(size(A, 1), size(A, 1)), D(size(A, 1)))
    call ldl_psd(A, L, D)
    X = B
    call solve_unit_lower(L, X)
    do i = 1, size(D)
      if (D(i) > 0.0_dp) then
        X(i, :) = X(i, :) / D(i)
      else
        X(i, :) = 0.0_dp
      end if
    end do
    call solve_unit_lower(L, X, transpose=.true.)
  end function psd_solve

  !> Lower triangular L with L L' = U U' for U of shape (k, n): the "tria"
  !> operation of square root filtering (DK 6.3), from the QR factorization
  !> U' = Q R, so that U = R' Q' and L = R'.
  function tria(U) result(L)
    real(dp), intent(in), contiguous :: U(:, :)
    real(dp), allocatable :: L(:, :)
    real(dp), allocatable :: A(:, :), tau(:), work(:)
    integer :: k, n, i, info

    k = size(U, 1)
    n = max(size(U, 2), k)          ! pad with zero columns if U is wide
    allocate (A(n, k), source=0.0_dp)
    A(1:size(U, 2), :) = transpose(U)
    allocate (tau(k), work(max(1, 64 * k)))
    call dgeqrf(n, k, A, n, tau, work, size(work), info)
    allocate (L(k, k), source=0.0_dp)
    do i = 1, k
      L(i:k, i) = A(i, i:k)
    end do
  end function tria

  !> True if all off-diagonal elements are below `tol` in absolute value.
  pure logical function is_diagonal(A, tol)
    real(dp), intent(in) :: A(:, :), tol
    integer :: i, j

    is_diagonal = .false.
    do j = 1, size(A, 2)
      do i = 1, size(A, 1)
        if (i /= j .and. abs(A(i, j)) >= tol) return
      end do
    end do
    is_diagonal = .true.
  end function is_diagonal

  !> A <- (A + A') / 2.
  subroutine symmetrize(A)
    real(dp), intent(inout), contiguous :: A(:, :)
    integer :: i, j

    do j = 2, size(A, 2)
      do i = 1, j - 1
        A(i, j) = 0.5_dp * (A(i, j) + A(j, i))
        A(j, i) = A(i, j)
      end do
    end do
  end subroutine symmetrize

  !> n x m identity-like matrix (ones on the leading diagonal).
  pure function eye(n, m) result(A)
    integer, intent(in) :: n
    integer, intent(in), optional :: m
    real(dp), allocatable :: A(:, :)
    integer :: i, mm

    mm = n
    if (present(m)) mm = m
    allocate (A(n, mm), source=0.0_dp)
    do i = 1, min(n, mm)
      A(i, i) = 1.0_dp
    end do
  end function eye
end module statespace_linalg
