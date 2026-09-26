``statespace_linalg``
=====================

Linear algebra over BLAS and LAPACK. Arguments are contiguous arrays.
``gemm``, ``gemv`` and ``chol_inv`` use inline loops for small operands,
where the per-call cost of BLAS and LAPACK dominates (see
:doc:`../../design/performance`).

.. code-block:: fortran

   subroutine gemm(transa, transb, alpha, A, B, beta, C)
   subroutine gemv(trans, alpha, A, x, beta, y)
   subroutine chol_inv(A, logdet, info)
   subroutine solve(A, B, info)
   subroutine solve_lyapunov(T, V, P, info)
   subroutine ldl_psd(A, L, D)
   subroutine solve_unit_lower(L, B, transpose)
   function psd_solve(A, B) result(X)
   function tria(U) result(L)
   pure logical function is_diagonal(A, tol)
   subroutine symmetrize(A)
   pure function eye(n, m) result(A)

``gemm``, ``gemv``
   :math:`C \leftarrow \alpha\, op(A)\, op(B) + \beta C` and
   :math:`y \leftarrow \alpha\, op(A) x + \beta y`, with op the identity
   ('N') or the transpose ('T'). As in BLAS, C and y are not read when
   :math:`\beta = 0`.
``chol_inv``
   Invert a symmetric positive definite matrix in place by Cholesky and
   return :math:`\log|A|`; ``info`` is ``SS_ERR_NOT_PD`` otherwise.
``solve``
   :math:`A X = B` for square A (LU); B is overwritten with X.
``solve_lyapunov``
   :math:`P = T P T' + V` by doubling: :math:`P_{k+1} = P_k + A_k P_k
   A_k'`, :math:`A_{k+1} = A_k^2`. ``SS_ERR_NOT_STATIONARY`` if T has an
   eigenvalue on or outside the unit circle.
``ldl_psd``
   :math:`A = L D L'` for symmetric positive semi-definite A; pivots below
   a relative tolerance are set to zero, so singular A is allowed.
``solve_unit_lower``
   :math:`B \leftarrow L^{-1} B` (or :math:`L'^{-1} B`) for unit lower
   triangular L.
``psd_solve``
   :math:`A^+ B` for symmetric positive semi-definite A, through
   ``ldl_psd``.
``tria``
   Lower triangular L with :math:`L L' = U U'`, the "tria" operation of
   square root filtering (DK §6.3), from the QR factorization of U'.
``is_diagonal``, ``symmetrize``, ``eye``
   Off-diagonal elements below ``tol``; :math:`A \leftarrow (A + A')/2`;
   the n × m matrix with ones on the leading diagonal.
