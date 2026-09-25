!> Generic example of Fortran's built-in matmul. statespace lib uses BLAS/LAPACK.
!> This example was constructed from LLM output plus this SO post:
!> https://stackoverflow.com/a/52601770.
!> This example is included just to demonstrate how to neatly print matrices,
!> which is tricky in Fortran due to the columnar memory mapping.

program intrinsic_matmul
    implicit none
    integer, parameter :: M = 3, N = 3, P = 2
    integer :: i, j
    real :: A(M, N), B(N, P), C(M, P)

    ! Initialize sample matrices
    A = reshape([1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 9.0], [M, N])
    B = reshape([9.0, 7.0, -1.0, 8.0, 6.0, -2.3], [N, P])
    
    write(*, *) "Matrix A"
    write(*, "(*(g0))") ((A(i,j), " ", j=1, N), new_line("A"), i=1, M)
    write(*, *) "Matrix B"
    write(*, "(*(g0))") ((B(i,j), " ", j=1, P), new_line("A"), i=1, N)
    
    ! Matrix Multiplication
    C = matmul(A, B)

    ! Print result
    print *, "Product AB"
    write(*, "(*(g0))") ((C(i,j), " ", j=1, P), new_line("A"), i=1, M)
end program intrinsic_matmul
