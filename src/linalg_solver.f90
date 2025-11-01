module linalg_solver
  use iso_fortran_env, only: real64
  implicit none
contains

  subroutine solve_linear_system(A, b, info)
    !! Solve A x = b with LAPACK DGESV (in-place on A and b)
    real(real64), intent(inout)         :: A(:,:)
    real(real64), intent(inout), target :: b(:)     ! target so we can make a 2-D view
    integer,      intent(out)           :: info

    ! LAPACK interface (pure Fortran). Interface has its own scope -> bring real64 in here too.
    interface
       subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         use iso_fortran_env, only: real64
         integer :: n, nrhs, lda, ldb, ipiv(*), info
         real(real64) :: a(lda,*), b(ldb,*)
       end subroutine dgesv
    end interface

    integer :: n, nrhs, lda, ldb
    integer, allocatable :: ipiv(:)
    real(real64), pointer :: B2D(:,:)

    n = size(A,1)
    if (size(A,2) /= n) then
       info = -1; return   ! not square
    end if
    if (size(b) /= n) then
       info = -2; return   ! rhs length mismatch
    end if

    nrhs = 1
    lda  = n
    ldb  = n

    ! 2-D view of b without copying, as LAPACK expects B(LDB,NRHS)
    B2D(1:n,1:1) => b

    allocate(ipiv(n))
    call dgesv(n, nrhs, A, lda, ipiv, B2D, ldb, info)
    deallocate(ipiv)
    nullify(B2D)
  end subroutine solve_linear_system

end module linalg_solver
