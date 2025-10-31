module linalg_solver
  use iso_fortran_env
  implicit none
contains

  subroutine solve_linear_system(A, b, info)
    !! Solves A x = b using LAPACK DGESV from OpenBLAS
    real(real64), intent(inout) :: A(:,:)
    real(real64), intent(inout) :: b(:)
    integer, intent(out) :: info

    interface
       subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info) bind(C, name="dgesv_")
         use iso_c_binding
         integer(c_int), value :: n, nrhs, lda, ldb
         integer(c_int) :: ipiv(*)
         real(c_double) :: a(lda,*), b(ldb,*)
         integer(c_int) :: info
       end subroutine dgesv
    end interface

    integer :: n, nrhs, lda, ldb
    integer, allocatable :: ipiv(:)

    n    = size(A,1)
    nrhs = 1
    lda  = n
    ldb  = n

    allocate(ipiv(n))
    call dgesv(n, nrhs, A, lda, ipiv, b, ldb, info)
    deallocate(ipiv)

  end subroutine solve_linear_system

end module linalg_solver
