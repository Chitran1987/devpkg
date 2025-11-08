module vec_utils
    use iso_fortran_env, only: real64
    implicit none
    public
    private :: rev_cmplx, rev_real, rev_int
    private :: print_mat_bool, print_mat_int, print_mat_real, print_mat_cmplx
    private :: lineval_r

    ! The rev(X) interface
    !!! call rev(X)  ---- reverses a vector whether real, int, cmplx
    interface rev
        module procedure rev_cmplx, rev_real, rev_int
    end interface rev
    
    ! The print_mat(X) interface
    !!! integer dmp !!create a dummy variable
    !! dmp = print_mat(X) !!works for X as int, bool, real, cmplx 
    interface print_mat
        module procedure print_mat_bool, print_mat_int, print_mat_real, print_mat_cmplx
    end interface print_mat
    
    ! The lineval(X) interface
    interface lineval
            module procedure lineval_r
    end interface lineval
contains

    !function produces a sequence with start, stop and length 
    function seq(st, en, len)
        real(real64), intent(in) :: st, en
        integer, intent(in) :: len
        integer :: i
        real(real64) :: del
        real(real64), dimension(len) :: seq 
        del = (en-st)/(len-1)
        do i = 1, len
            seq(i) = st + (i-1)*del 
        end do
    end function seq

    !function calculates the difference of subsequent nos in a vector
    function diff(X)
        real(real64), dimension(:), intent(in)::X
        real(real64), dimension(size(X)-1) :: diff
        integer :: i
        do i = 1, size(diff)
            diff(i) = X(i+1) - X(i)
        end do
    end function diff

    !function calculates the mean of a vector
    function mean(X)
        real(real64), intent(in), dimension(:) :: X
        real(real64) :: mean
        mean = sum(X)/size(X)
    end function mean

    !function calculates the var of a vector
    function var(X)
        real(real64), intent(in), dimension(:) :: X
        real(real64) :: var
        real(real64) :: dmp = 0
        integer :: i
        real :: a
        a = mean(X)
        do i = 1, size(X)
            dmp = dmp + (X(i) - a)**2
        end do
        !print *, dmp
        !print *, size(X)-1
        var = dmp/(size(X)-1)
        dmp = 0
    end function var

    !function calculates the sdev of a vector
    function sdev(X)
        real(real64), intent(in), dimension(:) :: X
        real(real64) :: sdev
        sdev = sqrt(var(X))
    end function sdev

    !function averages between successive elements of vector
    function mdpnt_vec(X)
        real(real64), intent(in), dimension(:) :: X
        real(real64), dimension(:), allocatable :: mdpnt_vec
        integer :: i, n
        n = size(X)
        allocate(mdpnt_vec(n-1))
        do i = 1, n-1
            mdpnt_vec(i) = (X(i) + X(i+1))/2
        end do
    end function mdpnt_vec

    !subroutine for reversing a real vector
    subroutine rev_real(X)
        real(real64) :: X(:)
        integer :: n, i
        real(real64) :: tmp
        n = size(X)
        do i = 1, n/2
            tmp = X(i)
            X(i) = X(n+1-i)
            X(n+1-i) = tmp
        end do

    end subroutine rev_real

    !subroutine for reversing a complex vector
    subroutine rev_cmplx(X)
        complex(real64) :: X(:)
        integer :: n, i
        complex(real64) :: tmp
        n = size(X)
        do i = 1, n/2
            tmp = X(i)
            X(i) = X(n+1-i)
            X(n+1-i) = tmp
        end do

    end subroutine rev_cmplx

    !Subroutine for reversing an integer vector
    subroutine rev_int(X)
        integer :: X(:)
        integer n, i
        integer :: tmp
        n = size(X)
        do i = 1, n/2
            tmp = X(i)
            X(i) = X(n+1-i)
            X(n+1-i) = tmp
        end do
    end subroutine rev_int

    !function for print_mat_real(M)
    function print_mat_real(M) result(dmp)
        real(real64), intent(in) :: M(:,:)
        integer i, dmp
        do i = 1, size(M, 1)
            print *, M(i, :)
        end do
        dmp = 0
    end function print_mat_real

    function print_mat_int(M) result(dmp)
        integer, intent(in) :: M(:,:)
        integer i, dmp
        do i = 1, size(M, 1)
            print *, M(i, :)
        end do
        dmp = 0
    end function print_mat_int

    function print_mat_bool(M) result(dmp)
        logical, intent(in) :: M(:,:)
        integer i, dmp
        do i = 1, size(M, 1)
            print *, M(i, :)
        end do
        dmp = 0
    end function print_mat_bool

    function print_mat_cmplx(M) result(dmp)
        complex(real64), intent(in) :: M(:,:)
        integer i, dmp
        do i = 1, size(M, 1)
            print *, M(i, :)
        end do
        dmp = 0
    end function print_mat_cmplx

    subroutine lineval_r(X, low, hi)
        real(real64) :: X(..)
        real(real64), intent(in) :: low, hi
        select rank(X)
        rank(0)
            X = X
        rank(1)
            X = ((hi - low)/(maxval(X) - minval(X)))*(X - minval(X)) + low
        rank(2)
            X = ((hi - low)/(maxval(X) - minval(X)))*(X - minval(X)) + low
        rank(3)
            X = ((hi - low)/(maxval(X) - minval(X)))*(X - minval(X)) + low
        rank(4)
            X = ((hi - low)/(maxval(X) - minval(X)))*(X - minval(X)) + low
        rank(5)
            X = ((hi - low)/(maxval(X) - minval(X)))*(X - minval(X)) + low
        rank default
            error stop "arrays over rank 5 are not supported"
        end select
    end subroutine lineval_r

    pure elemental real(real64) function gauss_1D(x, A, x0, sig) result(y)
        real(real64), intent(in) :: x, A, x0, sig
        y = A*exp((-1.0_real64*(x-x0)**2.0_real64/(2.0_real64*sig**2.0_real64)))
    end function gauss_1D

    
end module vec_utils
