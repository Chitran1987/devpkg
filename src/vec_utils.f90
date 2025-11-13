module vec_utils
    use iso_fortran_env, only: real64
    implicit none
    public
    private :: rev_cmplx, rev_real, rev_int
    private :: print_mat_bool, print_mat_int, print_mat_real, print_mat_cmplx
    private :: lineval_r

    !Define the constants------------------------------------------------------
    real(real64) :: pi = 4.0_real64*atan(1.0_real64)
    !End definition of constants-----------------------------------------------
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

    pure elemental real(real64) function gauss_1D_core(x, A, x0, sig) result(y)
        real(real64), intent(in) :: x
        real(real64), intent(in) :: A, x0, sig
        y = A*exp((-1.0_real64*(x-x0)**2.0_real64/(2.0_real64*sig**2.0_real64)))
    end function gauss_1D_core

    function gauss_1D(x, A, x0, sig) result(y)
        real(real64), intent(in) :: x(:)
        real(real64), optional :: A, x0, sig
        real(real64) :: A_alt, x0_alt, sig_alt
        real(real64) :: y(size(x))
        if ( .not. present(x0) ) then
            x0_alt = 0.0_real64
        else
            x0_alt = x0
        end if
        if ( .not. present(sig) ) then
            sig_alt = 1.0_real64
        else
            sig_alt = sig
        end if
        if ( .not. present(A) ) then
            A_alt = 1.0_real64/(sqrt(2*pi)*sig_alt)
        else
            A_alt = A
        end if
        y = gauss_1D_core(x, A_alt, x0_alt, sig_alt)
    end function gauss_1D

    function gauss_2D_nocorr(X, Y, Ax, Ay, x0, y0, sig_x, sig_y) result(tens)
        real(real64), intent(in) :: X(:), Y(:) !The main inputs
        real(real64), optional :: Ax, Ay, x0, y0, sig_x, sig_y !The optional inputs
        real(real64) :: tens(size(Y), size(X), 3) !The output
        real(real64) :: Xsp(size(Y), size(X)), Ysp(size(Y), size(X)), Gsp(size(Y), size(X)) !The dummys for each tensor position tens(:,:,1), tens(:,:,2) and tens(:,:,3)
        real(real64) :: Gx(size(X)), Gy(size(Y)), Gx_sp(size(Y), size(X)), Gy_sp(size(Y), size(X)) !The dummys for Gx*Gy at each position 
        real(real64) :: Ax_alt, Ay_alt, x0_alt, y0_alt, sig_x_alt, sig_y_alt !The alternate inputs
        !Error handling for Ax, Ay, x0, y0, sig_x, sig_y
        if ( .not. present(x0) ) then
            x0_alt = 0.0_real64
        else
            x0_alt = x0
        end if
        if ( .not. present(sig_x) ) then
            sig_x_alt = 1.0_real64
        else
            sig_x_alt = sig_x
        end if
        if ( .not. present(Ax) ) then
            Ax_alt = 1.0_real64/(sqrt(2*pi)*sig_x_alt)
        else
            Ax_alt = Ax
        end if
        if ( .not. present(y0) ) then
            y0_alt = 0.0_real64
        else
            y0_alt = y0
        end if
        if ( .not. present(sig_y) ) then
            sig_y_alt = 1.0_real64
        else
            sig_y_alt = sig_y
        end if
        if ( .not. present(Ay) ) then
            Ay_alt = 1.0_real64/(sqrt(2*pi)*sig_y_alt)
        else
            Ay_alt = Ay
        end if
        !
        !The core algorithm
        !First assume that everything is present
        Gx = gauss_1D_core(X, Ax_alt, x0_alt, sig_x_alt)
        Gy = gauss_1D_core(Y, Ay_alt, y0_alt, sig_y_alt)
        call rev(Gy) !Gy needs to be reversed
        Gx_sp = spread(Gx, dim=1, ncopies=size(Y))
        Gy_sp = spread(Gy, dim=2, ncopies=size(X))
        Gsp = Gx_sp*Gy_sp
        Xsp = spread(X, dim = 1, ncopies=size(Y))
        call rev(Y) !Y needs to be reversed for plotting as an X-Y plane
        Ysp = spread(Y, dim = 2, ncopies=size(X))
        tens(:,:,1) = Gsp
        tens(:,:,2) = Xsp
        tens(:,:,3) = Ysp

    end function gauss_2D_nocorr

    !Select a vertical line profile out of the Tm,n,p tensor  
    function lin_prof_v(M, v_val) result(res_mat)
        real(real64) :: M(:,:,:), v_val !inputs declaration
        real(real64) :: res_mat(size(M,1), 2) !output declaration
        real(real64) :: X_dumm(size(M,2)), dummy1(size(X_dumm)) !Dummy prior = 1
        integer :: imin(1) !Dummy prior = 2

        !core logic
        res_mat(:,1) = M(:,1,3) !The X-axis(distance vector) of the line profile is the Y-axis(3rd slice, all rows any single column) of the tensor 
        X_dumm = M(1,:,2) !Subset the X-axis of the tensor, to choose which X-value will be selected to draw the vertical line through 
        dummy1 = abs(X_dumm - v_val)
        !Get the index of the minimum value of dummy
        imin = minloc(dummy1)
        !Place the dataset from that index in the result matrix
        res_mat(:,2) = M(:,imin(1),1) 

    end function lin_prof_v

    !Select a horizontal line profile out of the Tm,n,p tensor
    function lin_prof_h(M, h_val) result(res_mat)
        real(real64) :: M(:,:,:), h_val !inputs declaration
        real(real64) :: res_mat(size(M,2), 2) !output declaration
        real(real64) :: Y_dumm(size(M,1)), dummy1(size(Y_dumm)) !Dummy prior = 1
        integer :: imin(1) !Dummy prior = 2

        !core logic
        res_mat(:,1) = M(1,:,2) !The X-axis(distance vector) of the line profile is the X-axis(2nd slice, all columns any single row) of the tensor
        Y_dumm = M(:,1,3) !Subset the Y-axis of the tensor, to choose which Y-value will be selected to draw a horizontal line through
        dummy1 = abs(Y_dumm - h_val)
        !Get the index of the minimum value of the dummy
        imin = minloc(dummy1)
        !Place the sataset from that index in the result matrix
        res_mat(:,2) = M(imin(1),:,1)
    end function lin_prof_h

    !Subroutine for creating a 2 slice X-Y (rank 2) tensor grid
    function grid_2(X, Y) result(tens)
        real(real64) :: X(:), Y(:) !Inputs 
        real(real64) :: tens(size(Y), size(X), 2) !Outputs - Tensor size 
        integer :: m, n
        m = size(X)
        n = size(Y)
        tens(:,:,1) = spread(X, 1, n)
        call rev(Y)
        tens(:,:,2) = spread(Y, 2, m)
    end function grid_2

end module vec_utils
