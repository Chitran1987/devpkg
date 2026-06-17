module stat_utils
    use iso_fortran_env, only: real64
    use vec_utils
    use plot_utils
    use calculus_utils
    use linalg_solver
    implicit none
    
contains
function lin_reg(dat) result(C) 
    real(real64) :: dat(:, :)
    real(real64), allocatable :: C(:),A(:),M(:,:), dmp(:)
    integer ::p, q, i, j, info
    p = size(dat, 1)
    q = size(dat, 2)
    allocate(C(q))
    allocate(A(q))
    allocate(dmp(q))
    allocate(M(q,q)) !This is a symmetric matrix, barring the last row and column
    do concurrent (j=1:q-1, i=1:q-1, i>=j) !Fill the first half
        A(j) = dot_product(dat(:,j),dat(:,q))
        M(i,j) = dot_product(dat(:,i), dat(:,j))
    end do
    do concurrent (j=2:q-1, i=1:q-1, i<j) !Write to the 2nd half
        M(j,i) = M(i,j)
    end do
    A(q) = sum(dat(:,q))
    do i = 1, q-1
        dmp(i) = sum(dat(:,i))
    end do
    dmp(q) = real(p, kind=real64)
    M(:,q) = dmp
    M(q,:) = dmp
    call solve_linear_system(A = M, b = A, info = info)
    if ( info == 0 ) then
        C = A
    else if ( info > 0 ) then
        error stop "Matrix M is singular — cannot invert."
    else
        error stop "Invalid arguments passed - please check arguments"
    end if
end function lin_reg

function lin_bg_sub_1D(dat, win, min_zero) result(ret_mat)
    real(real64) :: dat(:,:), win(:,:) !Input declaration
    real(real64), allocatable :: ret_mat(:,:) !Output declaration
    real(real64), allocatable :: new_dat(:,:) !matrix fed in linreg
    real(real64) :: X(size(dat, 1)), Y(size(dat, 1)) !The X and Y vectors built out of the dataset dat
    real(real64) :: coeff(2) !The regression co-efficients
    logical :: mask_int(size(win,1), size(dat,1)) !The mask matrix declaration
    logical :: mask(size(dat, 1)) !The mask vector declaration
    logical, intent(in), optional :: min_zero !The optional input
    logical :: minz !Need it because min_zero is optional
    integer :: n_win !The no of rows in the win matrix
    integer :: i !Internally usable integers

    !default values amongst arguments
    if ( .not. present(min_zero) ) then
        minz = .true.
    else
        minz = min_zero
    end if
    !global scope for entire function
    X = dat(:,1)
    Y = dat(:,2)

    !Create the actual mask matrix. 
    !See whether this matrix can be made cache friendly
    n_win = size(win, 1)
    do i = 1, n_win
        mask_int(i,:) = ( X >= win(i,1) ) .and. ( X <= win(i,2))
    end do

    !Create the mask vector
    mask = any(mask_int, dim=1)
    !Mask to get the new X and Y
    X = pack(X, mask)
    Y = pack(Y, mask)
    !Create the matrix which will be fed into linreg
    allocate(new_dat(size(X),2))
    new_dat(:,1) = X
    new_dat(:,2) = Y
    !get the regression co-efficients
    coeff = lin_reg(dat=new_dat)
    !Input the original matrix to the result
    ret_mat = dat
    ret_mat(:,2) = ret_mat(:,2) - (coeff(1)*ret_mat(:,1) + coeff(2))
    if ( minz .eqv. .true. ) then
        ret_mat(:,2) = ret_mat(:,2) - minval(ret_mat(:,2))
    end if
end function lin_bg_sub_1D

!function produces a uniformly distributed variable X ~ (min_val, max_val)
function unidistr(len, min_val, max_val) result(Y)
    integer :: len
    real(real64) :: min_val, max_val
    real(real64) :: Y(len)
    call random_number(Y)
    Y = Y*(max_val - min_val) + min_val
    return
end function unidistr

!Function for linear interpolation
function lin_interp(dat_X, dat_Y, X) result(Y)
    !Inputs 
    real(real64) :: dat_X(:), X(:)
    real(real64) :: dat_Y(size(dat_X))
    !Outputs 
    real(real64) :: Y(size(X))
    !Internal
    integer :: m, i, k
    real(real64) :: slope, c
    !logic code
    m = size(X)
    if (m > 1000) then
        !run a parallel loop
        !$omp parallel do
        do i=1,m
            !figure out which X(i) values it lies within which neighboring dat_X values
            k = maxloc(dat_X, dim=1, mask = dat_X <= X(i) ) !The index k | dat_X(k) <= X(i)
            if (dat_X(k) == X(i)) then
            Y(i) = dat_Y(k)
            else
            slope = ( dat_Y(k+1) - dat_Y(k) )/( dat_X(k+1) - dat_X(k) ) !slope
            c = ( dat_Y(k)*dat_X(k+1) - dat_Y(k+1)*dat_X(k) )/( dat_X(k+1) - dat_X(k) ) !Intercept
            Y(i) = slope*X(i) + c !draw a straight line between those X values
            endif
        end do
        !$omp parallel end do
    else
    !run a sequential loop
        do i=1,m
            !figure out which X(i) values it lies within which neighboring dat_X values
            k = maxloc(dat_X, dim=1, mask = dat_X <= X(i) ) !The index k | dat_X(k) <= X(i)
            if (dat_X(k) == X(i)) then
            Y(i) = dat_Y(k)
            else
            slope = ( dat_Y(k+1) - dat_Y(k) )/( dat_X(k+1) - dat_X(k) ) !slope
            c = ( dat_Y(k)*dat_X(k+1) - dat_Y(k+1)*dat_X(k) )/( dat_X(k+1) - dat_X(k) ) !Intercept
            Y(i) = slope*X(i) + c !draw a straight line between those X values
            endif
        end do
    endif
    return
end function lin_interp


!function produces a functioned f(X) distributed variable X ~ (min_val, max_val)
!f(X) could be a gaussian or could be anything else
!f(X) is denoted by the dataset dat_X, dat_Y
!Basically, data is generated by rejection sampling off of a uniform distribution
!The data does not count for the no. of remaining samples --> Use internally
!Uses linear interpolation internally --> slow --> use internally
function funcdistr(len, min_val, max_val, dat_X, dat_Y, scale) result(Y)
    !Inputs
    integer :: len
    real(real64) :: min_val, max_val, scale
    real(real64) :: dat_X(:)
    real(real64) :: dat_Y(size(dat_X))
    !Output
    real(real64), allocatable :: Y(:)
    !Internals
    real(real64), allocatable :: X1(:), Y1(:), Y1_predic(:)
    real(real64) :: min_dat_Y(1), max_dat_Y(1)
    real(real64), allocatable :: dat_X1(:), dat_Y1(:)
    real(real64) :: rat, span_ar
    logical, allocatable :: mask(:)
    !!!Real code----------------------------------------------!
    min_dat_Y = minval(dat_Y)
    max_dat_Y = maxval(dat_Y)
    span_ar = (max_val - min_val)*max_dat_Y(1) 
    !Get an idea of the acceptance area within the rejection block
    rat = integrate(X=dat_X, Y=dat_Y, xmin=min_val, xmax=max_val, Riemann = .true.)/span_ar
    !Build the first iteration of values
    X1 = unidistr(ceiling(scale*len/rat), min_val, max_val)
    Y1 = unidistr(len = ceiling(scale*len/rat), min_val = 0.0_real64, max_val = max_dat_Y(1))
    !Create the new interpolated distribution 
    dat_X1 = X1
    dat_Y1 = lin_interp(dat_X = dat_X, dat_Y = dat_Y, X = dat_X1)
    !Kick out values based on rejection sampling
    mask = Y1 <= dat_Y1
    Y = pack(X1, mask)
    return
end function funcdistr

!rejection sampling for ANY probability distr
!create a function for rejection sampling
!uses linear interp internally(through funcdistr) ---> slow
function rejection_sampling(len, min_val, max_val, dat_X, dat_Y)result(Y)
    !Inputs
    integer :: len
    real(real64) :: min_val, max_val
    real(real64) :: dat_X(:)
    real(real64) :: dat_Y(size(dat_X))
    !Output
    real(real64) :: Y(len)
    !Internals
    real(real64), allocatable :: Y1(:)
    real(real64) :: sc !scale
    !Logic
    sc = 2.0_real64
    Y1 = funcdistr(len, min_val, max_val, dat_X, dat_Y, sc)
    do while(size(Y1) < len)
        sc = sc + 1.0_real64
        Y1 = funcdistr(len, min_val, max_val, dat_X, dat_Y, sc)
    end do
    Y = Y1(1:len)
    return
end function rejection_sampling

!function produces a vector sampled from a simple gaussian distribution
!Use internally ---> can't control the length of output
function gaussdistr_int(len, x0, sig, scale) result(Y)
    !Inputs
    integer :: len
    real(real64) :: x0, sig, scale
    !Output
    real(real64), allocatable :: Y(:)
    !Internals
    real(real64), allocatable :: X1(:)
    real(real64), allocatable :: Y1(:), gauss_Y(:)
    real(real64) :: min_val, max_val, amp
    real(real64) :: rat, span_ar
    logical, allocatable :: mask(:)
    !!!Real code----------------------------------------------!
    !Build the first iteration of values
    min_val = x0 - 5.0_real64*sig
    max_val = x0 + 5.0_real64*sig
    amp = 1/sqrt(2*pi*(sig*sig))
    !Get an idea of the acceptance area within the rejection block
    span_ar = 10.0_real64/(sqrt(2.0_real64*pi)) !The spanning area
    rat = 1.0_real64/span_ar
    X1 = unidistr(len=ceiling(scale*len/rat), min_val = min_val, max_val = max_val)
    Y1 = unidistr(len=ceiling(scale*len/rat), min_val=0.0_real64, max_val=amp)
    gauss_Y = gauss_1D(X=X1, A=amp, x0=x0, sig=sig)
    mask = Y1 <= gauss_Y
    Y = pack(X1, mask)
    return
end function gaussdistr_int

!Quick, no linear interpolation used
!Expose
function sample_gauss(len, x0, sig) result(X)
    !Inputs
    integer :: len
    real(real64) :: x0, sig
    !Outputs
    real(real64) :: X(len)
    !Internal
    !real(real64) :: dat_X(10*len), dat_Y(10*len), amp
    real(real64), allocatable :: X1(:)
    real(real64) :: fact
    !Logic
    fact = 2.0_real64 !start point for factoring
    X1 = gaussdistr_int(len=len, x0=x0, sig=sig, scale=fact)
    do while(size(X1) < len)
        fact = fact + 1.0_real64
        X1 = gaussdistr_int(len=len, x0=x0, sig=sig, scale=fact)
    end do
    X = X1(1:len)
    return
end function sample_gauss

!Define a lorentzian function
pure elemental real(real64) function lorentzian(x, A, gamm, x0) result(y)
    !Inputs and Outputs
    real(real64), intent(in) :: x, A, gamm, x0
    y = A*(1/(pi*gamm))*(1/(1 + ((x - x0)/(gamm))**2))
    return
end function lorentzian

!function produces a vector sampled from a simple lorentzian distribution
!Use internally ---> can't control the length of output
!Length scale terminated to 99% of sampling probability range since lorentians have fat tails i.e. 2nd_moment -> +\infty
function lorentzdistr_int(len, x0, gamm, scale) result(Y)
    !Inputs
    integer :: len
    real(real64) :: x0, gamm, scale
    !Output
    real(real64), allocatable :: Y(:)
    !Internals
    real(real64), allocatable :: X1(:)
    real(real64), allocatable :: Y1(:), lor_Y(:)
    real(real64) :: min_val, max_val, amp
    real(real64) :: rat, span_ar
    logical, allocatable :: mask(:)
    !!!Real code----------------------------------------------!
    !Build the first iteration of values
    min_val = x0 - 63.66_real64*gamm
    max_val = x0 + 63.66_real64*gamm
    amp = 1.0_real64
    !Get an idea of the acceptance area within the rejection block
    span_ar = 127.32_real64/(pi*gamm) !The spanning area
    rat = 1.0_real64/span_ar
    X1 = unidistr(len=ceiling(scale*len/rat), min_val = min_val, max_val = max_val)
    Y1 = unidistr(len=ceiling(scale*len/rat), min_val=0.0_real64, max_val=amp)
    lor_Y = lorentzian(x=X1, A=1.0_real64, gamm=gamm, x0=x0)
    mask = Y1 <= lor_Y
    Y = pack(X1, mask)
    return
end function lorentzdistr_int

!Sample lorentz distribution
function sample_lorentz(len, x0, gamm) result(X)
    !Inputs
    real(real64) :: x0, gamm
    integer :: len
    !Outputs
    real(real64) :: X(len)
    !Internal
    real(real64), allocatable :: X1(:)
    real(real64) :: fact
    !Logic
    fact = 30.0_real64 !starting point for factoring
    X1 = lorentzdistr_int(len=len, x0=x0, gamm=gamm, scale=fact)
    !X = unidistr(len=len, min_val =0.0_real64, max_val=1.0_real64)
    do while(size(X1) < len)
        fact= 10.0*real64
        X1 = lorentzdistr_int(len=len, x0=x0, gamm=gamm, scale=fact)
    end do
    X = X1(1:len)
    return
end function sample_lorentz


end module stat_utils