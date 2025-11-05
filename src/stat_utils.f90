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
    
end module stat_utils