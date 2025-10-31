module stat_utils
    use iso_fortran_env, only: real64
    use vec_utils, plot_utils
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
    
end module stat_utils