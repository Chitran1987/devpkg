module calculus_utils
    use vec_utils
    use iso_fortran_env, only: real64
    implicit none
    
contains

    ! Integrates a 1D function between X and Y within limits xmin and xmax
    function integrate(X, Y, xmin, xmax, Riemann)
        real(real64), intent(in), dimension(:) :: X, Y
        real(real64) :: integrate
        real(real64), dimension(:), allocatable :: X_sub, Y_sub
        real(real64) :: xmin, xmax
        real(real64) :: sum
        logical, dimension(:), allocatable :: mask
        integer :: i
        logical, optional :: Riemann
        mask = ( (X >= xmin) .and. (X <= xmax) )
        X_sub = pack(X, mask) 
        Y_sub = pack(Y, mask)

        !Set Default Riemann argument to true
        if ( .not. present(Riemann) ) then
            Riemann = .true.
        end if

        !If Riemann Integratable
        if ( Riemann ) then
            sum = 0
            do i = 1, size(X_sub)-1
                sum = sum + Y_sub(i)*(X_sub(i+1)-X_sub(i))
            end do
            integrate = sum
            sum = 0
        else
            sum = 0
            do i = 1, size(X_sub)-1
                sum = sum + Y_sub(i)*(X_sub(i+1)-X_sub(i)) ! The Riemann rectangle
                sum = sum + 0.5*(Y_sub(i+1) - Y_sub(i))*(X_sub(i+1) - X_sub(i))
            end do
            integrate = sum
            sum = 0
        end if
        
    end function integrate


    !Differentiates a function Y = f(X) with an ord parameter
    recursive function differentiate(X, Y, ord) result(diffn)
        real(real64), dimension(:), intent(in) :: X, Y
        real(real64), dimension(:,:), allocatable :: diffn, dummy
        integer :: ord
        real(real64), dimension(:), allocatable :: X1, Y1, X2, Y2
        if ( ord == 1 ) then
            Y1 = diff(Y)/diff(X)
            X1 = mdpnt_vec(X)
            allocate(diffn(size(X1), 2))
            diffn(:,1) = X1
            diffn(:,2) = Y1
        else 
            dummy = differentiate(X, Y, 1)
            X2 = dummy(:,1)
            Y2 = dummy(:,2)
            diffn = differentiate(X2, Y2, ord-1)
        end if

    end function differentiate



end module calculus_utils