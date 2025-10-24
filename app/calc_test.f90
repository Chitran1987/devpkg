program calc_test
    !First write down the use statements
    use calculus_utils
    use vec_utils
    use plot_utils

    !Then write down the implicit none statements
    implicit none

    !Then declare the variables
    real(real64), dimension(:), allocatable :: X, Y
    real(real64) :: area
    real(real64), dimension(:,:), allocatable :: M

    !Assign and test
    allocate(X(10**4))
    allocate(Y(10**4))
    X = seq(-4.0_real64, 4.0_real64, 10**4)
    Y = sin(X)
    !
    call plot_gnu(X, Y)
    !Test the integrate function
    area = integrate(X, Y, 0.0_real64, 3.0_real64, Riemann = .false.)
    print *, "The area under the curve is ", area

    !Test the differentiate function
    M = differentiate(X, Y, 4)
    print *, "The sum of M, the differentiate matrix, is ", sum(M)
    !M = differentiate(M(:,1), M(:,2), 4)
    call plot_gnu(M(:,1), M(:,2))

    print *, "My fortran's default kind is", kind(1.0)

end program calc_test