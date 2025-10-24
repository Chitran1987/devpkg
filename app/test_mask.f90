program test_mask
    use vec_utils
    use plot_utils
    implicit none
    real, dimension(20) :: X, Y
    real, dimension(:), allocatable :: X_mask, Y_mask
    logical, dimension(:), allocatable :: mask
    X = seq(-10.0_real64, 10.0_real64, 20)
    Y = X**2
    print *, "The original array X"
    print *, X
    mask = ((X >= 5.0) .and. (X <= 10.0))
    X_mask = pack(X, mask)
    Y_mask = pack(Y, mask)
    print *, "The masked array X_mask"
    print *, X_mask
    print *, "The masked array Y_mask"
    print *, Y_mask
end program test_mask