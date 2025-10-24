program test_funcs
    use vec_utils
    implicit none
    real(real64), dimension(10) :: X
    X = seq(1.0_real64, 10.0_real64, 10)
    print *, X
    print *, mean(X)
    print *, sdev(X)
    print *, var(X)
    print *, sdev(X)
    print *, var(X)
end program test_funcs