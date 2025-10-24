program name
    use vec_utils
    use plot_utils
    implicit none
    real(real64), dimension(10**6) :: X, Y
    X = seq(-10.0_real64, 10.0_real64, 10**6)
    Y = sin(X)
    call plot_local(X, Y)
end program name