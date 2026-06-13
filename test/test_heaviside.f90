program test_heaviside
    use iso_fortran_env, only: real64
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1   
    implicit none
    real(real64), allocatable :: X(:), Y(:)
    real(real64) :: start
    start = 1.75_real64
    X = seq(st=-10.0_real64, en=10.0_real64, len=10000)
    Y = Heaviside(X=X, val=start, right=.false.)
    call plot_gnu(X=X, Y=Y)
end program test_heaviside