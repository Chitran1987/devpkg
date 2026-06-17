program test_sample_lorentz
    use iso_fortran_env, only: real64
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1 
    implicit none
    !code
    real(real64), allocatable :: X(:)
    X = sample_lorentz(len=1000, x0=-5.7_real64, gamm=0.05_real64)
    call writeXY(X=X, Y=X, nam='sample_lorentz')
    print*, size(x)
end program test_sample_lorentz