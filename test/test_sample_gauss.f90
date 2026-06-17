program test_sample_gauss
    use iso_fortran_env, only: real64
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1 
    implicit none
    !code
    real(real64) :: X(10000), Y(10000)
    X = sample_gauss(len=10000, mu=0.0_real64, sig=1.0_real64)
    Y = X
    call writeXY(X, Y, 'test_sample_gauss')
    
end program test_sample_gauss