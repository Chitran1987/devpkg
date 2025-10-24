program fft_func_attempt
    !use trialpkg, only: say_hello
    use, intrinsic :: iso_c_binding
    use fft_utils
    use plot_utils
    use vec_utils
    implicit none
    !include '../src/fftw3.f03'
    
    !First, you build the dataset
    real(real64), dimension(:), allocatable :: X, Y
    real(real64), allocatable :: M(:,:,:)
    X = seq(-10.0_real64, 10.0_real64, 10**4)
    !Y = sin(30*X) + cos(70*X) + sin(75*X) + sin(60*X + 2*asin(1.0)/3) + cos(110*X + 2*acos(-1.0)/3)
    Y = cos(30*X)
    M = fft_1D(X, Y)
    print *, "fft calculation over"
    call plot_gnu(M(:,1,1), M(:,2,1)) !The amplitude spectra
    call plot_gnu(M(:,1,2), M(:,2,2)) !The phase spectra
end program fft_func_attempt