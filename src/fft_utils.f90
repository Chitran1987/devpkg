module fft_utils
    use, intrinsic :: iso_c_binding
    use vec_utils
    implicit none
    include '../src/fftw3.f03'
    
    
contains

!function for calclating 1D ffts. Uses fftw3.f03
    function fft_1D(X, Y) 
        real(real64), intent(in), dimension(:) :: X, Y
        real(kind = c_float), dimension(size(X)) :: Y_c
        complex(c_float_complex), dimension(size(X)/2+1) :: fy 
        complex(c_double_complex), dimension(size(X)/2+1) :: fy64, fy64_1
        complex(c_double_complex), dimension(size(X)+1) ::fy64_tot
        real(real64), dimension(size(X)+1, 2, 2) :: fft_1D
        real(real64), dimension(size(X)+1) :: wf
        real(real64) :: del_X
        type(c_ptr) :: plan
        integer :: n
        n = size(X)
        del_X = mean(diff(X))

        !change type of Y to Y_c for input
        Y_c = real(Y, kind=c_double)
        
        !!! The three lines which actually do the fft work
        !call the plan ptr to do the dummy work
        plan = fftwf_plan_dft_r2c_1d(n, Y_c, fy, FFTW_ESTIMATE)
        !call the fftwf function to store the fft in fy
        call fftwf_execute_dft_r2c(plan, Y_c, fy)
        !destroy the fft plan ptr space that was allocated
        call fftwf_destroy_plan(plan)

        !Convert fy to fy64
        fy64 = cmplx(real(fy, kind=c_double), aimag(fy), kind=c_double_complex)
        fy64_1 = conjg(fy64)
        call rev(fy64_1)
        fy64_tot = [fy64_1(1:(size(fy64_1)-1)), fy64]
        !Define the wf
        wf = seq((1*2*asin(-1.0_real64))/del_X, (1*2*asin(1.0_real64))/del_X, n+1)
        !!The amplitude of the complex tensor
        fft_1D(:, 1, 1) = wf
        fft_1D(:, 2, 1) = abs(fy64_tot)
        !!The phase of the complex tensor
        fft_1D(:, 1, 2) = wf
        fft_1D(:, 2, 2) = atan2(aimag(fy64_tot), real(fy64_tot))
    end function fft_1D

    !Create a subroutine for calling 1D FFTs
    
end module fft_utils