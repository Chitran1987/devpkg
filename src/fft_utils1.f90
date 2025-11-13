module fft_utils1
    use, intrinsic :: iso_c_binding
    use iso_fortran_env, only: real64
    use vec_utils

    implicit none
    ! FFTW Fortran-2003 interface
    include '../src/fftw3.f03'

contains

    function fft_1D(X, Y) result(out)
        !! out(:,1,1) = angular frequency (centered around 0)
        !! out(:,2,1) = amplitude |F(ω)|
        !! out(:,1,2) = angular frequency (centered around 0)
        !! out(:,2,2) = phase arg(F(ω))

        real(real64), intent(in) :: X(:), Y(:)
        real(real64) :: out(size(X)+1, 2, 2)

        ! scalars
        integer :: n, nfreq
        integer :: i, idx
        real(real64) :: del_x, two_pi
        real(real64) :: omega_s, domega

        type(c_ptr) :: plan

        ! FFTW work arrays (double precision, C side)
        real(c_double)            :: y_c(size(X))
        complex(c_double_complex) :: fy(size(X)/2 + 1)

        ! Full spectrum & frequency axis in real64 (length n)
        complex(real64) :: fy_full(size(X))
        complex(real64) :: fy_centered(size(X))
        real(real64)    :: freq_centered(size(X))

        n      = size(X)              ! FFT length
        nfreq  = n/2 + 1              ! FFTW r2c length
        two_pi = 2.0_real64 * acos(-1.0_real64)

        ! step size Δx
        del_x = mean(diff(X))

        ! sampling angular frequency ω_s = 2π / Δx
        omega_s = two_pi / del_x
        domega  = omega_s / real(n, real64)     ! frequency bin width

        ! copy input to FFTW buffer (convert to c_double)
        y_c = real(Y, kind=c_double)

        ! --- FFT ---
        plan = fftw_plan_dft_r2c_1d(n, y_c, fy, FFTW_ESTIMATE)
        call fftw_execute_dft_r2c(plan, y_c, fy)
        call fftw_destroy_plan(plan)

        ! ------------------------------------------------------------------
        ! 1) Build full complex spectrum fy_full(k), k = 0..n-1 (1-based: 1..n)
        !    fy(1)       -> k = 0 (DC)
        !    fy(2:nfreq) -> k = 1..nfreq-1 (positive frequencies)
        !    fy_full(nfreq+1:n) = conjg of positive freqs in reverse
        !    This works for both even and odd n.
        ! ------------------------------------------------------------------

        ! lower half: k = 0..nfreq-1
        fy_full(1:nfreq) = cmplx(real(fy, kind=real64), aimag(fy), kind=real64)

        ! upper half: k = nfreq..n-1 (negative frequencies)
        if (n > 1) then
            fy_full(nfreq+1:n) = conjg( fy_full(nfreq-1:2:-1) )
        end if

        ! ------------------------------------------------------------------
        ! 2) Center the spectrum: order bins from negative to positive freq
        !
        !    We perform an "fftshift"-like circular shift.
        !
        !    For i = 0..n-1 (0-based index of centered array):
        !        original index k = (i + n/2) mod n
        !        freq_centered(i+1) = (i - n/2) * Δω
        !
        !    So freq_centered runs approximately from -ω_s/2 to +ω_s/2
        !    (exact endpoints depend on even/odd n, as usual for FFTs).
        ! ------------------------------------------------------------------
        do i = 0, n-1
            idx = mod(i + n/2, n)              ! 0-based index into fy_full
            fy_centered(i+1) = fy_full(idx+1)  ! convert 0-based -> 1-based
            freq_centered(i+1) = (real(i, real64) - real(n, real64)/2.0_real64) * domega
        end do

        ! ------------------------------------------------------------------
        ! 3) Fill output:
        !    - First n rows: centered spectrum & frequencies
        !    - Row n+1: periodic closure (same as row 1)
        ! ------------------------------------------------------------------

        ! amplitude
        out(1:n, 1, 1) = freq_centered
        out(1:n, 2, 1) = abs(fy_centered)

        ! phase
        out(1:n, 1, 2) = freq_centered
        out(1:n, 2, 2) = atan2(aimag(fy_centered), real(fy_centered))

        ! periodic closure so shape stays (n+1, 2, 2) as before
        out(n+1, :, :) = out(1, :, :)

    end function fft_1D

end module fft_utils1
