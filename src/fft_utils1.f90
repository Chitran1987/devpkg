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

    !fft_2D function core module
    function fft_2D_core_mod(tens, sampling_del) result(G)
        !! 2D FFT on a rank-3 tensor tens(m,n,p)
        !!
        !! Input layout:
        !!   tens(:,:,1) = data f(x,y)
        !!   tens(:,:,2) = X grid
        !!   tens(:,:,3) = Y grid
        !!
        !! Output layout G(m,n,4):
        !!   G(:,:,1) = |F(ω_x, ω_y)|     (amplitude spectrum, centered)
        !!   G(:,:,2) = arg(F(ω_x, ω_y)) (phase spectrum, centered)
        !!   G(:,:,3) = ω_x coordinates  (same shape as data)
        !!   G(:,:,4) = ω_y coordinates  (same shape as data)

        use, intrinsic :: iso_c_binding
        use iso_fortran_env, only: real64
        use omp_lib

        real(real64), intent(in) :: tens(:,:,:), sampling_del
        real(real64) :: G(size(tens,1), size(tens,2), 4)

        ! dimensions
        integer :: m, n
        integer :: ix, iy, jx, jy
        integer :: idx_x, idx_y

        ! grids and data
        real(real64) :: Z(size(tens,1), size(tens,2))
        real(real64) :: X(size(tens,1), size(tens,2))
        real(real64) :: Y(size(tens,1), size(tens,2))

        ! FFTW arrays (complex double on C side)
        complex(c_double_complex) :: z_c(size(tens,1), size(tens,2))
        complex(c_double_complex) :: Fy_c(size(tens,1), size(tens,2))

        ! Centered spectrum in real64
        complex(real64) :: Fy_centered(size(tens,1), size(tens,2))

        ! frequency grids
        real(real64) :: Wx(size(tens,1), size(tens,2))
        real(real64) :: Wy(size(tens,1), size(tens,2))

        ! scalars
        real(real64) :: two_pi
        real(real64) :: del_x, del_y
        real(real64) :: omega_sx, omega_sy
        real(real64) :: domega_x, domega_y

        type(c_ptr) :: plan

        ! -------------------------
        ! 0) Dimensions & unpack
        ! -------------------------
        m = size(tens,1)
        n = size(tens,2)

        Z = tens(:,:,1)
        X = tens(:,:,2)
        Y = tens(:,:,3)

        two_pi = 2.0_real64 * acos(-1.0_real64)

        ! -----------------------------------------------------
        ! 1) Sampling steps and angular sampling frequencies
        !    - assume X is constant along rows: use first row
        !    - assume Y is constant along columns: use first col
        ! -----------------------------------------------------
        del_x = mean(diff( X(1,:) ))
        del_y = mean(diff( Y(:,1) ))

        ! -------------------------------------------------------
        ! Error checking in comparison to sampling rate
        ! -------------------------------------------------------
        if ( abs(abs(del_x) - abs(del_y)) / min(abs(del_x), abs(del_y)) > sampling_del ) &
            error stop 'sampling rates for X and Y differ by more than sampling_del'

        omega_sx = two_pi / del_x
        omega_sy = two_pi / del_y

        domega_x = omega_sx / real(n, real64)  ! bin width in ω_x
        domega_y = omega_sy / real(m, real64)  ! bin width in ω_y

        ! -----------------------------------------------------
        ! 2) Prepare complex input for FFTW (imag part = 0)
        ! -----------------------------------------------------
        z_c = cmplx( real(Z, kind=c_double), 0.0_c_double )

        ! -----------------------------------------------------
        ! 3) 2D FFT (complex-to-complex) using FFTW
        !    Result Fy_c is in the standard (0..m-1, 0..n-1)
        !    "uncentered" ordering.
        !
        !    NOTE: FFTW planning is NOT thread-safe.
        !    We therefore serialize plan creation/destruction.
        ! -----------------------------------------------------
        !$omp critical(fftw_planner)
        plan = fftw_plan_dft_2d(m, n, z_c, Fy_c, FFTW_FORWARD, FFTW_ESTIMATE)
        !$omp end critical(fftw_planner)

        call fftw_execute_dft(plan, z_c, Fy_c)

        !$omp critical(fftw_planner)
        call fftw_destroy_plan(plan)
        !$omp end critical(fftw_planner)

        ! -----------------------------------------------------
        ! 4) Center the spectrum (2D fftshift)
        !
        !    For 0-based indices:
        !      jy = 0..m-1, jx = 0..n-1  (centered array indices)
        !      ky = (jy + m/2) mod m     (original indices)
        !      kx = (jx + n/2) mod n
        !
        !    Also build ω_x and ω_y grids:
        !      ω_x(jy,jx) = (jx - n/2) * Δω_x
        !      ω_y(jy,jx) = (jy - m/2) * Δω_y
        ! -----------------------------------------------------
        do jy = 0, m-1
            idx_y = mod(jy + m/2, m)   ! 0-based original y index
            do jx = 0, n-1
                idx_x = mod(jx + n/2, n)  ! 0-based original x index

                ! convert 0-based -> 1-based for Fortran arrays
                Fy_centered(jy+1, jx+1) = cmplx( &
                    real(Fy_c(idx_y+1, idx_x+1), kind=real64), &
                    aimag(Fy_c(idx_y+1, idx_x+1)), &
                    kind=real64 )

                Wx(jy+1, jx+1) = ( real(jx, real64) - real(n, real64)/2.0_real64 ) * domega_x
                Wy(jy+1, jx+1) = ( real(jy, real64) - real(m, real64)/2.0_real64 ) * domega_y
            end do
        end do

        ! -----------------------------------------------------
        ! 5) Fill output tensor G(m,n,4)
        ! -----------------------------------------------------
        ! amplitude
        G(:,:,1) = abs(Fy_centered)

        ! phase (atan2(Im, Re))
        G(:,:,2) = atan2( aimag(Fy_centered), real(Fy_centered) )

        ! ω_x and ω_y grids
        G(:,:,3) = Wx
        G(:,:,4) = Wy

    end function fft_2D_core_mod


    !fft_2D after zero padding
    function fft_2D(tens, sampling_del) result(G)
        real(real64) :: tens(:,:,:), sampling_del !Inputs
        real(real64), allocatable :: G(:,:,:) !Outputs
        real(real64), allocatable :: tens_padded(:,:,:) !Internal

        !Zero pad the tensor
        tens_padded = zero_pad_tens(tens)
        G = fft_2D_core_mod(tens_padded, sampling_del)
        return
    end function fft_2D

    function fft_2D_map(img_tens, Xspan, Yspan, k1st, k0 ) result(res_tens)
        use omp_lib

        real(real64) :: img_tens(:,:,:) !The image tensor on which the fft_2D_map() will be applied
        real(real64) :: Xspan, Yspan !The Xspan and the Yspan of the image that is needed to be subsetted
        !real(real64) :: k1st(:,:) !A structure of size k1st(:,4) where each row regards to a specific
        real(real64) :: k0(4) !The integration square around the (0,0) spot. Only kx_span and ky_span are needed
        real(real64) :: res_tens(size(img_tens,1), size(img_tens,2), size(img_tens,3))  !The res_tens showing the resulting tensor

        !--------------------Input argument k1st(n_spots,4)-------------------------------------------------------------!
        !!  The k1st(n_spots, 4) argument
        !!  n_spots = no. of first order/higher order spots around which integration is needed
        !!  each row has 4 column.
        !!  (1st column, 2nd coulmn) = (x,y) of bottom left corner of square
        !!  (3rd column, 4th coulmn) = (x,y) of top right corner of square
        real(real64) :: k1st(:,:) !A structure of size k1st(n_spots,4) where each row regards to a specific
        !--------------------Input argument k1st(n_spots,4)-------------------------------------------------------------!

        !--------------------Internal variables-------------------------------------------------------------------------!
        integer :: m, n !Lateral sizes of img_tens
        integer :: n_spots !The no. of spots over which to integrate
        integer :: p !Integrnal dummy
        integer :: i,j !Internal dummy
        real(real64) :: kx_lim(2), ky_lim(2) !The actual co-ordinates of the kx and ky over which to integrate for the higher order spots
        real(real64) :: k0x_lim(2), k0y_lim(2) !The actual co-ordinates of the kx and ky over which to integrate for the 0th order spots
        real(real64), allocatable :: tmp_img(:,:,:) !The temporary masked image zoomed into for each iteration
        real(real64), allocatable :: fft_tmp(:,:,:) !The fft of the temporary image created
        real(real64) :: center(2) !The center for the tmp_img created on each iteration
        logical, allocatable :: mask_recip_1st_ord_spots(:,:) !The mask for the reciprocal lattice for the first order spots
        logical, allocatable :: mask_recip_center(:,:) !The mask for the (0,0) spot
        real(real64) :: dmp_sum, dmp_sum_zero !The integrated values
        !--------------------Internal variables-------------------------------------------------------------------------!

        !!Transfer the X and Y matrixes to the result
        res_tens(:,:,2) = img_tens(:,:,2) !The X matrix transfer
        res_tens(:,:,3) = img_tens(:,:,3) !The Y matrix transfer

        !!Lateral sizes of the tensor img_tens
        m = size(img_tens, 1)
        n = size(img_tens, 2)

        !!No of spots to integrate over
        n_spots = size(k1st, 1)

        !!Now for the actual code
        k0x_lim = [k0(1), k0(3)]
        k0y_lim = [k0(2), k0(4)]

        !$omp parallel do collapse(2) default(none) &
        !$omp shared(img_tens, res_tens, m, n, n_spots, Xspan, Yspan, k1st, k0x_lim, k0y_lim) &
        !$omp private(i, j, p, center, tmp_img, fft_tmp, dmp_sum, dmp_sum_zero, kx_lim, ky_lim, mask_recip_1st_ord_spots, mask_recip_center)
        do j = 1, n, 1
            do i = 1, m, 1
                center = [img_tens(i,j,2), img_tens(i,j,3)]
                tmp_img = mask_tens_cent(tens=img_tens, cent=center, Xspan=Xspan, Yspan=Yspan) !The tmp_img has been created
                fft_tmp = fft_2D(tens = tmp_img, sampling_del = 0.1_real64) !The fft_tmp has been created
                !--------------Integration jobs over the fft_tmp image---------------------------------------------------!
                !--------Integrate over the 1st/higher order spots----------------!
                dmp_sum = 0.0_real64
                do p = 1, n_spots, 1
                    !!!!!Start the integrate2D_process
                    kx_lim = [k1st(p,1), k1st(p,3)]
                    ky_lim = [k1st(p,2), k1st(p,4)]
                    mask_recip_1st_ord_spots = ( (fft_tmp(:,:,3) >= kx_lim(1)) .and. (fft_tmp(:,:,3) <= kx_lim(2)) ) .and. ( (fft_tmp(:,:,4) >= ky_lim(1)) .and. (fft_tmp(:,:,4) <= ky_lim(2)) )
                    dmp_sum = dmp_sum + sum(fft_tmp(:,:,1), mask_recip_1st_ord_spots)
                end do
                !--------Integrate over the 1st/higher order spots----------------!
                !--------Integrate over the zero order spots----------------!
                mask_recip_center = ( (fft_tmp(:,:,3) >= k0x_lim(1)) .and. (fft_tmp(:,:,3) <= k0x_lim(2)) ) .and. ( (fft_tmp(:,:,4) >= k0y_lim(1)) .and. (fft_tmp(:,:,4) <= k0y_lim(2)) )
                dmp_sum_zero = sum(fft_tmp(:,:,1), mask_recip_center)             !--------Integrate over the zero order spots----------------!
                !--------------Integration jobs over the fft_tmp image---------------------------------------------------!
                res_tens(i,j,1) = dmp_sum/(dmp_sum_zero+10**(-5))
            end do
        end do
        !$omp end parallel do

    end function fft_2D_map


end module fft_utils1
