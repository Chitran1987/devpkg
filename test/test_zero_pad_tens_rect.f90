program test_zero_pad_rect
    use iso_fortran_env, only: real64
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1
    use lattice_utils
    implicit none
    real(real64), allocatable :: tens(:,:,:), tens_new(:,:,:)
    real(real64), allocatable :: f_tens(:,:,:), f_tens_new(:,:,:)
    real(real64), allocatable :: X(:), Y(:)
    X = seq(st=-10.0_real64, en = 10.0_real64, len = 1000)
    Y = seq(st = -5.0_real64, en = 5.0_real64, len = 500)

    ! Create the rectangular lattice and write it
    tens = rect_latt_sb(X = X, Y = Y, r_latt_x = 1.0_real64, r_latt_y = 0.5_real64, A = 1.0_real64, sig = 0.15_real64)
    call MatrixWrite(M = tens(:,:,1), nam = 'latt')

    !Create the zero padded square lattice and write it
    tens_new = zero_pad_tens(tens = tens)
    call MatrixWrite(M = tens_new(:,:,1), nam = 'latt_new')

    !Create the fft of the non-zero padded square lattice and write it
    f_tens = fft_2D(tens = tens, sampling_del = 0.1_real64)
    call MatrixWrite(M = f_tens(:,:,1), nam = 'f_latt')

    !Test the zeropad function X and Y matrixes
    call MatrixWrite(M = tens(:,:,2), nam = 'latt_old_X')
    call MatrixWrite(M = tens(:,:,3), nam = 'latt_old_Y')
    call MatrixWrite(M = tens_new(:,:,2), nam = 'latt_new_X' )
    call MatrixWrite(M = tens_new(:,:,3), nam = 'latt_new_Y')
    !Then come back and run this code

    !Create the fft of the zero padded squre lattice and write it
    f_tens_new = fft_2D(tens = tens_new, sampling_del = 0.1_real64)
    call MatrixWrite(M = f_tens_new(:,:,1), nam = 'f_latt_new')

end program test_zero_pad_rect