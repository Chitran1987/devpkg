program rect_lattice_XY_diff
    use iso_fortran_env, only: real64
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1
    use lattice_utils
    implicit none
    real(real64), allocatable :: X(:), Y(:), latt1(:,:,:), grid(:,:,:) !Lattice variables
    real(real64), allocatable :: f_latt(:,:,:)
    !build the square lattice
    X = seq(st=-10.0_real64, en=10.0_real64, len=1000)
    Y = seq(st = -8.0_real64, en = 8.0_real64, len = 800)
    !grid = grid_2(X, Y)
    allocate(latt1(size(Y), size(X), 3))
    !latt(:,:,2) = grid(:,:,1)
    !latt(:,:,3) = grid(:,:,2)    
    !latt(:,:,1) = sin(3*latt(:,:,2)+4*latt(:,:,3))*sin(4*latt(:,:,2)-3*latt(:,:,3))
    latt1 = rect_latt_sb(X = X, Y = Y, R_latt_x = 0.75_real64, R_latt_y = 0.5_real64, A = 1.0_real64, sig = 0.15_real64)
    call MatrixWrite(M=latt1(:,:,1), nam='rect_lattice_XY_diff')

    !build the Fourier Transform
    f_latt = fft_2D(tens = latt1, sampling_del=0.1_real64)
    call MatrixWrite(M=f_latt(:,:,1), nam = 'ft_rect_lattice_XY_diff')
    
end program rect_lattice_XY_diff