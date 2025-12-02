program rect_lattice_build
    use iso_fortran_env, only: real64
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1
    use lattice_utils
    implicit none
    real(real64), allocatable :: X(:), Y(:), latt(:,:,:), grid(:,:,:) !Lattice variables
    real(real64), allocatable :: f_latt(:,:,:)
    !build the square lattice
    X = seq(st=-10.0_real64, en=10.0_real64, len=500)
    Y = X
    grid = grid_2(X, Y)
    allocate(latt(size(Y), size(X), 3))
    !latt(:,:,2) = grid(:,:,1)
    !latt(:,:,3) = grid(:,:,2)    
    !latt(:,:,1) = sin(3*latt(:,:,2)+4*latt(:,:,3))*sin(4*latt(:,:,2)-3*latt(:,:,3))
    latt = rect_latt_sb(X = X, Y = Y, R_latt_x = 0.75_real64, R_latt_y = 0.5_real64, A = 1.0_real64, sig = 0.15_real64)
    call MatrixWrite(M=latt(:,:,1), nam='rect_lattice')

    !build the Fourier Transform
    f_latt = fft_2D(tens = latt, sampling_del = 0.1_real64)
    call MatrixWrite(M=f_latt(:,:,1), nam = 'ft_rect_lattice')
    
end program rect_lattice_build