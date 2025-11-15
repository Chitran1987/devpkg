program hex_lattice_build
    use iso_fortran_env, only: real64
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1
    implicit none
    real(real64), allocatable :: X(:), Y(:), latt(:,:,:), grid(:,:,:) !Lattice variables
    real(real64), allocatable :: f_latt(:,:,:)
    !build the square lattice
    X = seq(st=-10.0_real64, en=10.0_real64, len=500)
    Y = X
    grid = grid_2(X, Y)
    allocate(latt(size(Y), size(X), 3))
    latt(:,:,2) = grid(:,:,1)
    latt(:,:,3) = grid(:,:,2)    
    latt(:,:,1) = sin(3.0_real64*(1*latt(:,:,2)+0*latt(:,:,3)))
    latt(:,:,1) = latt(:,:,1) + sin(3.0_real64*(0.5_real64*latt(:,:,2) + (sqrt(3.0_real64)/2.0_real64)*latt(:,:,3)))
    latt(:,:,1) = latt(:,:,1) + sin(3.0_real64*(-0.5_real64*latt(:,:,2) + (sqrt(3.0_real64)/2.0_real64)*latt(:,:,3)))
    call MatrixWrite(M=latt(:,:,1), nam='lattice')

    !build the Fourier Transform
    f_latt = fft_2D(tens = latt)
    call MatrixWrite(M=f_latt(:,:,1), nam = 'ft_lattice')
    
end program hex_lattice_build