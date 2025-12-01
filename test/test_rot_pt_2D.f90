program rot_2D_test
    use iso_fortran_env, only: real64
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1
    use lattice_utils
    implicit none
    real(real64) :: vec(2), new_vec(2), ang
    vec = [1.0_real64, 1.0_real64]
    !ang = 2.0_real64*pi
    ang = 4.0_real64*pi
    new_vec = rot_pt_2D(pt=vec, alpha=ang)
    print*, new_vec
    
    !print*, pi
end program rot_2D_test