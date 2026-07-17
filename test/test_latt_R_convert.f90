program test_latt_conv
    use iso_fortran_env, only: real64
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1
    use lattice_utils 
    implicit none

    real(real64) :: R1(2), R2(2)
    real(real64) :: M(2,2)
    integer :: dmp
    R1 = [2.0_real64, 3.0_real64]
    R2 = [3.0_real64, 3.0_real64]
    M = latt_R_convert(R1, R2, 5)
    dmp = print_mat(M)
end program test_latt_conv