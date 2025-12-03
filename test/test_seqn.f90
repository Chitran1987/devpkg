program test_seqn
    use iso_fortran_env, only: real64
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1
    use lattice_utils
    implicit none
    real(real64) :: X(10)
    X = seqn(st=1.0_real64, len=10, del=-1.5_real64)
    print*, X
end program test_seqn