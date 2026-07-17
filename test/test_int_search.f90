program test_int_search
    use iso_fortran_env, only: real64
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1
    use lattice_utils 
    implicit none

    integer, allocatable :: M(:,:), dmp
    M = int_search(5)
    dmp = print_mat(M)
end program test_int_search