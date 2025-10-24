program mat_cmplx_prnt
    use vec_utils
    implicit none
    complex(real64) :: X(2, 2)
    integer dmp
    X = reshape([(1.0_real64, 1.0_real64), (2.0_real64, 2.0_real64), (3.0_real64, 3.0_real64), (4.0_real64, 4.0_real64)], shape=[2,2])
    dmp = print_mat(X)
end program mat_cmplx_prnt