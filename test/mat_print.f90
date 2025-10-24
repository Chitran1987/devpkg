program mat_print
    use vec_utils
    implicit none
    real(real64) :: X(2,2)
    integer dmp
    X = reshape([1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64], shape(X))
    dmp = print_mat(X)
end program mat_print