program line_eval
    use vec_utils
    implicit none
    real(real64) :: X(2,2)
    integer dmp
    X = reshape([1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64], shape = [2, 2])
    dmp = print_mat(X)
    call lineval(X, low = 0.0_real64, hi = 1.0_real64)
    dmp = print_mat(X)
end program line_eval