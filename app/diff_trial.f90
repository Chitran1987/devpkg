program diff_trial
    use vec_utils
    implicit none
    real(real64), dimension(10) :: X
    real(real64), dimension(9) :: Y
    X = seq(1.0_real64, 10.0_real64, 10)
    print *, X
    print *, diff(X)
end program diff_trial