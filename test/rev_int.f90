program chk_rev_int
    use vec_utils
    implicit none
    real(real64) :: X(5)
    integer :: Y(5)
    X = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
    print *, "X before rev", X
    Y = [1, 2, 3, 4, 5]
    print *, "Y before rev", Y
    call rev(X)
    print *, "X after rev", X
    call rev(Y)
    print *, "Y after rev", Y
end program chk_rev_int