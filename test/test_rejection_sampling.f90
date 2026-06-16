program test_rejection_sampling
    use iso_fortran_env, only: real64
    !use Statdistr
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1
    implicit none
    real(real64) :: X(10**4), Y(10**4), Y1(1000)
    integer :: i, unit
    X = seq(st=-10.0_real64, en=10.0_real64, len=10**4)
    Y = gauss_1D(x=X, A=1.0_real64, x0=0.0_real64, sig=0.5_real64)
    call plot_gnu(X, Y)
    Y1 = rejection_sampling(len=1000, min_val=-10.0_real64, max_val=10.0_real64, dat_X=X, dat_Y=Y)
    print*, Y1(1:10)
    !write the files to text
        open(newunit=unit, file="vector1.txt", status="replace", action="write")

    do i = 1, size(Y1)
        write(unit, *) Y1(i)
    end do
end program test_rejection_sampling
