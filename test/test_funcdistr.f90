program test_funcdistr
    use iso_fortran_env, only: real64
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1    
    implicit none
    real(real64) :: dat_X(10**3), dat_Y(10**3), scale
    real(real64), allocatable :: samp(:)
    integer :: unit, i
    dat_X = seq(st=-10.0_real64, en = +10.0_real64, len=10**3)
    dat_Y = dat_X**2.0_real64 + 3.0_real64
    call plot_gnu(dat_X, dat_Y)
    samp = funcdistr(len=10**3, min_val=-10.0_real64, max_val=10.0_real64, dat_X=dat_X, dat_Y=dat_Y, scale=1.5_real64)
    print*, samp(1:10)
    print*, size(samp)
    !write the files to text
        open(newunit=unit, file="vector.txt", status="replace", action="write")

    do i = 1, size(samp)
        write(unit, *) samp(i)
    end do

    close(unit)
end program test_funcdistr