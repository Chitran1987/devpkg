program test_unidistr
    use iso_fortran_env, only: real64
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1  
	implicit none
	real(real64) :: X(10)
	X = unidistr(len =size(X), min_val=10.0_real64, max_val = 20.0_real64)
	print*, X
end program test_unidistr