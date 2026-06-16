program test_lin_interp
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use fft_utils1 
  implicit none
  real(real64) :: X(10**1), Y(10**1), new_X(10**2), new_Y(10**2)
  X = seq(st=0.0_real64, en=10.0_real64, len=10**1)
  Y = X**2
  call plot_gnu(X, Y)
  new_X = seq(st=0.0_real64, en=10.0_real64, len=10**2)
  new_Y = lin_interp(dat_X = X, dat_Y = Y, X = new_X)
  call plot_gnu(new_X, new_Y)
end program test_lin_interp
