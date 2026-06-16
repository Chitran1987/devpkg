program test_interp_poly
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use fft_utils1 
  implicit none
  real(real64) :: X(100), Y(100), new_X(1000), new_Y(1000)
  X = seq(st = 0.0_real64, en=10.0_real64, len=100)
  Y = 3.0_real64*X**1 + 2.7_real64*X**2 + 12.0_real64*X + 1/((3.0_real64*X + 0.2)**(0.5))
  call plot_gnu(X, Y)
  new_X = seq(st=0.0_real64, en=10.0_real64, len=1000)
  new_Y = lin_interp(X, Y, new_X)
  call plot_gnu(new_X, new_Y)
end program test_interp_poly
