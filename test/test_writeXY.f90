program test_fft1_function
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use fft_utils1
  implicit none
  real(real64), allocatable :: X(:), Y(:)
  character(len=:), allocatable :: nam
  X = seq(-10.0_real64, 10.0_real64, 10**3)
  Y = X**2
  call writeXY(X = X, Y = Y, nam = 'test_writeXY')
end program test_fft1_function