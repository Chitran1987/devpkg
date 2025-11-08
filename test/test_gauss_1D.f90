program test_gauss_1D
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  implicit none
  real(real64), allocatable :: X(:), Y(:)
  X = seq(st=-10.0_real64, en=10.0_real64, len=10**3)
  Y = gauss_1D(x=X, A=1.0_real64, x0=1.5_real64, sig=1.5_real64)
  call plot_gnu(X, Y)
end program test_gauss_1D