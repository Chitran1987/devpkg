program test_fft1_function
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use fft_utils1
  implicit none
  real(real64), allocatable :: X(:), Y(:)
  real(real64), allocatable :: res(:,:,:)
  integer :: n
  n = 10000
  X = seq(st=-4.0_real64*pi, en=4.0_real64, len=n)
  Y = cos(7.5*X)
  call plot_local(X, Y)
  res = fft_1D(X,Y)
  call plot_oct(X = res(:,1,1), Y = res(:,2,1))
end program test_fft1_function