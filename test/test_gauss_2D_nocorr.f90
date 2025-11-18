program test_gauss_2D_nocorr
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  implicit none
  real(real64), allocatable :: X(:), Y(:), M(:,:,:)
  X = seq(st=-10.0_real64, en=10.0_real64, len=10**3)
  Y = seq(st=-10.0_real64, en=10.0_real64, len=1000)
  M = gauss_2D_nocorr(X = X, Y = Y, x0 = -1.5_real64, y0 = 7.7_real64)
  call MatrixWrite(M=M(:,:,1), nam='gauss_2D_nocorr_test_mat')
end program test_gauss_2D_nocorr