program test_gauss_2D_nocorr
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  implicit none
  real(real64), allocatable :: X(:), Y(:), M(:,:,:)
  X = seq(st=-8.0_real64, en=10.0_real64, len=10**3)
  Y = seq(st=-8.0_real64, en=10.0_real64, len=500)
  M = gauss_2D_nocorr(X, Y)
  print *, maxval(M(:,:,1))
end program test_gauss_2D_nocorr