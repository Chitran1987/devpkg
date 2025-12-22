program test_fft2D_map_function
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use lattice_utils
  use fft_utils1
  implicit none
  real(real64), allocatable :: X(:), Y(:)
  X = seq(st=-10.0_real64, en=10.0_real64, len=500)
  Y = sigmoid(X=X, k=10.0_real64, cutoff=0.0_real64)
  !call writeXY(X=X, Y=Y, nam='sigmoid')
  call plot_gnu(X, Y)
  Y = sigmoid_plat(X=X, k=10.0_real64, left_cut=-8.0_real64, right_cut=8.0_real64)
  call plot_gnu(X, Y)
end program test_fft2D_map_function