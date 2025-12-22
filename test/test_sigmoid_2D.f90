program test_sigmoid_2D
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use lattice_utils
  use fft_utils1
  implicit none
  real(real64), allocatable :: X(:), Y(:), XY(:,:,:), M(:,:)
  X = seq(st=-10.0_real64, en=10.0_real64, len=500)
  Y = X
  XY = grid_2(X,Y)
  M = sigmoid_2D(x=XY(:,:,1), y=XY(:,:,2), k=10.0_real64, x_lo=-15.0_real64, y_lo=-9.0_real64, x_hi=9.0_real64, y_hi=15.0_real64)
  call MatrixWrite(M=M, nam='test_sigmoid_2D')
end program test_sigmoid_2D