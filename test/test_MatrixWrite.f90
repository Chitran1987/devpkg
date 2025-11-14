program test_MatrixWrite
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use fft_utils1
  implicit none
  real(real64), allocatable :: X(:), Y(:), M(:,:,:)
  character(len=:), allocatable :: nam
  X = seq(-10.0_real64, 10.0_real64, 5)
  Y = X**2
  M = grid_2(X, Y)
  call MatrixWrite(M(:,:,1), 'test_MatrixWrite_X')
  call MatrixWrite(M(:,:,2), 'test_MatrixWrite_Y')
end program test_MatrixWrite