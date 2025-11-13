program test_fft2D_function
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use fft_utils1
  implicit none
  real(real64), allocatable :: X(:), Y(:), M(:,:,:)
  integer :: d1, d2
  X = seq(st=1.0_real64, en=5.0_real64, len=5)
  Y = seq(st=0.0_real64, en=4.0_real64, len=5)
  M = grid_2(X, Y)
  print *, "The X-spread matrix"
  d1 = print_mat(M(:,:,1))
  print *, "The Y-spread matrix"
  d2 = print_mat(M(:,:,2))
end program test_fft2D_function