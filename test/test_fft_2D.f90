program test_fft2D_function
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use fft_utils1
  implicit none
  real(real64), allocatable :: X(:), Y(:), Xsp(:,:), Ysp(:,:)
  real(real64), allocatable :: res(:,:,:), dumm(:,:), final(:,:,:), tens_XY(:,:,:)
  integer :: n, m
  n = 500
  m = 600
  X = seq(st=-4.0_real64*pi, en=4.0_real64, len=n)
  Y = seq(st=-3.0_real64*pi, en=-3.0_real64*pi, len=m)
  tens_XY = grid_2(X,Y)
  allocate(res(m, n, 3))
  res(:,:,2) = tens_XY(:,:,1)
  res(:,:,3) = tens_XY(:,:,2)
  res(:,:,1) = sin(-3.0_real64*res(:,:,2)+4.0_real64*res(:,:,3))
  final = fft_2D(res)
end program test_fft2D_function