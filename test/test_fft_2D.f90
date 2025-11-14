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
  m = 500
  X = seq(st=-5.0_real64*pi, en=5.0_real64, len=n)
  Y = seq(st=-5.0_real64*pi, en=5.0_real64*pi, len=m)
  tens_XY = grid_2(X,Y)
  allocate(res(m, n, 3))
  res(:,:,2) = tens_XY(:,:,1)
  res(:,:,3) = tens_XY(:,:,2)
  res(:,:,1) = sin(-1.0_real64*res(:,:,2)+ 0.5_real64*res(:,:,3))
  call MatrixWrite(res(:,:,1), 'lattice')
  final = fft_2D(res)
  call MatrixWrite(final(:,:,1), 'ft_lattice')
end program test_fft2D_function