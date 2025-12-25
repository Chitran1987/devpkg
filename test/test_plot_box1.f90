program test_plot_box1
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use lattice_utils
  use fft_utils1
  implicit none
  real(real64), allocatable :: tens(:,:,:), pl_tens(:,:,:), win(:,:,:), fwin(:,:,:)
  real(real64), allocatable :: X(:), Y(:)
  X = seq(st=-10.0_real64, en=10.0_real64, len=500)
  Y = X
  tens = square_latt_sb(X=X, Y=Y, R_latt=0.5_real64, A=1.0_real64, sig=0.1_real64)
  call MatrixWrite(M=tens(:,:,1), nam='test_plot_box1_0.0')
  pl_tens = plot_box1(img_tens=tens, box_vec= [-2.0_real64, -2.0_real64, 2.0_real64, 2.0_real64])
  call MatrixWrite(M=pl_tens(:,:,1), nam='test_plot_box1_0.1')
end program test_plot_box1