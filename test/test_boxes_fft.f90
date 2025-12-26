program test_plot_boxes_fft
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use lattice_utils
  use fft_utils1
  implicit none
  real(real64), allocatable :: tens(:,:,:), f_tens(:,:,:), f_tens_amp(:,:,:), pl_tens(:,:,:), win(:,:,:), fwin(:,:,:)
  real(real64), allocatable :: X(:), Y(:)
  real(real64), allocatable :: b_vec(:), b_mat(:,:), del
  X = seq(st=-10.0_real64, en=10.0_real64, len=500)
  Y = X
  del = mean(diff(X))
  tens = square_latt_sb(X=X, Y=Y, R_latt=0.5_real64, A=1.0_real64, sig=0.1_real64)
  f_tens = fft_2D(tens = tens, sampling_del = 0.1_real64)
  call MatrixWrite(M=f_tens(:,:,1), nam='test_plot_boxes_0.0')
  !Write the box_mat box
  b_vec = [11.5_real64, -1.0_real64, -1.0_real64, 11.5_real64, 13.5_real64, 1.0_real64, 1.0_real64, 13.5_real64]
  allocate(b_mat(2,4))
  b_mat = reshape(b_vec, shape(b_mat))
  allocate(f_tens_amp(size(f_tens, 1), size(f_tens, 2), 3))
  f_tens_amp(:,:,2:3) = f_tens(:,:,3:4)
  f_tens_amp(:,:,1) = f_tens(:,:,1)
  pl_tens = plot_boxes(img_tens = f_tens_amp, box_mat=b_mat, box_thick = del/2.0_real64, box_if = 0.1_real64)
  call MatrixWrite(M=pl_tens(:,:,1), nam='test_plot_boxes_0.1')
end program test_plot_boxes_fft