program test_mask_box
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use lattice_utils
  use fft_utils1
  implicit none
  real(real64), allocatable :: tens(:,:,:), f_tens(:,:,:), res_tens(:,:,:), f_tens_amp(:,:,:)
  real(real64), allocatable :: X(:), Y(:)
  real(real64) :: vec(4), del
  X = seq(st=-10.0_real64, en=10.0_real64, len=500)
  Y = X
  del = mean(diff(X))
  tens = square_latt_sb(X=X, Y=Y, R_latt=0.5_real64, A=1.0_real64, sig=0.1_real64)
  f_tens = fft_2D(tens=tens, sampling_del = 0.1_real64)
  call MatrixWrite(M=tens(:,:,1), nam='test_mask_box_0.0')
  call MatrixWrite(M=f_tens(:,:,1), nam='test_mask_box_0.1')
  !Transfer values to f_tens_amp
  allocate(f_tens_amp(size(tens,1), size(tens, 2), 3))
  f_tens_amp(:,:,2:3) = tens(:,:,2:3)
  f_tens_amp(:,:,1) = f_tens(:,:,1)
  vec = [-1.0_real64, -1.0_real64, 1.0_real64, 1.0_real64]
  res_tens = mask_box(tens=f_tens_amp, box_vec=vec)
  call MatrixWrite(M = res_tens(:,:,1), nam = 'test_mask_box_0.2')
end program test_mask_box