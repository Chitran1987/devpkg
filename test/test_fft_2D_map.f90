program test_fft2D_map_function
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use lattice_utils
  use fft_utils1
  implicit none
  real(real64), allocatable :: X(:), Y(:), latt(:,:,:), f_latt(:,:,:), zp_latt(:,:,:), f_latt1(:,:,:), map_img(:,:,:) !needed
  real(real64), allocatable :: mask_dummy(:,:,:) , k1st(:,:)
  real(real64) :: k0(4)
  X = seqn(st=-10.0_real64, len=200, del=0.1_real64) 
  Y = seqn(st=-10.0_real64, len=100, del=0.1_real64)
  print *, maxval(X)
  print *, maxval(Y)
  latt = square_latt_sb(X=X, Y=Y, R_latt=0.75_real64, A=1.0_real64, sig=0.25_real64)
  call MatrixWrite(M=latt(:,:,1), nam='test_fft2D_1')
  f_latt = fft_2D(tens=latt, sampling_del=0.1_real64)
  call MatrixWrite(M=f_latt(:,:,1), nam='test_fft2D_2')
  zp_latt = zero_pad_tens(tens=latt)
  call MatrixWrite(M=zp_latt(:,:,1), nam='test_fft_2D_3')
  f_latt1 = fft_2D(tens=zp_latt, sampling_del=0.1_real64)
  call MatrixWrite(M=f_latt1(:,:,1), nam='test_fft_2D_4')
  !Define the Xspan and Yspan
  mask_dummy = mask_tens_cent(tens=zp_latt, Xspan=2.0_real64, Yspan=2.0_real64, cent=[0.0_real64, 8.0_real64] )
  call MatrixWrite(M=mask_dummy(:,:,1), nam='test_fft_2D_5')
  f_latt1 = fft_2D(tens=mask_dummy, sampling_del=0.1_real64)
  call MatrixWrite(M=f_latt1(:,:,1), nam='test_fft2D_6')
  !Define k1st
  allocate(k1st(4,4))
  k1st = reshape([-2.0_real64, 6.5_real64, -2.0_real64, -10.5_real64, -6.5_real64, -2.0_real64, -10.5_real64, -2.0_real64, 2.0_real64, 10.5_real64, 2.0_real64, -6.5_real64, 10.5_real64, 2.0_real64, -6.5_real64, 2.0_real64], shape(k1st))
  k0 = [-4.0_real64, -4.0_real64, 4.0_real64, 4.0_real64]
  print *, 'mapping started'
  map_img = fft_2D_map(img_tens=zp_latt, Xspan = 2.0_real64, Yspan = 2.0_real64, k1st = k1st, k0 = k0 )
  call MatrixWrite(M = map_img(:,:,1), nam='fft_2D_7')
end program test_fft2D_map_function