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
  real(real64) :: k0(4), vec(16)
  logicaL, allocatable :: mask(:,:)

  X = seqn(st=-10.0_real64, len=300, del=0.05_real64) 
  Y = seqn(st=-10.0_real64, len=200, del=0.05_real64)
  print *, 'max value of X is =', maxval(X)
  print *, 'max value of Y is =', maxval(Y)

  !------------build the lattice-----------------------------------------------------!
  latt = rect_latt_sb(X=X, Y=Y, R_latt_x=0.5_real64, R_latt_y=0.8_real64, A = 1.0_real64, sig = 0.1_real64)

  !------------build the lattice-----------------------------------------------------!

  !-----------Disappear exactly one columns of atoms from the lattice----------------!
  mask = (latt(:,:,2) >= 0.8_real64) .and. (latt(:,:,2) <= 1.3_real64)
  where(mask)
    latt(:,:,1) = 0.0_real64
  end where
  call MatrixWrite(M = latt(:,:,1), nam = 'fft_2D_map_test2_0.0')
  !-----------Disappear exactly one columns of atoms from the lattice----------------!

  !----------call the fft_2D function------------------------------------------------!
  f_latt = fft_2D(tens=latt, sampling_del=0.1_real64)
  call MatrixWrite(M=f_latt(:,:,1), nam='fft_2D_map_test2_0.1')
  !----------call the fft_2D function------------------------------------------------!

  !-----------call the fft_2D_map_function-------------------------------------------!
  !-----build the k1st matrix
  allocate(k1st(4,4))
  vec = [-1.0_real64, 11.5_real64, -1.0_real64, -13.5_real64, &
          7.0_real64, -1.0_real64, -9.0_real64, -1.0_real64,  &
          1.0_real64, 13.5_real64, 1.0_real64, -11.5_real64,  &
          9.0_real64, 1.0_real64, -7.0_real64, 1.0_real64]
  k1st = reshape(vec, shape(k1st))
  k0 = [-2.0_real64, -2.0_real64, 2.0_real64, 2.0_real64]
  map_img = fft_2D_map(img_tens = latt, Xspan = 0.5_real64, Yspan = 0.8_real64, k1st = k1st, k0 = k0)
  call MatrixWrite(M=map_img(:,:,1), nam='fft_2D_map_test2_0.2')
  !-----------call the fft_2D_map_function-------------------------------------------!
end program test_fft2D_map_function