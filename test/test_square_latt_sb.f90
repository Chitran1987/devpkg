program test_fft2D_function
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use fft_utils1
  use lattice_utils
  implicit none
  real(real64), allocatable :: X(:), Y(:), Xsp(:,:), Ysp(:,:), tens(:,:,:), ftens(:,:,:)
  X = seq(st = -10.0_real64, en = 10.0_real64, len = 1000)
  Y = X
  tens = square_latt_sb(X = X, Y = Y, R_latt = 1.5_real64, A = 1.0_real64, sig = 0.25_real64)
  ftens = fft_2D(tens)
  call MatrixWrite(tens(:,:,1), 'sq_latt')
  call MatrixWrite(ftens(:,:,1), 'sq_latt_f')
end program test_fft2D_function