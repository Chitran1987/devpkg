program test_win_sigmoid
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use lattice_utils
  use fft_utils1
  implicit none
  real(real64), allocatable :: tens(:,:,:), ftens(:,:,:), win(:,:,:), fwin(:,:,:)
  real(real64), allocatable :: X(:), Y(:)
  X = seq(st=-10.0_real64, en=10.0_real64, len=500)
  Y = X
  tens = square_latt_sb(X=X, Y=Y, R_latt=0.5_real64, A=1.0_real64, sig=0.1_real64)
  call MatrixWrite(M=tens(:,:,1), nam='test_win_sig_0.0')
  ftens = fft_2D(tens=tens, sampling_del=0.1_real64)
  call MatrixWrite(M=ftens(:,:,1), nam='test_win_sig_0.1')
  win = window_sigmoid(tens=tens, cent=[5.0_real64, 5.0_real64], k = 10.0_real64, Xspan=0.75_real64, Yspan=0.75_real64)
  call MatrixWrite(M=win(:,:,1), nam='test_win_sig_1.0')
  fwin = fft_2D(tens=win, sampling_del=0.1_real64)
  call MatrixWrite(M=fwin(:,:,1), nam='test_win_sig_1.1')
end program test_win_sigmoid