program test_fft_latt_err
use iso_fortran_env, only: real64
use vec_utils
use stat_utils
use lattice_utils
use fft_utils1
implicit none
real(real64), allocatable :: X(:), Y1(:), Y2(:), Y(:)
real(real64), allocatable :: M(:,:,:), f_M(:,:,:)
real(real64) :: wind(1,2)
!Build the dataset
X = seqn(st=-10.0_real64, len = 1000, del = 0.01_real64)
Y = seqn(st=-5.0_real64, len = 1000, del = 0.01_real64)
M = square_latt_sb(X = X, Y = Y, R_latt = 0.5_real64, A = 1.0_real64, sig = 0.05_real64)
call MatrixWrite(M=M(:,:,1), nam='fft_2D_err_test')
f_M = fft_2D(tens = M, sampling_del = 0.01_real64)
call MatrixWrite(M=f_M(:,:,1), nam= 'fft_2D_err_test_fft')
end program test_fft_latt_err