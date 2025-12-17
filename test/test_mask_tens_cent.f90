program test_mask_tens_cent
use iso_fortran_env, only: real64
use vec_utils
use stat_utils
use lattice_utils
use fft_utils1
implicit none
real(real64), allocatable :: tens(:,:,:), tens_f(:,:,:), tens_mask(:,:,:), tens_mask_f(:,:,:), X(:), Y(:)
real(real64) :: center(2), spanX, spanY
!Build the sequences
X = seq(st = -10.0_real64, en = 10.0_real64, len = 500)
Y = X
!Build the lattice
tens = square_latt_sb(X=X, Y=Y, R_latt=1.0_real64, A=1.0_real64, sig = 0.2_real64)
call MatrixWrite(M=tens(:,:,1), nam='test_mask_tens')
tens_f = fft_2D(tens=tens, sampling_del=0.1_real64)
call MatrixWrite(M=tens_f(:,:,1), nam='test_mask_tens_f')
!xlim = [-5.0_real64, 5.0_real64]
!ylim = [-4.0_real64, 4.0_real64]
center = [1.0_real64, -1.0_real64]
tens_mask = mask_tens_cent(tens = tens,cent=center, Xspan = 2.0_real64, Yspan = 2.0_real64 )
call MatrixWrite(M=tens_mask(:,:,1), nam='test_mask_tens_masked')
tens_mask_f = fft_2D(tens=tens_mask, sampling_del=0.1_real64)
call MatrixWrite(M=tens_mask_f(:,:,1), nam='test_mask_tens_masked_f')
end program test_mask_tens_cent