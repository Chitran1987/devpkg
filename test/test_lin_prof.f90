program test_gauss_1D
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  implicit none
  real(real64), allocatable :: X(:), Y(:), G1(:,:,:), G2(:,:,:), lin_prof(:,:)
  X = seq(st=-10.0_real64, en=10.0_real64, len=10**3)
  Y = seq(st=-5.0_real64, en=5.0_real64, len=500)
  G1 = gauss_2D_nocorr(X, Y)
  G2 = gauss_2D_nocorr(x=X, y=Y, x0=1.5_real64, y0=-2.0_real64, sig_x = 2.3_real64, sig_y=1.7_real64)
  G1(:,:,1) = G1(:,:,1) + G2(:,:,1)
  lin_prof = lin_prof_v(G1, 1.5_real64)
  call plot_local(X=lin_prof(:,1), Y = lin_prof(:,2))
end program test_gauss_1D