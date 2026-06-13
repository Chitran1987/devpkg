program test_Heaviside_2D
    use iso_fortran_env, only: real64
    use vec_utils
    use stat_utils
    use plot_utils
    use calculus_utils
    use fft_utils1  
    implicit none
    real(real64), allocatable :: X(:), Y(:), XY(:,:,:)
    real(real64), allocatable :: M(:,:,:)
    logical :: r_t(2)
    X = seq(st=0.0_real64, en =10.0_real64, len=100)
    Y = seq(st=0.0_real64, en=10.0_real64, len=100)
    XY = grid_2(X=X, Y=Y)
    r_t = [.false. , .false.]
    M = Heaviside_2D(grid = XY, val = [5.0_real64, 5.0_real64], right_top = r_t)
    call MatrixWrite(M(:,:,1), nam = 'Heaviside_test')
end program test_Heaviside_2D
