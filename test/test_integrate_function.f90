program test_integrate_function
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  implicit none
  real(real64), allocatable :: X(:), Y(:)
  real(real64), allocatable :: res(:,:)
  integer :: n = 1*10**4
  real(real64) :: elap1, elap2
  integer :: start, finish, rate
  allocate(X(n))
  allocate(Y(n))
  X = seq(st=-4*pi, en=4*pi, len=n)
  Y = cos(X)
  call plot_local(X, Y)
  
  !function call 1
  res = integrate_function(X=X, Y=Y, y0=5.0_real64)
  call plot_local(X=res(:,1), Y=res(:,2))
  
  !function call 2
  res = integrate_func(X, Y, 1.0_real64)
  
  call plot_local(X = res(:,1), Y = res(:,2) )
end program test_integrate_function