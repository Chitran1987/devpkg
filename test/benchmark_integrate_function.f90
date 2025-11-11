program benchmark_integrate_function
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
  !call plot_local(X, Y)
  call system_clock(count = start, count_rate = rate) !system_call is a subroutine call which changes the values of start, rate and finish. It is not equivalent to functional assignment
  res = integrate_function(X=X, Y=Y, y0=5.0_real64)
  call system_clock(count = finish)
  elap1 = real(finish - start, real64)/real(rate, real64)
  print *, 'elap1 is = ', elap1
  call system_clock(count = start, count_rate = rate)
  res = integrate_func(X, Y, 1.0_real64)
  call system_clock(count = finish)
  elap2 = real(finish - start, real64)/real(rate, real64)
  print *, 'elap2 is = ', elap2
  !call plot_local(X = res(:,1), Y = res(:,2) )
end program benchmark_integrate_function