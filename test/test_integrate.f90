program test_integrate
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  implicit none
  real(real64), allocatable :: X(:), Y(:)
  real(real64) :: res
  integer :: n = 2*10**3
  allocate(X(n))
  allocate(Y(n))
  X = seq(st=-10.0_real64, en=10.0_real64, len=n)
  Y = X**2
  res = integrate(X=X, Y=Y, xmin=0.0_real64, xmax=0.1_real64, riemann = .false.)
  print *, 'The result of the integration is ', res
end program test_integrate