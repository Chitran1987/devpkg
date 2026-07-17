module devpkg
  use iso_fortran_env, only: real64
  use vec_utils
  use stat_utils
  use plot_utils
  use calculus_utils
  use fft_utils1
  use lattice_utils 
  implicit none
  private :: say_hello
contains
  subroutine say_hello
    print *, "Hello, devpkg!"
  end subroutine say_hello
end module devpkg
