module devpkg
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, devpkg!"
  end subroutine say_hello
end module devpkg
