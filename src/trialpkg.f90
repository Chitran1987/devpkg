module trialpkg
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, trialpkg!"
  end subroutine say_hello
end module trialpkg
