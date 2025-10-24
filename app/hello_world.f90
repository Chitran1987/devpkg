program name
    implicit none
    character(len = 50) :: nm1, nm2
    character(len = 50) :: age
    print *, "What is your first name ? "
    read *, nm1
    print *, "What is your last name ?"
    read *, nm2
    print *, "what is your age ?"
    read *, age
    print *, "Hello ",trim(nm1),". Your age is ",trim(age) 
end program name