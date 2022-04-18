program printName
    implicit none
    character (len=15) :: first_name
    call get_command_argument(1,first_name)
    !first_name = arg
    print "(1x,a)",first_name
end program printName
    
    