program printName
    implicit none
    character (len=15) :: first_name,the_name
    call get_command_argument(1,first_name)
    call get_command_argument(3,the_name)
    !first_name = arg
    print "(1x,a)",first_name
    print "(1x,a)",the_name
end program printName
    
    