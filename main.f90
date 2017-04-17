! This program reads observational data from Limoges et al. 2015 or 
! data from population synthesis code and calculates average 
! velocities, standart deviation and for observational data - 
! characteristics of bolometric magnitudes bins 
program velocities
    use observational, only: treatObservationalData
    use synthlimo, only: treatSynthDataWithCriterion
    use synthnolimo, only: treatSynthDataNoCriterion
  
    implicit none

    integer :: num_args  
    character(len = 12), dimension(:), allocatable :: args

    num_args = command_argument_count()
    if (num_args .eq. 0) then
        print*, "Please specify what data you want to process:"
        call printHelp
    else
        allocate(args(num_args))
        call get_command_argument(1, args(1))
        if (args(1) .eq. "-o") then
            call treatObservationalData
        else if (args(1) .eq. "-s") then
            if (num_args .eq. 1) then
                call treatSynthDataNoCriterion
            else
                call get_command_argument(2, args(2))
                if (args(2) .eq. "-l") then
                    call treatSynthDataWithCriterion
                else 
                    print*, "Wrong second argument"
                    call printHelp
                end if
            end if
        else 
            print*, "Wrong first argument"
            call printHelp
        end if
    end if

stop
end program

subroutine printHelp()
    implicit none
    print*, "   -o - for observational data"
    print*, "   -s - for data from population synthesis code without &
                     &Limoges criterion"
    print*, "   -s -l - for data from population synthesis code with &
                        &applied criterion of Limoges"
end subroutine printHelp