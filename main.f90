! This program reads observational data from Limoges et al. 2015 or 
! data from population synthesis code and calculates average 
! velocities, standart deviation and for observational data - 
! characteristics of bolometric magnitudes bins 
program velocities
    use observational, only: treatObservationalData
    use synthlimo, only: treatSynthDataWithCriterion
    use synthnolimo, only: treatSynthDataNoCriterion
  
    implicit none

    logical, parameter :: LIMOGES_CRIT_IS_USED = .false.
    logical, parameter :: OBSERV_DATA_IS_USED = .false.

    if (OBSERV_DATA_IS_USED) then
        call treatObservationalData
    else    
        if (LIMOGES_CRIT_IS_USED) then
            call treatSynthDataWithCriterion
        else
            call treatSynthDataNoCriterion
        end if
    end if
stop
end program
