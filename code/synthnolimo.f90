! This module treats data from population synthesis code
! No criterion is applied here
module synthnolimo

    use files, only: getNumberOfLines, &
                     getNewUnit
    use commons, only: INPUT_PATH, &
                       UW_PATH, &
                       VW_PATH, &
                       UV_PATH
    implicit none

    private

    public :: treatSynthDataNoCriterion

contains
    subroutine treatSynthDataNoCriterion
        integer :: numberOfWDs
        ! In this array we record velocity components for each WD
        ! in order to calculate standart deviation 
        real*8, dimension(:, :), allocatable :: velocityArray
        integer :: unitIn, &  ! Input data
                   unitVW, &  ! Outputs with velocities
                   unitUW, &
                   unitUV
        integer :: i
        real*8, dimension(22) :: inputData
        real*8, dimension(3) :: vel_hel, &  ! Heliocentric velocity
                                vel_sum, &
                                vel_avg, &
                                sigma
        real*8, dimension(3) :: sumOfRestsSquared = (/0.d0, 0.d0, 0.d0/)

        print*, "Data from population synthesis code"
        print*, "Limoges criterion is not used"

        numberOfWDs = getNumberOfLines(INPUT_PATH)
        allocate(velocityArray(3, numberOfWDs))
    
        open(getNewUnit(unitIn), file = INPUT_PATH, status='old')
        open(getNewUnit(unitVW), file = VW_PATH, status='old')
        open(getNewUnit(unitUW), file = UW_PATH, status='old')
        open(getNewUnit(unitUV), file = UV_PATH, status='old')

        do i = 1, numberOfWDs
            read(unitIn,*) inputData 
            vel_hel = inputData(20:22)
            ! Summing velocitiy components
            vel_sum = vel_sum + vel_hel
            ! Filling array for calculating SD
            velocityArray(:, i) = vel_hel     
        end do 
    
        vel_avg = vel_sum / dfloat(numberOfWDs)    
        ! Summing rests squared for each velocity component
        do i = 1, numberOfWDs
            sumOfRestsSquared = sumOfRestsSquared &
                                + (velocityArray(:,i)-vel_avg) ** 2
        end do    
        sigma = (sumOfRestsSquared / dfloat(numberOfWDs)) ** 0.5
    
        write(6,*) 'Average relative to Sun:',vel_avg
        write(6,*) 'Sigmas:                 ', sigma

        deallocate(velocityArray)
    end subroutine
end module
