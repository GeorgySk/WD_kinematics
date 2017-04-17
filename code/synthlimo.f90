! This module treats data from population synthesis code and applies
! criterion from Limoges et al. 2015
module synthlimo

    use files, only: getNumberOfLines, &
                     getNewUnit
    use commons, only: INPUT_PATH, &
                       VW_PATH, &
                       UW_PATH, &
                       UV_PATH, &
                       OUTPUT_FORMAT
    use astronomy, only: convertEquatorToGalact, &
                         convertGalacticToXYZ
    implicit none

    private

    public :: treatSynthDataWithCriterion

contains
    
    subroutine treatSynthDataWithCriterion

        integer :: numberOfWDs
        ! In this array we record velocity components for each WD
        ! in order to calculate standart deviation 
        real*8, dimension(:, :), allocatable :: velocityArray
        ! I/O units
        integer :: unitIn, &  ! Input data
                   unitVW, &  ! Outputs with velocities
                   unitUW, &
                   unitUV
        integer :: i
        real*8, dimension(22) :: inputData
        real*8 :: rightAscension, &  
                  declination, &
                  distanceInKpc
        real*8, dimension(3) :: vel_hel, &  ! Heliocentric velocity
                                coordOfWD
        real*8 :: distance
        real*8 :: l, &  ! Longitude
                  b  ! Lattitude
        real*8 :: highestCoord
        integer, dimension(3) :: counterOfWD=(/0, 0, 0/)
        real*8, dimension(3) :: vel_sum = (/0.d0, 0.d0, 0.d0/)
        real*8, dimension(3) :: vel_avg, &
                                sigma
        real*8, dimension(3) :: sumOfRestsSquared = (/0.d0, 0.d0, 0.d0/)



        print*, "Data from population synthesis code"
        print*, "Limoges criterion is used."

        numberOfWDs = getNumberOfLines(INPUT_PATH)
        allocate(velocityArray(3, numberOfWDs))
    
        open(getNewUnit(unitIn), file = INPUT_PATH, status='old')
        open(getNewUnit(unitVW), file = VW_PATH, status='old')
        open(getNewUnit(unitUW), file = UW_PATH, status='old')
        open(getNewUnit(unitUV), file = UV_PATH, status='old')

        
        do i = 1, numberOfWDs
            read(unitIn, *) inputData
            rightAscension = inputData(10); 
            declination = inputData(11); 
            distanceInKpc = inputData(12); 
            vel_hel = inputData(20:22)
            
            distance = distanceInKpc * 1.d3  ! Converting kpc to pc
          
            ! Converting Right Ascension and Declination to lattitude 
            ! and longitude
            call convertEquatorToGalact(rightAscension, declination, l, b) 
        
            ! Converting lattitude and longitude to X,Y,Z
            coordOfWD = convertGalacticToXYZ(distance, l, b)    
        
            highestCoord = maxval(abs(coordOfWD))
    
            ! Highest is X
            if (abs(highestCoord - abs(coordOfWD(1))) .lt. 1.0d-5) then
                
                ! Nº of WDs that will be taken into account 
                ! for plotting V and W
                counterOfWD(1) = counterOfWD(1) + 1
    
                ! Summing V and W components
                vel_sum(2) = vel_sum(2) + vel_hel(2)  
                vel_sum(3) = vel_sum(3) + vel_hel(3)
    
                ! Filling arrays for calculating SD
                velocityArray(2, counterOfWD(3)+counterOfWD(1)) = vel_hel(2)
                velocityArray(3, counterOfWD(2)+counterOfWD(1)) = vel_hel(3)
    
                ! outputs/vw.dat
                write(unitVW, OUTPUT_FORMAT) vel_hel(2), vel_hel(3)  
            
            ! Highest is Y
            else if (abs( highestCoord - abs(coordOfWD(2) )) .lt. 1.0d-5) then
                
                ! Nº of WDs that will be taken into account 
                ! for plotting U and W
                counterOfWD(2) = counterOfWD(2) + 1  
    
                ! Summing U and W components 
                vel_sum(1) = vel_sum(1) + vel_hel(1) ! U 
                vel_sum(3) = vel_sum(3) + vel_hel(3) ! W
    
                ! Filling arrays for calculating SD
                velocityArray(1, counterOfWD(2)+counterOfWD(3)) = vel_hel(1)
                velocityArray(3,counterOfWD(2)+counterOfWD(1)) = vel_hel(3)
    
                ! outputs/uw.dat
                write(unitUW, OUTPUT_FORMAT) vel_hel(1), vel_hel(3)   
            
            ! Highest is Z
            else if (abs( highestCoord - abs(coordOfWD(3) )) .lt. 1.0d-5) then
    
                ! Nº of WDs that will be taken into account 
                ! for plotting U and V
                counterOfWD(3) = counterOfWD(3) + 1  
    
                ! Summing U and V components 
                vel_sum(1) = vel_sum(1) + vel_hel(1) 
                vel_sum(2) = vel_sum(2) + vel_hel(2)
    
                ! Filling arrays for calculating SD             
                velocityArray(1, counterOfWD(2)+counterOfWD(3)) = vel_hel(1)
                velocityArray(2, counterOfWD(3)+counterOfWD(1)) = vel_hel(2)
    
                !outputs/uv.dat
                write(unitUV, OUTPUT_FORMAT) vel_hel(1), vel_hel(2)  
            
            else
                print*, "Error: couldn't determine highest coordinate."
            
            end if
        end do
    
        ! Calculating average U, V, W
        vel_avg(1) = vel_sum(1) / dfloat(counterOfWD(3)+counterOfWD(2)) 
        vel_avg(2) = vel_sum(2) / dfloat(counterOfWD(3)+counterOfWD(1)) 
        vel_avg(3) = vel_sum(3) / dfloat(counterOfWD(2)+counterOfWD(1)) 
        
        ! Summing rests squared for each velocity component
        do i = 1, counterOfWD(3)+counterOfWD(2)
            sumOfRestsSquared(1) = sumOfRestsSquared(1)&
                                   + (velocityArray(1, i)-vel_avg(1)) ** 2
        end do
        do i = 1, counterOfWD(3)+counterOfWD(1)
            sumOfRestsSquared(2) = sumOfRestsSquared(2)&
                                   + (velocityArray(2, i)-vel_avg(2)) ** 2
        end do
        do i = 1, counterOfWD(2)+counterOfWD(1)
            sumOfRestsSquared(3) = sumOfRestsSquared(3)&
                                   + (velocityArray(3, i)-vel_avg(3)) ** 2
        end do
    
        sigma(1) = (sumOfRestsSquared(1) / dfloat(counterOfWD(3) &
                   + counterOfWD(2))) ** 0.5d0
        sigma(2) = (sumOfRestsSquared(2) / dfloat(counterOfWD(3) &
                   + counterOfWD(1))) ** 0.5d0
        sigma(3) = (sumOfRestsSquared(3) / dfloat(counterOfWD(2) &
                   + counterOfWD(1))) ** 0.5d0

        write(6,*) 'Average relative to Sun:',vel_avg
        write(6,*) 'Sigmas:                 ', sigma

        deallocate(velocityArray)

    end subroutine

end module
