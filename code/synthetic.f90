! This module treats data from population synthesis code and applies
! criterion from Limoges et al. 2015
module synthetic
    

    use math, only: calculateStandartDeviation
    use files, only: getNumberOfLines, &
                     getNewUnit
    use commons, only: INPUT_PATH, &
                       VW_PATH, &
                       UW_PATH, &
                       UV_PATH, &
                       OUTPUT_FORMAT
    use astronomy, only: convertEquatorToGalact, &
                         convertGalacticToXYZ
    use criterion, only: applyLimogesCriterion
    
    implicit none

    public :: treatSynthData


contains

    
    subroutine treatSynthData(limogesCriterionIsUsed)

        logical, intent(in) :: limogesCriterionIsUsed
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
        character(len=1) :: highestCoord
        integer, dimension(3) :: counterOfWD=(/0, 0, 0/)
        real*8, dimension(3) :: vel_sum = (/0.d0, 0.d0, 0.d0/)
        real*8, dimension(3) :: vel_avg, &
                                sigma

        print*, "Data from population synthesis code"

        numberOfWDs = getNumberOfLines(INPUT_PATH)
        allocate(velocityArray(3, numberOfWDs))
        
        open(getNewUnit(unitIn), file = INPUT_PATH, status='old')
        open(getNewUnit(unitVW), file = VW_PATH, status='old')
        open(getNewUnit(unitUW), file = UW_PATH, status='old')
        open(getNewUnit(unitUV), file = UV_PATH, status='old')

        if (limogesCriterionIsUsed .eqv. .true.) then
            print*, "Limoges criterion is used."
            
            do i = 1, numberOfWDs
                read(unitIn, *) inputData
                rightAscension = inputData(10); 
                declination = inputData(11); 
                distanceInKpc = inputData(12); 
                vel_hel = inputData(20:22)
                
                distance = distanceInKpc * 1.d3  ! Converting kpc to pc
              
                ! Converting Right Ascension and Declination to lattitude 
                ! and longitude
                call convertEquatorToGalact(rightAscension, &
                                            declination, &
                                            l, &
                                            b) 
            
                ! Converting lattitude and longitude to X,Y,Z
                coordOfWD = convertGalacticToXYZ(distance, l, b)    
            
                call applyLimogesCriterion(coordOfWD,    vel_hel, &
                                           counterOfWD,  vel_sum, &
                                           velocityArray, highestCoord)
                select case (highestCoord)
                    case("x")
                        write(unitVW, OUTPUT_FORMAT) vel_hel(2), &
                                                     vel_hel(3)
                    case("y")
                        write(unitUW, OUTPUT_FORMAT) vel_hel(1), &
                                                     vel_hel(3)
                    case("z")
                        write(unitUV, OUTPUT_FORMAT) vel_hel(1), &
                                                     vel_hel(2)
                    case default
                        print*, "Error: couldn't determine highest coordinate"
                end select
            end do
        
            call calculateStandartDeviation(vel_sum, &
                                            counterOfWD, &
                                            velocityArray, &
                                            vel_avg, &
                                            sigma)
        else 
            print*, "Limoges criterion is not used"
   
            do i = 1, numberOfWDs
                read(unitIn,*) inputData 
                vel_hel = inputData(20:22)
                ! Summing velocitiy components
                vel_sum = vel_sum + vel_hel
                ! Filling array for calculating SD
                velocityArray(:, i) = vel_hel     
            end do 

            call calculateStandartDeviation(vel_sum, &
                                            numberOfWDs, &
                                            velocityArray, &
                                            vel_avg, &
                                            sigma)
        end if

        write(6,*) 'Average relative to Sun:', vel_avg
        write(6,*) 'Sigmas:                 ', sigma
    
        deallocate(velocityArray)

    end subroutine

    
end module
