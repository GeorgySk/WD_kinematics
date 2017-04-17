! This program reads observational data from Limoges et al. 2015 or 
! data from population synthesis code and calculates average 
! velocities, standart deviation and for observational data - 
! characteristics of bolometric magnitudes bins 
program velocities
use math, only: multiplyMatrixByVector
use astronomy, only: CTK, &
                     TRAN_MATR, &
                     convertGalacticToXYZ, &
                     convertEquatorToGalact, &
                     convertHoursToRad, &
                     convertDegreesToRad
use files, only: getNumberOfLines, &
                 getNewUnit
use observational, only: treatObservationalData
use synthlimo, only: treatSynthDataWithCriterion
  
implicit none

logical, parameter :: LIMOGES_CRIT_IS_USED = .true.
logical, parameter :: OBSERV_DATA_IS_USED = .false.
character(len = *), parameter :: INPUT_PATH = '/home/gemma/Documents/program&
                                              &/WD_population_40pc/output_data&
                                              &/boot_rowell_thin_1.out'
character(len = *), parameter :: VW_PATH = './outputs/vw.dat'
character(len = *), parameter :: UW_PATH = './outputs/uw.dat'
character(len = *), parameter :: UV_PATH = './outputs/uv.dat'
character(len = *), parameter :: MBOL_CLOUD_PATH = './outputs/mbol_cloud.dat'
character(len = *), parameter :: MBOL_AVG_PATH = './outputs/mbol_avg.dat'
character(len = *), parameter :: OBSERV_PATH = './inputs/observational.dat'
character(len = *), parameter :: OUTPUT_FORMAT = '(2(f12.6,3x))'
character(len = *), parameter :: CLOUD_FORMAT = '(4(f12.6,3x))'
character(len = *), parameter :: MBOL_AVG_FORMAT = '(7(f12.6,3x))'
! Binning parameters for bolometric magnitude
real*8, parameter :: MBOL_MIN = 5.75d0                                      
real*8, parameter :: MBOL_MAX = 20.75d0
real*8, parameter :: MBOL_INC = 0.5d0
integer, parameter :: NUM_OF_BINS = int((MBOL_MAX - MBOL_MIN) / MBOL_INC)
real*8, parameter :: VEL_RAD = 0.d0 ! Radial velocity  
! I/O units
integer unitIn, unitVW, unitUW, unitUV
real*8, dimension(3) :: vel_sum = (/0.d0, 0.d0, 0.d0/)
real*8, dimension(3) :: sumOfRestsSquared = (/0.d0, 0.d0, 0.d0/)
real*8, dimension(3) :: vel_hel, coordOfWD, vel_avg, sigma
integer, dimension(3) :: counterOfWD=(/0, 0, 0/)
real*8, dimension(22) :: inputData
real*8 distanceInKpc, distance, rightAscension, declination, l, b, &
       highestCoord
integer i, numberOfWDs
real*8, dimension(:, :), allocatable :: velocityArray
!_______________________________________________________________________________   

if (OBSERV_DATA_IS_USED) then
    call treatObservationalData
else    
    if (LIMOGES_CRIT_IS_USED) then
        call treatSynthDataWithCriterion
    else
        print*, "Data from population synthesis code"

        numberOfWDs = getNumberOfLines(INPUT_PATH)
        allocate(velocityArray(3, numberOfWDs))
    
        open(getNewUnit(unitIn), file = INPUT_PATH, status='old')
        open(getNewUnit(unitVW), file = VW_PATH, status='old')
        open(getNewUnit(unitUW), file = UW_PATH, status='old')
        open(getNewUnit(unitUV), file = UV_PATH, status='old')
        print*, "Limoges criterion is not used"
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
    
        sigma = (sumOfRestsSquared/dfloat(numberOfWDs)) ** 0.5
    
        write(6,*) 'Average relative to Sun:',vel_avg
        write(6,*) 'Sigmas:                 ', sigma

        deallocate(velocityArray)
    end if
    

end if

! deallocate(velocityArray)
    
stop
end program
